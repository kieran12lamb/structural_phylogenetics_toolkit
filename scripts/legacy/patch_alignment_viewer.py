#!/usr/bin/env python3
"""Patch build_dynamic_interactive_tree.py to add:
1. script tag for alignments_data.js in <head>
2. Tree container split into #treeCanvasWrapper + bottom #alignmentDrawer
3. Alignment viewer toolbar, controls, ruler, matrix canvas, consensus track, and taxa column
4. Full JavaScript alignment viewer engine with virtualized canvas rendering, ClustalX/Zappo/3Di color schemes,
   synchronized tree leaf sorting, metadata filter synchronization, and interactive cross-highlighting.
"""

import json
from pathlib import Path

target_path = Path("/Users/kieran.lamb/Github/Skills_Hackathon/scripts/build_dynamic_interactive_tree.py")
with open(target_path, "r", encoding="utf-8") as f:
    code = f.read()

# 1. ADD SCRIPT TAG FOR alignments_data.js in <head>
old_head_scripts = """  <script src="ca_500_structures.js"></script>
  <script src="ca_1193_structures.js"></script>
  <script src="ca_structures.js"></script>"""

new_head_scripts = """  <script src="ca_500_structures.js"></script>
  <script src="ca_1193_structures.js"></script>
  <script src="ca_structures.js"></script>
  <script src="alignments_data.js"></script>"""

assert old_head_scripts in code, "old_head_scripts not found"
code = code.replace(old_head_scripts, new_head_scripts, 1)

# 2. UPDATE treeContainer OPENING TAG
old_tree_open = """    <!-- MAIN SVG TREE CANVAS & FLOATING OVERLAYS -->
    <div id="treeContainer" class="flex-1 h-full relative overflow-hidden">"""

new_tree_open = """    <!-- MAIN WORKSPACE COLUMN (TREE CANVAS + BOTTOM ALIGNMENT DRAWER) -->
    <div id="treeContainer" class="flex-1 h-full flex flex-col relative overflow-hidden">
      
      <!-- TOP: TREE VIEWPORT AREA -->
      <div id="treeCanvasWrapper" class="flex-1 w-full h-full relative overflow-hidden">"""

assert old_tree_open in code, "old_tree_open not found"
code = code.replace(old_tree_open, new_tree_open, 1)

# 3. INSERT BOTTOM ALIGNMENT DRAWER AFTER FLOATING SEARCH BAR
old_tree_close = """      <!-- FLOATING SEARCH BAR (Top Right) -->
      <div class="absolute top-3 right-4 z-20 w-80">
        <div class="relative">
          <input type="text" id="searchInput" placeholder="Search metadata (Class, Binding, Taxon)..." oninput="handleSearch(this.value)" class="w-full pl-8 pr-3 py-1.5 rounded-xl bg-[var(--card-bg)] border border-[var(--border-color)] text-xs text-[var(--text-main)] placeholder-[var(--text-muted)] shadow-xl focus:outline-none focus:border-sky-400 backdrop-blur-md">
          <span class="absolute left-2.5 top-2 text-[var(--text-muted)] text-xs">🔍</span>
        </div>
      </div>

    </div>"""

alignment_drawer_html = """      <!-- FLOATING SEARCH BAR (Top Right) -->
      <div class="absolute top-3 right-4 z-20 w-80">
        <div class="relative">
          <input type="text" id="searchInput" placeholder="Search metadata (Class, Binding, Taxon)..." oninput="handleSearch(this.value)" class="w-full pl-8 pr-3 py-1.5 rounded-xl bg-[var(--card-bg)] border border-[var(--border-color)] text-xs text-[var(--text-main)] placeholder-[var(--text-muted)] shadow-xl focus:outline-none focus:border-sky-400 backdrop-blur-md">
          <span class="absolute left-2.5 top-2 text-[var(--text-muted)] text-xs">🔍</span>
        </div>
      </div>

      </div>

      <!-- BOTTOM: DOCKED ALIGNMENT VIEWER DRAWER -->
      <div id="alignmentDrawer" class="w-full border-t border-[var(--border-color)] bg-[var(--card-bg)] shrink-0 flex flex-col transition-all duration-200 z-30 shadow-2xl overflow-hidden" style="height: 240px;">
        <!-- DRAWER TOOLBAR -->
        <div class="h-9 px-3 border-b border-[var(--border-color)] bg-[var(--panel-bg)] flex items-center justify-between text-xs shrink-0 select-none">
          <!-- Left: Title & Summary -->
          <div class="flex items-center space-x-2 truncate">
            <span class="text-sm">🧬</span>
            <span class="font-bold text-[var(--text-main)] text-[11px] tracking-tight">Multiple Alignment Viewer</span>
            <span id="msaSummaryBadge" class="text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-sky-500/10 text-sky-400 border border-sky-500/30">1,193 seqs &bull; 533 cols</span>
          </div>

          <!-- Center: Mode & Color Schemes -->
          <div class="flex items-center space-x-2">
            <!-- Mode Toggle: AA vs 3Di -->
            <div class="flex bg-[var(--input-bg)] p-0.5 rounded border border-[var(--border-color)] text-[10px]">
              <button id="btnMsaModeAA" onclick="setMsaMode('aa')" class="px-2 py-0.5 rounded font-bold bg-sky-500 text-white transition shadow-sm cursor-pointer">
                🥩 Amino Acid (AA)
              </button>
              <button id="btnMsaMode3Di" onclick="setMsaMode('3di')" class="px-2 py-0.5 rounded font-medium text-[var(--text-muted)] hover:text-white transition cursor-pointer">
                🧊 3Di Structural
              </button>
            </div>

            <!-- Color Scheme Selector -->
            <select id="msaColorSelect" onchange="setMsaColorScheme(this.value)" class="bg-[var(--input-bg)] border border-[var(--border-color)] rounded px-2 py-0.5 text-[10.5px] text-[var(--text-main)] focus:outline-none focus:border-sky-400 font-medium">
              <option value="clustal">🎨 ClustalX Colors</option>
              <option value="zappo">🌈 Zappo (Physicochemical)</option>
              <option value="hydro">💧 Hydrophobicity</option>
              <option value="identity">🎯 Conservation Identity</option>
            </select>

            <!-- Row Ordering Selector -->
            <select id="msaSortSelect" onchange="setMsaSort(this.value)" class="bg-[var(--input-bg)] border border-[var(--border-color)] rounded px-2 py-0.5 text-[10.5px] text-[var(--text-main)] focus:outline-none focus:border-sky-400 font-medium">
              <option value="tree" selected>🌿 Tree Order (Synchronized)</option>
              <option value="selected">⭐ Selected First</option>
              <option value="alpha">🔤 Alphabetical</option>
            </select>
          </div>

          <!-- Right: Navigation, Zoom & Window Controls -->
          <div class="flex items-center space-x-2 text-[10px]">
            <!-- Jump to Position -->
            <div class="flex items-center space-x-1 bg-[var(--input-bg)] px-1.5 py-0.5 rounded border border-[var(--border-color)]">
              <span class="text-[var(--text-muted)]">Col:</span>
              <input type="number" id="msaPosInput" min="1" max="533" value="1" onchange="jumpMsaToPosition(this.value)" class="w-12 bg-transparent text-center text-[var(--text-main)] font-mono font-bold focus:outline-none text-[10px]">
              <span id="msaMaxColLabel" class="text-[var(--text-muted)] font-mono">/ 533</span>
            </div>

            <!-- Zoom Buttons -->
            <div class="flex items-center space-x-0.5 bg-[var(--input-bg)] p-0.5 rounded border border-[var(--border-color)]">
              <button onclick="zoomMsa(-2)" class="w-5 h-5 rounded hover:bg-slate-500/20 text-[var(--text-main)] font-bold flex items-center justify-center transition cursor-pointer" title="Zoom Out">&minus;</button>
              <span id="msaZoomBadge" class="px-1 text-[9.5px] font-mono text-[var(--text-muted)]">14px</span>
              <button onclick="zoomMsa(2)" class="w-5 h-5 rounded hover:bg-slate-500/20 text-[var(--text-main)] font-bold flex items-center justify-center transition cursor-pointer" title="Zoom In">+</button>
            </div>

            <div class="w-px h-3.5 bg-[var(--border-color)]"></div>

            <!-- Height Toggle -->
            <button id="btnMsaTall" onclick="toggleMsaHeight()" class="px-1.5 py-0.5 rounded hover:bg-slate-500/20 text-[var(--text-muted)] hover:text-white transition cursor-pointer" title="Toggle Height (Compact / Standard / Tall)">
              <span id="msaHeightIcon">↕ 240px</span>
            </button>

            <!-- Collapse / Expand Drawer -->
            <button id="btnMsaMinimize" onclick="toggleMsaDrawer()" class="px-2 py-0.5 rounded bg-slate-800 hover:bg-slate-700 text-[var(--text-main)] font-semibold transition flex items-center space-x-1 border border-[var(--border-color)] cursor-pointer" title="Minimize / Expand Alignment Viewer">
              <span id="msaDrawerIcon">▼</span>
              <span id="msaDrawerLabel">Min</span>
            </button>
          </div>
        </div>

        <!-- DRAWER CONTENT: RESIDUE MATRIX & TAXA LABELS -->
        <div id="msaViewport" class="flex-1 flex overflow-hidden relative select-none">
          <!-- Frozen Left Column: Taxa Names -->
          <div id="msaTaxaColumn" class="w-52 border-r border-[var(--border-color)] bg-[var(--card-bg)] shrink-0 flex flex-col overflow-hidden">
            <!-- Taxa Column Header -->
            <div class="h-6 px-3 border-b border-[var(--border-color)] bg-[var(--panel-bg)] flex items-center justify-between text-[9px] font-semibold text-[var(--text-muted)] uppercase tracking-wider shrink-0">
              <span>Taxon ID / Design</span>
              <span id="msaVisibleCount" class="font-mono">1,193</span>
            </div>
            <!-- Taxa Rows Container (synchronized scroll) -->
            <div id="msaTaxaList" class="flex-1 overflow-hidden relative custom-scroll" onwheel="onMsaWheel(event)">
            </div>
            <!-- Consensus Label Row -->
            <div class="h-7 px-3 border-t border-[var(--border-color)] bg-[var(--panel-bg)] flex items-center justify-between text-[10px] font-bold text-sky-400 shrink-0">
              <span>Consensus (Top 1)</span>
              <span id="msaConsensusScore" class="font-mono text-[9px] text-emerald-400">-</span>
            </div>
          </div>

          <!-- Right Center Column: Sequence Grid & Ruler -->
          <div id="msaGridContainer" class="flex-1 flex flex-col overflow-hidden relative bg-slate-950">
            <!-- Residue Coordinate Ruler -->
            <div id="msaRulerWrapper" class="h-6 border-b border-[var(--border-color)] bg-slate-900/90 overflow-hidden relative shrink-0">
              <canvas id="msaRulerCanvas" class="block w-full h-full"></canvas>
            </div>

            <!-- Virtualized Residue Matrix Canvas -->
            <div id="msaMatrixWrapper" class="flex-1 overflow-hidden relative cursor-crosshair" onwheel="onMsaWheel(event)">
              <canvas id="msaMatrixCanvas" class="block w-full h-full"></canvas>
              <!-- Interactive Cell Highlight Frame -->
              <div id="msaCellHover" class="absolute hidden border-2 border-sky-400 bg-sky-400/20 pointer-events-none rounded transition-all duration-75"></div>
            </div>

            <!-- Bottom Consensus & Conservation Bar Chart -->
            <div id="msaConsensusWrapper" class="h-7 border-t border-[var(--border-color)] bg-slate-900/90 overflow-hidden relative shrink-0">
              <canvas id="msaConsensusCanvas" class="block w-full h-full"></canvas>
            </div>
          </div>

          <!-- Tooltip for Residue Hover -->
          <div id="msaResidueTooltip" class="absolute hidden z-40 px-2.5 py-1.5 rounded-lg border border-[var(--border-color)] bg-[var(--tooltip-bg)] backdrop-blur-md shadow-2xl text-[10px] pointer-events-none space-y-0.5">
            <div class="flex items-center space-x-1 font-bold text-[var(--accent)]">
              <span id="msaTipResidue">-</span>
              <span id="msaTipPos" class="font-mono text-white/80">Pos -</span>
            </div>
            <div class="text-[9px] text-[var(--text-muted)]" id="msaTipTaxon">-</div>
            <div class="text-[9px] text-emerald-400 font-mono" id="msaTipConservation">-</div>
          </div>
        </div>
      </div>

    </div>"""

assert old_tree_close in code, "old_tree_close not found"
code = code.replace(old_tree_close, alignment_drawer_html, 1)

# 4. JAVASCRIPT: ADD ALIGNMENT VIEWER ENGINE
js_alignment_engine = """
    // =========================================================================
    // MULTIPLE SEQUENCE ALIGNMENT (MSA) VIEWER ENGINE
    // =========================================================================
    const msaState = {{
      mode: "aa", // "aa" or "3di"
      colorScheme: "clustal", // "clustal", "zappo", "hydro", "identity", "foldstate"
      sortOrder: "tree", // "tree", "selected", "alpha"
      scrollX: 0, // Column offset
      scrollY: 0, // Row offset
      cellWidth: 14,
      cellHeight: 18,
      heightMode: "standard", // "compact", "standard", "tall"
      isMinimized: false,
      hoveredCol: null,
      hoveredRow: null
    }};

    // Standard Amino Acid Color Matrices
    const AA_CLUSTAL_COLORS = {{
      'A': '#2563eb', 'I': '#2563eb', 'L': '#2563eb', 'M': '#2563eb', 'F': '#2563eb', 'W': '#2563eb', 'V': '#2563eb',
      'K': '#dc2626', 'R': '#dc2626',
      'D': '#db2777', 'E': '#db2777',
      'S': '#059669', 'T': '#059669', 'N': '#059669', 'Q': '#059669',
      'C': '#e11d48',
      'G': '#ea580c',
      'P': '#ca8a04',
      'H': '#0891b2', 'Y': '#0891b2',
      '-': '#1e293b'
    }};

    const AA_ZAPPO_COLORS = {{
      'I': '#ffafaf', 'L': '#ffafaf', 'V': '#ffafaf', 'A': '#ffafaf', 'M': '#ffafaf',
      'F': '#ffc800', 'W': '#ffc800', 'Y': '#ffc800',
      'K': '#4040ff', 'R': '#4040ff', 'H': '#4040ff',
      'D': '#ff0000', 'E': '#ff0000',
      'S': '#00ff00', 'T': '#00ff00', 'N': '#00ff00', 'Q': '#00ff00',
      'P': '#ffff00', 'G': '#ff8000',
      'C': '#ffc0cb',
      '-': '#1e293b'
    }};

    const AA_HYDRO_COLORS = {{
      'I': '#1e3a8a', 'V': '#1d4ed8', 'L': '#2563eb', 'F': '#3b82f6', 'C': '#60a5fa',
      'M': '#93c5fd', 'A': '#bfdbfe', 'W': '#dbeafe', 'G': '#64748b', 'T': '#e2e8f0',
      'S': '#f1f5f9', 'Y': '#f8fafc', 'P': '#fed7aa', 'H': '#fdba74', 'N': '#fb923c',
      'D': '#f97316', 'Q': '#ea580c', 'E': '#c2410c', 'K': '#9a3412', 'R': '#7c2d12',
      '-': '#1e293b'
    }};

    // FoldMason 3Di Structural Alphabet Colors (20 tertiary conformation states)
    const STRUCT_3DI_COLORS = {{
      'a': '#0284c7', 'b': '#0369a1', 'c': '#0ea5e9', 'd': '#38bdf8', // Alpha-helix core
      'e': '#d97706', 'f': '#b45309', 'g': '#f59e0b', 'h': '#fbbf24', 'i': '#fde68a', // Beta-sheet strands
      'k': '#16a34a', 'l': '#15803d', 'm': '#22c55e', 'n': '#4ade80', // Turns & sharp bends
      'p': '#9333ea', 'q': '#7e22ce', 'r': '#a855f7', 's': '#c084fc', 't': '#e9d5ff', 'v': '#6b21a8', // Loops & coils
      '-': '#1e293b'
    }};

    function getActiveAlignment() {{
      if (window.ALIGNMENTS_DATA && window.ALIGNMENTS_DATA[currentScale]) {{
        return window.ALIGNMENTS_DATA[currentScale];
      }}
      return null;
    }}

    function getMsaTaxaList() {{
      const align = getActiveAlignment();
      if (!align) return [];

      const seqs = align[msaState.mode] || {{}};
      let taxa = Object.keys(seqs);

      // Filter by active metadata filter
      if (filterState.isActive) {{
        const allowed = getFilteredTaxaSet();
        taxa = taxa.filter(t => allowed.has(t));
      }}

      if (msaState.sortOrder === "tree") {{
        // Follow top-to-bottom vertical leaf order in tree
        const leaves = getAllLeaves(activeTreeRoot);
        const order = {{}};
        leaves.forEach((l, idx) => {{ order[l.name] = idx; }});
        taxa.sort((a, b) => {{
          const idxA = order[a] !== undefined ? order[a] : 99999;
          const idxB = order[b] !== undefined ? order[b] : 99999;
          return idxA - idxB;
        }});
      }} else if (msaState.sortOrder === "selected") {{
        taxa.sort((a, b) => {{
          if (a === settings.selectedTaxon) return -1;
          if (b === settings.selectedTaxon) return 1;
          return a.localeCompare(b);
        }});
      }} else {{
        taxa.sort((a, b) => a.localeCompare(b));
      }}

      return taxa;
    }}

    function computeColumnConsensus(taxaList, alignLen) {{
      const align = getActiveAlignment();
      if (!align || taxaList.length === 0) return {{ consensus: "", conservation: [] }};

      const seqs = align[msaState.mode] || {{}};
      let consensus = "";
      const conservation = [];

      for (let c = 0; c < alignLen; c++) {{
        const counts = {{}};
        let totalNonGap = 0;

        for (let i = 0; i < taxaList.length; i++) {{
          const seq = seqs[taxaList[i]];
          if (seq && c < seq.length) {{
            const ch = seq[c];
            if (ch !== '-') {{
              counts[ch] = (counts[ch] || 0) + 1;
              totalNonGap++;
            }}
          }}
        }}

        let maxChar = "-";
        let maxCount = 0;
        for (let ch in counts) {{
          if (counts[ch] > maxCount) {{
            maxCount = counts[ch];
            maxChar = ch;
          }}
        }}

        consensus += maxChar;
        const score = taxaList.length > 0 ? (maxCount / taxaList.length) : 0;
        conservation.push(score);
      }}

      return {{ consensus, conservation }};
    }}

    function setMsaMode(mode) {{
      msaState.mode = mode;
      const btnAA = document.getElementById("btnMsaModeAA");
      const btn3Di = document.getElementById("btnMsaMode3Di");
      const selColor = document.getElementById("msaColorSelect");

      if (mode === "aa") {{
        if (btnAA) btnAA.className = "px-2 py-0.5 rounded font-bold bg-sky-500 text-white transition shadow-sm cursor-pointer";
        if (btn3Di) btn3Di.className = "px-2 py-0.5 rounded font-medium text-[var(--text-muted)] hover:text-white transition cursor-pointer";
        if (selColor) {{
          selColor.innerHTML = `
            <option value="clustal" selected>🎨 ClustalX Colors</option>
            <option value="zappo">🌈 Zappo (Physicochemical)</option>
            <option value="hydro">💧 Hydrophobicity</option>
            <option value="identity">🎯 Conservation Identity</option>
          `;
        }}
        msaState.colorScheme = "clustal";
      }} else {{
        if (btn3Di) btn3Di.className = "px-2 py-0.5 rounded font-bold bg-purple-500 text-white transition shadow-sm cursor-pointer";
        if (btnAA) btnAA.className = "px-2 py-0.5 rounded font-medium text-[var(--text-muted)] hover:text-white transition cursor-pointer";
        if (selColor) {{
          selColor.innerHTML = `
            <option value="foldstate" selected>🧊 3Di Fold Geometry (Helix/Strand/Loop)</option>
            <option value="identity">🎯 Conservation Identity</option>
          `;
        }}
        msaState.colorScheme = "foldstate";
      }}

      renderMsa();
    }}

    function setMsaColorScheme(cs) {{
      msaState.colorScheme = cs;
      renderMsa();
    }}

    function setMsaSort(so) {{
      msaState.sortOrder = so;
      renderMsa();
    }}

    function jumpMsaToPosition(col) {{
      const c = parseInt(col, 10);
      const align = getActiveAlignment();
      const maxCol = align ? align.length : 533;
      if (!isNaN(c) && c >= 1 && c <= maxCol) {{
        msaState.scrollX = Math.max(0, c - 1);
        renderMsa();
      }}
    }}

    function zoomMsa(delta) {{
      msaState.cellWidth = Math.max(8, Math.min(26, msaState.cellWidth + delta));
      const badge = document.getElementById("msaZoomBadge");
      if (badge) badge.textContent = `${{msaState.cellWidth}}px`;
      renderMsa();
    }}

    function toggleMsaHeight() {{
      const drawer = document.getElementById("alignmentDrawer");
      const icon = document.getElementById("msaHeightIcon");
      if (!drawer) return;

      if (msaState.heightMode === "standard") {{
        msaState.heightMode = "tall";
        drawer.style.height = "380px";
        if (icon) icon.textContent = "↕ 380px";
      }} else if (msaState.heightMode === "tall") {{
        msaState.heightMode = "compact";
        drawer.style.height = "160px";
        if (icon) icon.textContent = "↕ 160px";
      }} else {{
        msaState.heightMode = "standard";
        drawer.style.height = "240px";
        if (icon) icon.textContent = "↕ 240px";
      }}

      renderMsa();
      setTimeout(() => fitTreeToScreen(false), 220);
    }}

    function toggleMsaDrawer() {{
      const drawer = document.getElementById("alignmentDrawer");
      const icon = document.getElementById("msaDrawerIcon");
      const label = document.getElementById("msaDrawerLabel");
      if (!drawer) return;

      if (!msaState.isMinimized) {{
        msaState.prevHeight = drawer.style.height || "240px";
        drawer.style.height = "36px";
        msaState.isMinimized = true;
        if (icon) icon.textContent = "▲";
        if (label) label.textContent = "Exp";
      }} else {{
        drawer.style.height = msaState.prevHeight;
        msaState.isMinimized = false;
        if (icon) icon.textContent = "▼";
        if (label) label.textContent = "Min";
        renderMsa();
      }}

      setTimeout(() => fitTreeToScreen(false), 220);
    }}

    function onMsaWheel(e) {{
      e.preventDefault();
      const align = getActiveAlignment();
      if (!align) return;

      const taxa = getMsaTaxaList();
      const maxCols = align.length || 533;
      const maxRows = taxa.length;

      if (Math.abs(e.deltaX) > Math.abs(e.deltaY) || e.shiftKey) {{
        const delta = e.shiftKey ? e.deltaY : e.deltaX;
        msaState.scrollX = Math.max(0, Math.min(maxCols - 5, msaState.scrollX + Math.sign(delta) * 4));
      }} else {{
        msaState.scrollY = Math.max(0, Math.min(maxRows - 2, msaState.scrollY + Math.sign(e.deltaY) * 3));
      }}

      const posInput = document.getElementById("msaPosInput");
      if (posInput) posInput.value = Math.floor(msaState.scrollX) + 1;

      renderMsa();
    }}

    function renderMsa() {{
      if (msaState.isMinimized) return;
      const align = getActiveAlignment();
      if (!align) return;

      const taxa = getMsaTaxaList();
      const alignLen = align.length || 533;
      const seqs = align[msaState.mode] || {{}};

      // Update Summary Header Badge
      const badge = document.getElementById("msaSummaryBadge");
      const visCount = document.getElementById("msaVisibleCount");
      const maxColLbl = document.getElementById("msaMaxColLabel");
      const posInput = document.getElementById("msaPosInput");

      if (badge) badge.textContent = `${{taxa.length.toLocaleString()}} seqs • ${{alignLen}} cols (${{msaState.mode.toUpperCase()}})`;
      if (visCount) visCount.textContent = taxa.length.toLocaleString();
      if (maxColLbl) maxColLbl.textContent = `/ ${{alignLen}}`;
      if (posInput) posInput.max = alignLen;

      // Compute Column Consensus & Conservation
      const {{ consensus, conservation }} = computeColumnConsensus(taxa, alignLen);

      // Average Conservation Badge
      const avgScore = conservation.length > 0 ? (conservation.reduce((a, b) => a + b, 0) / conservation.length * 100).toFixed(1) : 0;
      const consBadge = document.getElementById("msaConsensusScore");
      if (consBadge) consBadge.textContent = `${{avgScore}}% Avg`;

      // 1. RENDER LEFT TAXA COLUMN
      const taxaListEl = document.getElementById("msaTaxaList");
      if (taxaListEl) {{
        taxaListEl.innerHTML = "";
        const visibleRowCount = Math.ceil(taxaListEl.clientHeight / msaState.cellHeight) + 1;
        const startR = Math.floor(msaState.scrollY);
        const endR = Math.min(taxa.length, startR + visibleRowCount);

        for (let i = startR; i < endR; i++) {{
          const tName = taxa[i];
          const m = TAXA_METADATA[tName] || {{}};
          const isSelected = (tName === settings.selectedTaxon);

          // Get category color dot
          const colDef = getActiveColorColumnDef();
          const catVal = (colDef && m[colDef.key] !== undefined) ? String(m[colDef.key]) : "Unclassified";
          const dotColor = (colDef && colDef.colors && colDef.colors[catVal]) ? colDef.colors[catVal] : "#38bdf8";

          const rowEl = document.createElement("div");
          rowEl.className = `flex items-center justify-between px-2.5 text-[10.5px] cursor-pointer truncate select-none border-b border-white/5 transition-all ${{
            isSelected ? "bg-sky-500/25 text-sky-300 font-bold border-l-2 border-sky-400" : "hover:bg-slate-500/15 text-[var(--text-muted)]"
          }}`;
          rowEl.style.height = `${{msaState.cellHeight}}px`;
          rowEl.style.lineHeight = `${{msaState.cellHeight}}px`;
          rowEl.title = `${{tName}} (${{catVal}}) - Click to focus in tree`;

          rowEl.innerHTML = `
            <div class="flex items-center space-x-1.5 truncate">
              <span class="w-2 h-2 rounded-full shrink-0" style="background-color: ${{dotColor}}"></span>
              <span class="truncate font-mono">${{tName}}</span>
            </div>
            <span class="text-[8.5px] font-mono text-slate-500 shrink-0">#${{i + 1}}</span>
          `;

          rowEl.onclick = () => {{
            selectTaxon(tName);
            highlightTangleTaxon(tName);
            centerOnSelection();
            renderMsa();
          }};
          rowEl.onmouseenter = () => highlightTangleTaxon(tName);
          rowEl.onmouseleave = () => clearTangleHighlight();

          taxaListEl.appendChild(rowEl);
        }}
      }}

      // 2. RENDER RULER CANVAS
      const rulerCanvas = document.getElementById("msaRulerCanvas");
      if (rulerCanvas) {{
        const rRect = rulerCanvas.parentElement.getBoundingClientRect();
        rulerCanvas.width = rRect.width;
        rulerCanvas.height = rRect.height;
        const rctx = rulerCanvas.getContext("2d");
        rctx.clearRect(0, 0, rulerCanvas.width, rulerCanvas.height);

        const cW = msaState.cellWidth;
        const startC = Math.floor(msaState.scrollX);
        const colCount = Math.ceil(rulerCanvas.width / cW) + 1;
        const endC = Math.min(alignLen, startC + colCount);

        rctx.fillStyle = "#94a3b8";
        rctx.font = "9px monospace";
        rctx.textBaseline = "middle";

        for (let c = startC; c < endC; c++) {{
          const x = (c - msaState.scrollX) * cW;
          const posNum = c + 1;

          if (posNum % 10 === 0 || posNum === 1) {{
            rctx.strokeStyle = "#475569";
            rctx.lineWidth = 1;
            rctx.beginPath();
            rctx.moveTo(x, 14);
            rctx.lineTo(x, 24);
            rctx.stroke();
            rctx.fillText(String(posNum), x + 2, 8);
          }} else if (posNum % 5 === 0) {{
            rctx.strokeStyle = "#334155";
            rctx.lineWidth = 1;
            rctx.beginPath();
            rctx.moveTo(x, 18);
            rctx.lineTo(x, 24);
            rctx.stroke();
          }}
        }}
      }}

      // 3. RENDER RESIDUE MATRIX CANVAS
      const matrixCanvas = document.getElementById("msaMatrixCanvas");
      if (matrixCanvas) {{
        const mRect = matrixCanvas.parentElement.getBoundingClientRect();
        matrixCanvas.width = mRect.width;
        matrixCanvas.height = mRect.height;
        const mctx = matrixCanvas.getContext("2d");
        mctx.clearRect(0, 0, matrixCanvas.width, matrixCanvas.height);

        const cW = msaState.cellWidth;
        const cH = msaState.cellHeight;

        const startC = Math.floor(msaState.scrollX);
        const colCount = Math.ceil(matrixCanvas.width / cW) + 1;
        const endC = Math.min(alignLen, startC + colCount);

        const startR = Math.floor(msaState.scrollY);
        const rowCount = Math.ceil(matrixCanvas.height / cH) + 1;
        const endR = Math.min(taxa.length, startR + rowCount);

        mctx.font = `bold ${{Math.max(9, Math.min(13, cW - 2))}}px 'SF Mono', Monaco, Menlo, Consolas, monospace`;
        mctx.textAlign = "center";
        mctx.textBaseline = "middle";

        for (let r = startR; r < endR; r++) {{
          const tName = taxa[r];
          const seq = seqs[tName] || "";
          const y = (r - msaState.scrollY) * cH;
          const isSelected = (tName === settings.selectedTaxon);

          for (let c = startC; c < endC; c++) {{
            const x = (c - msaState.scrollX) * cW;
            const ch = c < seq.length ? seq[c] : "-";

            // Determine Residue Background Color
            let bg = "#1e293b";
            if (ch === "-") {{
              bg = "#0f172a";
            }} else if (msaState.mode === "3di") {{
              if (msaState.colorScheme === "identity") {{
                const isCons = (ch === consensus[c]);
                bg = isCons ? "#10b981" : "#334155";
              }} else {{
                bg = STRUCT_3DI_COLORS[ch.toLowerCase()] || "#9333ea";
              }}
            }} else {{
              if (msaState.colorScheme === "zappo") {{
                bg = AA_ZAPPO_COLORS[ch.toUpperCase()] || "#1e293b";
              }} else if (msaState.colorScheme === "hydro") {{
                bg = AA_HYDRO_COLORS[ch.toUpperCase()] || "#1e293b";
              }} else if (msaState.colorScheme === "identity") {{
                const isCons = (ch === consensus[c]);
                const sc = conservation[c] || 0;
                bg = isCons ? (sc >= 0.8 ? "#10b981" : "#38bdf8") : "#334155";
              }} else {{
                bg = AA_CLUSTAL_COLORS[ch.toUpperCase()] || "#1e293b";
              }}
            }}

            mctx.fillStyle = bg;
            mctx.fillRect(x, y, cW - 0.5, cH - 0.5);

            // Text
            if (cW >= 9 && ch !== "-") {{
              mctx.fillStyle = "#ffffff";
              mctx.fillText(ch, x + cW / 2, y + cH / 2);
            }}
          }}

          if (isSelected) {{
            mctx.strokeStyle = "#38bdf8";
            mctx.lineWidth = 1.5;
            mctx.strokeRect(0, y, matrixCanvas.width, cH);
          }}
        }}
      }}

      // 4. RENDER CONSENSUS & CONSERVATION BAR CANVAS
      const consCanvas = document.getElementById("msaConsensusCanvas");
      if (consCanvas) {{
        const cRect = consCanvas.parentElement.getBoundingClientRect();
        consCanvas.width = cRect.width;
        consCanvas.height = cRect.height;
        const cctx = consCanvas.getContext("2d");
        cctx.clearRect(0, 0, consCanvas.width, consCanvas.height);

        const cW = msaState.cellWidth;
        const startC = Math.floor(msaState.scrollX);
        const colCount = Math.ceil(consCanvas.width / cW) + 1;
        const endC = Math.min(alignLen, startC + colCount);

        cctx.textAlign = "center";
        cctx.textBaseline = "middle";

        for (let c = startC; c < endC; c++) {{
          const x = (c - msaState.scrollX) * cW;
          const ch = consensus[c] || "-";
          const score = conservation[c] || 0;

          // Conservation Bar (bottom half)
          const barH = score * 14;
          cctx.fillStyle = score >= 0.8 ? "#10b981" : (score >= 0.5 ? "#38bdf8" : "#f59e0b");
          cctx.fillRect(x, 28 - barH, cW - 0.5, barH);

          // Consensus Character (top half)
          if (cW >= 9 && ch !== "-") {{
            cctx.font = `bold ${{Math.max(9, Math.min(12, cW - 2))}}px monospace`;
            cctx.fillStyle = "#f8fafc";
            cctx.fillText(ch, x + cW / 2, 7);
          }}
        }}
      }}
    }}
"""

cartesian_anchor = "    // CARTESIAN LAYOUT (Phylogram & Cladogram with Triangular Clades)"
assert cartesian_anchor in code, "cartesian_anchor not found"
code = code.replace(cartesian_anchor, js_alignment_engine + "\n    " + cartesian_anchor, 1)

# 5. HOOK renderMsa INTO renderTree AND switchDatasetScale
old_render_tree_end = """        updateMinimap();
      }} catch (err) {{
        console.error("Tree render error:", err);
      }}
    }}"""

new_render_tree_end = """        updateMinimap();
        renderMsa();
      }} catch (err) {{
        console.error("Tree render error:", err);
      }}
    }}"""

assert old_render_tree_end in code, "old_render_tree_end not found"
code = code.replace(old_render_tree_end, new_render_tree_end, 1)

# In switchDatasetScale: reset msaState.scrollX and scrollY, then renderMsa
old_scale_init = """      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      updateFilterUI();
      applyCurrentRooting();
      updateCladeManagementUI();
      updateCongruenceUI();"""

new_scale_init = """      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      updateFilterUI();
      applyCurrentRooting();
      updateCladeManagementUI();
      updateCongruenceUI();
      msaState.scrollX = 0;
      msaState.scrollY = 0;
      renderMsa();"""

assert old_scale_init in code, "old_scale_init not found"
code = code.replace(old_scale_init, new_scale_init, 1)

# In selectTaxon: scroll MSA to selected taxon
old_select_taxon = """    function selectTaxon(taxName) {{
      settings.selectedTaxon = taxName;"""

new_select_taxon = """    function selectTaxon(taxName) {{
      settings.selectedTaxon = taxName;
      if (typeof getMsaTaxaList === 'function') {{
        const msaTaxa = getMsaTaxaList();
        const idx = msaTaxa.indexOf(taxName);
        if (idx >= 0) {{
          msaState.scrollY = Math.max(0, idx - 2);
          renderMsa();
        }}
      }}"""

assert old_select_taxon in code, "old_select_taxon not found"
code = code.replace(old_select_taxon, new_select_taxon, 1)

# In DOMContentLoaded: add renderMsa and resize listener
old_init_block = """    // INITIALIZATION
    window.addEventListener("DOMContentLoaded", () => {{
      setTheme(settings.theme);
      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      applyCurrentRooting();
      updateCladeManagementUI();

      const firstTaxon = Object.keys(TAXA_METADATA)[0];
      if (firstTaxon) {{
        selectTaxon(firstTaxon);
      }}
    }});"""

new_init_block = """    // INITIALIZATION
    window.addEventListener("DOMContentLoaded", () => {{
      setTheme(settings.theme);
      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      updateFilterUI();
      applyCurrentRooting();
      updateCladeManagementUI();
      updateCongruenceUI();
      renderMsa();

      const firstTaxon = Object.keys(TAXA_METADATA)[0];
      if (firstTaxon) {{
        selectTaxon(firstTaxon);
      }}
    }});

    window.addEventListener("resize", () => {{
      renderTree();
      renderMsa();
    }});"""

assert old_init_block in code, "old_init_block not found"
code = code.replace(old_init_block, new_init_block, 1)

with open(target_path, "w", encoding="utf-8") as f:
    f.write(code)

print("Successfully applied Alignment Viewer patch to build_dynamic_interactive_tree.py!")
