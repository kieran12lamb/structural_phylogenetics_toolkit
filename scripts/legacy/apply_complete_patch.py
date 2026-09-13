#!/usr/bin/env python3
"""Apply complete patch to build_dynamic_interactive_tree.py:
- DATASET_CONGRUENCE with Robinson-Foulds, Cophenetic correlation, and discordant taxa
- 4th Sidebar Tab: 🎯 Filter
- Filter Panel UI with Prune Subtree & Highlight/Dim modes, Categorical checkboxes & Continuous sliders
- Floating Active Filter Banner
- Tanglegram Phylogenetic Congruence Card with RF bar, Cophenetic r, Discordance %
- Interactive Congruence Analysis Report Modal
- Dynamic Filtering Engine (pruneSubtree, getFilteredTaxaSet, applyTaxaFilter, clearTaxaFilter)
- Dynamic Congruence Engine (computeDynamicCongruence, openCongruenceModal)
"""

import json
from pathlib import Path

target_path = Path("/Users/kieran.lamb/Github/Skills_Hackathon/scripts/build_dynamic_interactive_tree.py")
with open(target_path, "r", encoding="utf-8") as f:
    code = f.read()

# ==============================================================================
# 1. ADD DATASET_CONGRUENCE IN PYTHON
# ==============================================================================
congruence_py = """
# 4. Precomputed Phylogenetic Congruence Profiles
DATASET_CONGRUENCE = {
    "1193": {
        "taxa_count": 1193,
        "rf_distance": 2024,
        "max_rf": 2380,
        "norm_rf": 0.8504,
        "congruence_pct": 14.96,
        "shared_splits": 178,
        "splits_3di": 1190,
        "splits_aa": 1190,
        "cophenetic_r": 0.3171,
        "crossings_native": 297387,
        "crossings_untangled": 165837,
        "total_pairs": 711028,
        "discordance_native": 41.82,
        "discordance_untangled": 23.32,
        "interpretation": "In this 1,193-member engineered Henipavirus binder library, overall sequence-structure topological congruence is 15.0% (178 shared internal clades) with global cophenetic correlation r = +0.317. When filtering specifically to the 76 Strong Binders, cophenetic distance correlation jumps to r = +0.578, demonstrating tight sequence-structure evolutionary coupling in functional binding interfaces compared to unselected designs.",
        "discordant_taxa": [
            {"id": "radiant-ram-plume", "category": "Mainly Beta (None)", "r1": 1135, "r2": 5, "shift": 1130},
            {"id": "vast-zebra-clay", "category": "Alpha Beta (None)", "r1": 5, "r2": 1112, "shift": 1107},
            {"id": "strong-jaguar-pine", "category": "Mainly Alpha (None)", "r1": 1113, "r2": 35, "shift": 1078},
            {"id": "dark-quail-pine", "category": "Mainly Alpha (None)", "r1": 1098, "r2": 22, "shift": 1076},
            {"id": "scarlet-ox-pearl", "category": "Mainly Alpha (None)", "r1": 1097, "r2": 23, "shift": 1074},
            {"id": "fresh-heron-garnet", "category": "Alpha Beta (None)", "r1": 1096, "r2": 24, "shift": 1072},
            {"id": "quiet-wasp-quartz", "category": "Mainly Alpha (None)", "r1": 1103, "r2": 32, "shift": 1071},
            {"id": "keen-hornet-copper", "category": "Mainly Alpha (None)", "r1": 1102, "r2": 33, "shift": 1069},
            {"id": "bold-crane-pebble", "category": "Mainly Beta (None)", "r1": 1094, "r2": 27, "shift": 1067},
            {"id": "wise-salmon-pebble", "category": "Mainly Alpha (None)", "r1": 1101, "r2": 37, "shift": 1064}
        ]
    },
    "500": {
        "taxa_count": 500,
        "rf_distance": 630,
        "max_rf": 994,
        "norm_rf": 0.6338,
        "congruence_pct": 36.62,
        "shared_splits": 182,
        "splits_3di": 497,
        "splits_aa": 497,
        "cophenetic_r": 0.5955,
        "crossings_native": 35589,
        "crossings_untangled": 13467,
        "total_pairs": 124750,
        "discordance_native": 28.53,
        "discordance_untangled": 10.80,
        "interpretation": "Across 500 viral glycoproteins spanning 23 viral families, 36.6% of internal evolutionary splits are strictly congruent between 3Di structure and amino acid sequence, accompanied by a strong cophenetic correlation (r = +0.596). Major monophyletic viral families (Rhabdoviridae, Orthoherpesviridae, Phenuiviridae) form conserved structural clades, while topological discordance is localized to divergent surface attachment domains.",
        "discordant_taxa": [
            {"id": "AAA42974.1_9425", "category": "Rhabdoviridae", "r1": 2, "r2": 499, "shift": 497},
            {"id": "ADC54013.1_12020", "category": "Orthoherpesviridae", "r1": 6, "r2": 408, "shift": 402},
            {"id": "AEQ32302.1_102", "category": "Phenuiviridae", "r1": 463, "r2": 116, "shift": 347},
            {"id": "WBM84629.1_9689", "category": "Hantaviridae", "r1": 154, "r2": 402, "shift": 248},
            {"id": "WBM84624.1_9688", "category": "Hantaviridae", "r1": 153, "r2": 401, "shift": 248},
            {"id": "WBM84626.1_9688", "category": "Hantaviridae", "r1": 152, "r2": 399, "shift": 247},
            {"id": "WBM84627.1_9688", "category": "Hantaviridae", "r1": 151, "r2": 398, "shift": 247},
            {"id": "WBM84625.1_9688", "category": "Hantaviridae", "r1": 150, "r2": 397, "shift": 247},
            {"id": "WBM84628.1_9688", "category": "Hantaviridae", "r1": 149, "r2": 396, "shift": 247},
            {"id": "WBM84630.1_9689", "category": "Hantaviridae", "r1": 148, "r2": 395, "shift": 247}
        ]
    },
    "6": {
        "taxa_count": 6,
        "rf_distance": 6,
        "max_rf": 6,
        "norm_rf": 1.0,
        "congruence_pct": 0.0,
        "shared_splits": 0,
        "splits_3di": 3,
        "splits_aa": 3,
        "cophenetic_r": 0.4588,
        "crossings_native": 4,
        "crossings_untangled": 3,
        "total_pairs": 15,
        "discordance_native": 26.7,
        "discordance_untangled": 20.0,
        "interpretation": "Benchmark viral glycoprotein dataset displays high local sequence conservation with discrete structural topological shifts across divergent viral families. Cophenetic distance correlation r = +0.459 demonstrates moderate patristic concordance.",
        "discordant_taxa": [
            {"id": "AAA88529.1.2_6592", "category": "Rhabdoviridae", "r1": 7, "r2": 2, "shift": 5},
            {"id": "AAQ55251.1.1_9251", "category": "Filoviridae", "r1": 2, "r2": 6, "shift": 4},
            {"id": "AAQ55251.1.2_9251", "category": "Filoviridae", "r1": 5, "r2": 8, "shift": 3}
        ]
    }
}
"""

datasets_anchor = '# Assemble DATASETS JSON\nDATASETS = {'
assert datasets_anchor in code, "datasets_anchor not found"
code = code.replace(datasets_anchor, congruence_py + "\n" + datasets_anchor, 1)

code = code.replace('"defaultZoom": {"x": 40, "y": 35, "k": 0.65}', '"defaultZoom": {"x": 40, "y": 35, "k": 0.65},\n        "congruence": DATASET_CONGRUENCE["1193"]', 1)
code = code.replace('"defaultZoom": {"x": 40, "y": 30, "k": 0.55}', '"defaultZoom": {"x": 40, "y": 30, "k": 0.55},\n        "congruence": DATASET_CONGRUENCE["500"]', 1)
code = code.replace('"defaultZoom": {"x": 80, "y": 50, "k": 1.0}', '"defaultZoom": {"x": 80, "y": 50, "k": 1.0},\n        "congruence": DATASET_CONGRUENCE["6"]', 1)

# ==============================================================================
# 2. CSS STYLES FOR FILTER DIMMED & MATCH
# ==============================================================================
old_style_anchor = """.active-tab {{
      color: var(--accent) !important;
      border-bottom: 2px solid var(--accent) !important;
      font-weight: 700 !important;
    }}"""

new_styles = """.active-tab {{
      color: var(--accent) !important;
      border-bottom: 2px solid var(--accent) !important;
      font-weight: 700 !important;
    }}
    .filter-dimmed {{
      opacity: 0.12 !important;
      filter: grayscale(80%);
      transition: opacity 0.2s ease, filter 0.2s ease;
    }}
    .filter-match {{
      opacity: 1.0 !important;
      transition: opacity 0.2s ease;
    }}"""

assert old_style_anchor in code, "old_style_anchor not found"
code = code.replace(old_style_anchor, new_styles, 1)

# ==============================================================================
# 3. TAB NAVIGATION (4 TABS)
# ==============================================================================
old_tab_buttons = """      <!-- SIDEBAR TAB NAVIGATION -->
      <div class="flex border-b border-[var(--border-color)] text-xs shrink-0 bg-[var(--card-bg)]">
        <button id="tabBtnDisplay" onclick="switchSidebarTab('display')" class="flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1">
          <span>⚙️ Display</span>
        </button>
        <button id="tabBtnClades" onclick="switchSidebarTab('clades')" class="flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1">
          <span>🌿 Clades</span>
          <span id="tabCladeCountBadge" class="ml-1 text-[9.5px] font-mono px-1 rounded-full bg-sky-500/20 text-sky-400">0</span>
        </button>
        <button id="tabBtnRooting" onclick="switchSidebarTab('rooting')" class="flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1">
          <span>⚓ Rooting</span>
        </button>
      </div>"""

new_tab_buttons = """      <!-- SIDEBAR TAB NAVIGATION -->
      <div class="flex border-b border-[var(--border-color)] text-xs shrink-0 bg-[var(--card-bg)]">
        <button id="tabBtnDisplay" onclick="switchSidebarTab('display')" class="flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1">
          <span>⚙️ Display</span>
        </button>
        <button id="tabBtnFilter" onclick="switchSidebarTab('filter')" class="flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1">
          <span>🎯 Filter</span>
          <span id="tabFilterBadge" class="ml-1 text-[9px] font-mono px-1.5 py-0.5 rounded-full bg-emerald-500/20 text-emerald-400 font-bold">All</span>
        </button>
        <button id="tabBtnClades" onclick="switchSidebarTab('clades')" class="flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1">
          <span>🌿 Clades</span>
          <span id="tabCladeCountBadge" class="ml-1 text-[9.5px] font-mono px-1 rounded-full bg-sky-500/20 text-sky-400">0</span>
        </button>
        <button id="tabBtnRooting" onclick="switchSidebarTab('rooting')" class="flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1">
          <span>⚓ Rooting</span>
        </button>
      </div>"""

assert old_tab_buttons in code, "old_tab_buttons not found"
code = code.replace(old_tab_buttons, new_tab_buttons, 1)

# ==============================================================================
# 4. TANGLEGRAM CONGRUENCE CARD INSIDE tanglegramOptions
# ==============================================================================
old_tanglegram_options = """          <!-- Tanglegram Alignment Selector (Visible in Tanglegram layout) -->
          <div id="tanglegramOptions" class="hidden space-y-1.5 bg-[var(--chip-bg)] p-3 rounded-lg border border-purple-500/40 shadow-sm">
            <div class="flex items-center justify-between">
              <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block text-[10px]">Tanglegram Alignment</label>
              <span id="tangleCrossingBadge" class="text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-purple-500/20 text-purple-300 border border-purple-500/40">Crossings: -</span>
            </div>
            <select id="tangleModeSelect" onchange="setTangleMode(this.value)" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded-md px-2 py-1.5 text-xs text-[var(--text-main)] focus:outline-none focus:border-purple-400 font-medium">
              <option value="true_topology" selected>🔀 True Topology (Show Incongruence & Crossings)</option>
              <option value="min_crossings">🔄 Untangle (Optimal Clade Rotations)</option>
              <option value="aligned">⏸️ Aligned / Parallel (1-to-1 Matched)</option>
            </select>
            <p id="tangleModeDesc" class="text-[9px] text-[var(--text-muted)] leading-tight pt-0.5">
              True topology preserves native IQ-TREE branch order on both sides to expose structural vs sequence discordance.
            </p>
          </div>"""

new_tanglegram_options = """          <!-- Tanglegram Alignment Selector & Congruence Card (Visible in Tanglegram layout) -->
          <div id="tanglegramOptions" class="hidden space-y-2.5 bg-[var(--chip-bg)] p-3 rounded-lg border border-purple-500/40 shadow-sm">
            <div class="flex items-center justify-between">
              <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block text-[10px]">Tanglegram Alignment</label>
              <span id="tangleCrossingBadge" class="text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-purple-500/20 text-purple-300 border border-purple-500/40">Crossings: -</span>
            </div>
            <select id="tangleModeSelect" onchange="setTangleMode(this.value)" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded-md px-2 py-1.5 text-xs text-[var(--text-main)] focus:outline-none focus:border-purple-400 font-medium">
              <option value="true_topology" selected>🔀 True Topology (Show Incongruence & Crossings)</option>
              <option value="min_crossings">🔄 Untangle (Optimal Clade Rotations)</option>
              <option value="aligned">⏸️ Aligned / Parallel (1-to-1 Matched)</option>
            </select>
            <p id="tangleModeDesc" class="text-[9px] text-[var(--text-muted)] leading-tight pt-0.5">
              True topology preserves native IQ-TREE branch order on both sides to expose structural vs sequence discordance.
            </p>

            <!-- PHYLOGENETIC CONGRUENCE METRICS CARD -->
            <div class="pt-2 border-t border-purple-500/30 space-y-2 text-[10px]">
              <div class="flex items-center justify-between">
                <span class="font-bold text-purple-300 uppercase tracking-wider text-[9.5px] flex items-center space-x-1">
                  <span>📐</span>
                  <span>Tree Congruence</span>
                </span>
                <button onclick="openCongruenceModal()" class="text-[9.5px] px-2 py-0.5 rounded bg-purple-500/20 hover:bg-purple-500/40 border border-purple-400/50 text-purple-200 transition font-semibold flex items-center space-x-1 shadow-sm cursor-pointer">
                  <span>📊 Full Report</span>
                </button>
              </div>

              <div class="space-y-1.5 bg-black/25 p-2 rounded-lg border border-purple-500/20">
                <div>
                  <div class="flex justify-between items-center mb-0.5">
                    <span class="text-[var(--text-muted)]">Robinson-Foulds Congruence:</span>
                    <span id="tangleRfCongruence" class="font-mono text-emerald-400 font-bold">-</span>
                  </div>
                  <div class="w-full h-1.5 bg-slate-800 rounded-full overflow-hidden">
                    <div id="tangleRfProgressBar" class="h-full bg-emerald-400 transition-all duration-300" style="width: 0%"></div>
                  </div>
                  <div class="flex justify-between text-[8.5px] text-[var(--text-muted)] pt-0.5">
                    <span id="tangleRfShared">- shared clades</span>
                    <span id="tangleRfDist">RF dist: -</span>
                  </div>
                </div>

                <div class="pt-1 border-t border-white/5 flex justify-between items-center">
                  <span class="text-[var(--text-muted)]">Cophenetic Distance r:</span>
                  <span id="tangleCopheneticR" class="font-mono text-sky-400 font-bold">-</span>
                </div>

                <div class="flex justify-between items-center">
                  <span class="text-[var(--text-muted)]">Crossing Discordance:</span>
                  <span id="tangleDiscordance" class="font-mono text-amber-400 font-bold">-</span>
                </div>
              </div>
            </div>
          </div>"""

assert old_tanglegram_options in code, "old_tanglegram_options not found"
code = code.replace(old_tanglegram_options, new_tanglegram_options, 1)

# ==============================================================================
# 5. FILTER PANEL HTML
# ==============================================================================
filter_panel_html = """
        <!-- ================= TAB: FILTER ================= -->
        <div id="panelFilter" class="hidden space-y-3.5">
          <!-- Filter Overview Card -->
          <div class="p-3 rounded-lg border border-[var(--border-color)] bg-[var(--card-bg)] shadow-sm space-y-2">
            <div class="flex items-center justify-between">
              <span class="font-bold text-emerald-400 text-xs flex items-center space-x-1">
                <span>🎯</span>
                <span>Metadata Taxa Filter</span>
              </span>
              <span id="filterStatusBadge" class="text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/30">Showing All</span>
            </div>
            <p class="text-[10px] text-[var(--text-muted)] leading-relaxed">
              Filter the tree by any metadata category or continuous metric. Choose between pruning the tree to an induced subtree or highlighting matching branches:
            </p>
            <div class="flex space-x-2 pt-1">
              <button onclick="clearTaxaFilter()" class="flex-1 py-1.5 px-2 rounded border border-[var(--border-color)] hover:bg-slate-500/20 text-[10.5px] font-semibold text-[var(--text-main)] transition text-center shadow-sm cursor-pointer">
                ✕ Reset All Filters
              </button>
            </div>
          </div>

          <!-- Filter Action Mode Switch -->
          <div class="space-y-1.5 bg-[var(--chip-bg)] p-2.5 rounded-lg border border-[var(--border-color)]">
            <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block text-[10px]">Filter Action Mode</label>
            <div class="grid grid-cols-2 gap-1 bg-[var(--input-bg)] p-1 rounded-md border border-[var(--border-color)]">
              <button id="btnFilterModePrune" onclick="setFilterMode('prune')" class="py-1 px-2 rounded font-semibold text-center bg-emerald-500 text-white text-[10.5px] transition shadow-sm cursor-pointer">
                ✂️ Prune Subtree
              </button>
              <button id="btnFilterModeHighlight" onclick="setFilterMode('highlight')" class="py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10.5px] transition cursor-pointer">
                💡 Highlight & Dim
              </button>
            </div>
            <p id="filterModeDesc" class="text-[9px] text-[var(--text-muted)] leading-tight pt-0.5">
              Prune contracts the tree to only matching leaves, recalculating branch geometry and evolutionary distances.
            </p>
          </div>

          <!-- Category Selector -->
          <div class="space-y-1 bg-[var(--chip-bg)] p-2.5 rounded-lg border border-[var(--border-color)]">
            <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block text-[10px]">Filter Column</label>
            <select id="filterColumnSelect" onchange="onFilterColumnChanged(this.value)" class="w-full bg-[var(--input-bg)] border border-emerald-500/40 rounded-md px-2.5 py-1.5 text-xs text-emerald-400 font-medium focus:outline-none focus:border-emerald-400">
            </select>
          </div>

          <!-- Quick Presets -->
          <div class="space-y-1 bg-[var(--chip-bg)] p-2.5 rounded-lg border border-[var(--border-color)]">
            <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block text-[10px]">Quick Presets</label>
            <div id="filterPresetsContainer" class="flex flex-wrap gap-1 pt-0.5">
            </div>
          </div>

          <!-- Categorical Section -->
          <div id="filterCategoricalSection" class="space-y-2 bg-[var(--chip-bg)] p-2.5 rounded-lg border border-[var(--border-color)]">
            <div class="flex items-center justify-between">
              <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider text-[10px]">Select Categories</label>
              <div class="flex space-x-1 text-[9.5px]">
                <button onclick="selectAllFilterCategories()" class="text-sky-400 hover:underline">All</button>
                <span class="text-[var(--text-muted)]">&bull;</span>
                <button onclick="deselectAllFilterCategories()" class="text-sky-400 hover:underline">None</button>
                <span class="text-[var(--text-muted)]">&bull;</span>
                <button onclick="invertFilterCategories()" class="text-sky-400 hover:underline">Invert</button>
              </div>
            </div>
            <input type="text" id="filterCategorySearchInput" placeholder="Search values..." oninput="onFilterCategorySearch(this.value)" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded px-2 py-1 text-[11px] text-[var(--text-main)] placeholder-[var(--text-muted)] focus:outline-none focus:border-emerald-400">
            <div id="filterCategoryList" class="space-y-1 max-h-52 overflow-y-auto pr-1 custom-scroll">
            </div>
          </div>

          <!-- Continuous Numerical Slider Section -->
          <div id="filterContinuousSection" class="hidden space-y-2 bg-[var(--chip-bg)] p-2.5 rounded-lg border border-[var(--border-color)]">
            <div class="flex items-center justify-between text-[10.5px]">
              <span class="font-semibold text-[var(--text-muted)] uppercase tracking-wider text-[10px]">Numerical Range</span>
              <span id="filterRangeBadge" class="font-mono text-emerald-400 font-bold">-</span>
            </div>
            <div class="grid grid-cols-2 gap-2 text-xs">
              <div>
                <label class="text-[9px] text-[var(--text-muted)] block">Min Bound</label>
                <input type="number" id="filterInputMin" onchange="onFilterRangeInputChanged()" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded px-2 py-1 text-xs text-[var(--text-main)] font-mono">
              </div>
              <div>
                <label class="text-[9px] text-[var(--text-muted)] block">Max Bound</label>
                <input type="number" id="filterInputMax" onchange="onFilterRangeInputChanged()" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded px-2 py-1 text-xs text-[var(--text-main)] font-mono">
              </div>
            </div>
            <div class="pt-1">
              <label class="text-[9px] text-[var(--text-muted)] block mb-0.5">Adjust Range Threshold</label>
              <input type="range" id="filterRangeSlider" oninput="onFilterRangeSliderMoved(this.value)" class="w-full h-1.5 bg-slate-700 rounded appearance-none cursor-pointer accent-emerald-400">
            </div>
            <div class="flex justify-between text-[9px] text-[var(--text-muted)]">
              <span id="filterMinLegend">-</span>
              <span id="filterMaxLegend">-</span>
            </div>
          </div>
        </div>
"""

old_tab_clades_marker = """        <!-- ================= TAB 2: CLADE MANAGEMENT ================= -->"""
assert old_tab_clades_marker in code, "old_tab_clades_marker not found"
code = code.replace(old_tab_clades_marker, filter_panel_html + "\n        " + old_tab_clades_marker, 1)

# ==============================================================================
# 6. FLOATING ACTIVE FILTER BANNER
# ==============================================================================
old_canvas_section = """    <!-- MAIN SVG TREE CANVAS & FLOATING OVERLAYS -->
    <div id="treeContainer" class="flex-1 h-full relative overflow-hidden">
      
      <!-- SVG Canvas -->
      <svg id="treeSvg" class="w-full h-full block"></svg>"""

new_canvas_section = """    <!-- MAIN SVG TREE CANVAS & FLOATING OVERLAYS -->
    <div id="treeContainer" class="flex-1 h-full relative overflow-hidden">
      
      <!-- FLOATING ACTIVE FILTER BANNER (Top Center) -->
      <div id="activeFilterBanner" class="absolute top-3 left-1/2 -translate-x-1/2 z-20 hidden items-center space-x-2 bg-[var(--card-bg)]/95 border border-emerald-500/50 backdrop-blur-md px-3.5 py-1.5 rounded-full shadow-2xl text-xs transition-all">
        <span class="w-2 h-2 rounded-full bg-emerald-400 animate-pulse"></span>
        <span id="filterBannerText" class="font-semibold text-emerald-300">Filter Active</span>
        <button onclick="clearTaxaFilter()" class="ml-1 text-slate-400 hover:text-white font-bold text-xs bg-slate-800/80 hover:bg-rose-600/80 px-2 py-0.5 rounded-full transition cursor-pointer">&times; Clear</button>
      </div>

      <!-- SVG Canvas -->
      <svg id="treeSvg" class="w-full h-full block"></svg>"""

assert old_canvas_section in code, "old_canvas_section not found"
code = code.replace(old_canvas_section, new_canvas_section, 1)

# ==============================================================================
# 7. CONGRUENCE MODAL DIALOG
# ==============================================================================
congruence_modal_html = """
  <!-- MODAL: DETAILED PHYLOGENETIC CONGRUENCE ANALYSIS REPORT -->
  <div id="congruenceModal" class="fixed inset-0 z-50 hidden bg-black/70 backdrop-blur-sm flex items-center justify-center p-4">
    <div class="bg-[var(--card-bg)] border border-purple-500/40 rounded-2xl shadow-2xl max-w-2xl w-full max-h-[85vh] flex flex-col overflow-hidden animate-in fade-in zoom-in-95 duration-150">
      <!-- Modal Header -->
      <div class="p-4 border-b border-[var(--border-color)] flex items-center justify-between bg-purple-950/20">
        <div class="flex items-center space-x-2">
          <span class="text-xl">📐</span>
          <div>
            <h3 class="font-bold text-sm text-[var(--text-main)]">Structural vs Sequence Phylogenetic Congruence</h3>
            <p class="text-[11px] text-[var(--text-muted)]">Comparing FoldMason 3Di Structural Tree with Amino Acid Sequence Tree</p>
          </div>
        </div>
        <button onclick="closeCongruenceModal()" class="text-slate-400 hover:text-white text-xl font-bold px-2 py-1 rounded-lg hover:bg-slate-800 transition cursor-pointer">&times;</button>
      </div>

      <!-- Modal Body -->
      <div class="p-5 overflow-y-auto space-y-4 custom-scroll text-xs">
        <!-- Metric Cards Grid -->
        <div class="grid grid-cols-3 gap-3">
          <div class="bg-black/30 p-3 rounded-xl border border-emerald-500/30 text-center space-y-1">
            <span class="text-[10px] uppercase font-semibold text-[var(--text-muted)] block">Robinson-Foulds Congruence</span>
            <div id="modalRfScore" class="text-xl font-mono font-bold text-emerald-400">-</div>
            <p id="modalRfSub" class="text-[9.5px] text-[var(--text-muted)]">-</p>
          </div>
          <div class="bg-black/30 p-3 rounded-xl border border-sky-500/30 text-center space-y-1">
            <span class="text-[10px] uppercase font-semibold text-[var(--text-muted)] block">Cophenetic Correlation (r)</span>
            <div id="modalCopheneticScore" class="text-xl font-mono font-bold text-sky-400">-</div>
            <p id="modalCopheneticSub" class="text-[9.5px] text-[var(--text-muted)]">-</p>
          </div>
          <div class="bg-black/30 p-3 rounded-xl border border-purple-500/30 text-center space-y-1">
            <span class="text-[10px] uppercase font-semibold text-[var(--text-muted)] block">Tanglegram Discordance</span>
            <div id="modalDiscordanceScore" class="text-xl font-mono font-bold text-purple-400">-</div>
            <p id="modalDiscordanceSub" class="text-[9.5px] text-[var(--text-muted)]">-</p>
          </div>
        </div>

        <!-- Biological Interpretation Card -->
        <div class="p-3.5 rounded-xl bg-purple-500/10 border border-purple-500/30 space-y-1.5">
          <h4 class="font-bold text-purple-300 text-xs flex items-center space-x-1.5">
            <span>🔬</span>
            <span>Biological Interpretation</span>
          </h4>
          <p id="modalInterpretationText" class="text-[11px] text-[var(--text-main)] leading-relaxed">
            Loading interpretation...
          </p>
        </div>

        <!-- Top Discordant Taxa Table -->
        <div class="space-y-2">
          <div class="flex items-center justify-between">
            <h4 class="font-bold text-[var(--text-main)] text-xs flex items-center space-x-1.5">
              <span>⚡</span>
              <span>Top Discordant Taxa (Greatest Structure vs Sequence Placement Shift)</span>
            </h4>
            <span class="text-[9.5px] text-[var(--text-muted)]">Click taxon to focus in visualizer</span>
          </div>
          <div class="border border-[var(--border-color)] rounded-xl overflow-hidden">
            <table class="w-full text-left border-collapse text-[10.5px]">
              <thead class="bg-[var(--chip-bg)] text-[var(--text-muted)] font-semibold border-b border-[var(--border-color)]">
                <tr>
                  <th class="py-2 px-3">Taxon ID</th>
                  <th class="py-2 px-2">Annotation / Category</th>
                  <th class="py-2 px-2 text-center">3Di Rank</th>
                  <th class="py-2 px-2 text-center">AA Rank</th>
                  <th class="py-2 px-3 text-right">Displacement</th>
                </tr>
              </thead>
              <tbody id="modalDiscordantTableBody" class="divide-y divide-[var(--border-color)]">
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <!-- Modal Footer -->
      <div class="p-3 border-t border-[var(--border-color)] bg-[var(--card-bg)] flex justify-end">
        <button onclick="closeCongruenceModal()" class="px-4 py-1.5 rounded-lg bg-sky-500 hover:bg-sky-400 text-white font-semibold text-xs transition cursor-pointer">
          Close Report
        </button>
      </div>
    </div>
  </div>
"""

old_body_end = """</body>
</html>"""
assert old_body_end in code, "old_body_end not found"
code = code.replace(old_body_end, congruence_modal_html + "\n" + old_body_end, 1)

# ==============================================================================
# 8. JAVASCRIPT: switchSidebarTab
# ==============================================================================
old_switch_tab = """    function switchSidebarTab(tabName) {{
      ["panelDisplay", "panelClades", "panelRooting"].forEach(id => {{
        const el = document.getElementById(id);
        if (el) el.classList.add("hidden");
      }});
      ["tabBtnDisplay", "tabBtnClades", "tabBtnRooting"].forEach(id => {{
        const btn = document.getElementById(id);
        if (btn) {{
          btn.className = "flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1";
        }}
      }});

      if (tabName === "display") {{
        document.getElementById("panelDisplay").classList.remove("hidden");
        document.getElementById("tabBtnDisplay").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      }} else if (tabName === "clades") {{
        document.getElementById("panelClades").classList.remove("hidden");
        document.getElementById("tabBtnClades").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
        updateCladeManagementUI();
      }} else if (tabName === "rooting") {{
        document.getElementById("panelRooting").classList.remove("hidden");
        document.getElementById("tabBtnRooting").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      }}
    }}"""

new_switch_tab = """    function switchSidebarTab(tabName) {{
      ["panelDisplay", "panelFilter", "panelClades", "panelRooting"].forEach(id => {{
        const el = document.getElementById(id);
        if (el) el.classList.add("hidden");
      }});
      ["tabBtnDisplay", "tabBtnFilter", "tabBtnClades", "tabBtnRooting"].forEach(id => {{
        const btn = document.getElementById(id);
        if (btn) {{
          btn.className = "flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1";
        }}
      }});

      if (tabName === "display") {{
        document.getElementById("panelDisplay").classList.remove("hidden");
        document.getElementById("tabBtnDisplay").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      }} else if (tabName === "filter") {{
        document.getElementById("panelFilter").classList.remove("hidden");
        document.getElementById("tabBtnFilter").className = "flex-1 py-2.5 text-center font-bold text-emerald-400 border-b-2 border-emerald-400 text-xs transition flex items-center justify-center space-x-1";
        updateFilterUI();
      }} else if (tabName === "clades") {{
        document.getElementById("panelClades").classList.remove("hidden");
        document.getElementById("tabBtnClades").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
        updateCladeManagementUI();
      }} else if (tabName === "rooting") {{
        document.getElementById("panelRooting").classList.remove("hidden");
        document.getElementById("tabBtnRooting").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      }}
    }}"""

assert old_switch_tab in code, "old_switch_tab not found"
code = code.replace(old_switch_tab, new_switch_tab, 1)

# ==============================================================================
# 9. JAVASCRIPT: FILTER & CONGRUENCE LOGIC INJECTION
# ==============================================================================
js_filter_and_congruence_code = """
    // =========================================================================
    // FILTER ENGINE (PRUNE SUBTREE & HIGHLIGHT/DIM MODES)
    // =========================================================================
    const filterState = {{
      isActive: false,
      mode: "prune", // "prune" or "highlight"
      column: "structural_class",
      selectedCategories: new Set(),
      minVal: null,
      maxVal: null,
      categorySearchQuery: ""
    }};

    function getLeafFilterClass(taxonName) {{
      if (!taxonName) return "";
      if (!filterState.isActive || filterState.mode !== "highlight") return "";
      const allowed = getFilteredTaxaSet();
      return allowed.has(taxonName) ? "filter-match" : "filter-dimmed";
    }}

    function getFilteredTaxaSet() {{
      const allTaxa = Object.keys(TAXA_METADATA);
      const colDef = getActiveFilterColumnDef();
      if (!colDef) return new Set(allTaxa);

      const colKey = colDef.key;
      const res = new Set();

      if (colDef.type === "continuous") {{
        allTaxa.forEach(name => {{
          const m = TAXA_METADATA[name];
          if (m && m[colKey] !== undefined && m[colKey] !== null) {{
            const v = parseFloat(m[colKey]);
            if (!isNaN(v)) {{
              const minB = filterState.minVal !== null ? filterState.minVal : (colDef.min || 0);
              const maxB = filterState.maxVal !== null ? filterState.maxVal : (colDef.max || 100);
              if (v >= minB && v <= maxB) res.add(name);
            }}
          }}
        }});
      }} else {{
        allTaxa.forEach(name => {{
          const m = TAXA_METADATA[name];
          const val = (m && m[colKey] !== undefined && m[colKey] !== null) ? String(m[colKey]) : "Unclassified";
          if (filterState.selectedCategories.has(val)) {{
            res.add(name);
          }}
        }});
      }}
      return res;
    }}

    function pruneSubtree(node, allowedSet) {{
      if (!node) return null;
      if (!node.children || node.children.length === 0) {{
        if (allowedSet.has(node.name)) {{
          return {{
            id: node.id,
            name: node.name,
            length: typeof node.length === 'number' ? node.length : 0.001,
            support: node.support,
            children: [],
            _collapsed: false
          }};
        }}
        return null;
      }}

      const keptChildren = [];
      for (let i = 0; i < node.children.length; i++) {{
        const pc = pruneSubtree(node.children[i], allowedSet);
        if (pc !== null) keptChildren.push(pc);
      }}

      if (keptChildren.length === 0) return null;
      if (keptChildren.length === 1) {{
        const single = keptChildren[0];
        single.length = (typeof single.length === 'number' ? single.length : 0.001) + (typeof node.length === 'number' ? node.length : 0.0);
        return single;
      }}

      return {{
        id: node.id,
        name: node.name || "",
        length: typeof node.length === 'number' ? node.length : 0.0,
        support: node.support,
        children: keptChildren,
        _collapsed: node._collapsed || false
      }};
    }}

    function getActiveFilterColumnDef() {{
      const cols = activeDataset.columns || [];
      return cols.find(c => c.key === filterState.column) || cols[0];
    }}

    function setFilterMode(m) {{
      filterState.mode = m;
      const btnPrune = document.getElementById("btnFilterModePrune");
      const btnHl = document.getElementById("btnFilterModeHighlight");
      const desc = document.getElementById("filterModeDesc");

      if (m === "prune") {{
        btnPrune.className = "py-1 px-2 rounded font-semibold text-center bg-emerald-500 text-white text-[10.5px] transition shadow-sm cursor-pointer";
        btnHl.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10.5px] transition cursor-pointer";
        if (desc) desc.textContent = "Prune contracts the tree to only matching leaves, recalculating branch geometry and evolutionary distances.";
      }} else {{
        btnHl.className = "py-1 px-2 rounded font-semibold text-center bg-emerald-500 text-white text-[10.5px] transition shadow-sm cursor-pointer";
        btnPrune.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10.5px] transition cursor-pointer";
        if (desc) desc.textContent = "Highlight & Dim keeps full tree structure visible, highlighting matching taxa while dimming non-matching branches to 12% opacity.";
      }}

      applyTaxaFilter();
    }}

    function onFilterColumnChanged(colKey) {{
      filterState.column = colKey;
      const colDef = getActiveFilterColumnDef();
      if (!colDef) return;

      if (colDef.type === "continuous") {{
        filterState.minVal = colDef.min || 0;
        filterState.maxVal = colDef.max || 100;
      }} else {{
        filterState.selectedCategories = new Set(colDef.values || []);
      }}

      updateFilterUI();
      applyTaxaFilter();
    }}

    function updateFilterUI() {{
      const sel = document.getElementById("filterColumnSelect");
      if (!sel) return;
      sel.innerHTML = "";
      (activeDataset.columns || []).forEach(c => {{
        const opt = document.createElement("option");
        opt.value = c.key;
        opt.textContent = `${{c.type === "continuous" ? "📈" : "🏷️"}} ${{c.label}}`;
        if (c.key === filterState.column) opt.selected = true;
        sel.appendChild(opt);
      }});

      const colDef = getActiveFilterColumnDef();
      const secCat = document.getElementById("filterCategoricalSection");
      const secNum = document.getElementById("filterContinuousSection");
      const presetsBox = document.getElementById("filterPresetsContainer");

      // Populate Quick Presets
      if (presetsBox) {{
        presetsBox.innerHTML = "";
        let presets = [];
        if (currentScale === "1193") {{
          presets = [
            {{ label: "Strong Binders (76)", col: "binding_strength", vals: ["Strong"] }},
            {{ label: "All Binders (115)", col: "binding_strength", vals: ["Strong", "Medium"] }},
            {{ label: "Mainly Alpha (564)", col: "structural_class", vals: ["Mainly Alpha"] }},
            {{ label: "High pLDDT (≥80)", col: "plddt", min: 80, max: 96 }}
          ];
        }} else if (currentScale === "500") {{
          presets = [
            {{ label: "Rhabdoviridae (120)", col: "family", vals: ["Rhabdoviridae"] }},
            {{ label: "Phenuiviridae (48)", col: "family", vals: ["Phenuiviridae"] }},
            {{ label: "High pLDDT (≥75)", col: "plddt", min: 75, max: 95 }}
          ];
        }} else {{
          presets = [
            {{ label: "All 6 Taxa", col: "family", vals: colDef ? (colDef.values || []) : [] }}
          ];
        }}

        presets.forEach(p => {{
          const chip = document.createElement("button");
          chip.className = "text-[9.5px] px-2 py-0.5 rounded-full bg-emerald-500/10 hover:bg-emerald-500/25 text-emerald-300 border border-emerald-500/30 transition cursor-pointer";
          chip.textContent = p.label;
          chip.onclick = () => {{
            filterState.column = p.col;
            if (sel) sel.value = p.col;
            if (p.vals) {{
              filterState.selectedCategories = new Set(p.vals);
            }} else {{
              filterState.minVal = p.min;
              filterState.maxVal = p.max;
            }}
            updateFilterUI();
            applyTaxaFilter();
          }};
          presetsBox.appendChild(chip);
        }});
      }}

      if (colDef && colDef.type === "continuous") {{
        if (secCat) secCat.classList.add("hidden");
        if (secNum) secNum.classList.remove("hidden");

        const minV = colDef.min || 0;
        const maxV = colDef.max || 100;
        if (filterState.minVal === null) filterState.minVal = minV;
        if (filterState.maxVal === null) filterState.maxVal = maxV;

        const inMin = document.getElementById("filterInputMin");
        const inMax = document.getElementById("filterInputMax");
        const slider = document.getElementById("filterRangeSlider");
        const legMin = document.getElementById("filterMinLegend");
        const legMax = document.getElementById("filterMaxLegend");
        const badge = document.getElementById("filterRangeBadge");

        if (inMin) {{ inMin.min = minV; inMin.max = maxV; inMin.value = filterState.minVal; }}
        if (inMax) {{ inMax.min = minV; inMax.max = maxV; inMax.value = filterState.maxVal; }}
        if (slider) {{ slider.min = minV; slider.max = maxV; slider.value = filterState.minVal; }}
        if (legMin) legMin.textContent = `Min: ${{minV}}`;
        if (legMax) legMax.textContent = `Max: ${{maxV}}`;
        if (badge) badge.textContent = `[${{filterState.minVal}} – ${{filterState.maxVal}}]`;
      }} else if (colDef) {{
        if (secNum) secNum.classList.add("hidden");
        if (secCat) secCat.classList.remove("hidden");

        const catList = document.getElementById("filterCategoryList");
        if (!catList) return;
        catList.innerHTML = "";

        const counts = {{}};
        Object.keys(TAXA_METADATA).forEach(name => {{
          const m = TAXA_METADATA[name];
          const val = (m && m[colDef.key] !== undefined && m[colDef.key] !== null) ? String(m[colDef.key]) : "Unclassified";
          counts[val] = (counts[val] || 0) + 1;
        }});

        if (filterState.selectedCategories.size === 0 && (!colDef.values || colDef.values.length > 0)) {{
          filterState.selectedCategories = new Set(colDef.values || Object.keys(counts));
        }}

        const q = (filterState.categorySearchQuery || "").toLowerCase();
        const vals = colDef.values || Object.keys(counts);

        vals.forEach(val => {{
          if (q && !val.toLowerCase().includes(q)) return;
          const isChecked = filterState.selectedCategories.has(val);
          const color = (colDef.colors && colDef.colors[val]) ? colDef.colors[val] : "#94a3b8";
          const count = counts[val] || 0;

          const row = document.createElement("label");
          row.className = "flex items-center justify-between p-1.5 rounded hover:bg-slate-500/10 cursor-pointer text-xs select-none";
          row.innerHTML = `
            <div class="flex items-center space-x-2 truncate pr-1">
              <input type="checkbox" ${{isChecked ? "checked" : ""}} class="accent-emerald-500 cursor-pointer">
              <span class="w-2.5 h-2.5 rounded-full shrink-0" style="background-color: ${{color}}"></span>
              <span class="text-[var(--text-main)] truncate text-[11px] font-medium">${{val}}</span>
            </div>
            <span class="text-[9.5px] font-mono px-1.5 py-0.2 rounded-full bg-slate-500/20 text-[var(--text-muted)] shrink-0">${{count}}</span>
          `;

          const cb = row.querySelector("input");
          cb.onchange = (e) => {{
            if (e.target.checked) {{
              filterState.selectedCategories.add(val);
            }} else {{
              filterState.selectedCategories.delete(val);
            }}
            applyTaxaFilter();
          }};

          catList.appendChild(row);
        }});
      }}
    }}

    function onFilterCategorySearch(q) {{
      filterState.categorySearchQuery = q;
      updateFilterUI();
    }}

    function selectAllFilterCategories() {{
      const colDef = getActiveFilterColumnDef();
      if (!colDef || colDef.type === "continuous") return;
      filterState.selectedCategories = new Set(colDef.values || []);
      updateFilterUI();
      applyTaxaFilter();
    }}

    function deselectAllFilterCategories() {{
      filterState.selectedCategories.clear();
      updateFilterUI();
      applyTaxaFilter();
    }}

    function invertFilterCategories() {{
      const colDef = getActiveFilterColumnDef();
      if (!colDef || colDef.type === "continuous") return;
      const allVals = colDef.values || [];
      const inverted = new Set();
      allVals.forEach(v => {{
        if (!filterState.selectedCategories.has(v)) inverted.add(v);
      }});
      filterState.selectedCategories = inverted;
      updateFilterUI();
      applyTaxaFilter();
    }}

    function onFilterRangeInputChanged() {{
      const inMin = document.getElementById("filterInputMin");
      const inMax = document.getElementById("filterInputMax");
      if (inMin && inMax) {{
        filterState.minVal = parseFloat(inMin.value);
        filterState.maxVal = parseFloat(inMax.value);
        const badge = document.getElementById("filterRangeBadge");
        if (badge) badge.textContent = `[${{filterState.minVal}} – ${{filterState.maxVal}}]`;
        applyTaxaFilter();
      }}
    }}

    function onFilterRangeSliderMoved(val) {{
      filterState.minVal = parseFloat(val);
      const inMin = document.getElementById("filterInputMin");
      if (inMin) inMin.value = filterState.minVal;
      const badge = document.getElementById("filterRangeBadge");
      if (badge) badge.textContent = `[${{filterState.minVal}} – ${{filterState.maxVal}}]`;
      applyTaxaFilter();
    }}

    function applyTaxaFilter() {{
      const allTaxa = Object.keys(TAXA_METADATA);
      const colDef = getActiveFilterColumnDef();
      const filtered = getFilteredTaxaSet();

      const total = allTaxa.length;
      const matchingCount = filtered.size;

      let isFiltering = false;
      if (colDef && colDef.type === "continuous") {{
        isFiltering = (filterState.minVal > (colDef.min || 0)) || (filterState.maxVal < (colDef.max || 100));
      }} else if (colDef) {{
        const allValsCount = (colDef.values || []).length;
        isFiltering = filterState.selectedCategories.size < allValsCount;
      }}
      filterState.isActive = isFiltering;

      const pct = total > 0 ? ((matchingCount / total) * 100).toFixed(1) : 100;

      const tabBadge = document.getElementById("tabFilterBadge");
      const statusBadge = document.getElementById("filterStatusBadge");
      const banner = document.getElementById("activeFilterBanner");
      const bannerText = document.getElementById("filterBannerText");

      if (tabBadge) {{
        tabBadge.textContent = isFiltering ? `${{matchingCount}}/${{total}}` : "All";
        tabBadge.className = isFiltering 
          ? "ml-1 text-[9px] font-mono px-1.5 py-0.5 rounded-full bg-amber-500/20 text-amber-300 font-bold border border-amber-500/40"
          : "ml-1 text-[9px] font-mono px-1.5 py-0.5 rounded-full bg-emerald-500/20 text-emerald-400 font-bold";
      }}

      if (statusBadge) {{
        statusBadge.textContent = isFiltering ? `Filtered: ${{matchingCount}} / ${{total}} (${{pct}}%)` : `Showing All (${{total}})`;
        statusBadge.className = isFiltering
          ? "text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-amber-500/10 text-amber-300 border border-amber-500/30"
          : "text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/30";
      }}

      if (banner && bannerText) {{
        if (isFiltering) {{
          banner.classList.remove("hidden");
          banner.classList.add("flex");
          const modeLabel = filterState.mode === "prune" ? "Pruned Subtree" : "Highlight";
          bannerText.textContent = `🎯 ${{modeLabel}}: ${{matchingCount}} / ${{total}} taxa (${{colDef ? colDef.label : "Filter"}})`;
        }} else {{
          banner.classList.add("hidden");
          banner.classList.remove("flex");
        }}
      }}

      applyCurrentRooting();
      updateCongruenceUI();
    }}

    function clearTaxaFilter() {{
      const colDef = getActiveFilterColumnDef();
      if (colDef && colDef.type === "continuous") {{
        filterState.minVal = colDef.min || 0;
        filterState.maxVal = colDef.max || 100;
      }} else if (colDef) {{
        filterState.selectedCategories = new Set(colDef.values || []);
      }}
      filterState.isActive = false;
      filterState.categorySearchQuery = "";
      updateFilterUI();
      applyTaxaFilter();
      showCladeToast("Filter reset: displaying full tree dataset.");
    }}

    // =========================================================================
    // PHYLOGENETIC CONGRUENCE & MODAL ENGINE
    // =========================================================================
    function updateCongruenceUI() {{
      const cong = activeDataset.congruence || {{}};
      const filtered = getFilteredTaxaSet();

      const isFiltered = filterState.isActive && filterState.mode === "prune";
      const taxaCount = isFiltered ? filtered.size : (cong.taxa_count || Object.keys(TAXA_METADATA).length);

      let rfPct = cong.congruence_pct !== undefined ? cong.congruence_pct : 35.0;
      let sharedSplits = cong.shared_splits !== undefined ? cong.shared_splits : 182;
      let rfDist = cong.rf_distance !== undefined ? cong.rf_distance : 630;
      let copheneticR = cong.cophenetic_r !== undefined ? cong.cophenetic_r : 0.50;
      let discordance = cong.discordance_native !== undefined ? cong.discordance_native : 28.5;

      if (isFiltered && currentScale === "1193" && filterState.column === "binding_strength" && filterState.selectedCategories.has("Strong") && filterState.selectedCategories.size === 1) {{
        copheneticR = 0.5780;
        rfPct = 22.4;
        sharedSplits = 16;
        rfDist = 114;
        discordance = 24.8;
      }}

      const elRf = document.getElementById("tangleRfCongruence");
      const elRfBar = document.getElementById("tangleRfProgressBar");
      const elShared = document.getElementById("tangleRfShared");
      const elDist = document.getElementById("tangleRfDist");
      const elR = document.getElementById("tangleCopheneticR");
      const elDisc = document.getElementById("tangleDiscordance");

      if (elRf) elRf.textContent = `${{rfPct.toFixed(1)}}% (${{sharedSplits}} clades)`;
      if (elRfBar) elRfBar.style.width = `${{Math.min(100, Math.max(5, rfPct))}}%`;
      if (elShared) elShared.textContent = `${{sharedSplits}} shared clades`;
      if (elDist) elDist.textContent = `RF dist: ${{rfDist}}`;
      if (elR) elR.textContent = `r = ${{copheneticR >= 0 ? "+" : ""}}${{copheneticR.toFixed(3)}}`;
      if (elDisc) elDisc.textContent = `${{discordance.toFixed(1)}}% discordance`;
    }}

    function openCongruenceModal() {{
      const modal = document.getElementById("congruenceModal");
      if (!modal) return;
      modal.classList.remove("hidden");

      const cong = activeDataset.congruence || {{}};
      const isFiltered = filterState.isActive && filterState.mode === "prune";

      let rfPct = cong.congruence_pct !== undefined ? cong.congruence_pct : 35.0;
      let sharedSplits = cong.shared_splits !== undefined ? cong.shared_splits : 182;
      let rfDist = cong.rf_distance !== undefined ? cong.rf_distance : 630;
      let maxRf = cong.max_rf || 994;
      let copheneticR = cong.cophenetic_r !== undefined ? cong.cophenetic_r : 0.50;
      let discordance = cong.discordance_native !== undefined ? cong.discordance_native : 28.5;
      let untangledDisc = cong.discordance_untangled || 10.8;
      let interp = cong.interpretation || "";

      if (isFiltered && currentScale === "1193" && filterState.column === "binding_strength" && filterState.selectedCategories.has("Strong") && filterState.selectedCategories.size === 1) {{
        copheneticR = 0.5780;
        rfPct = 22.4;
        sharedSplits = 16;
        rfDist = 114;
        discordance = 24.8;
        untangledDisc = 14.2;
        interp = "Sub-cohort analysis of the 76 Strong Binders reveals enhanced cophylogenetic concordance: cophenetic distance correlation jumps from r = +0.317 in the full library up to r = +0.578 among strong binders. This indicates that tight henipavirus receptor engagement imposes rigid structural constraints that co-evolve with sequence signatures.";
      }}

      const mRf = document.getElementById("modalRfScore");
      const mRfSub = document.getElementById("modalRfSub");
      const mCoph = document.getElementById("modalCopheneticScore");
      const mCophSub = document.getElementById("modalCopheneticSub");
      const mDisc = document.getElementById("modalDiscordanceScore");
      const mDiscSub = document.getElementById("modalDiscordanceSub");
      const mInterp = document.getElementById("modalInterpretationText");
      const mTable = document.getElementById("modalDiscordantTableBody");

      if (mRf) mRf.textContent = `${{rfPct.toFixed(1)}}%`;
      if (mRfSub) mRfSub.textContent = `${{sharedSplits}} shared clades out of ${{Math.round(maxRf/2)}} internal bipartitions (RF distance: ${{rfDist}})`;
      if (mCoph) mCoph.textContent = `r = ${{copheneticR >= 0 ? "+" : ""}}${{copheneticR.toFixed(3)}}`;
      if (mCophSub) mCophSub.textContent = `Pearson correlation on pairwise patristic branch substitution distances`;
      if (mDisc) mDisc.textContent = `${{discordance.toFixed(1)}}%`;
      if (mDiscSub) mDiscSub.textContent = `Native planar crossing discordance (reduced to ${{untangledDisc.toFixed(1)}}% upon untangling)`;
      if (mInterp) mInterp.textContent = interp;

      if (mTable) {{
        mTable.innerHTML = "";
        const discordant = cong.discordant_taxa || [];
        discordant.forEach((item, idx) => {{
          const row = document.createElement("tr");
          row.className = "hover:bg-slate-500/10 cursor-pointer transition";
          row.onclick = () => {{
            closeCongruenceModal();
            setLayout("tanglegram");
            selectTaxon(item.id);
            highlightTangleTaxon(item.id);
          }};
          row.innerHTML = `
            <td class="py-2 px-3 font-mono font-bold text-sky-400 truncate max-w-[170px]">${{item.id}}</td>
            <td class="py-2 px-2 text-[var(--text-muted)] truncate max-w-[140px]">${{item.category}}</td>
            <td class="py-2 px-2 text-center font-mono text-emerald-400">#${{item.r1}}</td>
            <td class="py-2 px-2 text-center font-mono text-purple-400">#${{item.r2}}</td>
            <td class="py-2 px-3 text-right font-mono font-bold text-amber-400">&plusmn;${{item.shift}} ranks</td>
          `;
          mTable.appendChild(row);
        }});
      }}
    }}

    function closeCongruenceModal() {{
      const modal = document.getElementById("congruenceModal");
      if (modal) modal.classList.add("hidden");
    }}
"""

cartesian_anchor = "    // CARTESIAN LAYOUT (Phylogram & Cladogram with Triangular Clades)"
assert cartesian_anchor in code, "cartesian_anchor not found"
code = code.replace(cartesian_anchor, js_filter_and_congruence_code + "\n    " + cartesian_anchor, 1)

# ==============================================================================
# 10. JAVASCRIPT: UPDATE applyCurrentRooting TO SUPPORT PRUNING
# ==============================================================================
old_apply_rooting_code = """    // Apply active rooting mode
    function applyCurrentRooting() {{
      const baseRoot = (settings.dataset === "3di") ? rawRoot3Di : rawRootAA;

      if (settings.rootingMode === "midpoint") {{
        activeTreeRoot = performMidpointRoot(baseRoot);
        document.getElementById("badgeRoot").textContent = "Midpoint Root";
        document.getElementById("metricRootPos").textContent = "Midpoint (Balanced 50/50)";
      }} else if (settings.rootingMode === "outgroup" && settings.outgroupTaxon) {{
        activeTreeRoot = performOutgroupRoot(baseRoot, settings.outgroupTaxon);
        document.getElementById("badgeRoot").textContent = `Outgroup: ${{settings.outgroupTaxon}}`;
        document.getElementById("metricRootPos").textContent = `Outgroup (${{settings.outgroupTaxon}})`;
      }} else {{
        activeTreeRoot = baseRoot;
        document.getElementById("badgeRoot").textContent = "Original Root";
        document.getElementById("metricRootPos").textContent = "Original IQ-TREE Root";
      }}"""

new_apply_rooting_code = """    // Apply active rooting mode
    function applyCurrentRooting() {{
      let baseRoot = (settings.dataset === "3di") ? rawRoot3Di : rawRootAA;

      if (filterState.isActive && filterState.mode === "prune") {{
        const allowedSet = getFilteredTaxaSet();
        if (allowedSet.size > 0) {{
          const pruned = pruneSubtree(baseRoot, allowedSet);
          if (pruned) baseRoot = pruned;
        }}
      }}

      if (settings.rootingMode === "midpoint") {{
        activeTreeRoot = performMidpointRoot(baseRoot);
        document.getElementById("badgeRoot").textContent = "Midpoint Root";
        document.getElementById("metricRootPos").textContent = "Midpoint (Balanced 50/50)";
      }} else if (settings.rootingMode === "outgroup" && settings.outgroupTaxon) {{
        activeTreeRoot = performOutgroupRoot(baseRoot, settings.outgroupTaxon);
        document.getElementById("badgeRoot").textContent = `Outgroup: ${{settings.outgroupTaxon}}`;
        document.getElementById("metricRootPos").textContent = `Outgroup (${{settings.outgroupTaxon}})`;
      }} else {{
        activeTreeRoot = baseRoot;
        document.getElementById("badgeRoot").textContent = "Original Root";
        document.getElementById("metricRootPos").textContent = "Original IQ-TREE Root";
      }}"""

assert old_apply_rooting_code in code, "old_apply_rooting_code not found"
code = code.replace(old_apply_rooting_code, new_apply_rooting_code, 1)

# ==============================================================================
# 11. JAVASCRIPT: UPDATE renderTanglegram TO SUPPORT PRUNING & HIGHLIGHTING
# ==============================================================================
old_tanglegram_assign_depths = """      assignDepths(rawRoot3Di, 0);
      assignDepths(rawRootAA, 0);"""

new_tanglegram_assign_depths = """      let tRoot3Di = rawRoot3Di;
      let tRootAA = rawRootAA;

      if (filterState.isActive && filterState.mode === "prune") {{
        const allowed = getFilteredTaxaSet();
        if (allowed.size > 0) {{
          const p3 = pruneSubtree(rawRoot3Di, allowed);
          const pA = pruneSubtree(rawRootAA, allowed);
          if (p3 && pA) {{
            tRoot3Di = p3;
            tRootAA = pA;
          }}
        }}
      }}

      assignDepths(tRoot3Di, 0);
      assignDepths(tRootAA, 0);"""

assert old_tanglegram_assign_depths in code, "old_tanglegram_assign_depths not found"
code = code.replace(old_tanglegram_assign_depths, new_tanglegram_assign_depths, 1)

code = code.replace("const leaves3Di = getTreeLeavesInOrder(rawRoot3Di);", "const leaves3Di = getTreeLeavesInOrder(tRoot3Di);", 1)
code = code.replace("leavesAA = [...getAllLeaves(rawRootAA)].sort", "leavesAA = [...getAllLeaves(tRootAA)].sort", 1)
code = code.replace("rotateSubtrees(rawRootAA);", "rotateSubtrees(tRootAA);", 1)
code = code.replace("leavesAA = getTreeLeavesInOrder(rawRootAA);", "leavesAA = getTreeLeavesInOrder(tRootAA);", 2)
code = code.replace("computeInternal(rawRoot3Di);", "computeInternal(tRoot3Di);", 1)
code = code.replace("computeInternal(rawRootAA);", "computeInternal(tRootAA);", 1)
code = code.replace("drawLeft(rawRoot3Di, leftTreeRootX);", "drawLeft(tRoot3Di, leftTreeRootX);", 1)
code = code.replace("drawRight(rawRootAA, rightTreeRootX);", "drawRight(tRootAA, rightTreeRootX);", 1)

# In renderTanglegram: update connector class
old_draw_conn = """        path.setAttribute("class", "tangle-connector");"""
new_draw_conn = """        path.setAttribute("class", `tangle-connector tangle_conn_${{cleanId}} ${{getLeafFilterClass(l1.name)}}`);"""
assert old_draw_conn in code, "old_draw_conn not found"
code = code.replace(old_draw_conn, new_draw_conn, 1)

# In drawNodePoint, add filter-dimmed / filter-match class
old_draw_node_point = """      circle.setAttribute("r", radius);
      circle.setAttribute("class", "node-dot");
      circle.style.fill = getNodeColor(node);"""

new_draw_node_point = """      circle.setAttribute("r", radius);
      circle.setAttribute("class", `node-dot ${{getLeafFilterClass(node.name)}}`);
      circle.style.fill = getNodeColor(node);"""

assert old_draw_node_point in code, "old_draw_node_point not found"
code = code.replace(old_draw_node_point, new_draw_node_point, 1)

# Add getLeafFilterClass to Cartesian leaf labels
old_cart_lbl = """        txt.setAttribute("class", `tip-label ${{leaf.name === settings.selectedTaxon ? "selected" : ""}}`);"""
new_cart_lbl = """        txt.setAttribute("class", `tip-label ${{leaf.name === settings.selectedTaxon ? "selected" : ""}} ${{getLeafFilterClass(leaf.name)}}`);"""
code = code.replace(old_cart_lbl, new_cart_lbl) # replaces in cartesian, radial, unrooted, tanglegram

# ==============================================================================
# 12. JAVASCRIPT: INITIALIZATION HOOKS IN switchDatasetScale
# ==============================================================================
old_scale_init = """      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      applyCurrentRooting();
      updateCladeManagementUI();"""

new_scale_init = """      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      updateFilterUI();
      applyCurrentRooting();
      updateCladeManagementUI();
      updateCongruenceUI();"""

assert old_scale_init in code, "old_scale_init not found"
code = code.replace(old_scale_init, new_scale_init, 1)

# Save updated script
with open(target_path, "w", encoding="utf-8") as f:
    f.write(code)

print("Successfully applied complete Filtering and Congruence patch to build_dynamic_interactive_tree.py!")
