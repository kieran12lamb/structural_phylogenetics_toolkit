#!/usr/bin/env python3
"""Patch build_dynamic_interactive_tree.py to add:
1. Dynamic Metadata-Based Tree Filtering (Prune Subtree & Highlight/Dim modes, categorical & continuous)
2. Tanglegram Quantitative Phylogenetic Congruence Engine (Robinson-Foulds distance, Cophenetic correlation, Entanglement)
3. Congruence Analysis Report modal with top discordant taxa and biological interpretation
"""

import json
from pathlib import Path

target_file = Path("/Users/kieran.lamb/Github/Skills_Hackathon/scripts/build_dynamic_interactive_tree.py")
with open(target_file, "r") as f:
    text = f.read()

# 1. Update the TAB BUTTONS: Add 4th tab 'tabBtnFilter'
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

assert old_tab_buttons in text, "old_tab_buttons not found"
text = text.replace(old_tab_buttons, new_tab_buttons, 1)

# 2. Add Congruence metrics card inside tanglegramOptions
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

              <div class="space-y-1.5 bg-black/25 p-2.5 rounded-lg border border-purple-500/20">
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

assert old_tanglegram_options in text, "old_tanglegram_options not found"
text = text.replace(old_tanglegram_options, new_tanglegram_options, 1)

# 3. Add Filter Panel right before TAB 2 (Clades)
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
              <button onclick="clearTaxaFilter()" class="flex-1 py-1.5 px-2 rounded border border-[var(--border-color)] hover:bg-slate-500/20 text-[10.5px] font-semibold text-[var(--text-main)] transition text-center shadow-sm">
                ✕ Reset All Filters
              </button>
            </div>
          </div>

          <!-- Filter Action Mode Switch -->
          <div class="space-y-1.5 bg-[var(--chip-bg)] p-2.5 rounded-lg border border-[var(--border-color)]">
            <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block text-[10px]">Filter Action Mode</label>
            <div class="grid grid-cols-2 gap-1 bg-[var(--input-bg)] p-1 rounded-md border border-[var(--border-color)]">
              <button id="btnFilterModePrune" onclick="setFilterMode('prune')" class="py-1 px-2 rounded font-semibold text-center bg-emerald-500 text-white text-[10.5px] transition shadow-sm">
                ✂️ Prune Subtree
              </button>
              <button id="btnFilterModeHighlight" onclick="setFilterMode('highlight')" class="py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10.5px] transition">
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
assert old_tab_clades_marker in text, "old_tab_clades_marker not found"
text = text.replace(old_tab_clades_marker, filter_panel_html + "\n        " + old_tab_clades_marker, 1)

# 4. Add Floating Filter Banner
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
        <button onclick="clearTaxaFilter()" class="ml-1 text-slate-400 hover:text-white font-bold text-xs bg-slate-800/80 hover:bg-rose-600/80 px-2 py-0.5 rounded-full transition">&times; Clear</button>
      </div>

      <!-- SVG Canvas -->
      <svg id="treeSvg" class="w-full h-full block"></svg>"""

assert old_canvas_section in text, "old_canvas_section not found"
text = text.replace(old_canvas_section, new_canvas_section, 1)

# 5. Add Congruence Modal dialog before </body>
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
assert old_body_end in text, "old_body_end not found"
text = text.replace(old_body_end, congruence_modal_html + "\n" + old_body_end, 1)

# 6. Update switchSidebarTab JavaScript function to handle 'filter' tab
old_switch_tab = """    function switchSidebarTab(tabName) {
      ["panelDisplay", "panelClades", "panelRooting"].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.classList.add("hidden");
      });
      ["tabBtnDisplay", "tabBtnClades", "tabBtnRooting"].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
          el.className = "flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1";
        }
      });
      if (tabName === "display") {
        document.getElementById("panelDisplay").classList.remove("hidden");
        document.getElementById("tabBtnDisplay").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      } else if (tabName === "clades") {
        document.getElementById("panelClades").classList.remove("hidden");
        document.getElementById("tabBtnClades").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      } else if (tabName === "rooting") {
        document.getElementById("panelRooting").classList.remove("hidden");
        document.getElementById("tabBtnRooting").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      }
    }"""

new_switch_tab = """    function switchSidebarTab(tabName) {
      ["panelDisplay", "panelFilter", "panelClades", "panelRooting"].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.classList.add("hidden");
      });
      ["tabBtnDisplay", "tabBtnFilter", "tabBtnClades", "tabBtnRooting"].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
          el.className = "flex-1 py-2.5 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-xs transition flex items-center justify-center space-x-1";
        }
      });
      if (tabName === "display") {
        document.getElementById("panelDisplay").classList.remove("hidden");
        document.getElementById("tabBtnDisplay").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      } else if (tabName === "filter") {
        document.getElementById("panelFilter").classList.remove("hidden");
        document.getElementById("tabBtnFilter").className = "flex-1 py-2.5 text-center font-bold text-emerald-400 border-b-2 border-emerald-400 text-xs transition flex items-center justify-center space-x-1";
      } else if (tabName === "clades") {
        document.getElementById("panelClades").classList.remove("hidden");
        document.getElementById("tabBtnClades").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      } else if (tabName === "rooting") {
        document.getElementById("panelRooting").classList.remove("hidden");
        document.getElementById("tabBtnRooting").className = "flex-1 py-2.5 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-xs transition flex items-center justify-center space-x-1";
      }
    }"""

assert old_switch_tab in text, "old_switch_tab not found"
text = text.replace(old_switch_tab, new_switch_tab, 1)

# Write modified file
with open(target_file, "w") as f:
    f.write(text)

print("HTML structure successfully patched with Filter Tab, Banner, and Congruence Modal!")
