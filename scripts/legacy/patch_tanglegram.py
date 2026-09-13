#!/usr/bin/env python3
"""Patch build_dynamic_interactive_tree.py with multi-mode Tanglegram support:
- True Topology (independent branch orders on both sides, real discordance & crossings)
- Untangle (Barycenter subtree rotation optimization without altering clades)
- Aligned / Parallel (1-to-1 matching)
"""

from pathlib import Path

target_file = Path("/Users/kieran.lamb/Github/Skills_Hackathon/scripts/build_dynamic_interactive_tree.py")
with open(target_file, "r") as f:
    text = f.read()

# 1. HTML Markup: Add tanglegramOptions card right after treeDatasetSection
old_markup = """          <!-- Tree Dataset Selector -->
          <div id="treeDatasetSection">
            <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block mb-1 text-[10.5px]">Tree Dataset (Single Tree)</label>
            <select id="datasetSelect" onchange="switchDataset(this.value)" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded-md px-2.5 py-1.5 text-xs text-[var(--text-main)] focus:outline-none focus:border-sky-400">
              <option value="3di">3Di Structural Tree (FoldMason + Q.3Di.LLM)</option>
              <option value="aa">Amino Acid Sequence Tree (IQ-TREE LG+G4)</option>
            </select>
          </div>"""

new_markup = """          <!-- Tree Dataset Selector -->
          <div id="treeDatasetSection">
            <label class="font-semibold text-[var(--text-muted)] uppercase tracking-wider block mb-1 text-[10.5px]">Tree Dataset (Single Tree)</label>
            <select id="datasetSelect" onchange="switchDataset(this.value)" class="w-full bg-[var(--input-bg)] border border-[var(--border-color)] rounded-md px-2.5 py-1.5 text-xs text-[var(--text-main)] focus:outline-none focus:border-sky-400">
              <option value="3di">3Di Structural Tree (FoldMason + Q.3Di.LLM)</option>
              <option value="aa">Amino Acid Sequence Tree (IQ-TREE LG+G4)</option>
            </select>
          </div>

          <!-- Tanglegram Alignment Selector (Visible in Tanglegram layout) -->
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

assert old_markup in text, "old_markup not found"
text = text.replace(old_markup, new_markup, 1)

# 2. Add tangleMode to settings
old_settings = """      rootingMode: "midpoint",
      outgroupTaxon: null,
      zoom: {{ x: 40, y: 35, k: 0.65 }}
    }};"""

new_settings = """      rootingMode: "midpoint",
      outgroupTaxon: null,
      tangleMode: "true_topology",
      zoom: {{ x: 40, y: 35, k: 0.65 }}
    }};"""

assert old_settings in text, "old_settings not found"
text = text.replace(old_settings, new_settings, 1)

# 3. Update setLayout to toggle single tree dataset vs tanglegram options, and add setTangleMode
old_set_layout = """    function setLayout(layoutName) {{
      settings.layout = layoutName;
      ["btnRect", "btnClado", "btnRadial", "btnUnrooted", "btnTangle"].forEach(id => {{
        const btn = document.getElementById(id);
        if (btn) btn.className = "py-1.5 rounded font-medium text-center hover:bg-slate-500/20 text-[var(--text-muted)] text-[10.5px] transition";
      }});
      const activeBtnMap = {{
        "rectangular": "btnRect",
        "cladogram": "btnClado",
        "radial": "btnRadial",
        "unrooted": "btnUnrooted",
        "tanglegram": "btnTangle"
      }};
      if (activeBtnMap[layoutName]) {{
        document.getElementById(activeBtnMap[layoutName]).className = "py-1.5 rounded font-medium text-center bg-sky-500 text-white text-[10.5px] transition";
      }}
      renderTree();
      fitTreeToScreen(true);
    }}"""

new_set_layout = """    function setLayout(layoutName) {{
      settings.layout = layoutName;
      ["btnRect", "btnClado", "btnRadial", "btnUnrooted", "btnTangle"].forEach(id => {{
        const btn = document.getElementById(id);
        if (btn) btn.className = "py-1.5 rounded font-medium text-center hover:bg-slate-500/20 text-[var(--text-muted)] text-[10.5px] transition";
      }});
      const activeBtnMap = {{
        "rectangular": "btnRect",
        "cladogram": "btnClado",
        "radial": "btnRadial",
        "unrooted": "btnUnrooted",
        "tanglegram": "btnTangle"
      }};
      if (activeBtnMap[layoutName]) {{
        document.getElementById(activeBtnMap[layoutName]).className = "py-1.5 rounded font-medium text-center bg-sky-500 text-white text-[10.5px] transition";
      }}

      // Toggle single-tree dataset selector vs tanglegram alignment options
      const dsSec = document.getElementById("treeDatasetSection");
      const tangleOpt = document.getElementById("tanglegramOptions");
      if (layoutName === "tanglegram") {{
        if (dsSec) dsSec.classList.add("hidden");
        if (tangleOpt) tangleOpt.classList.remove("hidden");
      }} else {{
        if (dsSec) dsSec.classList.remove("hidden");
        if (tangleOpt) tangleOpt.classList.add("hidden");
      }}

      renderTree();
      fitTreeToScreen(true);
    }}

    function setTangleMode(mode) {{
      settings.tangleMode = mode;
      const descMap = {{
        "true_topology": "True topology preserves native IQ-TREE branch order on both sides to expose structural vs sequence discordance.",
        "min_crossings": "Rotates internal clades of the sequence tree to minimize crossing tangles while preserving 100% of phylogenetic clades.",
        "aligned": "Orders right-hand leaves directly alongside matching left-hand leaves for parallel 1-to-1 visual comparison."
      }};
      const descEl = document.getElementById("tangleModeDesc");
      if (descEl) descEl.textContent = descMap[mode] || "";
      renderTree();
      updateMinimap();
    }}"""

assert old_set_layout in text, "old_set_layout not found"
text = text.replace(old_set_layout, new_set_layout, 1)

# 4. Replace renderTanglegram with complete multi-mode implementation
old_render_tangle = text[text.find("    // DUAL NON-OVERLAPPING TANGLEGRAM LAYOUT") : text.find("    // 11. TOOLTIPS & ACTIONS")]

new_render_tangle = """    // DUAL TANGLEGRAM LAYOUT (With True Topology, Min-Crossing Untangle, and Aligned Modes)
    function renderTanglegram(g) {{
      const spacing = settings.verticalSpacing;
      const leftTreeRootX = 50;
      const leftTreeSpan = 260;
      const leftLabelX = leftTreeRootX + leftTreeSpan + 15;
      const leftConnectorX = leftLabelX + 135;
      const connectorSpan = 280;
      const rightConnectorX = leftConnectorX + connectorSpan;
      const rightLabelX = rightConnectorX + 135;
      const rightTreeLeavesX = rightLabelX + 15;
      const rightTreeSpan = 260;
      const rightTreeRootX = rightTreeLeavesX + rightTreeSpan;

      assignDepths(rawRoot3Di, 0);
      assignDepths(rawRootAA, 0);

      // 1. Natural DFS post-order leaf traversal for 3Di
      function getTreeLeavesInOrder(node) {{
        if (!node.children || node.children.length === 0) return [node];
        let res = [];
        node.children.forEach(c => {{
          res = res.concat(getTreeLeavesInOrder(c));
        }});
        return res;
      }}

      const leaves3Di = getTreeLeavesInOrder(rawRoot3Di);
      leaves3Di.forEach((leaf, idx) => {{
        leaf.y = 50 + idx * spacing;
      }});

      const map3DiIndex = {{}};
      leaves3Di.forEach((l, idx) => {{
        map3DiIndex[l.name] = idx;
      }});

      // 2. Order AA Tree Leaves based on settings.tangleMode
      let leavesAA = [];

      if (settings.tangleMode === "aligned") {{
        // Aligned / Parallel: order right-hand leaves to match left-hand leaves
        leavesAA = [...getAllLeaves(rawRootAA)].sort((a, b) => {{
          const idxA = map3DiIndex[a.name] !== undefined ? map3DiIndex[a.name] : 9999;
          const idxB = map3DiIndex[b.name] !== undefined ? map3DiIndex[b.name] : 9999;
          return idxA - idxB;
        }});
        leavesAA.forEach((leaf, idx) => {{
          leaf.y = 50 + idx * spacing;
        }});
      }} else if (settings.tangleMode === "min_crossings") {{
        // Optimal Subtree Rotations (Barycenter Heuristic): rotate children around internal nodes
        // to minimize crossings while strictly preserving phylogenetic clades
        function getAvg3Di(node) {{
          const lf = getTreeLeavesInOrder(node);
          const idxs = lf.map(l => map3DiIndex[l.name] !== undefined ? map3DiIndex[l.name] : 0);
          return idxs.reduce((a, b) => a + b, 0) / (idxs.length || 1);
        }}
        function rotateSubtrees(node) {{
          if (!node.children || node.children.length <= 1) return;
          node.children.forEach(rotateSubtrees);
          node.children.sort((a, b) => getAvg3Di(a) - getAvg3Di(b));
        }}
        rotateSubtrees(rawRootAA);
        leavesAA = getTreeLeavesInOrder(rawRootAA);
        leavesAA.forEach((leaf, idx) => {{
          leaf.y = 50 + idx * spacing;
        }});
      }} else {{
        // "true_topology": Independent DFS post-order traversal of the AA tree (Native IQ-TREE order)
        // Reveals true topological discordance and crossing tangles!
        leavesAA = getTreeLeavesInOrder(rawRootAA);
        leavesAA.forEach((leaf, idx) => {{
          leaf.y = 50 + idx * spacing;
        }});
      }}

      // Compute internal Y coordinates for both trees
      function computeInternal(node) {{
        if (!node.children || node.children.length === 0) return;
        node.children.forEach(computeInternal);
        const childY = node.children.map(c => c.y).filter(y => typeof y === 'number' && !isNaN(y));
        if (childY.length > 0) {{
          node.y = childY.reduce((a, b) => a + b, 0) / childY.length;
        }}
      }}
      computeInternal(rawRoot3Di);
      computeInternal(rawRootAA);

      const maxDepth3Di = Math.max(...leaves3Di.map(l => settings.branchLengths ? l.depth : l.cladoDepth)) || 1.0;
      const maxDepthAA = Math.max(...leavesAA.map(l => settings.branchLengths ? l.depth : l.cladoDepth)) || 1.0;

      // Draw Left Tree (3Di Structural, branching rightwards towards center)
      function drawLeft(node, currentX) {{
        node.x = currentX;
        if (!node.children || node.children.length === 0) return;

        const childY = node.children.map(c => c.y).filter(y => typeof y === 'number');
        const minY = Math.min(...childY);
        const maxY = Math.max(...childY);

        const vLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
        vLine.setAttribute("x1", node.x);
        vLine.setAttribute("y1", minY);
        vLine.setAttribute("x2", node.x);
        vLine.setAttribute("y2", maxY);
        vLine.setAttribute("class", "branch-path");
        g.appendChild(vLine);

        node.children.forEach(child => {{
          const bLen = settings.branchLengths ? (typeof child.length === 'number' ? child.length : 1.0) : 1.0;
          const childX = node.x + (bLen / maxDepth3Di) * leftTreeSpan;
          child.x = childX;

          const hLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          hLine.setAttribute("x1", node.x);
          hLine.setAttribute("y1", child.y);
          hLine.setAttribute("x2", childX);
          hLine.setAttribute("y2", child.y);
          hLine.setAttribute("class", "branch-path");
          g.appendChild(hLine);

          drawLeft(child, childX);
        }});
        drawNodePoint(g, node);
      }}

      // Draw Right Tree (AA Sequence, branching leftwards towards center)
      function drawRight(node, currentX) {{
        node.x = currentX;
        if (!node.children || node.children.length === 0) return;

        const childY = node.children.map(c => c.y).filter(y => typeof y === 'number');
        const minY = Math.min(...childY);
        const maxY = Math.max(...childY);

        const vLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
        vLine.setAttribute("x1", node.x);
        vLine.setAttribute("y1", minY);
        vLine.setAttribute("x2", node.x);
        vLine.setAttribute("y2", maxY);
        vLine.setAttribute("class", "branch-path");
        g.appendChild(vLine);

        node.children.forEach(child => {{
          const bLen = settings.branchLengths ? (typeof child.length === 'number' ? child.length : 1.0) : 1.0;
          const childX = node.x - (bLen / maxDepthAA) * rightTreeSpan;
          child.x = childX;

          const hLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          hLine.setAttribute("x1", node.x);
          hLine.setAttribute("y1", child.y);
          hLine.setAttribute("x2", childX);
          hLine.setAttribute("y2", child.y);
          hLine.setAttribute("class", "branch-path");
          g.appendChild(hLine);

          drawRight(child, childX);
        }});
        drawNodePoint(g, node);
      }}

      drawLeft(rawRoot3Di, leftTreeRootX);
      drawRight(rawRootAA, rightTreeRootX);

      // Section Header Badges
      const leftHeader = document.createElementNS("http://www.w3.org/2000/svg", "text");
      leftHeader.setAttribute("x", leftTreeRootX);
      leftHeader.setAttribute("y", 25);
      leftHeader.setAttribute("fill", "var(--accent)");
      leftHeader.setAttribute("font-weight", "bold");
      leftHeader.setAttribute("font-size", "12px");
      leftHeader.textContent = "3Di Structural Tree (FoldMason + Q.3Di.LLM)";
      g.appendChild(leftHeader);

      const rightHeader = document.createElementNS("http://www.w3.org/2000/svg", "text");
      rightHeader.setAttribute("x", rightTreeRootX);
      rightHeader.setAttribute("y", 25);
      rightHeader.setAttribute("text-anchor", "end");
      rightHeader.setAttribute("fill", "#a855f7");
      rightHeader.setAttribute("font-weight", "bold");
      rightHeader.setAttribute("font-size", "12px");
      rightHeader.textContent = "Amino Acid Tree (IQ-TREE LG+G4)";
      g.appendChild(rightHeader);

      // Count Inversions / Crossings
      const mapAA = {{}};
      leavesAA.forEach(l => {{ mapAA[l.name] = l; }});

      const arrCross = [];
      leaves3Di.forEach(l1 => {{
        const l2 = mapAA[l1.name];
        if (l2 && typeof l2.y === 'number') arrCross.push(l2.y);
      }});
      let crossings = 0;
      for (let i = 0; i < arrCross.length; i++) {{
        for (let j = i + 1; j < arrCross.length; j++) {{
          if (arrCross[i] > arrCross[j]) crossings++;
        }}
      }}
      const badgeCross = document.getElementById("tangleCrossingBadge");
      if (badgeCross) {{
        badgeCross.textContent = `Crossings: ${{crossings.toLocaleString()}}`;
      }}

      // Draw Tanglegram Connecting Curves and Labels
      leaves3Di.forEach(l1 => {{
        drawNodePoint(g, l1);
        const l2 = mapAA[l1.name];
        if (!l2) return;
        drawNodePoint(g, l2);

        const cleanId = l1.name.replace(/[^a-zA-Z0-9]/g, "_");

        // Hairline extension from 3Di leaf to left label
        if (leftLabelX - 5 > l1.x) {{
          const extL = document.createElementNS("http://www.w3.org/2000/svg", "line");
          extL.setAttribute("x1", l1.x);
          extL.setAttribute("y1", l1.y);
          extL.setAttribute("x2", leftLabelX - 5);
          extL.setAttribute("y2", l1.y);
          extL.setAttribute("stroke", "var(--grid-line)");
          extL.setAttribute("stroke-width", "0.6px");
          extL.setAttribute("stroke-dasharray", "2,2");
          g.appendChild(extL);
        }}

        // Left Label
        const txt1 = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt1.setAttribute("id", "tangle_label_left_" + cleanId);
        txt1.setAttribute("x", leftLabelX);
        txt1.setAttribute("y", l1.y + 3.5);
        txt1.setAttribute("class", `tip-label ${{l1.name === settings.selectedTaxon ? "selected" : ""}}`);
        txt1.textContent = l1.name;
        txt1.onclick = () => selectTaxon(l1.name);
        txt1.onmouseenter = () => highlightTangleTaxon(l1.name);
        txt1.onmouseleave = () => clearTangleHighlight();
        g.appendChild(txt1);

        // Hairline extension from AA leaf to right label
        if (l2.x > rightLabelX + 5) {{
          const extR = document.createElementNS("http://www.w3.org/2000/svg", "line");
          extR.setAttribute("x1", rightLabelX + 5);
          extR.setAttribute("y1", l2.y);
          extR.setAttribute("x2", l2.x);
          extR.setAttribute("y2", l2.y);
          extR.setAttribute("stroke", "var(--grid-line)");
          extR.setAttribute("stroke-width", "0.6px");
          extR.setAttribute("stroke-dasharray", "2,2");
          g.appendChild(extR);
        }}

        // Right Label
        const txt2 = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt2.setAttribute("id", "tangle_label_right_" + cleanId);
        txt2.setAttribute("x", rightLabelX);
        txt2.setAttribute("y", l2.y + 3.5);
        txt2.setAttribute("text-anchor", "end");
        txt2.setAttribute("class", `tip-label ${{l1.name === settings.selectedTaxon ? "selected" : ""}}`);
        txt2.textContent = l2.name;
        txt2.onclick = () => selectTaxon(l1.name);
        txt2.onmouseenter = () => highlightTangleTaxon(l1.name);
        txt2.onmouseleave = () => clearTangleHighlight();
        g.appendChild(txt2);

        // Cubic Bézier Curve across Middle Channel
        const x1 = leftConnectorX;
        const y1 = l1.y;
        const x2 = rightConnectorX;
        const y2 = l2.y;
        const cp1x = x1 + (x2 - x1) * 0.5;
        const cp2x = x2 - (x2 - x1) * 0.5;

        const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
        path.setAttribute("id", "tangle_" + cleanId);
        path.setAttribute("d", `M ${{x1}} ${{y1}} C ${{cp1x}} ${{y1}}, ${{cp2x}} ${{y2}}, ${{x2}} ${{y2}}`);
        path.setAttribute("class", "tangle-connector");
        path.style.stroke = getNodeColor(l1);

        path.onmouseenter = () => highlightTangleTaxon(l1.name);
        path.onmouseleave = () => clearTangleHighlight();
        path.onclick = () => selectTaxon(l1.name);
        g.appendChild(path);
      }});
    }}

"""

assert len(old_render_tangle) > 500, "old_render_tangle slice invalid"
text = text.replace(old_render_tangle, new_render_tangle, 1)

with open(target_file, "w") as f:
    f.write(text)

print("Successfully patched build_dynamic_interactive_tree.py with multi-mode Tanglegram!")
