    const DATASETS = /*__DATASETS_JSON__*/;

    let currentScale = "1193";
    let activeDataset = DATASETS[currentScale];
    let NEWICK_3DI = activeDataset.newick_3di;
    let NEWICK_AA = activeDataset.newick_aa;
    let NEWICK_ESM2 = activeDataset.newick_esm2 || activeDataset.newick_esm2_cosine;
    let NEWICK_ESM2_COSINE = activeDataset.newick_esm2_cosine;
    let NEWICK_ESM2_EUCLIDEAN = activeDataset.newick_esm2_euclidean;
    let NEWICK_ESM2_L1 = activeDataset.newick_esm2_l1;
    function syncDatasetAlignmentCoverage(scale) {
      if (!window.ALIGNMENTS_DATA || !window.ALIGNMENTS_DATA[scale]) return;
      const alnData = window.ALIGNMENTS_DATA[scale];
      const covMap = alnData.coverage || {};
      const ds = DATASETS[scale];
      if (!ds) return;

      if (ds.taxa) {
        Object.keys(ds.taxa).forEach(tid => {
          if (covMap[tid] !== undefined) {
            ds.taxa[tid].alignment_coverage = Math.round(covMap[tid] * 100);
          }
        });
      }

      if (ds.columns && !ds.columns.some(c => c.key === "alignment_coverage")) {
        ds.columns.push({
          key: "alignment_coverage",
          label: "Alignment Coverage (%)",
          type: "continuous",
          min: 0,
          max: 100
        });
      }
    }

    Object.keys(DATASETS).forEach(sc => {
      syncDatasetAlignmentCoverage(sc);
    });
    let TAXA_METADATA = activeDataset.taxa;

    let rawTrees = {
      "3di": null,
      "aa": null,
      "esm2_cosine": null,
      "esm2_euclidean": null,
      "esm2_l1": null
    };
    let rawRoot3Di = null;
    let rawRootAA = null;
    let rawRootESM2 = null;

    let umapPointCoords = {};

    const tangleState = {
      activeLeftRoot: null,
      activeRightRoot: null,
      leftLeaves: [],
      rightLeaves: [],
      leftTreeRootX: 50,
      rightTreeRootX: 1150,
      leftConnectorX: 460,
      rightConnectorX: 740,
      leftColor: "var(--accent)",
      rightColor: "#a855f7"
    };

    const settings = {
      dataset: "3di",
      embedMetric: "cosine",
      tangleCompare: "3di_vs_aa",
      layout: "rectangular",
      alignLabels: false,
      branchLengths: true,
      showMeta: false,
      colorColumn: activeDataset.defaultColorCol || "structural_class",
      cladeGroupColumn: activeDataset.defaultCladeCol || "structural_class",
      tipLabelColumn: "taxon_id",
      verticalSpacing: activeDataset.defaultSpacing || 14,
      nodeRadius: activeDataset.defaultRadius || 2.8,
      labelSize: 10,
      zoomAdaptiveLabels: true,
      curvedConnectors: true,
      plddtGlow: true,
      selectedTaxon: null,
      theme: localStorage.getItem("phylo_theme") || "dark",
      cladeHomogeneity: 75,
      minCladeSize: 3,
      rootingMode: "midpoint",
      outgroupTaxon: null,
      tangleMode: "true_topology",
      scopedClade: null,
      zoom: { x: 40, y: 35, k: 0.65 },
      treeRotation: 0,
      radialArc: 360,
      radialRadiusScale: 1.0,
      unrootedScale: 1.0,
      labelOrientation: "radial",
      concentricRings: true,
      cladeSectors: true,
      branchWidth: 1.4,
      colorBranchesBySupport: false,
      supportMetric: "ufboot",
      supportPalette: "traffic",
      showGrid: true,
      unrootedLengthMode: "sqrt",
      unrootedDaylight: true,
      unrootedLabelFilter: "smart",
      staggerLabels: true
    };

    // 1. ACTIVE DATASET & METADATA ACCESSORS
    function getActiveDataset() {
      return DATASETS[currentScale];
    }

    function getActiveColumns() {
      return (getActiveDataset() && getActiveDataset().columns) || [];
    }

    function getActiveColorColumnDef() {
      const cols = getActiveColumns();
      return cols.find(c => c.key === settings.colorColumn) || cols[0];
    }

    function getActiveCladeColumnDef() {
      const cols = getActiveColumns().filter(c => c.type === "categorical");
      return cols.find(c => c.key === settings.cladeGroupColumn) || cols[0];
    }

    function getLeafLabelText(leafName) {
      if (!leafName) return "";
      if (!settings.tipLabelColumn || settings.tipLabelColumn === "taxon_id" || settings.tipLabelColumn === "taxon_name") {
        return leafName;
      }

      const meta = (typeof TAXA_METADATA !== 'undefined' && TAXA_METADATA) ? TAXA_METADATA[leafName] : null;
      if (!meta) return leafName;

      if (settings.tipLabelColumn === "id_and_color") {
        const colDef = getActiveColorColumnDef();
        const colorVal = (colDef && meta[colDef.key] !== undefined) ? meta[colDef.key] : "";
        return colorVal ? `${leafName} (${colorVal})` : leafName;
      }

      const val = meta[settings.tipLabelColumn];
      if (val !== undefined && val !== null && String(val).trim() !== "") {
        return String(val);
      }
      return leafName;
    }

    // 2. PARSE NEWICK
    let nodeIdCounter = 0;
    function parseNewick(str) {
      let s = str.trim().replace(/;$/, "");
      let cursor = 0;

      function parseNode() {
        let node = { 
          id: "node_" + (++nodeIdCounter),
          children: [], 
          name: "", 
          length: 0.0, 
          support: null,
          ufboot: null,
          alrt: null,
          supportRaw: "",
          _collapsed: false
        };

        if (s[cursor] === '(') {
          cursor++;
          while (cursor < s.length) {
            node.children.push(parseNode());
            if (s[cursor] === ',') {
              cursor++;
            } else if (s[cursor] === ')') {
              cursor++;
              break;
            }
          }
        }

        let token = "";
        while (cursor < s.length && s[cursor] !== ':' && s[cursor] !== ',' && s[cursor] !== ')' && s[cursor] !== ';') {
          token += s[cursor];
          cursor++;
        }
        token = token.trim();
        if (token) {
          if (node.children.length > 0) {
            node.supportRaw = token;
            if (token.includes('/')) {
              const parts = token.split('/');
              const alrt = parseFloat(parts[0]);
              const ufboot = parseFloat(parts[1]);
              if (!isNaN(alrt)) node.alrt = alrt;
              if (!isNaN(ufboot)) node.ufboot = ufboot;
              node.support = !isNaN(ufboot) ? ufboot : (!isNaN(alrt) ? alrt : null);
            } else {
              let sup = parseFloat(token);
              if (!isNaN(sup)) {
                if (sup <= 1.0 && sup > 0.0) sup = sup * 100.0;
                node.support = sup;
                node.ufboot = sup;
                node.alrt = sup;
              } else {
                node.name = token;
              }
            }
          } else {
            node.name = token;
          }
        }

        if (cursor < s.length && s[cursor] === ':') {
          cursor++;
          let lenStr = "";
          while (cursor < s.length && s[cursor] !== ',' && s[cursor] !== ')' && s[cursor] !== ';') {
            lenStr += s[cursor];
            cursor++;
          }
          let l = parseFloat(lenStr.trim());
          if (!isNaN(l)) node.length = Math.max(0.0001, l);
        }

        return node;
      }

      return parseNode();
    }

    function parseAllActiveTrees() {
      rawTrees["3di"] = parseNewick(NEWICK_3DI);
      rawTrees["aa"] = parseNewick(NEWICK_AA);
      delete rawTrees["esm2_cosine"];
      delete rawTrees["esm2_euclidean"];
      delete rawTrees["esm2_l1"];
      delete rawTrees["esm2"];

      if (NEWICK_ESM2_COSINE) rawTrees["esm2_cosine"] = parseNewick(NEWICK_ESM2_COSINE);
      if (NEWICK_ESM2_EUCLIDEAN) rawTrees["esm2_euclidean"] = parseNewick(NEWICK_ESM2_EUCLIDEAN);
      if (NEWICK_ESM2_L1) rawTrees["esm2_l1"] = parseNewick(NEWICK_ESM2_L1);

      rawRoot3Di = rawTrees["3di"];
      rawRootAA = rawTrees["aa"];
      const metricKey = "esm2_" + (settings.embedMetric || "cosine");
      rawRootESM2 = rawTrees[metricKey] || rawTrees["esm2_cosine"] || null;
    }

    function updateModalityOptions() {
      const ds = activeDataset || (typeof DATASETS !== "undefined" && DATASETS[currentScale]) || {};
      const hasEsm = Boolean(ds.has_esm && (ds.newick_esm2_cosine || ds.newick_esm2));

      // 1. Single Tree Dataset Selector (#datasetSelect)
      const optEsm = document.querySelector("#datasetSelect option[value='esm2']");
      if (optEsm) {
        optEsm.disabled = !hasEsm;
        if (!hasEsm) {
          optEsm.textContent = "🤖 ESM-2 PLM Tree (Not Run in Pipeline)";
          optEsm.className = "text-slate-500 bg-slate-900/80 italic cursor-not-allowed";
        } else {
          optEsm.textContent = "🤖 ESM-2 PLM Tree (Hierarchical Clustering)";
          optEsm.className = "text-[var(--text-main)]";
        }
      }
      if (!hasEsm && settings.dataset === "esm2") {
        settings.dataset = "3di";
        const dsSel = document.getElementById("datasetSelect");
        if (dsSel) dsSel.value = "3di";
        switchDataset("3di");
      }

      // 2. Embedding Metric Sub-Section (#embedMetricSubSection)
      const embedSub = document.getElementById("embedMetricSubSection");
      if (embedSub) {
        if (!hasEsm || settings.dataset !== "esm2") {
          embedSub.classList.add("hidden");
        } else {
          embedSub.classList.remove("hidden");
        }
      }

      // 3. Tanglegram Comparison Pairs (#tangleCompareSelect)
      const tangleSelect = document.getElementById("tangleCompareSelect");
      if (tangleSelect) {
        tangleSelect.querySelectorAll("option").forEach(opt => {
          if (opt.value.includes("esm2")) {
            opt.disabled = !hasEsm;
            if (!hasEsm) {
              if (!opt.dataset.origText) opt.dataset.origText = opt.textContent;
              if (!opt.textContent.includes("(Not Run)")) {
                opt.textContent = opt.dataset.origText + " (Not Run)";
              }
              opt.className = "text-slate-500 bg-slate-900/80 italic cursor-not-allowed";
            } else if (opt.dataset.origText) {
              opt.textContent = opt.dataset.origText;
              opt.className = "text-[var(--text-main)]";
            }
          }
        });
        tangleSelect.querySelectorAll("optgroup").forEach(og => {
          if (og.label && (og.label.includes("PLM") || og.label.includes("Embeddings"))) {
            og.disabled = !hasEsm;
            og.className = !hasEsm ? "text-slate-600 italic" : "";
          }
        });
        if (!hasEsm && settings.tangleCompare && settings.tangleCompare.includes("esm2")) {
          settings.tangleCompare = "3di_vs_aa";
          tangleSelect.value = "3di_vs_aa";
        }
      }

      // 4. Clade Partition Source (#cladePartitionSourceSelect)
      const cladeOptEsm = document.querySelector("#cladePartitionSourceSelect option[value='esm2']");
      if (cladeOptEsm) {
        cladeOptEsm.disabled = !hasEsm;
        if (!hasEsm) {
          if (!cladeOptEsm.dataset.origText) cladeOptEsm.dataset.origText = cladeOptEsm.textContent;
          if (!cladeOptEsm.textContent.includes("(Not Run)")) {
            cladeOptEsm.textContent = cladeOptEsm.dataset.origText + " (Not Run)";
          }
          cladeOptEsm.className = "text-slate-500 bg-slate-900/80 italic cursor-not-allowed";
        } else if (cladeOptEsm.dataset.origText) {
          cladeOptEsm.textContent = cladeOptEsm.dataset.origText;
          cladeOptEsm.className = "text-[var(--text-main)]";
        }
      }
      if (!hasEsm && cladePartitionState.source === "esm2") {
        cladePartitionState.source = "3di";
        const cladeSrcSel = document.getElementById("cladePartitionSourceSelect");
        if (cladeSrcSel) cladeSrcSel.value = "3di";
        setCladePartitionSource("3di");
      }

      // 5. UMAP Layout Button (#btnUmap) & ESM View Mode Section (#esmViewModeSection)
      const btnUmap = document.getElementById("btnUmap");
      const hasUmap = Boolean(hasEsm && ds.esm2_umap && Object.keys(ds.esm2_umap).length > 0);
      if (btnUmap) {
        if (!hasUmap) {
          btnUmap.disabled = true;
          btnUmap.title = "ESM-2 UMAP (Not Available for this Dataset)";
          btnUmap.className = "py-1.5 rounded font-medium text-center opacity-40 cursor-not-allowed text-slate-500 text-[10.5px] transition";
        } else {
          btnUmap.disabled = false;
          btnUmap.title = "ESM-2 PLM 2D UMAP Scatter Projection";
          btnUmap.className = (settings.layout === "umap")
            ? "py-1.5 rounded font-medium text-center bg-purple-600 text-white text-[10.5px] transition shadow-sm cursor-pointer"
            : "py-1.5 rounded font-medium text-center hover:bg-purple-500/20 text-purple-400 border border-purple-500/30 text-[10.5px] transition flex items-center justify-center space-x-0.5 cursor-pointer";
        }
      }

      if (!hasUmap && settings.layout === "umap") {
        setLayout("rectangular");
      }

      const esmViewSec = document.getElementById("esmViewModeSection");
      if (esmViewSec) {
        if (hasUmap && settings.dataset === "esm2") {
          esmViewSec.classList.remove("hidden");
          const btnEsmTree = document.getElementById("btnEsmViewTree");
          const btnEsmUmap = document.getElementById("btnEsmViewUmap");
          if (btnEsmTree && btnEsmUmap) {
            if (settings.layout === "umap") {
              btnEsmUmap.className = "py-1 px-2 rounded font-semibold text-center bg-purple-600 text-white text-[10.5px] transition shadow-sm cursor-pointer flex items-center justify-center space-x-1";
              btnEsmTree.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-purple-500/20 text-[10.5px] transition cursor-pointer flex items-center justify-center space-x-1";
            } else {
              btnEsmTree.className = "py-1 px-2 rounded font-semibold text-center bg-purple-600 text-white text-[10.5px] transition shadow-sm cursor-pointer flex items-center justify-center space-x-1";
              btnEsmUmap.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-purple-500/20 text-[10.5px] transition cursor-pointer flex items-center justify-center space-x-1";
            }
          }
        } else {
          esmViewSec.classList.add("hidden");
        }
      }
    }

    parseAllActiveTrees();
    let activeTreeRoot = rawRoot3Di;

    // 3. TREE DIAMETER & REROOTING ENGINE
    function buildGraphFromTree(rootNode) {
      const adj = {};
      const nodesMap = {};

      function traverse(n, p = null, w = 0) {
        nodesMap[n.id] = n;
        if (!adj[n.id]) adj[n.id] = [];
        if (p !== null) {
          adj[n.id].push({ neighbor: p.id, w: w });
          adj[p.id].push({ neighbor: n.id, w: w });
        }
        if (n.children) {
          n.children.forEach(c => traverse(c, n, c.length || 0.001));
        }
      }
      traverse(rootNode);
      return { adj, nodesMap };
    }

    function computeTreeDiameter(graph) {
      const leafIds = Object.keys(graph.nodesMap).filter(id => !graph.nodesMap[id].children || graph.nodesMap[id].children.length === 0);
      if (leafIds.length < 2) return null;
      const leafSet = new Set(leafIds);

      function getFarthest(startId) {
        const dist = {};
        const parent = {};
        dist[startId] = 0;
        parent[startId] = { p: null, w: 0 };
        const queue = [startId];
        let farthestId = startId;
        let maxDist = 0;

        while (queue.length > 0) {
          const curr = queue.shift();
          const d = dist[curr];
          if (d > maxDist && leafSet.has(curr)) {
            maxDist = d;
            farthestId = curr;
          }

          const neighbors = graph.adj[curr] || [];
          for (let i = 0; i < neighbors.length; i++) {
            const edge = neighbors[i];
            const nb = edge.neighbor;
            if (dist[nb] === undefined) {
              dist[nb] = d + edge.w;
              parent[nb] = { p: curr, w: edge.w };
              queue.push(nb);
            }
          }
        }
        return { farthestId, maxDist, parent };
      }

      const r1 = getFarthest(leafIds[0]);
      const r2 = getFarthest(r1.farthestId);

      const path = [];
      let curr = r2.farthestId;
      while (curr !== null) {
        const info = r2.parent[curr];
        path.push({ id: curr, w: info.w });
        curr = info.p;
      }
      path.reverse();

      return {
        l1: r1.farthestId,
        l2: r2.farthestId,
        diameter: r2.maxDist,
        path: path
      };
    }

    function rerootAtEdge(graph, uId, vId, distU, distV) {
      const newRoot = {
        id: "reroot_" + (++nodeIdCounter),
        children: [],
        name: "",
        length: 0.0,
        support: null,
        _collapsed: false
      };

      function buildSubtree(currId, parentId) {
        const orig = graph.nodesMap[currId];
        const copy = {
          id: orig.id,
          name: orig.name,
          length: 0.0,
          support: orig.support,
          ufboot: orig.ufboot,
          alrt: orig.alrt,
          supportRaw: orig.supportRaw || "",
          children: [],
          _collapsed: orig._collapsed || false
        };

        const neighbors = graph.adj[currId] || [];
        neighbors.forEach(edge => {
          const nb = edge.neighbor;
          if (nb !== parentId) {
            const childNode = buildSubtree(nb, currId);
            childNode.length = edge.w;
            copy.children.push(childNode);
          }
        });
        return copy;
      }

      const childU = buildSubtree(uId, vId);
      childU.length = distU;
      const childV = buildSubtree(vId, uId);
      childV.length = distV;

      newRoot.children = [childU, childV];
      return newRoot;
    }

    function performMidpointRoot(rootNode) {
      const graph = buildGraphFromTree(rootNode);
      const diamInfo = computeTreeDiameter(graph);
      if (!diamInfo || !diamInfo.path || diamInfo.path.length < 2) return rootNode;

      const targetMid = diamInfo.diameter / 2.0;
      let cum = 0;
      let midU = null, midV = null, dU = 0, dV = 0;

      for (let i = 0; i < diamInfo.path.length - 1; i++) {
        const u = diamInfo.path[i].id;
        const v = diamInfo.path[i + 1].id;
        const w = diamInfo.path[i + 1].w;
        if (cum + w >= targetMid) {
          midU = u;
          midV = v;
          dU = targetMid - cum;
          dV = w - dU;
          break;
        }
        cum += w;
      }

      if (!midU) return rootNode;
      const rooted = rerootAtEdge(graph, midU, midV, dU, dV);
      rooted._diamInfo = diamInfo;
      return rooted;
    }

    function performOutgroupRoot(rootNode, taxonName) {
      const graph = buildGraphFromTree(rootNode);
      const targetId = Object.keys(graph.nodesMap).find(id => graph.nodesMap[id].name === taxonName);
      if (!targetId) return rootNode;

      const edges = graph.adj[targetId];
      if (!edges || edges.length === 0) return rootNode;
      const parentEdge = edges[0];
      const halfW = parentEdge.w / 2.0;
      return rerootAtEdge(graph, targetId, parentEdge.neighbor, halfW, halfW);
    }

    function performNodeReroot(rootNode, nodeId) {
      const graph = buildGraphFromTree(rootNode);
      const targetNode = graph.nodesMap[nodeId];
      if (!targetNode) return rootNode;

      const edges = graph.adj[nodeId];
      if (!edges || edges.length === 0) return rootNode;
      const parentEdge = edges[0];
      const halfW = parentEdge.w / 2.0;
      return rerootAtEdge(graph, nodeId, parentEdge.neighbor, halfW, halfW);
    }

    // Apply active rooting mode
    function applyCurrentRooting() {
      let baseRoot = rawRoot3Di;
      if (settings.dataset === "aa") {
        baseRoot = rawRootAA;
      } else if (settings.dataset === "esm2") {
        const metricKey = "esm2_" + (settings.embedMetric || "cosine");
        baseRoot = rawTrees[metricKey] || rawTrees["esm2_cosine"] || rawRoot3Di;
      }

      if (filterState.isActive && filterState.mode === "prune") {
        const allowedSet = getFilteredTaxaSet();
        if (allowedSet.size > 0) {
          const pruned = pruneSubtree(baseRoot, allowedSet);
          if (pruned) baseRoot = pruned;
        }
      }

      if (settings.scopedClade && settings.scopedClade.taxa && settings.scopedClade.taxa.size > 0) {
        const pruned = pruneSubtree(baseRoot, settings.scopedClade.taxa);
        if (pruned) baseRoot = pruned;
      }

      if (settings.rootingMode === "midpoint") {
        activeTreeRoot = performMidpointRoot(baseRoot);
        document.getElementById("badgeRoot").textContent = "Midpoint Root";
        document.getElementById("metricRootPos").textContent = "Midpoint (Balanced 50/50)";
      } else if (settings.rootingMode === "outgroup" && settings.outgroupTaxon) {
        activeTreeRoot = performOutgroupRoot(baseRoot, settings.outgroupTaxon);
        document.getElementById("badgeRoot").textContent = `Outgroup: ${settings.outgroupTaxon}`;
        document.getElementById("metricRootPos").textContent = `Outgroup (${settings.outgroupTaxon})`;
      } else {
        activeTreeRoot = baseRoot;
        document.getElementById("badgeRoot").textContent = "Original Root";
        document.getElementById("metricRootPos").textContent = "Original IQ-TREE Root";
      }

      // Bypass single-child roots on trimmed / pruned trees so root branches cleanly
      while (activeTreeRoot && activeTreeRoot.children && activeTreeRoot.children.length === 1) {
        activeTreeRoot = activeTreeRoot.children[0];
      }

      // Update diameter metrics
      const diam = activeTreeRoot._diamInfo || computeTreeDiameter(buildGraphFromTree(baseRoot));
      if (diam) {
        document.getElementById("metricDiameter").textContent = diam.diameter.toFixed(4) + " subs/site";
        const g = buildGraphFromTree(baseRoot);
        const n1 = g.nodesMap[diam.l1] ? (g.nodesMap[diam.l1].name || diam.l1) : diam.l1;
        const n2 = g.nodesMap[diam.l2] ? (g.nodesMap[diam.l2].name || diam.l2) : diam.l2;
        document.getElementById("metricFarthest").textContent = `${n1} ↔ ${n2}`;
      }

      renderTree();
      fitTreeToScreen(true);
      updateCladeManagementUI();
      if (typeof renderMsa === 'function') {
        renderMsa();
      }
    }

    function setRootingMode(mode) {
      settings.rootingMode = mode;
      ["btnRootMidpoint", "btnRootOriginal", "btnRootOutgroup"].forEach(id => {
        const btn = document.getElementById(id);
        if (btn) {
          btn.className = "py-1.5 rounded font-medium text-center text-[10px] transition text-[var(--text-muted)] hover:bg-slate-500/20";
        }
      });
      if (mode === "midpoint") {
        document.getElementById("btnRootMidpoint").className = "py-1.5 rounded font-medium text-center text-[10px] transition bg-sky-500 text-white";
      } else if (mode === "original") {
        document.getElementById("btnRootOriginal").className = "py-1.5 rounded font-medium text-center text-[10px] transition bg-sky-500 text-white";
      } else if (mode === "outgroup") {
        document.getElementById("btnRootOutgroup").className = "py-1.5 rounded font-medium text-center text-[10px] transition bg-sky-500 text-white";
      }
      applyCurrentRooting();
    }

    function setOutgroupTaxon(taxName) {
      settings.outgroupTaxon = taxName;
      document.getElementById("selectedOutgroupBadge").textContent = taxName;
      if (settings.rootingMode === "outgroup") {
        applyCurrentRooting();
      }
    }

    let pendingHoverNode = null;
    function triggerHoverReroot() {
      if (!pendingHoverNode) return;
      activeTreeRoot = performNodeReroot(activeTreeRoot, pendingHoverNode.id);
      settings.rootingMode = "custom";
      document.getElementById("badgeRoot").textContent = `Rerooted at ${pendingHoverNode.name || "Internal Node"}`;
      document.getElementById("metricRootPos").textContent = `Manual Edge (${pendingHoverNode.name || pendingHoverNode.id})`;
      ["btnRootMidpoint", "btnRootOriginal", "btnRootOutgroup"].forEach(id => {
        const btn = document.getElementById(id);
        if (btn) btn.className = "py-1.5 rounded font-medium text-center text-[10px] transition text-[var(--text-muted)] hover:bg-slate-500/20";
      });
      hideTooltipImmediate();
      renderTree();
      fitTreeToScreen(true);
      showCladeToast(`Tree rerooted at ${pendingHoverNode.name || "selected edge"}!`);
    }

    // 4. CLADE COLLAPSE & CLASSIFICATION ENGINE
    function getAllLeaves(node) {
      if (!node) return [];
      if (!node.children || node.children.length === 0) {
        return [node];
      }
      let leaves = [];
      node.children.forEach(c => {
        leaves = leaves.concat(getAllLeaves(c));
      });
      return leaves;
    }

    function getVisibleLeaves(node) {
      if (!node) return [];
      if (node._collapsed) {
        return [node];
      }
      if (!node.children || node.children.length === 0) {
        return [node];
      }
      let leaves = [];
      node.children.forEach(c => {
        leaves = leaves.concat(getVisibleLeaves(c));
      });
      return leaves;
    }

    function getCladeInfo(node) {
      const leaves = getAllLeaves(node);
      const cladeColDef = getActiveCladeColumnDef();
      const colKey = cladeColDef ? cladeColDef.key : "structural_class";

      const counts = {};
      let maxCount = 0;
      let dominantCategory = "Mixed";

      leaves.forEach(l => {
        const m = TAXA_METADATA[l.name];
        if (m && m[colKey] !== undefined && m[colKey] !== null) {
          const cat = String(m[colKey]);
          counts[cat] = (counts[cat] || 0) + 1;
          if (counts[cat] > maxCount) {
            maxCount = counts[cat];
            dominantCategory = cat;
          }
        }
      });

      const homogeneity = leaves.length > 0 ? (maxCount / leaves.length) * 100 : 0;
      let cladeColor = "#38bdf8";
      if (cladeColDef && cladeColDef.colors && cladeColDef.colors[dominantCategory]) {
        cladeColor = cladeColDef.colors[dominantCategory];
      }

      return {
        count: leaves.length,
        dominantCategory,
        homogeneity,
        cladeColor,
        colKey,
        colLabel: cladeColDef ? cladeColDef.label : "Category"
      };
    }

    function toggleCladeCollapse(node) {
      if (!node.children || node.children.length === 0) return;
      node._collapsed = !node._collapsed;
      renderTree();
      updateCladeManagementUI();
      const info = getCladeInfo(node);
      showCladeToast(`${node._collapsed ? "Collapsed" : "Expanded"} clade (${info.count} taxa &bull; ${info.dominantCategory})`);
    }

    function toggleCategoryClade(catVal, forceCollapse) {
      const cladeColDef = getActiveCladeColumnDef();
      const colKey = cladeColDef ? cladeColDef.key : "structural_class";
      const threshold = settings.cladeHomogeneity || 75;

      let affected = 0;
      function scan(node) {
        if (!node.children || node.children.length === 0) return;
        const info = getCladeInfo(node);
        if (info.dominantCategory === catVal && info.homogeneity >= threshold && info.count >= settings.minCladeSize) {
          node._collapsed = forceCollapse;
          affected++;
          if (forceCollapse) return;
        }
        node.children.forEach(scan);
      }
      scan(activeTreeRoot);
      renderTree();
      fitTreeToScreen(true);
      updateCladeManagementUI();
      showCladeToast(`${forceCollapse ? "Collapsed" : "Expanded"} clades for ${catVal}`);
    }

    function collapseByCurrentCategory() {
      const cladeColDef = getActiveCladeColumnDef();
      const colKey = cladeColDef ? cladeColDef.key : "structural_class";
      const threshold = settings.cladeHomogeneity || 75;

      let count = 0;
      function scan(node) {
        if (!node.children || node.children.length === 0) return;
        const info = getCladeInfo(node);
        if (info.homogeneity >= threshold && info.count >= settings.minCladeSize) {
          node._collapsed = true;
          count++;
          return;
        }
        node.children.forEach(scan);
      }
      scan(activeTreeRoot);
      renderTree();
      fitTreeToScreen(true);
      updateCladeManagementUI();
      showCladeToast(`Collapsed ${count} pure ${cladeColDef.label} clades (${threshold}%+ purity)`);
    }

    function collapseSubclades(minDepth = 4) {
      assignDepths(activeTreeRoot, 0);
      let count = 0;
      function scan(node) {
        if (!node.children || node.children.length === 0) return;
        if (node.cladoDepth >= minDepth && node.children.length > 0) {
          node._collapsed = true;
          count++;
          return;
        }
        node.children.forEach(scan);
      }
      scan(activeTreeRoot);
      renderTree();
      fitTreeToScreen(true);
      updateCladeManagementUI();
      showCladeToast(`Collapsed ${count} sub-clades at depth ≥ ${minDepth}`);
    }

    function collapseTopLineages(targetLeaves = 10) {
      function uncollapseAll(n) {
        n._collapsed = false;
        if (n.children) n.children.forEach(uncollapseAll);
      }
      uncollapseAll(activeTreeRoot);

      function collapseRecursive(n) {
        if (!n.children || n.children.length === 0) return;
        if (n.cladoDepth >= 2) {
          n._collapsed = true;
          return;
        }
        n.children.forEach(collapseRecursive);
      }
      collapseRecursive(activeTreeRoot);

      renderTree();
      fitTreeToScreen(true);
      updateCladeManagementUI();
      showCladeToast(`Summarized into top major structural lineages`);
    }

    function expandAllClades() {
      function scan(node) {
        node._collapsed = false;
        if (node.children) node.children.forEach(scan);
      }
      scan(activeTreeRoot);
      renderTree();
      fitTreeToScreen(true);
      updateCladeManagementUI();
      showCladeToast("Expanded all clades across active tree.");
    }

    function setHomogeneity(val) {
      settings.cladeHomogeneity = parseInt(val);
      document.getElementById("homogeneityVal").textContent = val + "%";
      collapseByCurrentCategory();
    }

    function updateCladeManagementUI() {
      const allLeaves = getAllLeaves(activeTreeRoot);
      const visibleLeaves = getVisibleLeaves(activeTreeRoot);

      const statBadge = document.getElementById("cladeStatBadge");
      if (statBadge) {
        if (visibleLeaves.length === allLeaves.length) {
          statBadge.textContent = `${allLeaves.length} / ${allLeaves.length} visible`;
          statBadge.className = "text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/30";
        } else {
          statBadge.textContent = `${visibleLeaves.length} / ${allLeaves.length} visible`;
          statBadge.className = "text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-amber-500/10 text-amber-400 border border-amber-500/30";
        }
      }

      let collapsedNodesCount = 0;
      function countCollapsed(n) {
        if (n._collapsed) collapsedNodesCount++;
        if (n.children) n.children.forEach(countCollapsed);
      }
      countCollapsed(activeTreeRoot);

      const tabCount = document.getElementById("tabCladeCountBadge");
      if (tabCount) tabCount.textContent = collapsedNodesCount;

      const container = document.getElementById("cladeFamilyList");
      if (!container) return;
      container.innerHTML = "";

      const cladeColDef = getActiveCladeColumnDef();
      if (!cladeColDef) return;

      const listHeader = document.getElementById("cladeListHeader");
      if (listHeader) listHeader.textContent = `Categories in: ${cladeColDef.label}`;

      const catCounts = {};
      const catPlddtSum = {};
      allLeaves.forEach(l => {
        const m = TAXA_METADATA[l.name];
        if (m && m[cladeColDef.key] !== undefined && m[cladeColDef.key] !== null) {
          const val = String(m[cladeColDef.key]);
          catCounts[val] = (catCounts[val] || 0) + 1;
          const p = parseFloat(m.plddt) || 0;
          catPlddtSum[val] = (catPlddtSum[val] || 0) + p;
        }
      });

      const uniqueVals = cladeColDef.values || Object.keys(catCounts);
      uniqueVals.forEach(val => {
        const count = catCounts[val] || 0;
        if (count === 0) return;
        const avgPlddt = (count > 0 && catPlddtSum[val]) ? (catPlddtSum[val] / count).toFixed(1) + " pLDDT" : "";
        const col = (cladeColDef.colors && cladeColDef.colors[val]) ? cladeColDef.colors[val] : "#94a3b8";

        const row = document.createElement("div");
        row.className = "flex items-center justify-between p-2 rounded-lg bg-[var(--card-bg)] border border-[var(--border-color)] hover:border-slate-500 transition";
        row.innerHTML = `
          <div class="flex items-center space-x-2 truncate">
            <span class="w-2.5 h-2.5 rounded-full shrink-0" style="background-color: ${col};"></span>
            <div class="truncate">
              <span class="font-semibold text-[var(--text-main)] truncate text-[11px]">${val}</span>
              <span class="text-[9.5px] text-[var(--text-muted)] ml-1 font-mono">${count} taxa &bull; ${avgPlddt}</span>
            </div>
          </div>
          <div class="flex items-center space-x-1 shrink-0 ml-2">
            <button onclick="toggleCategoryClade('${val.replace(/'/g, "\\\\\\'")}', true)" class="px-1.5 py-0.5 rounded bg-sky-500/10 hover:bg-sky-500/20 text-sky-400 font-medium text-[9px] transition" title="Collapse ${val} clades">▾ Collapse</button>
            <button onclick="toggleCategoryClade('${val.replace(/'/g, "\\\\\\'")}', false)" class="px-1.5 py-0.5 rounded bg-slate-500/10 hover:bg-slate-500/20 text-[var(--text-muted)] hover:text-[var(--text-main)] font-medium text-[9px] transition" title="Expand ${val} clades">▴ Expand</button>
          </div>
        `;
        container.appendChild(row);
      });

      updateCladePartitionUI();
    }
    // =========================================================================
    // SILHOUETTE & CLADE-BASED SCOPED PARTITIONING ENGINE
    // =========================================================================
    const CLADE_PALETTE = [
      "#38bdf8", "#ec4899", "#10b981", "#f59e0b", "#8b5cf6",
      "#06b6d4", "#f43f5e", "#84cc16", "#eab308", "#6366f1",
      "#14b8a6", "#d946ef", "#22c55e", "#f97316", "#a855f7"
    ];

    const cladePartitionState = {
      source: "esm2", // "esm2", "3di", "aa"
      k: 3,
      hideMinorClades: true, // Hide singletons and subclades with <5 sequences by default
      minCladeThreshold: 5,
      highlightedCladeId: null
    };

    function getActiveSilhouetteData() {
      const src = cladePartitionState.source || "esm2";
      const ds = activeDataset || getActiveDataset();
      if (!ds) return null;
      return ds["silhouette_" + src] || (src === "esm2" ? ds.silhouette : null) || null;
    }

    function toggleMinCladeFilter() {
      cladePartitionState.hideMinorClades = !cladePartitionState.hideMinorClades;
      const btn = document.getElementById("btnToggleMinCladeFilter");
      if (btn) {
        if (cladePartitionState.hideMinorClades) {
          btn.className = "badge-sky text-[8.5px] px-1.5 py-0.2 rounded font-semibold cursor-pointer transition hover:opacity-90";
          btn.innerHTML = "🛡️ &ge;5 Seqs: ON";
        } else {
          btn.className = "text-[8.5px] px-1.5 py-0.2 rounded bg-slate-800 text-slate-400 border border-[var(--border-color)] font-semibold cursor-pointer transition hover:bg-slate-700";
          btn.innerHTML = "🛡️ &ge;5 Seqs: OFF";
        }
      }
      updateCladePartitionUI();
    }

    function setCladePartitionSource(source) {
      cladePartitionState.source = source;
      const sil = getActiveSilhouetteData();
      if (sil && sil.best_k) {
        cladePartitionState.k = sil.best_k;
      }
      const slider = document.getElementById("cladeKSlider");
      if (slider) {
        const totalTaxa = Object.keys(TAXA_METADATA).length || 10;
        const maxK = sil && sil.profile ? Math.max(...sil.profile.map(p => p.k)) : Math.min(40, totalTaxa - 1);
        slider.max = Math.max(5, maxK);
        slider.value = cladePartitionState.k;
      }
      const kValEl = document.getElementById("cladeKVal");
      if (kValEl) kValEl.textContent = `k = ${cladePartitionState.k}`;

      updateCladePartitionUI();
    }

    function setCladePartitionK(kVal) {
      const k = parseInt(kVal);
      cladePartitionState.k = k;
      const slider = document.getElementById("cladeKSlider");
      if (slider && parseInt(slider.value) !== k) slider.value = k;
      const kValEl = document.getElementById("cladeKVal");
      if (kValEl) kValEl.textContent = `k = ${k}`;
      updateCladePartitionUI();
    }

    function renderSilhouetteSparkline() {
      const svgEl = document.getElementById("silhouetteSparklineSvg");
      if (!svgEl) return;
      const sil = getActiveSilhouetteData();
      const silBox = document.getElementById("silhouetteProfileBox");

      if (!sil || !sil.profile || sil.profile.length === 0) {
        if (silBox) silBox.classList.add("hidden");
        return;
      }
      if (silBox) silBox.classList.remove("hidden");

      // Update Title Label
      const titleMap = { "esm2": "ESM-2 PLM", "3di": "3Di Structure", "aa": "Amino Acid Sequence" };
      const titleEl = document.getElementById("silProfileTitle");
      if (titleEl) {
        titleEl.innerHTML = `<span>📊</span> ${titleMap[cladePartitionState.source] || "Phylogenetic"} Silhouette Profile S(k)`;
      }

      // Update Peak Badge and Slider Max
      const peakBadge = document.getElementById("silPeakBadge");
      if (peakBadge) {
        peakBadge.textContent = `Peak k=${sil.best_k} (S=${sil.best_score.toFixed(3)})`;
      }
      const slider = document.getElementById("cladeKSlider");
      if (slider) {
        const maxK = Math.max(...sil.profile.map(p => p.k));
        slider.max = maxK;
      }

      // Quick Peak Chips (Mathematical suggestions)
      const chipsEl = document.getElementById("silPeakChips");
      if (chipsEl) {
        chipsEl.innerHTML = "";
        const peaks = sil.peaks || [];
        peaks.forEach(p => {
          const btn = document.createElement("button");
          const isCurrent = (p.k === cladePartitionState.k);
          btn.className = `px-2 py-0.5 rounded font-mono font-bold text-[9.5px] transition cursor-pointer flex items-center gap-1 ${
            isCurrent ? "bg-emerald-500 text-white shadow-sm ring-1 ring-emerald-300" : "bg-emerald-500/20 hover:bg-emerald-500/35 text-emerald-300 border border-emerald-500/40"
          }`;
          btn.innerHTML = `<span>⭐</span> k=${p.k} (${p.score.toFixed(2)})`;
          btn.onclick = () => setCladePartitionK(p.k);
          chipsEl.appendChild(btn);
        });
      }

      // Render SVG Sparkline
      const width = 280;
      const height = 75;
      svgEl.setAttribute("viewBox", `0 0 ${width} ${height}`);
      svgEl.innerHTML = "";

      const profile = sil.profile;
      const minK = 2;
      const maxK = Math.max(...profile.map(p => p.k)) || 35;
      const maxScore = Math.max(...profile.map(p => p.score), 0.75);

      const padX = 20;
      const padY = 12;
      const scaleX = (k) => padX + ((k - minK) / (maxK - minK || 1)) * (width - 2 * padX);
      const scaleY = (s) => height - padY - (Math.max(0, s) / maxScore) * (height - 2 * padY);

      // SVG Definitions (Gradient)
      const defs = document.createElementNS("http://www.w3.org/2000/svg", "defs");
      defs.innerHTML = `
        <linearGradient id="silGrad" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0%" stop-color="#10b981" stop-opacity="0.35"/>
          <stop offset="100%" stop-color="#10b981" stop-opacity="0.0"/>
        </linearGradient>
      `;
      svgEl.appendChild(defs);

      // Baseline grid
      const grid = document.createElementNS("http://www.w3.org/2000/svg", "g");
      const y05 = scaleY(0.5);
      grid.innerHTML = `
        <line x1="${padX}" y1="${scaleY(0)}" x2="${width - padX}" y2="${scaleY(0)}" stroke="var(--border-color)" stroke-width="0.75"/>
        <line x1="${padX}" y1="${y05}" x2="${width - padX}" y2="${y05}" stroke="#10b981" stroke-dasharray="2 2" stroke-width="0.75" opacity="0.4"/>
        <text x="${width - padX - 2}" y="${y05 - 2}" fill="#10b981" opacity="0.6" font-size="7" font-family="monospace" text-anchor="end">S=0.5</text>
      `;
      svgEl.appendChild(grid);

      // Sparkline Path & Fill Area
      let pathD = "";
      let areaD = `M ${scaleX(profile[0].k)} ${scaleY(0)} `;

      profile.forEach((p, i) => {
        const x = scaleX(p.k);
        const y = scaleY(p.score);
        if (i === 0) {
          pathD += `M ${x} ${y}`;
          areaD += `L ${x} ${y}`;
        } else {
          pathD += ` L ${x} ${y}`;
          areaD += ` L ${x} ${y}`;
        }
      });
      areaD += ` L ${scaleX(profile[profile.length - 1].k)} ${scaleY(0)} Z`;

      const strokeCol = isDarkTheme() ? "#10b981" : "#047857";
      const areaPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
      areaPath.setAttribute("d", areaD);
      areaPath.setAttribute("fill", "url(#silGrad)");
      svgEl.appendChild(areaPath);

      const linePath = document.createElementNS("http://www.w3.org/2000/svg", "path");
      linePath.setAttribute("d", pathD);
      linePath.setAttribute("fill", "none");
      linePath.setAttribute("stroke", strokeCol);
      linePath.setAttribute("stroke-width", "2");
      linePath.setAttribute("stroke-linecap", "round");
      linePath.setAttribute("stroke-linejoin", "round");
      svgEl.appendChild(linePath);

      // Data Points
      profile.forEach(p => {
        const x = scaleX(p.k);
        const y = scaleY(p.score);
        const isCurrent = (p.k === cladePartitionState.k);
        const isPeak = p.is_peak;

        const gPoint = document.createElementNS("http://www.w3.org/2000/svg", "g");
        gPoint.style.cursor = "pointer";
        gPoint.onclick = () => setCladePartitionK(p.k);

        const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
        circle.setAttribute("cx", x);
        circle.setAttribute("cy", y);
        circle.setAttribute("r", isCurrent ? 5.5 : (isPeak ? 4.5 : 2.5));
        circle.setAttribute("fill", isCurrent ? "#38bdf8" : (isPeak ? "#fbbf24" : "#10b981"));
        circle.setAttribute("stroke", isCurrent ? "#ffffff" : (isPeak ? "#78350f" : "#064e3b"));
        circle.setAttribute("stroke-width", isCurrent ? "2" : "1");

        const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
        title.textContent = `k=${p.k}: Silhouette S=${p.score.toFixed(4)}${isPeak ? " ⭐ PEAK" : ""} (Click to select)`;
        gPoint.appendChild(title);
        gPoint.appendChild(circle);

        if (isPeak) {
          const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
          txt.setAttribute("x", x);
          txt.setAttribute("y", y - 6);
          txt.setAttribute("font-size", "7.5");
          txt.setAttribute("fill", "#fbbf24");
          txt.setAttribute("text-anchor", "middle");
          txt.setAttribute("font-weight", "bold");
          txt.textContent = "⭐";
          gPoint.appendChild(txt);
        }

        svgEl.appendChild(gPoint);
      });
    }

    function getCladesForCurrentCut() {
      const k = cladePartitionState.k;
      const source = cladePartitionState.source;
      const sil = getActiveSilhouetteData();

      let rawClades = [];

      if (sil && sil.profile && sil.profile.length > 0) {
        const prof = sil.profile.find(p => p.k === k) || sil.profile[0];
        const clusters = prof.clusters || {};
        Object.keys(clusters).forEach((cid, idx) => {
          rawClades.push({
            id: idx + 1,
            name: `Clade ${idx + 1}`,
            taxa: clusters[cid] || [],
            score: prof.score
          });
        });
      } else {
        // Tree-based bipartition / greedy cut for 3Di or AA
        const srcRoot = (source === "aa") ? rawRootAA : rawRoot3Di;
        if (srcRoot) {
          let frontier = [srcRoot];
          while (frontier.length < k) {
            let bestIdx = -1;
            let maxLeaves = -1;
            for (let i = 0; i < frontier.length; i++) {
              const n = frontier[i];
              if (n.children && n.children.length > 1) {
                const leafCount = getAllLeaves(n).length;
                if (leafCount > maxLeaves) {
                  maxLeaves = leafCount;
                  bestIdx = i;
                }
              }
            }
            if (bestIdx === -1) break;
            const toExpand = frontier.splice(bestIdx, 1)[0];
            frontier.push(...toExpand.children);
          }
          frontier.forEach((node, idx) => {
            rawClades.push({
              id: idx + 1,
              name: `Lineage ${idx + 1}`,
              taxa: getAllLeaves(node).map(l => l.name),
              score: null
            });
          });
        }
      }

      // Compute dominant metadata annotations for each clade
      const colDef = getActiveCladeColumnDef() || getActiveColorColumnDef();
      const colKey = colDef ? colDef.key : null;

      rawClades.forEach(c => {
        const counts = {};
        c.taxa.forEach(t => {
          const m = TAXA_METADATA[t];
          if (m && colKey && m[colKey] !== undefined) {
            const v = String(m[colKey]);
            counts[v] = (counts[v] || 0) + 1;
          }
        });
        let domVal = "";
        let domCount = 0;
        Object.entries(counts).forEach(([v, cnt]) => {
          if (cnt > domCount) {
            domCount = cnt;
            domVal = v;
          }
        });
        c.dominantVal = domVal;
        c.dominantPct = c.taxa.length > 0 ? Math.round((domCount / c.taxa.length) * 100) : 0;
        c.dominantCount = domCount;
      });

      // Sort by size descending
      rawClades.sort((a, b) => b.taxa.length - a.taxa.length);
      // Re-index cleanly
      rawClades.forEach((c, idx) => {
        c.displayId = idx + 1;
      });

      return rawClades;
    }

    function updateCladePartitionUI() {
      const source = cladePartitionState.source;
      const k = cladePartitionState.k;
      const sil = getActiveSilhouetteData();

      // Update K Score Badge
      const kScoreEl = document.getElementById("cladeKScore");
      if (kScoreEl) {
        if (sil && sil.profile) {
          const prof = sil.profile.find(p => p.k === k);
          kScoreEl.textContent = prof ? `S = ${prof.score.toFixed(3)}` : "";
          kScoreEl.classList.remove("hidden");
        } else {
          kScoreEl.textContent = "Tree Cut";
          kScoreEl.classList.remove("hidden");
        }
      }

      renderSilhouetteSparkline();

      // Render Clade Roster with Singleton / Minor Subclade Handling (<5 sequences)
      const clades = getCladesForCurrentCut();
      const rosterEl = document.getElementById("cladePartitionRoster");
      const countEl = document.getElementById("cladeRosterCount");

      const majorClades = clades.filter(c => c.taxa.length >= cladePartitionState.minCladeThreshold);
      const minorClades = clades.filter(c => c.taxa.length < cladePartitionState.minCladeThreshold);

      if (countEl) {
        if (cladePartitionState.hideMinorClades) {
          countEl.textContent = `${majorClades.length} Major (${cladePartitionState.minCladeThreshold}+ seqs)`;
        } else {
          countEl.textContent = `${clades.length} Clades`;
        }
      }

      const btnExitRoster = document.getElementById("btnExitScopeRoster");
      if (btnExitRoster) {
        if (settings.scopedClade) {
          btnExitRoster.classList.remove("hidden");
        } else {
          btnExitRoster.classList.add("hidden");
        }
      }

      if (!rosterEl) return;
      rosterEl.innerHTML = "";

      const totalCohortTaxa = Object.keys(TAXA_METADATA).length || 1;

      // Helper to render a clade card
      function createCladeCard(c, isCompact = false) {
        const color = (customPaletteState.isActive && customPaletteState.applyToClades && customPaletteState.palette.length > 0)
          ? customPaletteState.palette[(c.displayId - 1) % customPaletteState.palette.length]
          : CLADE_PALETTE[(c.displayId - 1) % CLADE_PALETTE.length];
        const isScoped = settings.scopedClade && (settings.scopedClade.id === c.id || settings.scopedClade.name === c.name);
        const pctOfCohort = ((c.taxa.length / totalCohortTaxa) * 100).toFixed(1);

        const card = document.createElement("div");
        card.className = `p-2.5 rounded-lg border transition space-y-1.5 ${
          isScoped 
            ? "border-sky-400 bg-sky-500/15 ring-2 ring-sky-500/30 shadow-md" 
            : (isCompact ? "border-[var(--border-color)]/60 bg-[var(--card-bg)]/80 hover:border-slate-500" : "border-[var(--border-color)] bg-[var(--chip-bg)] hover:border-slate-500")
        }`;

        card.innerHTML = `
          <div class="flex items-center justify-between">
            <div class="flex items-center space-x-2 truncate">
              <span class="w-2.5 h-2.5 rounded-full shrink-0 shadow-sm" style="background-color: ${color};"></span>
              <span class="font-bold text-xs text-[var(--text-main)] truncate">${c.name}</span>
              ${isScoped ? '<span class="px-1.5 py-0.2 rounded text-[8.5px] font-bold bg-sky-500 text-white animate-pulse">ACTIVE SCOPE</span>' : ''}
            </div>
            <span class="text-[10px] font-mono font-bold" style="color: var(--badge-sky-text);">${c.taxa.length.toLocaleString()} taxa (${pctOfCohort}%)</span>
          </div>

          <div class="flex items-center justify-between text-[9.5px] text-[var(--text-muted)]">
            <div class="truncate mr-2">
              ${c.dominantVal ? `<span>Majority: <strong class="text-[var(--text-main)]">${c.dominantVal}</strong> (${c.dominantPct}%)</span>` : '<span>Diverse lineage</span>'}
            </div>
            ${c.score !== null && c.score !== undefined ? `<span class="font-mono font-bold shrink-0" style="color: var(--badge-emerald-text);">Cohesion S=${c.score.toFixed(3)}</span>` : ''}
          </div>

          <div class="flex items-center space-x-1.5 pt-1 border-t border-[var(--border-color)]/60">
            ${isScoped ? `
              <button onclick="exitCladeScope()" class="flex-1 py-1 px-2 rounded badge-rose font-bold text-[10px] transition text-center cursor-pointer">
                ✕ Reset to Full Cohort
              </button>
            ` : `
              <button onclick="scopeToCladeByIndex(${c.id})" class="flex-1 py-1 px-2 rounded badge-sky font-bold text-[10px] transition flex items-center justify-center space-x-1 shadow-sm cursor-pointer">
                <span>🔍</span>
                <span>Scope Entire Analysis</span>
              </button>
            `}
            <button onclick="highlightCladeTaxaByIndex(${c.id})" class="py-1 px-2.5 rounded bg-[var(--chip-bg)] hover:bg-slate-500/20 text-[var(--text-main)] border border-[var(--border-color)] text-[10px] font-medium transition cursor-pointer" title="Highlight this clade in current view without filtering others">
              👁️ View
            </button>
          </div>
        `;
        return card;
      }

      // 1. Render Major Clades (taxa >= 5)
      majorClades.forEach(c => {
        rosterEl.appendChild(createCladeCard(c, false));
      });

      // 2. Render Minor Clades / Singletons (taxa < 5)
      if (minorClades.length > 0) {
        if (cladePartitionState.hideMinorClades) {
          const minorTotalTaxa = minorClades.reduce((sum, c) => sum + c.taxa.length, 0);
          const detailsEl = document.createElement("details");
          detailsEl.className = "group bg-slate-900/60 rounded-lg border border-[var(--border-color)] overflow-hidden shadow-sm";
          
          detailsEl.innerHTML = `
            <summary class="p-2 text-[10px] font-medium text-[var(--text-muted)] cursor-pointer flex items-center justify-between hover:text-slate-300 transition select-none">
              <span class="flex items-center gap-1.5">
                <span>🔍</span>
                <span>${minorClades.length} Minor Lineages &amp; Singletons (&lt;5 seqs, ${minorTotalTaxa} taxa total)</span>
              </span>
              <span class="text-[9px] group-open:rotate-180 transition-transform">▼</span>
            </summary>
            <div id="minorCladesContainer" class="p-2 pt-0 space-y-1.5 border-t border-[var(--border-color)]/40 mt-1 max-h-48 overflow-y-auto custom-scroll">
            </div>
          `;

          const container = detailsEl.querySelector("#minorCladesContainer");
          minorClades.forEach(c => {
            container.appendChild(createCladeCard(c, true));
          });

          rosterEl.appendChild(detailsEl);
        } else {
          // If hide filter is OFF, render minor clades inline
          minorClades.forEach(c => {
            rosterEl.appendChild(createCladeCard(c, false));
          });
        }
      }
    }

    function scopeToCladeByIndex(cladeId) {
      const clades = getCladesForCurrentCut();
      const clade = clades.find(c => c.id === cladeId);
      if (!clade || clade.taxa.length === 0) return;

      const cladeTitle = `${clade.name}${clade.dominantVal ? " (" + clade.dominantVal + ")" : ""}`;
      settings.scopedClade = {
        id: clade.id,
        name: cladeTitle,
        taxa: new Set(clade.taxa),
        taxaArray: clade.taxa,
        score: clade.score,
        source: cladePartitionState.source,
        k: cladePartitionState.k
      };

      // Update Floating Banner
      const banner = document.getElementById("scopedCladeBanner");
      const nameEl = document.getElementById("scopedCladeName");
      const countEl = document.getElementById("scopedCladeCount");
      const scoreBadge = document.getElementById("scopedCladeScoreBadge");
      if (nameEl) nameEl.textContent = cladeTitle;
      if (countEl) countEl.textContent = clade.taxa.length.toLocaleString();
      if (scoreBadge) {
        scoreBadge.textContent = (clade.score !== null && clade.score !== undefined) 
          ? `Silhouette S = ${clade.score.toFixed(3)}` 
          : `Divergence Cut k=${cladePartitionState.k}`;
      }
      if (banner) {
        banner.classList.remove("hidden");
        banner.classList.add("flex");
      }

      // Shift filter banner down if both are active to prevent overlap
      const fBanner = document.getElementById("activeFilterBanner");
      if (fBanner && !fBanner.classList.contains("hidden")) {
        fBanner.classList.add("top-11");
        fBanner.classList.remove("top-3");
      }

      // Update Sidebar Scope Badge
      const scopeBadge = document.getElementById("cladeScopeActiveBadge");
      if (scopeBadge) {
        scopeBadge.textContent = `Scoped: ${clade.name}`;
        scopeBadge.className = "text-[9px] font-mono px-2 py-0.5 rounded-full bg-sky-500/20 text-sky-300 border border-sky-500/40 font-bold";
      }

      // Update Badge Taxa
      const bTaxa = document.getElementById("badgeTaxa");
      if (bTaxa) {
        bTaxa.textContent = `${clade.taxa.length.toLocaleString()} Taxa (Scoped ${clade.name})`;
      }

      // Auto-select a representative taxon in this clade if selected is outside
      if (!settings.selectedTaxon || !settings.scopedClade.taxa.has(settings.selectedTaxon)) {
        const firstTaxon = clade.taxa[0];
        if (firstTaxon) selectTaxon(firstTaxon);
      }

      // Apply Pruning and Rerender
      applyCurrentRooting();
      showCladeToast(`Isolated analysis to ${cladeTitle} (${clade.taxa.length} taxa).`);
    }

    function exitCladeScope() {
      settings.scopedClade = null;

      // Hide Floating Banner
      const banner = document.getElementById("scopedCladeBanner");
      if (banner) {
        banner.classList.add("hidden");
        banner.classList.remove("flex");
      }

      // Restore filter banner position
      const fBanner = document.getElementById("activeFilterBanner");
      if (fBanner) {
        fBanner.classList.remove("top-11");
        fBanner.classList.remove("top-14");
        fBanner.classList.add("top-3");
      }

      // Update Sidebar Scope Badge
      const scopeBadge = document.getElementById("cladeScopeActiveBadge");
      if (scopeBadge) {
        scopeBadge.textContent = "Full Cohort";
        scopeBadge.className = "text-[9px] font-mono px-2 py-0.5 rounded-full bg-slate-700/60 text-slate-400 border border-slate-600/40";
      }

      // Restore Cohort Badge Taxa
      const bTaxa = document.getElementById("badgeTaxa");
      if (bTaxa) {
        if (currentScale === "1193") bTaxa.textContent = "1,193 ESMFold Designs";
        else if (currentScale === "500") bTaxa.textContent = "500 Viral Structures";
        else bTaxa.textContent = "6 Benchmark Taxa";
      }

      // Reapply Rooting to restore full tree
      applyCurrentRooting();
      if (typeof renderMsa === 'function') {
        if (typeof msaState !== 'undefined') {
          msaState._minimapCacheKey = null;
          msaState.scrollY = 0;
        }
        renderMsa();
      }
      showCladeToast("Reset scope to full cohort view.");
    }

    function highlightCladeTaxaByIndex(cladeId) {
      const clades = getCladesForCurrentCut();
      const clade = clades.find(c => c.id === cladeId);
      if (!clade || clade.taxa.length === 0) return;

      // Focus on first taxon or highlight
      const firstTaxon = clade.taxa[0];
      if (firstTaxon) {
        selectTaxon(firstTaxon);
        centerOnSelection();
      }
      showCladeToast(`Focused on ${clade.name} (${clade.taxa.length} taxa).`);
    }


    function showCladeToast(msg) {
      const t = document.getElementById("cladeToast");
      if (t) {
        t.textContent = msg;
        t.classList.remove("hidden");
        clearTimeout(t._timer);
        t._timer = setTimeout(() => {
          t.classList.add("hidden");
        }, 3200);
      }
      if (typeof showToastNotification === 'function') {
        showToastNotification(msg);
      }
    }

    // =========================================================================
    // PIPELINE STUDIO & VIRO3D DATASET REQUEST ENGINE
    // =========================================================================
    const pipelineStudioState = {
      source: "viro3d", // "viro3d", "alphafold", or "local"
      isRunning: false,
      pollInterval: null
    };

    function setPipelineSource(src) {
      pipelineStudioState.source = src;
      const btnViro = document.getElementById("btnPipeSourceViro");
      const btnAfdb = document.getElementById("btnPipeSourceAlphaFold");
      const btnLocal = document.getElementById("btnPipeSourceLocal");
      const secViro = document.getElementById("pipeViroSection");
      const secAfdb = document.getElementById("pipeAlphaFoldSection");
      const secLocal = document.getElementById("pipeLocalSection");

      const activeClasses = "py-1 px-1 rounded font-semibold text-center bg-sky-500 text-slate-950 text-[10px] transition shadow-sm cursor-pointer";
      const viroActiveClasses = "py-1 px-1 rounded font-semibold text-center bg-amber-500 text-slate-950 text-[10px] transition shadow-sm cursor-pointer";
      const inactiveClasses = "py-1 px-1 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10px] transition cursor-pointer";

      if (btnViro) btnViro.className = src === "viro3d" ? viroActiveClasses : inactiveClasses;
      if (btnAfdb) btnAfdb.className = src === "alphafold" ? activeClasses : inactiveClasses;
      if (btnLocal) btnLocal.className = src === "local" ? viroActiveClasses : inactiveClasses;

      if (secViro) secViro.classList.toggle("hidden", src !== "viro3d");
      if (secAfdb) secAfdb.classList.toggle("hidden", src !== "alphafold");
      if (secLocal) secLocal.classList.toggle("hidden", src !== "local");

      updatePipelineCommandPreview();
    }

    function setPipelinePreset(preset) {
      const input = document.getElementById("pipeQualifierInput");
      if (input) {
        input.value = preset;
      }
      const outDirInput = document.getElementById("pipeOutputDirInput");
      if (outDirInput) {
        outDirInput.value = `${preset}_workflow`;
      }
      updatePipelineCommandPreview();
    }

    function setAlphaFoldPreset(preset) {
      const input = document.getElementById("pipeAfdbInput");
      if (input) {
        input.value = preset;
      }
      const outDirInput = document.getElementById("pipeOutputDirInput");
      if (outDirInput) {
        outDirInput.value = "afdb_custom_workflow";
      }
      updatePipelineCommandPreview();
    }

    function updatePipelineCommandPreview() {
      const disp = document.getElementById("pipeCommandDisplay");
      if (!disp) return;

      const src = pipelineStudioState.source;
      const outDir = (document.getElementById("pipeOutputDirInput")?.value || "viral_custom_workflow").trim();
      const treeType = document.getElementById("pipeTreeTypeSelect")?.value || "both";
      const aligner = document.getElementById("pipeAlignerSelect")?.value || "foldmason";
      const matrix = document.getElementById("pipeMatrixSelect")?.value || "alphafold";
      const bootstrap = document.getElementById("pipeBootstrapToggle")?.checked ?? true;
      const threads = document.getElementById("pipeThreadsSelect")?.value || "AUTO";

      let parts = ["python scripts/viral_phylogenetics.py pipeline"];

      if (src === "viro3d") {
        const qual = (document.getElementById("pipeQualifierInput")?.value || "glycoprotein").trim();
        const count = document.getElementById("pipeCountSlider")?.value || "50";
        parts.push(`--source viro3d`);
        parts.push(`--qualifier ${qual}`);
        parts.push(`--count ${count}`);
      } else if (src === "alphafold") {
        const afInput = (document.getElementById("pipeAfdbInput")?.value || "Q99720, P87666, P0DTC2").trim();
        const afFmt = document.getElementById("pipeAfdbFormatSelect")?.value || "pdb";
        const afPae = document.getElementById("pipeAfdbPaeToggle")?.checked;
        parts.push(`--source alphafold`);
        parts.push(`-u "${afInput}"`);
        if (afFmt !== "pdb") parts.push(`--format ${afFmt}`);
        if (afPae) parts.push(`--download-pae`);
      } else {
        const folder = (document.getElementById("pipeLocalDirInput")?.value || "300_rdrp/structures").trim();
        const meta = (document.getElementById("pipeLocalMetaInput")?.value || "").trim();
        parts.push(`--input-folder ${folder}`);
        if (meta) parts.push(`--metadata ${meta}`);
      }

      if (aligner === "mafft") {
        parts.push("--aligner mafft");
      }

      const minCov = document.getElementById("pipeCoverageInput")?.value;
      const multiAln = document.getElementById("pipeMultiAlnToggle")?.checked;
      if (minCov && minCov !== "0" && minCov !== "100") {
        const covVal = (parseFloat(minCov) / 100).toFixed(2);
        parts.push(`--min-coverage ${covVal}`);
      }
      if (multiAln) {
        parts.push("--multi-alignment");
      }
      parts.push(`--tree-type ${treeType}`);
      parts.push(`--matrix ${matrix}`);
      if (bootstrap) {
        parts.push("--bootstrap 1000");
      } else {
        parts.push("--bootstrap 0");
      }
      parts.push(`--threads ${threads}`);

      const embedToggle = document.getElementById("pipeEmbedToggle")?.checked;
      if (embedToggle) {
        const embModel = document.getElementById("pipeEmbedModelSelect")?.value || "esm2";
        const embClust = document.getElementById("pipeEmbedClusteringSelect")?.value || "upgma";
        parts.push("--embed");
        parts.push(`--embed-model ${embModel}`);
        parts.push(`--embed-clustering ${embClust}`);
      }

      parts.push(`--output-dir ${outDir}`);

      disp.textContent = parts.join(" \\\n  ");
    }

    function copyPipelineCommand() {
      const disp = document.getElementById("pipeCommandDisplay");
      const badge = document.getElementById("pipeCommandCopiedBadge");
      if (!disp) return;

      const singleLineCmd = disp.textContent.split("\\n").map(s => s.trim().replace(/\\$/, "")).join(" ").replace(/\\s+/g, " ").trim();
      navigator.clipboard.writeText(singleLineCmd).then(() => {
        if (badge) {
          badge.classList.remove("hidden");
          setTimeout(() => badge.classList.add("hidden"), 2200);
        }
      }).catch(() => {
        // Fallback
        const ta = document.createElement("textarea");
        ta.value = singleLineCmd;
        document.body.appendChild(ta);
        ta.select();
        document.execCommand("copy");
        document.body.removeChild(ta);
        if (badge) {
          badge.classList.remove("hidden");
          setTimeout(() => badge.classList.add("hidden"), 2200);
        }
      });
    }

    function downloadPipelineScript() {
      const disp = document.getElementById("pipeCommandDisplay");
      if (!disp) return;

      const scriptContent = `#!/usr/bin/env bash\nset -e\n\n# Viral Structural Phylogenetics Automated Pipeline Script\n# Generated by Interactive Visualization Suite\n\n${disp.textContent}\n`;
      const blob = new Blob([scriptContent], { type: "application/x-sh" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "run_viral_pipeline.sh";
      a.click();
      URL.revokeObjectURL(url);
    }

    async function queryViro3dApi() {
      const qual = document.getElementById("pipeQualifierInput")?.value?.trim() || "glycoprotein";
      const count = document.getElementById("pipeCountSlider")?.value || 50;
      const statusDiv = document.getElementById("viro3dCheckStatus");

      if (!statusDiv) return;
      statusDiv.classList.remove("hidden");
      statusDiv.innerHTML = `<div class="text-amber-400 animate-pulse flex items-center space-x-1.5 font-mono text-[10px]"><span>⏳</span><span>Querying Viro3D database for "${qual}"...</span></div>`;

      try {
        let res = null;
        try {
          const isLocal = window.location.protocol === "http:" && (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1");
          const bridgeUrl = isLocal ? `/api/viro3d_query?qualifier=${encodeURIComponent(qual)}&count=10` : `http://localhost:8000/api/viro3d_query?qualifier=${encodeURIComponent(qual)}&count=10`;
          const bridgeRes = await fetch(bridgeUrl, { signal: AbortSignal.timeout(3000) });
          if (bridgeRes.ok) {
            res = await bridgeRes.json();
          }
        } catch (_) {}

        if (!res) {
          const directUrl = `https://viro3d.cvr.gla.ac.uk/api/proteins/protein_name/?qualifier=${encodeURIComponent(qual)}&page_size=10&page_num=1`;
          const directRes = await fetch(directUrl, { mode: "cors", signal: AbortSignal.timeout(4000) });
          if (directRes.ok) {
            res = await directRes.json();
          }
        }

        if (res && res.protein_structures && res.protein_structures.length > 0) {
          const items = res.protein_structures;
          const totalHits = res.total_records || items.length;
          const families = [...new Set(items.map(it => it.family).filter(Boolean))].slice(0, 4);

          statusDiv.innerHTML = `
            <div class="space-y-1">
              <div class="flex items-center justify-between text-emerald-400 font-bold">
                <span>✅ Found in Viro3D</span>
                <span class="font-mono text-[9px] px-1 rounded bg-emerald-500/20">${totalHits} hits</span>
              </div>
              <div class="text-[9.5px] text-slate-300">
                <span class="text-[var(--text-muted)]">Families:</span> ${families.join(", ") || "Diverse lineages"}
              </div>
              <div class="text-[9px] text-slate-400 max-h-16 overflow-y-auto font-mono space-y-0.5 custom-scroll pt-0.5">
                ${items.slice(0, 4).map(it => `<div>&bull; <strong>${it.record_id || it.accession}</strong>: ${it.protein_name || it.product || qual}</div>`).join("")}
              </div>
              <div class="pt-1 flex justify-between items-center text-[9px]">
                <a href="https://viro3d.cvr.gla.ac.uk" target="_blank" class="text-sky-400 hover:underline">Open Viro3D Portal ↗</a>
                <span class="text-emerald-400 font-semibold">Ready to fetch!</span>
              </div>
            </div>
          `;
        } else {
          statusDiv.innerHTML = `
            <div class="space-y-1 text-slate-300">
              <div class="flex items-center justify-between text-amber-400 font-bold">
                <span>🌐 Viro3D Target Configured</span>
                <a href="https://viro3d.cvr.gla.ac.uk" target="_blank" class="text-sky-400 hover:underline font-mono text-[9px]">Portal ↗</a>
              </div>
              <div class="text-[9.5px]">Target: <strong>${qual}</strong> (${count} requested)</div>
              <div class="text-[9px] text-[var(--text-muted)]">Execute command below to download structures via the Viro3D REST client.</div>
            </div>
          `;
        }
      } catch (err) {
        statusDiv.innerHTML = `
          <div class="space-y-1 text-slate-300">
            <div class="flex items-center justify-between text-amber-400 font-bold">
              <span>🌐 Viro3D Target Configured</span>
              <a href="https://viro3d.cvr.gla.ac.uk" target="_blank" class="text-sky-400 hover:underline font-mono text-[9px]">Portal ↗</a>
            </div>
            <div class="text-[9.5px]">Target: <strong>${qual}</strong> (${count} requested)</div>
            <div class="text-[9px] text-[var(--text-muted)]">Execute command below to download structures via the Viro3D REST client.</div>
          </div>
        `;
      }
    }

    async function queryAlphaFoldApi() {
      const rawInput = document.getElementById("pipeAfdbInput")?.value?.trim() || "Q99720";
      const statusDiv = document.getElementById("afdbCheckStatus");
      if (!statusDiv) return;

      statusDiv.classList.remove("hidden");
      statusDiv.innerHTML = `<div class="text-sky-400 animate-pulse flex items-center space-x-1.5 font-mono text-[10px]"><span>⏳</span><span>Querying AlphaFold Database for "${rawInput}"...</span></div>`;

      // Extract accessions or keyword
      let ids = rawInput.split(",").map(s => s.trim().toUpperCase()).filter(s => s.length > 0);

      // If text query rather than accessions, try searching UniProt API
      if (ids.length === 1 && !/^[A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2}$/i.test(ids[0]) && !ids[0].includes("_")) {
        try {
          const uSearchUrl = `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(ids[0])}&size=4&fields=accession,id,gene_names,organism_name,length`;
          const uRes = await fetch(uSearchUrl, { mode: "cors", signal: AbortSignal.timeout(4000) });
          if (uRes.ok) {
            const uData = await uRes.json();
            const resolved = (uData.results || []).map(r => r.primaryAccession).filter(Boolean);
            if (resolved.length > 0) {
              ids = resolved;
            }
          }
        } catch (_) {}
      }

      ids = ids.slice(0, 6);

      try {
        const entries = [];
        for (const uid of ids) {
          try {
            const afUrl = `https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(uid)}`;
            const resp = await fetch(afUrl, { mode: "cors", signal: AbortSignal.timeout(3500) });
            if (resp.ok) {
              const data = await resp.json();
              if (data && data.length > 0) {
                const best = data.find(e => e.uniprotAccession === uid) || data[0];
                entries.push(best);
              }
            }
          } catch (_) {}
        }

        if (entries.length > 0) {
          statusDiv.innerHTML = `
            <div class="space-y-2">
              <div class="flex items-center justify-between text-sky-400 font-bold border-b border-slate-800 pb-1">
                <span class="flex items-center space-x-1"><span>🧬</span><span>AlphaFold Database (${entries.length} structures)</span></span>
                <a href="https://alphafold.ebi.ac.uk" target="_blank" class="text-sky-400 hover:underline font-mono text-[9px]">AFDB Portal ↗</a>
              </div>
              <div class="grid grid-cols-1 gap-1.5 max-h-48 overflow-y-auto custom-scroll pr-1">
                ${entries.map(e => {
                  const plddt = e.globalMetricValue ? Number(e.globalMetricValue).toFixed(1) : "N/A";
                  let plddtBadge = "bg-slate-700 text-slate-300";
                  if (e.globalMetricValue) {
                    const val = Number(e.globalMetricValue);
                    if (val >= 90) plddtBadge = "bg-emerald-500/20 text-emerald-300 border border-emerald-500/40";
                    else if (val >= 70) plddtBadge = "bg-sky-500/20 text-sky-300 border border-sky-500/40";
                    else if (val >= 50) plddtBadge = "bg-amber-500/20 text-amber-300 border border-amber-500/40";
                    else plddtBadge = "bg-rose-500/20 text-rose-300 border border-rose-500/40";
                  }
                  return `
                    <div class="p-1.5 rounded bg-slate-800/80 border border-slate-700 text-[9.5px] space-y-0.5">
                      <div class="flex items-center justify-between">
                        <span class="font-mono font-bold text-sky-300">${e.uniprotAccession} <span class="text-slate-400 font-normal">(${e.uniprotId || e.gene || ""})</span></span>
                        <span class="font-mono text-[9px] px-1.5 py-0.2 rounded-full font-bold ${plddtBadge}" title="Mean pLDDT Confidence">${plddt} pLDDT</span>
                      </div>
                      <div class="text-[9px] text-slate-300 truncate">${e.uniprotDescription || e.organismScientificName || "Viral protein"}</div>
                      <div class="flex items-center justify-between text-[8.5px] text-slate-400 pt-0.5">
                        <span>Organism: <em class="text-slate-300">${e.organismScientificName || "Unknown"}</em></span>
                        <span class="font-mono">${e.sequenceEnd || "?"} aa</span>
                      </div>
                      <div class="flex items-center justify-end space-x-2 pt-0.5 text-[8.5px] font-mono">
                        ${e.pdbUrl ? `<a href="${e.pdbUrl}" download class="text-amber-400 hover:underline">⬇ PDB</a>` : ""}
                        ${e.cifUrl ? `<a href="${e.cifUrl}" download class="text-sky-400 hover:underline">⬇ mmCIF</a>` : ""}
                        ${e.paeDocUrl ? `<a href="${e.paeDocUrl}" target="_blank" class="text-purple-400 hover:underline">PAE JSON ↗</a>` : ""}
                      </div>
                    </div>
                  `;
                }).join("")}
              </div>
              <div class="text-[9px] text-emerald-400 font-medium pt-0.5 flex justify-between items-center">
                <span>✅ All models validated in AFDB</span>
                <span class="font-mono text-slate-400">Ready for FoldMason!</span>
              </div>
            </div>
          `;
        } else {
          statusDiv.innerHTML = `
            <div class="space-y-1 text-slate-300">
              <div class="flex items-center justify-between text-sky-400 font-bold">
                <span>🧬 AlphaFold DB Configured</span>
                <a href="https://alphafold.ebi.ac.uk" target="_blank" class="text-sky-400 hover:underline font-mono text-[9px]">Portal ↗</a>
              </div>
              <div class="text-[9.5px]">Target IDs: <strong class="font-mono text-sky-300">${rawInput}</strong></div>
              <div class="text-[9px] text-[var(--text-muted)]">Run the generated CLI command below to download structures directly from EBI AlphaFold DB.</div>
            </div>
          `;
        }
      } catch (err) {
        statusDiv.innerHTML = `
          <div class="space-y-1 text-slate-300">
            <div class="flex items-center justify-between text-sky-400 font-bold">
              <span>🧬 AlphaFold DB Configured</span>
              <a href="https://alphafold.ebi.ac.uk" target="_blank" class="text-sky-400 hover:underline font-mono text-[9px]">Portal ↗</a>
            </div>
            <div class="text-[9.5px]">Target IDs: <strong class="font-mono text-sky-300">${rawInput}</strong></div>
            <div class="text-[9px] text-[var(--text-muted)]">Run the generated CLI command below to download structures directly from EBI AlphaFold DB.</div>
          </div>
        `;
      }
    }

    let pipelinePollTimer = null;


    async function runPipelineInBrowser() {
      const consoleWrap = document.getElementById("pipeConsoleWrapper");
      const consoleOut = document.getElementById("pipeConsoleOutput");
      const statusLbl = document.getElementById("pipeConsoleStatus");

      if (consoleWrap) consoleWrap.classList.remove("hidden");
      if (statusLbl) statusLbl.textContent = "Connecting to Local Bridge...";

      try {
        const isLocal = window.location.protocol === "http:" && (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1");
        const bridgeUrl = isLocal ? "/api/run_pipeline" : "http://localhost:8000/api/run_pipeline";

        const payload = {
          source: pipelineStudioState.source,
          qualifier: document.getElementById("pipeQualifierInput")?.value?.trim() || "glycoprotein",
          count: parseInt(document.getElementById("pipeCountSlider")?.value || 50),
          input_folder: document.getElementById("pipeLocalDirInput")?.value?.trim() || "300_rdrp/structures",
          metadata: document.getElementById("pipeLocalMetaInput")?.value?.trim() || "",
          output_dir: document.getElementById("pipeOutputDirInput")?.value?.trim() || "viral_custom_workflow",
          aligner: document.getElementById("pipeAlignerSelect")?.value || "foldmason",
          tree_type: document.getElementById("pipeTreeTypeSelect")?.value || "both",
          matrix: document.getElementById("pipeMatrixSelect")?.value || "alphafold",
          bootstrap: document.getElementById("pipeBootstrapToggle")?.checked ?? true,
          threads: document.getElementById("pipeThreadsSelect")?.value || "AUTO",
          embed: document.getElementById("pipeEmbedToggle")?.checked ?? false,
          embed_model: document.getElementById("pipeEmbedModelSelect")?.value || "esm2",
          embed_clustering: document.getElementById("pipeEmbedClusteringSelect")?.value || "upgma"
        };

        const res = await fetch(bridgeUrl, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(payload)
        });

        if (res.ok) {
          if (statusLbl) statusLbl.textContent = "Pipeline Running in Background...";
          if (consoleOut) {
            consoleOut.innerHTML = `
              <div class="text-emerald-400 font-bold">🚀 Job started successfully!</div>
              <div class="text-slate-400">${payload.source === 'viro3d' ? ('Target: Viro3D ' + payload.qualifier + ' (' + payload.count + ' seqs)') : ('Folder: ' + payload.input_folder)}</div>
              <div class="text-slate-500 pt-1">// Streaming live execution logs...</div>
            `;
          }
          startPipelinePolling();
        } else {
          const errData = await res.json().catch(() => ({}));
          if (consoleOut) {
            consoleOut.innerHTML = `
              <div class="text-amber-400 font-bold">⚠️ Bridge Notice: ${errData.error || "Could not launch job"}</div>
              <div class="text-slate-400 pt-1">You can run this directly in your terminal:</div>
              <div class="text-emerald-300 font-mono mt-1 p-1 bg-black/40 rounded">${document.getElementById("pipeCommandDisplay")?.textContent || ""}</div>
            `;
          }
        }
      } catch (e) {
        if (statusLbl) statusLbl.textContent = "Local Server Bridge Inactive";
        if (consoleOut) {
          consoleOut.innerHTML = `
            <div class="text-sky-300 font-bold">💡 How to Run Directly from this Page:</div>
            <div class="text-slate-300 text-[9px] leading-relaxed pt-1">
              To trigger background execution directly from the web browser, start the local execution server:
            </div>
            <div class="text-amber-300 font-mono text-[9px] bg-black/60 p-1.5 rounded my-1 select-all">
              python scripts/serve_interactive.py
            </div>
            <div class="text-slate-400 text-[8.5px]">
              Or click <strong>"Copy Command"</strong> above to run in your current terminal.
            </div>
          `;
        }
      }
    }

    function startPipelinePolling() {
      if (pipelinePollTimer) clearInterval(pipelinePollTimer);
      const isLocal = window.location.protocol === "http:" && (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1");
      const statusUrl = isLocal ? "/api/status" : "http://localhost:8000/api/status";

      pipelinePollTimer = setInterval(async () => {
        try {
          const res = await fetch(statusUrl);
          if (!res.ok) return;
          const data = await res.json();
          const consoleOut = document.getElementById("pipeConsoleOutput");
          const statusLbl = document.getElementById("pipeConsoleStatus");
          const spinner = document.getElementById("pipeRunningSpinner");

          if (consoleOut && data.logs && data.logs.length > 0) {
            consoleOut.innerHTML = data.logs.map(line => {
              let clr = "text-slate-300";
              if (line.includes("Error") || line.includes("ERROR")) clr = "text-rose-400 font-bold";
              else if (line.includes("Successfully") || line.includes("complete") || line.includes("Saved")) clr = "text-emerald-400";
              else if (line.includes("[IQ-TREE") || line.includes("[FoldMason")) clr = "text-sky-300";
              else if (line.includes("[Viro3D]")) clr = "text-amber-300";
              return `<div class="${clr}">${line.replace(/</g, "&lt;")}</div>`;
            }).join("");
            consoleOut.scrollTop = consoleOut.scrollHeight;
          }

          if (!data.is_running) {
            clearInterval(pipelinePollTimer);
            pipelinePollTimer = null;
            if (spinner) spinner.className = "hidden";
            if (statusLbl) {
              statusLbl.textContent = (data.returncode === 0) ? "✅ Pipeline Complete!" : `Exited with status ${data.returncode}`;
            }
          }
        } catch (_) {}
      }, 1200);
    }

    function clearPipelineConsole() {
      const consoleOut = document.getElementById("pipeConsoleOutput");
      if (consoleOut) consoleOut.innerHTML = `<div class="text-slate-500">// Console cleared.</div>`;
    }

    // 5. SIDEBAR TAB NAVIGATION (5 Tabs: Display, Filter, Clades, Rooting, Pipeline Studio)
    function switchSidebarTab(tabName) {
      ["panelDisplay", "panelFilter", "panelClades", "panelRooting", "panelPipeline"].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.classList.add("hidden");
      });
      ["tabBtnDisplay", "tabBtnFilter", "tabBtnClades", "tabBtnRooting", "tabBtnPipeline"].forEach(id => {
        const btn = document.getElementById(id);
        if (btn) {
          btn.className = "flex-1 py-2 text-center font-medium text-[var(--text-muted)] hover:text-[var(--text-main)] border-b-2 border-transparent text-[11px] transition flex items-center justify-center space-x-0.5";
        }
      });

      if (tabName === "display") {
        document.getElementById("panelDisplay").classList.remove("hidden");
        document.getElementById("tabBtnDisplay").className = "flex-1 py-2 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-[11px] transition flex items-center justify-center space-x-0.5";
      } else if (tabName === "filter") {
        document.getElementById("panelFilter").classList.remove("hidden");
        document.getElementById("tabBtnFilter").className = "flex-1 py-2 text-center font-bold text-emerald-400 border-b-2 border-emerald-400 text-[11px] transition flex items-center justify-center space-x-0.5";
        updateFilterUI();
      } else if (tabName === "clades") {
        document.getElementById("panelClades").classList.remove("hidden");
        document.getElementById("tabBtnClades").className = "flex-1 py-2 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-[11px] transition flex items-center justify-center space-x-0.5";
        updateCladeManagementUI();
      } else if (tabName === "rooting") {
        document.getElementById("panelRooting").classList.remove("hidden");
        document.getElementById("tabBtnRooting").className = "flex-1 py-2 text-center font-bold text-sky-400 border-b-2 border-sky-400 text-[11px] transition flex items-center justify-center space-x-0.5";
      } else if (tabName === "pipeline") {
        document.getElementById("panelPipeline").classList.remove("hidden");
        document.getElementById("tabBtnPipeline").className = "flex-1 py-2 text-center font-bold text-amber-400 border-b-2 border-amber-400 text-[11px] transition flex items-center justify-center space-x-0.5";
        updatePipelineCommandPreview();
      }
    }

    // 6. THEME TOGGLING & HIGH-CONTRAST THEMES
    const THEMES_META = {
      custom: { name: "Custom", icon: "✨", isDark: false },
      dark: { name: "Midnight", icon: "🌙", isDark: true },
      obsidian: { name: "Obsidian", icon: "🖤", isDark: true },
      forest: { name: "Forest", icon: "🌲", isDark: true },
      dracula: { name: "Dracula", icon: "🧛", isDark: true },
      cyberpunk: { name: "Cyberpunk", icon: "⚡", isDark: true },
      steel: { name: "Steel", icon: "🛡️", isDark: true },
      espresso: { name: "Espresso", icon: "☕", isDark: true },
      light: { name: "Pure White", icon: "☀️", isDark: false },
      publication: { name: "Pure White", icon: "☀️", isDark: false },
      solarized: { name: "Solarized", icon: "📜", isDark: false },
      nordic: { name: "Nordic", icon: "❄️", isDark: false },
      parchment: { name: "Parchment", icon: "📖", isDark: false },
      mint: { name: "Mint", icon: "🌿", isDark: false },
      high_contrast: { name: "Monochrome", icon: "👁️", isDark: false }
    };

    function isDarkTheme(t = settings.theme) {
      if (t === "custom") {
        return typeof customThemeState !== "undefined" ? !customThemeState.isLight : false;
      }
      const lightThemes = new Set(["light", "publication", "solarized", "nordic", "parchment", "mint", "high_contrast"]);
      return !lightThemes.has(t);
    }
    // ==========================================
    // CUSTOM THEME STUDIO ENGINE & PRESETS
    // ==========================================
    let customThemeState = {
      name: "Custom Light",
      isLight: true,
      vars: {
        "--bg-main": "#ffffff",
        "--panel-bg": "rgba(248, 250, 252, 0.98)",
        "--card-bg": "#ffffff",
        "--chip-bg": "#f8fafc",
        "--border-color": "#e2e8f0",
        "--text-main": "#0f172a",
        "--text-muted": "#64748b",
        "--branch-stroke": "#1e293b",
        "--accent": "#0284c7"
      }
    };

    let tempTheme = JSON.parse(JSON.stringify(customThemeState));

    const THEME_PARAM_DEFS = [
      { key: "--bg-main", label: "Workspace Background", desc: "Tree canvas background" },
      { key: "--panel-bg", label: "Panel & Toolbar", desc: "Top header and toolbar backgrounds" },
      { key: "--card-bg", label: "Card & Drawer Surface", desc: "Sidebar cards, modals, and MSA" },
      { key: "--border-color", label: "Border & Dividers", desc: "Subtle separation lines" },
      { key: "--text-main", label: "Primary Text Ink", desc: "Main typography and titles" },
      { key: "--text-muted", label: "Muted Subtext", desc: "Sub-labels and metadata text" },
      { key: "--branch-stroke", label: "Phylogenetic Branches", desc: "Branch stroke color" },
      { key: "--accent", label: "Accent & Highlights", desc: "Selected nodes and highlights" }
    ];

    const THEME_PRESETS = [
      {
        name: "Linear Crisp Light",
        icon: "💎",
        isLight: true,
        vars: {
          "--bg-main": "#ffffff",
          "--panel-bg": "#f8fafc",
          "--card-bg": "#ffffff",
          "--chip-bg": "#f8fafc",
          "--border-color": "#e2e8f0",
          "--text-main": "#0f172a",
          "--text-muted": "#64748b",
          "--branch-stroke": "#1e293b",
          "--accent": "#0284c7"
        }
      },
      {
        name: "GitHub Minimal Light",
        icon: "🌿",
        isLight: true,
        vars: {
          "--bg-main": "#ffffff",
          "--panel-bg": "#f6f8fa",
          "--card-bg": "#ffffff",
          "--chip-bg": "#f6f8fa",
          "--border-color": "#d0d7de",
          "--text-main": "#1f2328",
          "--text-muted": "#656d76",
          "--branch-stroke": "#24292f",
          "--accent": "#0969da"
        }
      },
      {
        name: "Editorial Parchment",
        icon: "📖",
        isLight: true,
        vars: {
          "--bg-main": "#fbf9f4",
          "--panel-bg": "#f5f1e8",
          "--card-bg": "#ffffff",
          "--chip-bg": "#fcfbf9",
          "--border-color": "#e7dfcf",
          "--text-main": "#1c1917",
          "--text-muted": "#57534e",
          "--branch-stroke": "#292524",
          "--accent": "#b45309"
        }
      },
      {
        name: "Cafe Latte Cream",
        icon: "☕",
        isLight: true,
        vars: {
          "--bg-main": "#fbf7f4",
          "--panel-bg": "#f3eae3",
          "--card-bg": "#ffffff",
          "--chip-bg": "#f7f1ec",
          "--border-color": "#e6d7cc",
          "--text-main": "#2b1810",
          "--text-muted": "#785a4a",
          "--branch-stroke": "#45281b",
          "--accent": "#b45309"
        }
      },
      {
        name: "Pastel Rose Blossom",
        icon: "🌸",
        isLight: true,
        vars: {
          "--bg-main": "#fff5f7",
          "--panel-bg": "#fde8ed",
          "--card-bg": "#ffffff",
          "--chip-bg": "#fdf2f4",
          "--border-color": "#fecdd3",
          "--text-main": "#4c0519",
          "--text-muted": "#9f1239",
          "--branch-stroke": "#881337",
          "--accent": "#e11d48"
        }
      },
      {
        name: "Cyber Neon Glow",
        icon: "🌌",
        isLight: false,
        vars: {
          "--bg-main": "#0d0221",
          "--panel-bg": "#150535",
          "--card-bg": "#1f0b4a",
          "--chip-bg": "#2d1b69",
          "--border-color": "#3b2069",
          "--text-main": "#fef08a",
          "--text-muted": "#f472b6",
          "--branch-stroke": "#06b6d4",
          "--accent": "#f43f5e"
        }
      },
      {
        name: "Obsidian Velvet",
        icon: "🖤",
        isLight: false,
        vars: {
          "--bg-main": "#030712",
          "--panel-bg": "#0a0f1d",
          "--card-bg": "#111827",
          "--chip-bg": "#1f2937",
          "--border-color": "#374151",
          "--text-main": "#f9fafb",
          "--text-muted": "#9ca3af",
          "--branch-stroke": "#38bdf8",
          "--accent": "#2dd4bf"
        }
      },
      {
        name: "Monochrome AAA",
        icon: "👁️",
        isLight: true,
        vars: {
          "--bg-main": "#ffffff",
          "--panel-bg": "#ffffff",
          "--card-bg": "#ffffff",
          "--chip-bg": "#f4f4f5",
          "--border-color": "#71717a",
          "--text-main": "#000000",
          "--text-muted": "#27272a",
          "--branch-stroke": "#000000",
          "--accent": "#000000"
        }
      }
    ];

    function initCustomThemeFromStorage() {
      try {
        const stored = localStorage.getItem("phylo_custom_theme");
        if (stored) {
          const parsed = JSON.parse(stored);
          if (parsed && parsed.vars) {
            customThemeState = Object.assign(customThemeState, parsed);
            tempTheme = JSON.parse(JSON.stringify(customThemeState));
          }
        }
      } catch (e) {
        console.warn("Could not load custom theme from localStorage", e);
      }
      updateThemeCustomBadgeInMenu();
    }

    function updateThemeCustomBadgeInMenu() {
      const lbl = document.getElementById("themeCustomNameLabel");
      if (lbl && customThemeState.name) lbl.textContent = customThemeState.name;
    }

    function openThemeModal() {
      tempTheme = JSON.parse(JSON.stringify(customThemeState));
      renderThemePresets();
      renderThemeInputs();
      updateThemeModalPreview();
      const modal = document.getElementById("themeModal");
      if (modal) modal.classList.remove("hidden");
      const menu = document.getElementById("themeMenu");
      if (menu) menu.classList.add("hidden");
    }

    function closeThemeModal() {
      const modal = document.getElementById("themeModal");
      if (modal) modal.classList.add("hidden");
      // If active theme is not custom, restore current theme styling
      if (settings.theme !== "custom") {
        clearCustomCssVariables();
      } else {
        applyCustomThemeProperties();
      }
    }

    function setCustomThemeFoundation(isLight) {
      tempTheme.isLight = isLight;
      const btnL = document.getElementById("btnThemeLightMode");
      const btnD = document.getElementById("btnThemeDarkMode");
      if (btnL && btnD) {
        if (isLight) {
          btnL.className = "px-3 py-1 rounded font-bold transition cursor-pointer bg-sky-500 text-white shadow-sm";
          btnD.className = "px-3 py-1 rounded font-medium transition cursor-pointer text-[var(--text-muted)] hover:text-[var(--text-main)]";
        } else {
          btnD.className = "px-3 py-1 rounded font-bold transition cursor-pointer bg-purple-500 text-white shadow-sm";
          btnL.className = "px-3 py-1 rounded font-medium transition cursor-pointer text-[var(--text-muted)] hover:text-[var(--text-main)]";
        }
      }
      updateThemeModalPreview();
    }

    function renderThemePresets() {
      const container = document.getElementById("themePresetsGrid");
      if (!container) return;
      container.innerHTML = "";

      THEME_PRESETS.forEach(p => {
        const btn = document.createElement("button");
        btn.className = "p-1.5 rounded-lg border border-[var(--border-color)] bg-[var(--card-bg)] hover:bg-slate-500/15 text-left text-[10.5px] transition cursor-pointer flex items-center space-x-1.5 shadow-sm";
        btn.onclick = () => loadThemePreset(p);
        btn.innerHTML = `
          <span>${p.icon}</span>
          <span class="truncate font-medium text-[var(--text-main)]">${p.name}</span>
        `;
        container.appendChild(btn);
      });
    }

    function loadThemePreset(preset) {
      tempTheme.name = preset.name;
      tempTheme.isLight = preset.isLight;
      tempTheme.vars = Object.assign({}, preset.vars);
      setCustomThemeFoundation(preset.isLight);
      renderThemeInputs();
      updateThemeModalPreview();
      // Live reflect on document
      previewCustomCssVariables(tempTheme);
    }

    function renderThemeInputs() {
      const container = document.getElementById("themeInputsGrid");
      if (!container) return;
      container.innerHTML = "";

      THEME_PARAM_DEFS.forEach(param => {
        const val = tempTheme.vars[param.key] || "#ffffff";
        // Convert rgba or complex color to hex for input[type=color]
        let hexVal = val;
        if (!hexVal.startsWith("#")) {
          hexVal = tempTheme.isLight ? "#ffffff" : "#0b1120";
        } else if (hexVal.length === 4) {
          hexVal = "#" + hexVal[1] + hexVal[1] + hexVal[2] + hexVal[2] + hexVal[3] + hexVal[3];
        }

        const row = document.createElement("div");
        row.className = "flex items-center justify-between p-2 rounded-lg bg-[var(--input-bg)] border border-[var(--border-color)] space-x-2";
        row.innerHTML = `
          <div class="min-w-0 flex-1">
            <span class="text-[11px] font-semibold text-[var(--text-main)] block truncate">${param.label}</span>
            <span class="text-[9px] text-[var(--text-muted)] block truncate">${param.desc}</span>
          </div>
          <div class="flex items-center space-x-1.5 shrink-0">
            <input type="color" value="${hexVal}" class="w-7 h-7 rounded border border-[var(--border-color)] bg-transparent cursor-pointer" oninput="onThemeColorWheelChange('${param.key}', this.value)">
            <input type="text" value="${val}" class="w-16 px-1.5 py-1 text-[10px] font-mono text-[var(--text-main)] bg-[var(--card-bg)] border border-[var(--border-color)] rounded focus:outline-none focus:border-sky-400 uppercase" onchange="onThemeHexChange('${param.key}', this.value)">
          </div>
        `;
        container.appendChild(row);
      });

      // Update JSON textarea
      const ta = document.getElementById("themeJsonTextarea");
      if (ta) ta.value = JSON.stringify(tempTheme, null, 2);
    }

    function onThemeColorWheelChange(varKey, hexValue) {
      tempTheme.vars[varKey] = hexValue;
      renderThemeInputs();
      updateThemeModalPreview();
      previewCustomCssVariables(tempTheme);
    }

    function onThemeHexChange(varKey, hexValue) {
      let val = hexValue.trim();
      if (!val.startsWith("#") && /^[0-9a-fA-F]{3,6}$/.test(val)) val = "#" + val;
      tempTheme.vars[varKey] = val;
      renderThemeInputs();
      updateThemeModalPreview();
      previewCustomCssVariables(tempTheme);
    }

    function updateThemeModalPreview() {
      const pBox = document.getElementById("themePreviewBox");
      if (!pBox) return;
      pBox.style.backgroundColor = tempTheme.vars["--bg-main"] || "#ffffff";
      pBox.style.borderColor = tempTheme.vars["--border-color"] || "#e2e8f0";
    }

    function previewCustomCssVariables(theme) {
      const isL = theme.isLight;
      for (const [k, v] of Object.entries(theme.vars)) {
        document.documentElement.style.setProperty(k, v);
      }
      document.documentElement.style.setProperty("--tip-label", theme.vars["--text-main"]);
      document.documentElement.style.setProperty("--input-bg", theme.vars["--card-bg"]);
      document.documentElement.style.setProperty("--grid-line", isL ? "rgba(0, 0, 0, 0.03)" : "rgba(51, 65, 85, 0.25)");

      // High-contrast adaptive badges
      document.documentElement.style.setProperty("--badge-sky-text", isL ? "#0369a1" : "#38bdf8");
      document.documentElement.style.setProperty("--badge-sky-bg", isL ? "rgba(3, 105, 161, 0.1)" : "rgba(56, 189, 248, 0.15)");
      document.documentElement.style.setProperty("--badge-sky-border", isL ? "rgba(3, 105, 161, 0.3)" : "rgba(56, 189, 248, 0.35)");

      document.documentElement.style.setProperty("--badge-emerald-text", isL ? "#047857" : "#34d399");
      document.documentElement.style.setProperty("--badge-emerald-bg", isL ? "rgba(4, 120, 87, 0.1)" : "rgba(16, 185, 129, 0.15)");
      document.documentElement.style.setProperty("--badge-emerald-border", isL ? "rgba(4, 120, 87, 0.3)" : "rgba(16, 185, 129, 0.35)");

      document.documentElement.style.setProperty("--badge-rose-text", isL ? "#be123c" : "#fb7185");
      document.documentElement.style.setProperty("--badge-rose-bg", isL ? "rgba(190, 18, 60, 0.1)" : "rgba(244, 63, 94, 0.15)");
      document.documentElement.style.setProperty("--badge-rose-border", isL ? "rgba(190, 18, 60, 0.3)" : "rgba(244, 63, 94, 0.35)");

      document.documentElement.style.setProperty("--badge-amber-text", isL ? "#b45309" : "#fbbf24");
      document.documentElement.style.setProperty("--badge-amber-bg", isL ? "rgba(180, 83, 9, 0.1)" : "rgba(245, 158, 11, 0.15)");
      document.documentElement.style.setProperty("--badge-amber-border", isL ? "rgba(180, 83, 9, 0.3)" : "rgba(245, 158, 11, 0.35)");

      document.documentElement.style.setProperty("--badge-purple-text", isL ? "#7e22ce" : "#c084fc");
      document.documentElement.style.setProperty("--badge-purple-bg", isL ? "rgba(126, 34, 206, 0.1)" : "rgba(168, 85, 247, 0.15)");
      document.documentElement.style.setProperty("--badge-purple-border", isL ? "rgba(126, 34, 206, 0.3)" : "rgba(168, 85, 247, 0.35)");

      renderTree();
      renderMsa();
    }

    function applyCustomThemeProperties() {
      previewCustomCssVariables(customThemeState);
    }

    function clearCustomCssVariables() {
      for (const k of Object.keys(customThemeState.vars)) {
        document.documentElement.style.removeProperty(k);
      }
      document.documentElement.style.removeProperty("--tip-label");
      document.documentElement.style.removeProperty("--input-bg");
      document.documentElement.style.removeProperty("--grid-line");
      document.documentElement.style.removeProperty("--badge-sky-text");
      document.documentElement.style.removeProperty("--badge-sky-bg");
      document.documentElement.style.removeProperty("--badge-sky-border");
      document.documentElement.style.removeProperty("--badge-emerald-text");
      document.documentElement.style.removeProperty("--badge-emerald-bg");
      document.documentElement.style.removeProperty("--badge-emerald-border");
      document.documentElement.style.removeProperty("--badge-rose-text");
      document.documentElement.style.removeProperty("--badge-rose-bg");
      document.documentElement.style.removeProperty("--badge-rose-border");
      document.documentElement.style.removeProperty("--badge-amber-text");
      document.documentElement.style.removeProperty("--badge-amber-bg");
      document.documentElement.style.removeProperty("--badge-amber-border");
      document.documentElement.style.removeProperty("--badge-purple-text");
      document.documentElement.style.removeProperty("--badge-purple-bg");
      document.documentElement.style.removeProperty("--badge-purple-border");
    }

    function applyCustomThemeStudio() {
      customThemeState = JSON.parse(JSON.stringify(tempTheme));
      localStorage.setItem("phylo_custom_theme", JSON.stringify(customThemeState));
      updateThemeCustomBadgeInMenu();
      setTheme("custom");
      closeThemeModal();
      showCladeToast(`Applied custom theme "${customThemeState.name}".`);
    }

    function resetCustomThemeToDefault() {
      localStorage.removeItem("phylo_custom_theme");
      tempTheme = JSON.parse(JSON.stringify(THEME_PRESETS[0]));
      setCustomThemeFoundation(tempTheme.isLight);
      renderThemeInputs();
      updateThemeModalPreview();
      previewCustomCssVariables(tempTheme);
    }

    function copyThemeJson() {
      const ta = document.getElementById("themeJsonTextarea");
      if (ta) {
        navigator.clipboard.writeText(ta.value).then(() => {
          showCladeToast("Theme JSON copied to clipboard!");
        });
      }
    }

    function importThemeJson() {
      const ta = document.getElementById("themeJsonTextarea");
      if (!ta) return;
      try {
        const parsed = JSON.parse(ta.value);
        if (parsed && parsed.vars) {
          tempTheme = Object.assign(tempTheme, parsed);
          setCustomThemeFoundation(Boolean(tempTheme.isLight));
          renderThemeInputs();
          updateThemeModalPreview();
          previewCustomCssVariables(tempTheme);
          showCladeToast("Loaded theme JSON!");
        } else {
          alert("Invalid theme JSON format. Must contain 'vars' object.");
        }
      } catch (e) {
        alert("Invalid JSON: " + e.message);
      }
    }


    function toggleThemeMenu(e) {
      if (e) e.stopPropagation();
      const menu = document.getElementById("themeMenu");
      if (menu) menu.classList.toggle("hidden");
    }

    document.addEventListener("click", (e) => {
      const wrapper = document.getElementById("themeMenuWrapper");
      if (wrapper && !wrapper.contains(e.target)) {
        const menu = document.getElementById("themeMenu");
        if (menu && !menu.classList.contains("hidden")) {
          menu.classList.add("hidden");
        }
      }
    });

    function setTheme(t) {
      settings.theme = t;
      document.documentElement.setAttribute("data-theme", t);
      localStorage.setItem("phylo_theme", t);

      if (t === "custom") {
        applyCustomThemeProperties();
        if (THEMES_META.custom) {
          THEMES_META.custom.name = customThemeState.name || "Custom";
          THEMES_META.custom.isDark = !customThemeState.isLight;
        }
      } else {
        clearCustomCssVariables();
      }

      const meta = THEMES_META[t] || { name: t, icon: "🎨" };
      const iconEl = document.getElementById("themeActiveIcon");
      const lblEl = document.getElementById("themeActiveLabel");
      if (iconEl) iconEl.textContent = meta.icon;
      if (lblEl) lblEl.textContent = meta.name;
      
      const menu = document.getElementById("themeMenu");
      if (menu) menu.classList.add("hidden");

      updateLegend();
      renderTree();
      updateMinimap();
      if (typeof sidebarViewer !== "undefined" && sidebarViewer && sidebarViewer.draw) {
        sidebarViewer.draw();
      }
      if (typeof headerLogoViewer !== "undefined" && headerLogoViewer && headerLogoViewer.draw) {
        headerLogoViewer.draw();
      }
      if (typeof renderSilhouetteSparkline === "function") {
        renderSilhouetteSparkline();
      }
      if (typeof renderMsa === "function") {
        if (typeof msaState !== "undefined") msaState._minimapCacheKey = null;
        renderMsa();
      }
    }

    function toggleTheme() {
      const current = document.documentElement.getAttribute("data-theme") || "dark";
      setTheme(isDarkTheme(current) ? "light" : "dark");
    }

    // 7. 3D PROTEIN VIEWER (C-ALPHA BACKBONE ENGINE)
    class ProteinViewer3D {
      constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) return;
        this.ctx = this.canvas.getContext("2d");
        this.atoms = [];
        this.currentTaxon = null;
        this.rotX = 0.25;
        this.rotY = 0.0;
        this.baseScale = 1.0;
        this.zoom = 1.0;
        this.autoRotate = true;
        this.animId = null;
        this.isDragging = false;
        this.lastMouse = { x: 0, y: 0 };
        this.initEvents();
      }

      initEvents() {
        this.canvas.addEventListener("mousedown", (e) => {
          this.isDragging = true;
          this.autoRotate = false;
          this.lastMouse = { x: e.clientX, y: e.clientY };
        });

        window.addEventListener("mousemove", (e) => {
          if (!this.isDragging) return;
          const dx = e.clientX - this.lastMouse.x;
          const dy = e.clientY - this.lastMouse.y;
          this.rotY += dx * 0.015;
          this.rotX += dy * 0.015;
          this.lastMouse = { x: e.clientX, y: e.clientY };
          this.draw();
        });

        window.addEventListener("mouseup", () => {
          if (this.isDragging) {
            this.isDragging = false;
            setTimeout(() => { this.autoRotate = true; }, 1800);
          }
        });

        this.canvas.addEventListener("wheel", (e) => {
          e.preventDefault();
          const factor = e.deltaY < 0 ? 1.1 : 0.9;
          this.zoom = Math.min(Math.max(0.4, this.zoom * factor), 4.0);
          this.draw();
        }, { passive: false });
      }

      loadStructure(taxonId) {
        const db = window.CA_STRUCTURES || {};
        if (!db) {
          this.atoms = [];
          this.draw();
          return;
        }

        // Resilient structure key lookup (exact, stripped, genbank, or stem prefix)
        let raw = db[taxonId];
        if (!raw) {
          const cleanId = taxonId.trim();
          raw = db[cleanId];
          if (!raw) {
            const meta = (typeof TAXA_METADATA !== 'undefined' && TAXA_METADATA[taxonId]) ? TAXA_METADATA[taxonId] : {};
            if (meta.genbank && db[meta.genbank]) {
              raw = db[meta.genbank];
            } else {
              const stem = taxonId.split("_")[0];
              const matchKey = Object.keys(db).find(k => k === stem || k.startsWith(stem));
              if (matchKey) raw = db[matchKey];
            }
          }
        }

        if (!raw || raw.length === 0) {
          this.atoms = [];
          this.draw();
          return;
        }
        this.currentTaxon = taxonId;

        // Retina High-DPI Canvas Buffer Setup
        const dpr = Math.max(1, window.devicePixelRatio || 1);
        const cw = this.canvas.clientWidth || parseInt(this.canvas.getAttribute("width")) || 280;
        const ch = this.canvas.clientHeight || parseInt(this.canvas.getAttribute("height")) || 176;

        this.canvas.width = Math.round(cw * dpr);
        this.canvas.height = Math.round(ch * dpr);

        // Compute Centroid & Normalize Coordinates
        let cx = 0, cy = 0, cz = 0;
        for (let i = 0; i < raw.length; i++) {
          cx += raw[i][0];
          cy += raw[i][1];
          cz += raw[i][2];
        }
        cx /= raw.length;
        cy /= raw.length;
        cz /= raw.length;

        let maxR = 0.001;
        this.atoms = raw.map(a => {
          const x = a[0] - cx;
          const y = a[1] - cy;
          const z = a[2] - cz;
          const d = Math.sqrt(x * x + y * y + z * z);
          if (d > maxR) maxR = d;
          return {
            x, y, z,
            plddt: a[3],
            resnum: a[4],
            resname: a[5] || "GLY",
            color: this.getPlddtColor(a[3])
          };
        });

        const minDim = Math.min(this.canvas.width, this.canvas.height) / 2.0;
        this.baseScale = (minDim * 0.72) / maxR;
        this.zoom = 1.0;
        this.draw();
        this.startAnimation();
      }

      getPlddtColor(val) {
        if (val >= 90) return "#2563eb";
        if (val >= 70) return "#38bdf8";
        if (val >= 50) return "#facc15";
        return "#f97316";
      }

      startAnimation() {
        if (this.animId) cancelAnimationFrame(this.animId);
        const loop = () => {
          if (this.autoRotate) {
            this.rotY += 0.006;
            this.draw();
          }
          this.animId = requestAnimationFrame(loop);
        };
        this.animId = requestAnimationFrame(loop);
      }

      stopAnimation() {
        if (this.animId) cancelAnimationFrame(this.animId);
        this.animId = null;
      }

      draw() {
        const w = this.canvas.width;
        const h = this.canvas.height;
        const ctx = this.ctx;
        ctx.clearRect(0, 0, w, h);

        const isDark = isDarkTheme();
        ctx.fillStyle = isDark ? "#020617" : (settings.theme === "solarized" ? "#fdf6e3" : (settings.theme === "nordic" ? "#eceff4" : "#f8fafc"));
        ctx.fillRect(0, 0, w, h);

        if (!this.atoms || this.atoms.length === 0) {
          ctx.fillStyle = isDark ? "#64748b" : "#94a3b8";
          ctx.font = "20px system-ui";
          ctx.textAlign = "center";
          ctx.fillText("No structure loaded", w / 2, h / 2);
          return;
        }

        const cosX = Math.cos(this.rotX), sinX = Math.sin(this.rotX);
        const cosY = Math.cos(this.rotY), sinY = Math.sin(this.rotY);
        const scale = this.baseScale * this.zoom;
        const centerX = w / 2;
        const centerY = h / 2;

        const projected = this.atoms.map(a => {
          const x1 = a.x * cosY + a.z * sinY;
          const z1 = -a.x * sinY + a.z * cosY;
          const y2 = a.y * cosX - z1 * sinX;
          const z2 = a.y * sinX + z1 * cosX;
          return {
            px: centerX + x1 * scale,
            py: centerY + y2 * scale,
            pz: z2,
            color: a.color
          };
        });

        // Draw C-alpha Backbone Tube Segments
        for (let i = 0; i < projected.length - 1; i++) {
          const p1 = projected[i];
          const p2 = projected[i + 1];

          ctx.beginPath();
          ctx.moveTo(p1.px, p1.py);
          ctx.lineTo(p2.px, p2.py);
          ctx.strokeStyle = p1.color;
          ctx.lineWidth = 3.6;
          ctx.lineCap = "round";
          ctx.stroke();
        }

        // Draw atom nodes sorted by depth
        projected.sort((a, b) => a.pz - b.pz);
        for (let i = 0; i < projected.length; i++) {
          const p = projected[i];
          ctx.beginPath();
          ctx.arc(p.px, p.py, 2.6, 0, 2 * Math.PI);
          ctx.fillStyle = p.color;
          ctx.fill();
        }
      }
    }

    const tooltipViewer = new ProteinViewer3D("tooltipCanvas");
    const sidebarViewer = new ProteinViewer3D("sidebarCanvas");

    // 7b. MINI HEADER LOGO VIEWER (ROTATING C-ALPHA BACKBONE)
    class MiniHeaderLogoViewer {
      constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) return;
        this.ctx = this.canvas.getContext("2d");
        this.rotY = 0.0;
        this.rotX = 0.35;
        this.atoms = [];
        this.animId = null;
        this.init();
      }

      init() {
        const dpr = Math.max(1, window.devicePixelRatio || 1);
        this.dpr = dpr;
        this.canvas.width = Math.round(20 * dpr);
        this.canvas.height = Math.round(20 * dpr);

        const db = window.CA_STRUCTURES || {};
        const sampleTaxon = Object.keys(db)[0];
        let raw = sampleTaxon ? db[sampleTaxon] : null;
        if (raw && raw.length > 25) {
          const start = Math.floor(raw.length * 0.2);
          raw = raw.slice(start, start + 35);
        } else {
          raw = [];
          for (let i = 0; i < 30; i++) {
            const angle = i * 1.7;
            const r = i < 16 ? 4.8 : 3.2 + (i - 16) * 0.4;
            const x = Math.cos(angle) * r;
            const y = (i - 15) * 1.5;
            const z = Math.sin(angle) * r;
            const plddt = 80 + Math.sin(i * 0.4) * 16;
            raw.push([x, y, z, plddt]);
          }
        }

        let cx = 0, cy = 0, cz = 0;
        for (let i = 0; i < raw.length; i++) {
          cx += raw[i][0]; cy += raw[i][1]; cz += raw[i][2];
        }
        cx /= raw.length; cy /= raw.length; cz /= raw.length;

        let maxR = 0.001;
        this.atoms = raw.map(a => {
          const x = a[0] - cx, y = a[1] - cy, z = a[2] - cz;
          const d = Math.sqrt(x*x + y*y + z*z);
          if (d > maxR) maxR = d;
          const plddt = a[3] !== undefined ? a[3] : 88;
          let colorDark = "#38bdf8", colorLight = "#0284c7";
          if (plddt >= 90) { colorDark = "#38bdf8"; colorLight = "#0284c7"; }
          else if (plddt >= 75) { colorDark = "#34d399"; colorLight = "#059669"; }
          else if (plddt >= 60) { colorDark = "#fbbf24"; colorLight = "#d97706"; }
          else { colorDark = "#f43f5e"; colorLight = "#e11d48"; }
          return { x, y, z, colorDark, colorLight };
        });

        this.baseScale = (this.canvas.width * 0.38) / maxR;
        this.startLoop();
      }

      startLoop() {
        const loop = () => {
          this.rotY += 0.006;
          this.draw();
          this.animId = requestAnimationFrame(loop);
        };
        this.animId = requestAnimationFrame(loop);
      }

      draw() {
        const w = this.canvas.width;
        const h = this.canvas.height;
        const ctx = this.ctx;
        ctx.clearRect(0, 0, w, h);

        const isDark = typeof isDarkTheme === "function" ? isDarkTheme() : true;
        const cosX = Math.cos(this.rotX), sinX = Math.sin(this.rotX);
        const cosY = Math.cos(this.rotY), sinY = Math.sin(this.rotY);
        const scale = this.baseScale;
        const centerX = w / 2;
        const centerY = h / 2;

        const projected = this.atoms.map(a => {
          const x1 = a.x * cosY + a.z * sinY;
          const z1 = -a.x * sinY + a.z * cosY;
          const y2 = a.y * cosX - z1 * sinX;
          const z2 = a.y * sinX + z1 * cosX;
          return {
            px: centerX + x1 * scale,
            py: centerY + y2 * scale,
            pz: z2,
            color: isDark ? a.colorDark : a.colorLight
          };
        });

        for (let i = 0; i < projected.length - 1; i++) {
          const p1 = projected[i];
          const p2 = projected[i + 1];
          ctx.beginPath();
          ctx.moveTo(p1.px, p1.py);
          ctx.lineTo(p2.px, p2.py);
          ctx.strokeStyle = p1.color;
          ctx.lineWidth = 1.6 * (this.dpr || 1);
          ctx.lineCap = "round";
          ctx.stroke();
        }

        for (let i = 0; i < projected.length; i += 2) {
          const p = projected[i];
          ctx.beginPath();
          ctx.arc(p.px, p.py, 1.1 * (this.dpr || 1), 0, 2 * Math.PI);
          ctx.fillStyle = p.color;
          ctx.fill();
        }
      }
    }

    let headerLogoViewer = null;
    function initHeaderLogo() {
      if (!headerLogoViewer && document.getElementById("headerLogoCanvas")) {
        headerLogoViewer = new MiniHeaderLogoViewer("headerLogoCanvas");
      }
    }

    // 8. COLOR LOGIC & VIRIDIS INTERPOLATOR
    function interpolateViridis(t) {
      const stops = [
        [68, 1, 84],
        [59, 82, 139],
        [33, 145, 140],
        [94, 201, 98],
        [253, 231, 37]
      ];
      const p = Math.max(0, Math.min(t, 1)) * (stops.length - 1);
      const i = Math.floor(p);
      const frac = p - i;
      if (i >= stops.length - 1) {
        const [r, g, b] = stops[stops.length - 1];
        return `rgb(${r},${g},${b})`;
      }
      const [r1, g1, b1] = stops[i];
      const [r2, g2, b2] = stops[i + 1];
      const r = Math.round(r1 + (r2 - r1) * frac);
      const g = Math.round(g1 + (g2 - g1) * frac);
      const b = Math.round(b1 + (b2 - b1) * frac);
      return `rgb(${r},${g},${b})`;
    }

    // 8B. BRANCH BOOTSTRAP SUPPORT COLOR RESOLUTION
    function getBranchStrokeColor(node) {
      if (!settings.colorBranchesBySupport || !node) return "";
      
      let val = null;
      if (settings.supportMetric === "alrt") {
        if (node.alrt !== undefined && node.alrt !== null && !isNaN(node.alrt)) {
          val = node.alrt;
        } else if (node.support !== undefined && node.support !== null && !isNaN(node.support)) {
          val = node.support;
        }
      } else {
        if (node.ufboot !== undefined && node.ufboot !== null && !isNaN(node.ufboot)) {
          val = node.ufboot;
        } else if (node.support !== undefined && node.support !== null && !isNaN(node.support)) {
          val = node.support;
        }
      }

      if (val === null || isNaN(val)) return "";
      if (val <= 1.0 && val > 0.0) val = val * 100.0;

      if (settings.supportPalette === "traffic") {
        if (val >= 95) return "#10b981"; // Strong support (Emerald green)
        if (val >= 70) return "#f59e0b"; // Moderate support (Amber)
        return "#ef4444"; // Weak support (Rose red)
      } else if (settings.supportPalette === "viridis") {
        return interpolateViridis(Math.max(0, Math.min(100, val)) / 100.0);
      } else if (settings.supportPalette === "grayscale") {
        const alpha = Math.max(0.18, Math.min(1.0, val / 100.0));
        return `rgba(56, 189, 248, ${alpha.toFixed(2)})`;
      }
      return "";
    }

    // 8C. VERTICAL BRANCH GRADIENT RESOLUTION
    let verticalGradCounter = 0;
    function applyVerticalBranchGradient(vLine, node, minY, maxY, defs) {
      if (!settings.colorBranchesBySupport || !node || !defs) return;
      if (typeof minY !== 'number' || typeof maxY !== 'number' || isNaN(minY) || isNaN(maxY)) return;

      const span = maxY - minY;
      const realPoints = [];

      // 1. Collect all child horizontal branches with real support values/colors
      if (node.children && node.children.length > 0) {
        node.children.forEach(child => {
          const col = getBranchStrokeColor(child);
          if (col && typeof child.y === 'number' && !isNaN(child.y)) {
            realPoints.push({ y: child.y, color: col });
          }
        });
      }

      // 2. Collect incoming branch at node.y if it has a real support value/color
      const parentCol = getBranchStrokeColor(node);
      if (parentCol && typeof node.y === 'number' && !isNaN(node.y)) {
        realPoints.push({ y: node.y, color: parentCol });
      }

      if (realPoints.length === 0) return;

      // Deduplicate points with nearly identical Y coordinates (within 0.5px)
      realPoints.sort((a, b) => a.y - b.y);
      const uniquePoints = [];
      realPoints.forEach(pt => {
        if (uniquePoints.length === 0 || Math.abs(pt.y - uniquePoints[uniquePoints.length - 1].y) > 0.5) {
          uniquePoints.push(pt);
        }
      });

      if (uniquePoints.length === 1 || span <= 0.5) {
        vLine.style.stroke = uniquePoints[0].color;
        return;
      }

      const firstCol = uniquePoints[0].color;
      const allSame = uniquePoints.every(p => p.color === firstCol);
      if (allSame) {
        vLine.style.stroke = firstCol;
        return;
      }

      // 3. Construct vertical linearGradient with userSpaceOnUse for robust line rendering
      const gradId = "vgrad_" + (++verticalGradCounter);
      const grad = document.createElementNS("http://www.w3.org/2000/svg", "linearGradient");
      grad.setAttribute("id", gradId);
      grad.setAttribute("gradientUnits", "userSpaceOnUse");
      grad.setAttribute("x1", node.x);
      grad.setAttribute("y1", minY);
      grad.setAttribute("x2", node.x);
      grad.setAttribute("y2", maxY);

      // Pad top if first real value is below minY
      if (uniquePoints[0].y > minY + 0.1) {
        const stopStart = document.createElementNS("http://www.w3.org/2000/svg", "stop");
        stopStart.setAttribute("offset", "0%");
        stopStart.setAttribute("stop-color", uniquePoints[0].color);
        grad.appendChild(stopStart);
      }

      uniquePoints.forEach(pt => {
        const frac = Math.max(0, Math.min(1, (pt.y - minY) / span));
        const stop = document.createElementNS("http://www.w3.org/2000/svg", "stop");
        stop.setAttribute("offset", (frac * 100).toFixed(2) + "%");
        stop.setAttribute("stop-color", pt.color);
        grad.appendChild(stop);
      });

      // Pad bottom if last real value is above maxY
      const lastPt = uniquePoints[uniquePoints.length - 1];
      if (lastPt.y < maxY - 0.1) {
        const stopEnd = document.createElementNS("http://www.w3.org/2000/svg", "stop");
        stopEnd.setAttribute("offset", "100%");
        stopEnd.setAttribute("stop-color", lastPt.color);
        grad.appendChild(stopEnd);
      }

      defs.appendChild(grad);
      vLine.style.stroke = `url(#${gradId})`;
    }

    // =========================================================================
    // COLOR PALETTE STUDIO ENGINE (CUSTOM HEX LIST & INTERACTIVE COLOR WHEEL)
    // =========================================================================
    const DEFAULT_EXPANDED_PALETTE = [
      "#38bdf8", "#f97316", "#a855f7", "#10b981", "#ec4899",
      "#eab308", "#06b6d4", "#8b5cf6", "#f43f5e", "#14b8a6",
      "#6366f1", "#84cc16", "#e11d48", "#0284c7", "#ca8a04",
      "#d946ef", "#059669", "#b45309", "#4f46e5", "#f59e0b",
      "#22c55e", "#0ea5e9", "#d97706", "#9333ea", "#2dd4bf",
      "#fb7185", "#3b82f6", "#16a34a", "#c026d3", "#64748b",
      "#e2e8f0", "#78716c", "#a3e635", "#34d399", "#818cf8",
      "#c084fc", "#f472b6", "#fb923c", "#facc15", "#4ade80"
    ];

    const PALETTE_PRESETS = [
      { id: "botanical", name: "Botanical Earth", colors: ["#B9554E", "#627B08", "#267567", "#294719", "#72A183"] },
      { id: "ocean", name: "Ocean Depth", colors: ["#0284c7", "#06b6d4", "#14b8a6", "#10b981", "#6366f1"] },
      { id: "sunset", name: "Sunset Ember", colors: ["#e11d48", "#f43f5e", "#f97316", "#f59e0b", "#eab308"] },
      { id: "cyber", name: "Neon Synth", colors: ["#a855f7", "#ec4899", "#06b6d4", "#10b981", "#facc15"] },
      { id: "viridis", name: "Viridis Bio", colors: ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"] },
      { id: "ictv", name: "ICTV Classic", colors: ["#38bdf8", "#ec4899", "#10b981", "#f97316", "#8b5cf6", "#06b6d4", "#eab308", "#a855f7"] }
    ];

    let customPaletteState = {
      isActive: false,
      palette: ["#B9554E", "#627B08", "#267567", "#294719", "#72A183"],
      applyToCategorical: true,
      applyToContinuous: true,
      applyToClades: true
    };

    let tempPalette = [...customPaletteState.palette];

    function hexToRgb(hex) {
      let c = (hex || "#38bdf8").replace(/^#/, '');
      if (c.length === 3) {
        c = c.split('').map(x => x + x).join('');
      }
      const num = parseInt(c, 16);
      if (isNaN(num)) return [56, 189, 248];
      return [(num >> 16) & 255, (num >> 8) & 255, num & 255];
    }

    function interpolatePalette(t, paletteArray) {
      if (!paletteArray || paletteArray.length === 0) return "#38bdf8";
      if (paletteArray.length === 1) return paletteArray[0];
      const stops = paletteArray.map(hexToRgb);
      const p = Math.max(0, Math.min(t, 1)) * (stops.length - 1);
      const i = Math.floor(p);
      const frac = p - i;
      if (i >= stops.length - 1) {
        const [r, g, b] = stops[stops.length - 1];
        return `rgb(${r},${g},${b})`;
      }
      const [r1, g1, b1] = stops[i];
      const [r2, g2, b2] = stops[i + 1];
      const r = Math.round(r1 + (r2 - r1) * frac);
      const g = Math.round(g1 + (g2 - g1) * frac);
      const b = Math.round(b1 + (b2 - b1) * frac);
      return `rgb(${r},${g},${b})`;
    }

    function parseHexList(text) {
      const matches = (text || "").match(/#?[0-9a-fA-F]{6}\b|#?[0-9a-fA-F]{3}\b/g) || [];
      return matches.map(m => m.startsWith('#') ? m.toUpperCase() : ('#' + m.toUpperCase())).map(m => {
        if (m.length === 4) {
          return '#' + m[1] + m[1] + m[2] + m[2] + m[3] + m[3];
        }
        return m;
      });
    }

    function getCategoryColorFromCustomPalette(colDef, val) {
      if (!customPaletteState.isActive || !customPaletteState.palette || customPaletteState.palette.length === 0) {
        return (colDef.colors && colDef.colors[val]) ? colDef.colors[val] : "#94a3b8";
      }
      const vals = colDef.values || (colDef.colors ? Object.keys(colDef.colors) : []);
      const idx = vals.indexOf(val);
      if (idx >= 0) {
        return customPaletteState.palette[idx % customPaletteState.palette.length];
      }
      // Deterministic string hash fallback
      let hash = 0;
      for (let i = 0; i < val.length; i++) {
        hash = (hash << 5) - hash + val.charCodeAt(i);
        hash |= 0;
      }
      const pos = Math.abs(hash) % customPaletteState.palette.length;
      return customPaletteState.palette[pos];
    }

    function openPaletteModal() {
      tempPalette = [...customPaletteState.palette];
      const modal = document.getElementById("paletteModal");
      if (modal) modal.classList.remove("hidden");
      renderPalettePresets();
      renderPaletteSwatches();
      updatePalettePreviews();
    }

    function closePaletteModal() {
      const modal = document.getElementById("paletteModal");
      if (modal) modal.classList.add("hidden");
    }

    function renderPalettePresets() {
      const container = document.getElementById("palettePresetsContainer");
      if (!container) return;
      container.innerHTML = "";
      PALETTE_PRESETS.forEach(preset => {
        const btn = document.createElement("button");
        btn.className = "flex items-center space-x-1.5 px-2.5 py-1 rounded-lg border border-[var(--border-color)] bg-[var(--card-bg)] hover:bg-slate-500/15 text-[11px] text-[var(--text-main)] font-medium transition cursor-pointer shadow-sm";
        btn.title = `Preset: ${preset.name} (${preset.colors.length} colors)`;
        
        const swatchesHtml = preset.colors.map(c => `<span class="w-2 h-2 rounded-full inline-block" style="background-color: ${c};"></span>`).join('');
        btn.innerHTML = `<span>${preset.name}</span><span class="flex space-x-0.5 ml-1">${swatchesHtml}</span>`;
        btn.onclick = () => {
          tempPalette = [...preset.colors];
          renderPaletteSwatches();
          const hexInp = document.getElementById("paletteHexInput");
          if (hexInp) hexInp.value = tempPalette.join(", ");
        };
        container.appendChild(btn);
      });

      // Default 40-color preset button
      const defBtn = document.createElement("button");
      defBtn.className = "flex items-center space-x-1 px-2.5 py-1 rounded-lg border border-slate-600/40 bg-slate-700/30 hover:bg-slate-700/60 text-[11px] text-[var(--text-muted)] hover:text-white transition cursor-pointer";
      defBtn.innerHTML = `<span>🔄 Default 40-Color</span>`;
      defBtn.onclick = () => {
        tempPalette = [...DEFAULT_EXPANDED_PALETTE.slice(0, 10)];
        renderPaletteSwatches();
        const hexInp = document.getElementById("paletteHexInput");
        if (hexInp) hexInp.value = tempPalette.join(", ");
      };
      container.appendChild(defBtn);
    }

    function renderPaletteSwatches() {
      const listEl = document.getElementById("paletteSwatchesList");
      if (!listEl) return;
      listEl.innerHTML = "";

      tempPalette.forEach((hex, index) => {
        const card = document.createElement("div");
        card.className = "flex items-center space-x-2 p-1.5 rounded-lg bg-[var(--card-bg)] border border-[var(--border-color)] shadow-sm group";
        
        card.innerHTML = `
          <div class="relative w-7 h-7 rounded-md border border-white/20 shadow-inner shrink-0 cursor-pointer overflow-hidden flex items-center justify-center" style="background-color: ${hex};" title="Click to open Colour Wheel">
            <input type="color" value="${hex}" oninput="onColorWheelChange(${index}, this.value)" class="opacity-0 absolute inset-0 w-full h-full cursor-pointer">
            <span class="text-[9px] opacity-0 group-hover:opacity-100 font-bold drop-shadow text-white pointer-events-none">🎨</span>
          </div>
          <input type="text" value="${hex}" onchange="onHexTextChange(${index}, this.value)" class="w-20 bg-[var(--input-bg)] border border-[var(--border-color)] rounded px-1.5 py-0.5 text-[11px] font-mono text-[var(--text-main)] uppercase text-center focus:outline-none focus:border-purple-400">
          <div class="flex items-center space-x-0.5 ml-auto">
            ${index > 0 ? `<button onclick="movePaletteColor(${index}, -1)" class="text-[10px] px-1 py-0.5 text-[var(--text-muted)] hover:text-[var(--text-main)] transition cursor-pointer" title="Move Left">&larr;</button>` : ''}
            ${index < tempPalette.length - 1 ? `<button onclick="movePaletteColor(${index}, 1)" class="text-[10px] px-1 py-0.5 text-[var(--text-muted)] hover:text-[var(--text-main)] transition cursor-pointer" title="Move Right">&rarr;</button>` : ''}
            <button onclick="removePaletteColor(${index})" class="text-[12px] px-1 py-0.5 text-[var(--text-muted)] hover:text-rose-400 transition cursor-pointer ${tempPalette.length <= 2 ? 'opacity-30 cursor-not-allowed' : ''}" title="Remove Color">&times;</button>
          </div>
        `;
        listEl.appendChild(card);
      });

      // Update hex input string if not currently focused
      const hexInput = document.getElementById("paletteHexInput");
      if (hexInput && document.activeElement !== hexInput) {
        hexInput.value = tempPalette.join(", ");
      }

      const badge = document.getElementById("paletteColorCountBadge");
      if (badge) badge.textContent = `${tempPalette.length} Colors`;

      updatePalettePreviews();
    }

    function onColorWheelChange(index, value) {
      tempPalette[index] = value.toUpperCase();
      const hexInput = document.getElementById("paletteHexInput");
      if (hexInput) hexInput.value = tempPalette.join(", ");
      renderPaletteSwatches();
    }

    function onHexTextChange(index, value) {
      const parsed = parseHexList(value);
      if (parsed.length > 0) {
        tempPalette[index] = parsed[0];
      }
      const hexInput = document.getElementById("paletteHexInput");
      if (hexInput) hexInput.value = tempPalette.join(", ");
      renderPaletteSwatches();
    }

    function onPaletteHexInput(text) {
      const parsed = parseHexList(text);
      if (parsed.length >= 2) {
        tempPalette = parsed;
        renderPaletteSwatches();
      }
    }

    function addPaletteColor() {
      const defaultNextColors = ["#B9554E", "#627B08", "#267567", "#294719", "#72A183", "#38bdf8", "#f97316", "#a855f7", "#10b981", "#eab308"];
      const nextColor = defaultNextColors[tempPalette.length % defaultNextColors.length];
      tempPalette.push(nextColor);
      renderPaletteSwatches();
    }

    function removePaletteColor(index) {
      if (tempPalette.length <= 2) {
        alert("A palette must contain at least 2 colors.");
        return;
      }
      tempPalette.splice(index, 1);
      renderPaletteSwatches();
    }

    function movePaletteColor(index, delta) {
      const target = index + delta;
      if (target < 0 || target >= tempPalette.length) return;
      const tmp = tempPalette[index];
      tempPalette[index] = tempPalette[target];
      tempPalette[target] = tmp;
      renderPaletteSwatches();
    }

    function shufflePalette() {
      for (let i = tempPalette.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [tempPalette[i], tempPalette[j]] = [tempPalette[j], tempPalette[i]];
      }
      renderPaletteSwatches();
    }

    function reversePalette() {
      tempPalette.reverse();
      renderPaletteSwatches();
    }

    function copyPaletteHexList() {
      const text = tempPalette.join(", ");
      navigator.clipboard.writeText(text).then(() => {
        showCladeToast("Copied palette hex list to clipboard!");
      }).catch(() => {
        prompt("Copy palette hex list:", text);
      });
    }

    function updatePalettePreviews() {
      // Continuous preview
      const gradPreview = document.getElementById("paletteGradientPreview");
      if (gradPreview) {
        gradPreview.style.background = `linear-gradient(90deg, ${tempPalette.join(", ")})`;
      }

      // Categorical preview
      const catContainer = document.getElementById("paletteCategoriesPreview");
      const labelEl = document.getElementById("paletteCategoryColLabel");
      const colDef = getActiveColorColumnDef();

      if (labelEl) {
        labelEl.textContent = colDef ? colDef.label : "Active Categories";
      }

      if (catContainer) {
        catContainer.innerHTML = "";
        const sampleCategories = (colDef && colDef.values && colDef.values.length > 0)
          ? colDef.values.slice(0, 12)
          : ["Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta"];

        sampleCategories.forEach((cat, idx) => {
          const color = tempPalette[idx % tempPalette.length];
          const pill = document.createElement("span");
          pill.className = "inline-flex items-center space-x-1 px-2 py-0.5 rounded-full text-[10px] font-medium border border-white/10 shadow-sm";
          pill.style.backgroundColor = color + "26"; // 15% opacity tint
          pill.style.borderColor = color + "66";
          pill.style.color = color;
          pill.innerHTML = `<span class="w-1.5 h-1.5 rounded-full" style="background-color: ${color};"></span><span>${cat}</span>`;
          catContainer.appendChild(pill);
        });
      }
    }

    function applyPaletteStudio() {
      if (tempPalette.length < 2) {
        alert("Please provide at least 2 colors.");
        return;
      }
      customPaletteState.isActive = true;
      customPaletteState.palette = [...tempPalette];
      localStorage.setItem("phylo_custom_palette", JSON.stringify(customPaletteState));

      updatePaletteMiniStrip();
      updateLegend();
      renderTree();
      updateMinimap();
      if (typeof renderMsa === "function") {
        if (typeof msaState !== "undefined") msaState._minimapCacheKey = null;
        renderMsa();
      }
      updateCladeManagementUI();
      if (typeof renderMsa === "function") {
        if (typeof msaState !== "undefined") msaState._minimapCacheKey = null;
        renderMsa();
      }
      closePaletteModal();
      showCladeToast(`Applied custom palette (${customPaletteState.palette.length} colors).`);
    }

    function resetToDefaultPalette() {
      customPaletteState.isActive = false;
      localStorage.removeItem("phylo_custom_palette");
      tempPalette = ["#B9554E", "#627B08", "#267567", "#294719", "#72A183"];

      updatePaletteMiniStrip();
      updateLegend();
      renderTree();
      updateMinimap();
      updateCladeManagementUI();
      closePaletteModal();
      showCladeToast("Reset colors to default dataset palettes.");
    }

    function updatePaletteMiniStrip() {
      const stripEl = document.getElementById("paletteMiniStrip");
      const badgeEl = document.getElementById("paletteActiveBadge");
      if (!stripEl) return;
      stripEl.innerHTML = "";

      const activeColors = (customPaletteState.isActive && customPaletteState.palette.length > 0)
        ? customPaletteState.palette
        : DEFAULT_EXPANDED_PALETTE.slice(0, 8);

      activeColors.forEach(c => {
        const span = document.createElement("span");
        span.className = "flex-1 h-full";
        span.style.backgroundColor = c;
        stripEl.appendChild(span);
      });

      if (badgeEl) {
        if (customPaletteState.isActive) {
          badgeEl.textContent = `Custom (${customPaletteState.palette.length})`;
          badgeEl.className = "badge-purple text-[9px] font-mono px-1.5 py-0.5 rounded-full font-bold";
        } else {
          badgeEl.textContent = "Default";
          badgeEl.className = "text-[9px] font-mono px-1.5 py-0.5 rounded-full bg-slate-700/60 text-slate-400 border border-slate-600/40";
        }
      }
    }

    function getNodeColor(node) {
      if (node._collapsed) {
        const info = getCladeInfo(node);
        return info.cladeColor || "#38bdf8";
      }

      if (settings.colorColumn === "solid") {
        return (customPaletteState.isActive && customPaletteState.palette.length > 0) 
          ? customPaletteState.palette[0] 
          : "#38bdf8";
      }

      const colDef = getActiveColorColumnDef();
      if (!colDef) return "#38bdf8";

      const meta = TAXA_METADATA[node.name];
      if (!meta) return "#94a3b8";

      const val = meta[colDef.key];
      if (val === undefined || val === null || val === "") return "#94a3b8";

      if (colDef.type === "categorical") {
        if (customPaletteState.isActive && customPaletteState.applyToCategorical && customPaletteState.palette.length > 0) {
          return getCategoryColorFromCustomPalette(colDef, String(val));
        }
        return (colDef.colors && colDef.colors[val]) ? colDef.colors[val] : "#94a3b8";
      } else if (colDef.type === "continuous") {
        const num = parseFloat(val);
        if (isNaN(num)) return "#94a3b8";

        if (customPaletteState.isActive && customPaletteState.applyToContinuous && customPaletteState.palette.length > 0) {
          const min = colDef.min || 0;
          const max = colDef.max || 100;
          const t = Math.max(0, Math.min((num - min) / (max - min || 1), 1.0));
          return interpolatePalette(t, customPaletteState.palette);
        }

        if (colDef.key.toLowerCase().includes("plddt")) {
          if (num >= 90) return "#2563eb";
          if (num >= 70) return "#38bdf8";
          if (num >= 50) return "#facc15";
          return "#f97316";
        }

        const min = colDef.min || 0;
        const max = colDef.max || 100;
        const t = Math.max(0, Math.min((num - min) / (max - min || 1), 1.0));
        return interpolateViridis(t);
      }

      return "#38bdf8";
    }

    // 9. DYNAMIC METADATA SELECTORS & LEGEND
    function populateMetadataSelectors() {
      const colorSel = document.getElementById("colorColumnSelect");
      const cladeSel = document.getElementById("cladeColumnSelect");
      const tipLabelSel = document.getElementById("tipLabelColumnSelect");
      if (!colorSel || !cladeSel) return;

      const cols = getActiveColumns();
      const catCols = cols.filter(c => c.type === "categorical");

      colorSel.innerHTML = "";
      cols.forEach(c => {
        const opt = document.createElement("option");
        opt.value = c.key;
        const icon = (c.type === "categorical") ? "🏷️" : "📈";
        opt.textContent = `${icon} ${c.label}`;
        if (c.key === settings.colorColumn) opt.selected = true;
        colorSel.appendChild(opt);
      });
      const solidOpt = document.createElement("option");
      solidOpt.value = "solid";
      solidOpt.textContent = "⚪ Solid (Uniform Sky Blue)";
      if (settings.colorColumn === "solid") solidOpt.selected = true;
      colorSel.appendChild(solidOpt);

      cladeSel.innerHTML = "";
      catCols.forEach(c => {
        const opt = document.createElement("option");
        opt.value = c.key;
        opt.textContent = `🌿 ${c.label}`;
        if (c.key === settings.cladeGroupColumn) opt.selected = true;
        cladeSel.appendChild(opt);
      });

      if (tipLabelSel) {
        tipLabelSel.innerHTML = "";
        const availableKeys = new Set(["taxon_id", "id_and_color"]);

        // 1. Default Taxon ID
        const idOpt = document.createElement("option");
        idOpt.value = "taxon_id";
        idOpt.textContent = "🆔 Taxon ID / Accession";
        tipLabelSel.appendChild(idOpt);

        // 2. Active Columns (Categorical & Continuous)
        cols.forEach(c => {
          const opt = document.createElement("option");
          opt.value = c.key;
          const icon = (c.type === "categorical") ? "🏷️" : "📈";
          opt.textContent = `${icon} ${c.label}`;
          tipLabelSel.appendChild(opt);
          availableKeys.add(c.key);
        });

        // 3. Extra Taxon Metadata Fields (e.g. virus, organism, etc.)
        if (typeof TAXA_METADATA !== 'undefined' && TAXA_METADATA) {
          const firstTaxon = Object.keys(TAXA_METADATA)[0];
          const firstMeta = firstTaxon ? (TAXA_METADATA[firstTaxon] || {}) : {};
          const skipKeys = new Set(["id", "color", "taxon_id", "taxon_name", ...cols.map(c => c.key)]);
          Object.keys(firstMeta).forEach(k => {
            if (!skipKeys.has(k)) {
              const opt = document.createElement("option");
              opt.value = k;
              const cleanK = k.replace(/_/g, " ").replace(/\b\w/g, l => l.toUpperCase());
              opt.textContent = `📋 ${cleanK}`;
              tipLabelSel.appendChild(opt);
              availableKeys.add(k);
            }
          });
        }

        // 4. Combined: Taxon ID + Color Attribute
        const comboOpt = document.createElement("option");
        comboOpt.value = "id_and_color";
        comboOpt.textContent = "✨ Taxon ID + Active Color";
        tipLabelSel.appendChild(comboOpt);

        // Check if existing setting is still valid for this dataset; fallback if not
        if (!settings.tipLabelColumn || !availableKeys.has(settings.tipLabelColumn)) {
          settings.tipLabelColumn = "taxon_id";
        }
        tipLabelSel.value = settings.tipLabelColumn;

        const badge = document.getElementById("tipLabelActiveBadge");
        if (badge) {
          const activeOpt = tipLabelSel.selectedOptions[0];
          const activeText = activeOpt ? activeOpt.textContent.replace(/^[^a-zA-Z0-9]+/, "") : "Taxon ID";
          badge.textContent = activeText.split(" (")[0].slice(0, 16) || "Taxon ID";
        }
      }
    }

    function setTipLabelColumn(colKey) {
      settings.tipLabelColumn = colKey;
      const badge = document.getElementById("tipLabelActiveBadge");
      const sel = document.getElementById("tipLabelColumnSelect");
      if (badge && sel) {
        const activeOpt = sel.selectedOptions[0];
        const activeText = activeOpt ? activeOpt.textContent.replace(/^[^a-zA-Z0-9]+/, "") : "Taxon ID";
        badge.textContent = activeText.split(" (")[0].slice(0, 16) || "Taxon ID";
      }
      renderTree();
      updateMinimap();
      if (typeof renderMsa === 'function') renderMsa();
      if (settings.selectedTaxon) selectTaxon(settings.selectedTaxon);
    }

    function setColorColumn(colKey) {
      settings.colorColumn = colKey;
      updateLegend();
      renderTree();
      updateMinimap();
      if (typeof renderMsa === 'function') renderMsa();
      if (settings.selectedTaxon) selectTaxon(settings.selectedTaxon);
    }

    function setCladeGroupColumn(colKey) {
      settings.cladeGroupColumn = colKey;
      updateCladeManagementUI();
      renderTree();
    }

    function updateLegend() {
      // Branch Support Legend Section
      const bsLegend = document.getElementById("branchSupportLegend");
      const bsLabel = document.getElementById("branchSupportMetricLabel");
      const bsBody = document.getElementById("branchSupportLegendBody");
      if (bsLegend && bsBody) {
        if (!settings.colorBranchesBySupport) {
          bsLegend.classList.add("hidden");
        } else {
          bsLegend.classList.remove("hidden");
          if (bsLabel) bsLabel.textContent = settings.supportMetric === "alrt" ? "SH-aLRT" : "UFboot";
          if (settings.supportPalette === "traffic") {
            bsBody.innerHTML = `
              <div class="grid grid-cols-3 gap-1 pt-0.5">
                <div class="flex items-center space-x-1 text-[9px] text-[var(--text-muted)]">
                  <span class="w-2 h-2 rounded-full bg-emerald-500 shrink-0"></span>
                  <span class="truncate">≥95% Strong</span>
                </div>
                <div class="flex items-center space-x-1 text-[9px] text-[var(--text-muted)]">
                  <span class="w-2 h-2 rounded-full bg-amber-400 shrink-0"></span>
                  <span class="truncate">70-94% Mod</span>
                </div>
                <div class="flex items-center space-x-1 text-[9px] text-[var(--text-muted)]">
                  <span class="w-2 h-2 rounded-full bg-rose-500 shrink-0"></span>
                  <span class="truncate">&lt;70% Weak</span>
                </div>
              </div>
            `;
          } else if (settings.supportPalette === "viridis") {
            bsBody.innerHTML = `
              <div class="pt-0.5 space-y-1">
                <div class="h-2 w-full rounded-full border border-[var(--border-color)]/50" style="background: linear-gradient(90deg, rgb(68,1,84), rgb(59,82,139), rgb(33,145,140), rgb(94,201,98), rgb(253,231,37));"></div>
                <div class="flex justify-between text-[8.5px] font-mono text-[var(--text-muted)]">
                  <span>0%</span>
                  <span>50%</span>
                  <span>100%</span>
                </div>
              </div>
            `;
          } else if (settings.supportPalette === "grayscale") {
            bsBody.innerHTML = `
              <div class="pt-0.5 space-y-1">
                <div class="h-2 w-full rounded-full border border-[var(--border-color)]/50" style="background: linear-gradient(90deg, rgba(56,189,248,0.18), rgba(56,189,248,1.0));"></div>
                <div class="flex justify-between text-[8.5px] font-mono text-[var(--text-muted)]">
                  <span>Faded (0%)</span>
                  <span>Vivid (100%)</span>
                </div>
              </div>
            `;
          }
        }
      }

      const container = document.getElementById("legendGrid");
      const titleEl = document.getElementById("legendTitle");
      const badgeEl = document.getElementById("legendCountBadge");
      if (!container) return;
      container.innerHTML = "";

      if (settings.colorColumn === "solid") {
        if (titleEl) titleEl.textContent = "Uniform Accent Color";
        if (badgeEl) badgeEl.textContent = "Solid";
        container.innerHTML = `
          <div class="flex items-center space-x-2 text-[10px] text-[var(--text-muted)]">
            <span class="w-2.5 h-2.5 rounded-full bg-sky-400"></span>
            <span>All nodes colored uniformly</span>
          </div>`;
        return;
      }

      const colDef = getActiveColorColumnDef();
      if (!colDef) return;

      if (titleEl) titleEl.textContent = `Legend: ${colDef.label}`;

      if (colDef.type === "continuous") {
        if (badgeEl) badgeEl.textContent = `${colDef.min} – ${colDef.max}`;
        const isPlddt = colDef.key.toLowerCase().includes("plddt");
        const gradDiv = document.createElement("div");
        gradDiv.className = "w-full space-y-1.5 pt-1";

        let gradStyle;
        if (customPaletteState.isActive && customPaletteState.applyToContinuous && customPaletteState.palette.length > 0) {
          gradStyle = `background: linear-gradient(90deg, ${customPaletteState.palette.join(", ")});`;
        } else if (isPlddt) {
          gradStyle = "background: linear-gradient(90deg, #f97316 0%, #facc15 35%, #38bdf8 70%, #2563eb 100%);";
        } else {
          gradStyle = "background: linear-gradient(90deg, rgb(68,1,84) 0%, rgb(59,82,139) 25%, rgb(33,145,140) 50%, rgb(94,201,98) 75%, rgb(253,231,37) 100%);";
        }

        gradDiv.innerHTML = `
          <div class="h-2.5 w-full rounded-full border border-[var(--border-color)] shadow-inner" style="${gradStyle}"></div>
          <div class="flex justify-between text-[9px] font-mono text-[var(--text-muted)]">
            <span>${colDef.min}</span>
            <span class="font-semibold text-[var(--text-main)]">${colDef.label}</span>
            <span>${colDef.max}</span>
          </div>
        `;
        container.appendChild(gradDiv);
        return;
      }

      // Categorical Legend
      const counts = {};
      Object.values(TAXA_METADATA).forEach(m => {
        if (m && m[colDef.key] !== undefined && m[colDef.key] !== null) {
          const v = String(m[colDef.key]);
          counts[v] = (counts[v] || 0) + 1;
        }
      });

      const uniqueVals = colDef.values || Object.keys(counts);
      if (badgeEl) badgeEl.textContent = `${uniqueVals.length} categories`;

      uniqueVals.forEach(val => {
        const count = counts[val] || 0;
        if (count === 0 && colDef.cardinality > 25) return;
        let col;
        if (customPaletteState.isActive && customPaletteState.applyToCategorical && customPaletteState.palette.length > 0) {
          col = getCategoryColorFromCustomPalette(colDef, String(val));
        } else {
          col = (colDef.colors && colDef.colors[val]) ? colDef.colors[val] : "#94a3b8";
        }

        const div = document.createElement("div");
        div.className = "flex items-center justify-between text-[10px] text-[var(--text-muted)] cursor-pointer hover:text-[var(--text-main)] group py-0.5";
        div.innerHTML = `
          <div class="flex items-center space-x-1.5 truncate">
            <span class="w-2.5 h-2.5 rounded-full shrink-0" style="background-color: ${col};"></span>
            <span class="truncate group-hover:underline">${val}</span>
          </div>
          <span class="text-[9px] font-mono text-[var(--text-muted)] ml-1 shrink-0">${count}</span>
        `;
        div.onclick = () => {
          toggleCategoryClade(val, true);
        };
        container.appendChild(div);
      });
    }

    // 10. TREE RENDERING PIPELINE
    const svg = document.getElementById("treeSvg");

    function assignDepths(node, currentDepth = 0, currentClado = 0) {
      node.depth = currentDepth;
      node.cladoDepth = currentClado;
      if (node.children) {
        node.children.forEach(c => {
          const bLen = (typeof c.length === 'number' && !isNaN(c.length)) ? c.length : 1.0;
          assignDepths(c, currentDepth + bLen, currentClado + 1);
        });
      }
    }

    function setEsmViewMode(mode) {
      if (mode === "umap") {
        if (settings.dataset !== "esm2") {
          settings.dataset = "esm2";
          const dsSel = document.getElementById("datasetSelect");
          if (dsSel) dsSel.value = "esm2";
          switchDataset("esm2");
        }
        setLayout("umap");
      } else {
        setLayout("rectangular");
      }
    }

    function renderUmapScatter(g, defs) {
      umapPointCoords = {};
      const ds = activeDataset || (typeof DATASETS !== "undefined" && DATASETS[currentScale]) || {};
      const umapData = ds.esm2_umap;

      if (!umapData || Object.keys(umapData).length === 0) {
        const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
        text.setAttribute("x", "600");
        text.setAttribute("y", "350");
        text.setAttribute("text-anchor", "middle");
        text.setAttribute("fill", "var(--text-muted)");
        text.setAttribute("font-size", "14px");
        text.setAttribute("font-family", "system-ui, sans-serif");
        text.textContent = "🪐 No ESM-2 UMAP embeddings available for this dataset.";
        g.appendChild(text);

        const sub = document.createElementNS("http://www.w3.org/2000/svg", "text");
        sub.setAttribute("x", "600");
        sub.setAttribute("y", "380");
        sub.setAttribute("text-anchor", "middle");
        sub.setAttribute("fill", "var(--text-muted)");
        sub.setAttribute("font-size", "11px");
        sub.setAttribute("font-family", "system-ui, sans-serif");
        sub.textContent = "Run the PLM embedding pipeline with --model esm2 to extract latent representations.";
        g.appendChild(sub);
        return;
      }

      // 1. Filtered taxa check
      const allowedSet = (filterState.isActive && filterState.mode === "prune") ? getFilteredTaxaSet() : null;
      const cladeSet = (settings.scopedClade && settings.scopedClade.taxa && settings.scopedClade.taxa.size > 0) ? settings.scopedClade.taxa : null;

      let allTaxa = Object.keys(umapData);
      if (allowedSet) {
        allTaxa = allTaxa.filter(t => allowedSet.has(t));
      }
      if (cladeSet) {
        allTaxa = allTaxa.filter(t => cladeSet.has(t));
      }

      if (allTaxa.length === 0) {
        const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
        text.setAttribute("x", "600");
        text.setAttribute("y", "350");
        text.setAttribute("text-anchor", "middle");
        text.setAttribute("fill", "var(--text-muted)");
        text.setAttribute("font-size", "13px");
        text.textContent = "⚠️ All taxa filtered out in active UMAP view.";
        g.appendChild(text);
        return;
      }

      // 2. Coordinate Extents
      let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
      allTaxa.forEach(t => {
        const pt = umapData[t];
        if (pt && pt.length >= 2) {
          if (pt[0] < minX) minX = pt[0];
          if (pt[0] > maxX) maxX = pt[0];
          if (pt[1] < minY) minY = pt[1];
          if (pt[1] > maxY) maxY = pt[1];
        }
      });

      if (minX === maxX) { minX -= 1; maxX += 1; }
      if (minY === maxY) { minY -= 1; maxY += 1; }

      const width = 1000;
      const height = 650;
      const padX = 120;
      const padY = 70;

      function mapX(x) {
        return padX + ((x - minX) / (maxX - minX)) * width;
      }
      function mapY(y) {
        return padY + height - ((y - minY) / (maxY - minY)) * height;
      }

      // 3. Background Grid & Axis Crosshairs
      const gridG = document.createElementNS("http://www.w3.org/2000/svg", "g");
      gridG.setAttribute("class", "umap-grid");

      const frame = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      frame.setAttribute("x", padX);
      frame.setAttribute("y", padY);
      frame.setAttribute("width", width);
      frame.setAttribute("height", height);
      frame.setAttribute("fill", "none");
      frame.setAttribute("stroke", "var(--border-color)");
      frame.setAttribute("stroke-width", "1");
      frame.setAttribute("stroke-dasharray", "4,4");
      frame.setAttribute("opacity", "0.6");
      gridG.appendChild(frame);

      for (let step = 1; step < 5; step++) {
        const gx = padX + (width / 5) * step;
        const gy = padY + (height / 5) * step;

        const vx = document.createElementNS("http://www.w3.org/2000/svg", "line");
        vx.setAttribute("x1", gx);
        vx.setAttribute("x2", gx);
        vx.setAttribute("y1", padY);
        vx.setAttribute("y2", padY + height);
        vx.setAttribute("stroke", "var(--border-color)");
        vx.setAttribute("stroke-width", "0.75");
        vx.setAttribute("stroke-dasharray", "2,4");
        vx.setAttribute("opacity", "0.35");
        gridG.appendChild(vx);

        const hy = document.createElementNS("http://www.w3.org/2000/svg", "line");
        hy.setAttribute("x1", padX);
        hy.setAttribute("x2", padX + width);
        hy.setAttribute("y1", gy);
        hy.setAttribute("y2", gy);
        hy.setAttribute("stroke", "var(--border-color)");
        hy.setAttribute("stroke-width", "0.75");
        hy.setAttribute("stroke-dasharray", "2,4");
        hy.setAttribute("opacity", "0.35");
        gridG.appendChild(hy);
      }

      if (minX <= 0 && maxX >= 0) {
        const zx = mapX(0);
        const zeroV = document.createElementNS("http://www.w3.org/2000/svg", "line");
        zeroV.setAttribute("x1", zx);
        zeroV.setAttribute("x2", zx);
        zeroV.setAttribute("y1", padY);
        zeroV.setAttribute("y2", padY + height);
        zeroV.setAttribute("stroke", "var(--text-muted)");
        zeroV.setAttribute("stroke-width", "1");
        zeroV.setAttribute("opacity", "0.3");
        gridG.appendChild(zeroV);
      }
      if (minY <= 0 && maxY >= 0) {
        const zy = mapY(0);
        const zeroH = document.createElementNS("http://www.w3.org/2000/svg", "line");
        zeroH.setAttribute("x1", padX);
        zeroH.setAttribute("x2", padX + width);
        zeroH.setAttribute("y1", zy);
        zeroH.setAttribute("y2", zy);
        zeroH.setAttribute("stroke", "var(--text-muted)");
        zeroH.setAttribute("stroke-width", "1");
        zeroH.setAttribute("opacity", "0.3");
        gridG.appendChild(zeroH);
      }

      const lblX = document.createElementNS("http://www.w3.org/2000/svg", "text");
      lblX.setAttribute("x", padX + width / 2);
      lblX.setAttribute("y", padY + height + 32);
      lblX.setAttribute("text-anchor", "middle");
      lblX.setAttribute("fill", "var(--text-muted)");
      lblX.setAttribute("font-size", "11.5px");
      lblX.setAttribute("font-weight", "600");
      lblX.setAttribute("letter-spacing", "0.05em");
      lblX.textContent = "UMAP-1 (ESM-2 Latent Dimension 1)";
      gridG.appendChild(lblX);

      const lblY = document.createElementNS("http://www.w3.org/2000/svg", "text");
      lblY.setAttribute("transform", `translate(${padX - 25}, ${padY + height / 2}) rotate(-90)`);
      lblY.setAttribute("text-anchor", "middle");
      lblY.setAttribute("fill", "var(--text-muted)");
      lblY.setAttribute("font-size", "11.5px");
      lblY.setAttribute("font-weight", "600");
      lblY.setAttribute("letter-spacing", "0.05em");
      lblY.textContent = "UMAP-2 (ESM-2 Latent Dimension 2)";
      gridG.appendChild(lblY);

      g.appendChild(gridG);

      // 4. Header Info in Viewport
      const infoG = document.createElementNS("http://www.w3.org/2000/svg", "g");
      const infoTitle = document.createElementNS("http://www.w3.org/2000/svg", "text");
      infoTitle.setAttribute("x", padX + 8);
      infoTitle.setAttribute("y", padY - 24);
      infoTitle.setAttribute("fill", "var(--text-main)");
      infoTitle.setAttribute("font-size", "14px");
      infoTitle.setAttribute("font-weight", "700");
      infoTitle.textContent = "🪐 ESM-2 (1,280-dim) UMAP Projection Space";
      infoG.appendChild(infoTitle);

      const colDef = getActiveColorColumnDef();
      const infoSub = document.createElementNS("http://www.w3.org/2000/svg", "text");
      infoSub.setAttribute("x", padX + 8);
      infoSub.setAttribute("y", padY - 10);
      infoSub.setAttribute("fill", "var(--text-muted)");
      infoSub.setAttribute("font-size", "10.5px");
      infoSub.textContent = `${allTaxa.length.toLocaleString()} taxa • Metric: Cosine • Colored by ${colDef ? colDef.label : "Category"}`;
      infoG.appendChild(infoSub);
      g.appendChild(infoG);

      // 5. Draw Taxa Scatter Points
      const baseRadius = Math.max(3.0, (settings.nodeRadius || 3.0) * 1.6);
      const pointsG = document.createElementNS("http://www.w3.org/2000/svg", "g");
      pointsG.setAttribute("class", "umap-points");

      const labelsG = document.createElementNS("http://www.w3.org/2000/svg", "g");
      labelsG.setAttribute("class", "umap-labels");

      allTaxa.forEach(taxonName => {
        const pt = umapData[taxonName];
        if (!pt || pt.length < 2) return;

        const cx = mapX(pt[0]);
        const cy = mapY(pt[1]);
        umapPointCoords[taxonName] = { x: cx, y: cy };
        const color = getNodeColor({ name: taxonName });
        const isSelected = (taxonName === settings.selectedTaxon);
        const filterClass = getLeafFilterClass(taxonName);
        const isDimmed = (filterClass === "filter-dimmed");

        const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
        circle.setAttribute("cx", cx);
        circle.setAttribute("cy", cy);
        circle.setAttribute("r", isSelected ? baseRadius * 1.8 : (isDimmed ? baseRadius * 0.75 : baseRadius));
        circle.setAttribute("fill", color);
        circle.setAttribute("stroke", isSelected ? "#ffffff" : "rgba(0,0,0,0.35)");
        circle.setAttribute("stroke-width", isSelected ? "2.5" : "0.75");
        circle.setAttribute("opacity", isDimmed ? "0.15" : "0.88");
        circle.setAttribute("class", `umap-point cursor-pointer transition hover:opacity-100 ${isSelected ? "selected" : ""}`);
        circle.setAttribute("data-taxon", taxonName);

        circle.onmouseenter = (e) => {
          circle.setAttribute("r", isSelected ? baseRadius * 2.2 : baseRadius * 1.4);
          circle.setAttribute("stroke", "#ffffff");
          circle.setAttribute("stroke-width", "2");
          showNodeTooltip(e, { name: taxonName });
        };
        circle.onmouseleave = () => {
          circle.setAttribute("r", isSelected ? baseRadius * 1.8 : (isDimmed ? baseRadius * 0.75 : baseRadius));
          circle.setAttribute("stroke", isSelected ? "#ffffff" : "rgba(0,0,0,0.35)");
          circle.setAttribute("stroke-width", isSelected ? "2.5" : "0.75");
          scheduleHideTooltip();
        };
        circle.onclick = (e) => {
          e.stopPropagation();
          selectTaxon(taxonName, false);
        };

        pointsG.appendChild(circle);

        if (isSelected) {
          const ring = document.createElementNS("http://www.w3.org/2000/svg", "circle");
          ring.setAttribute("cx", cx);
          ring.setAttribute("cy", cy);
          ring.setAttribute("r", baseRadius * 2.6);
          ring.setAttribute("fill", "none");
          ring.setAttribute("stroke", color);
          ring.setAttribute("stroke-width", "2");
          ring.setAttribute("stroke-dasharray", "3,3");
          pointsG.appendChild(ring);
        }

        if (settings.showTipLabels && (!isDimmed || isSelected)) {
          const labelText = getLeafLabelText(taxonName);
          const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
          txt.setAttribute("x", cx + baseRadius + 4);
          txt.setAttribute("y", cy + 3.5);
          txt.setAttribute("fill", isSelected ? "var(--accent)" : "var(--text-main)");
          txt.setAttribute("font-size", `${Math.max(7.5, Math.min(11, (settings.labelSize || 9) * 0.95))}px`);
          txt.setAttribute("font-family", "system-ui, sans-serif");
          txt.setAttribute("font-weight", isSelected ? "700" : "500");
          txt.setAttribute("class", `tip-label pointer-events-none ${isSelected ? "selected" : ""}`);
          txt.textContent = labelText;
          labelsG.appendChild(txt);
        }
      });

      g.appendChild(pointsG);
      if (settings.showTipLabels) {
        g.appendChild(labelsG);
      }
    }

    function renderTree() {
      try {
        svg.innerHTML = "";

        const defs = document.createElementNS("http://www.w3.org/2000/svg", "defs");
        svg.appendChild(defs);

        const g = document.createElementNS("http://www.w3.org/2000/svg", "g");
        g.setAttribute("id", "mainViewport");
        g.setAttribute("transform", `translate(${settings.zoom.x}, ${settings.zoom.y}) scale(${settings.zoom.k})`);
        svg.appendChild(g);

        if (settings.layout === "tanglegram") {
          renderTanglegram(g, defs);
        } else if (settings.layout === "umap") {
          renderUmapScatter(g, defs);
        } else {
          assignDepths(activeTreeRoot, 0);

          const allLeaves = getAllLeaves(activeTreeRoot);
          const visibleLeaves = getVisibleLeaves(activeTreeRoot);

          const maxDepth = Math.max(...allLeaves.map(l => settings.branchLengths ? l.depth : l.cladoDepth)) || 1.0;

          if (settings.layout === "radial") {
            renderRadialTree(g, visibleLeaves, maxDepth);
          } else if (settings.layout === "unrooted") {
            renderUnrootedTree(g, visibleLeaves, maxDepth);
          } else {
            renderCartesianTree(g, visibleLeaves, maxDepth, defs);
          }
        }

        updateMinimap();
        renderMsa();
      } catch (err) {
        console.error("Tree render error:", err);
      }
    }


    // =========================================================================
    // FILTER ENGINE (PRUNE SUBTREE & HIGHLIGHT/DIM MODES)
    // =========================================================================
    const filterState = {
      isActive: false,
      mode: "prune", // "prune" or "highlight"
      column: "structural_class",
      selectedCategories: new Set(),
      minVal: null,
      maxVal: null,
      categorySearchQuery: ""
    };

    function getLeafFilterClass(taxonName) {
      if (!taxonName) return "";
      if (!filterState.isActive || filterState.mode !== "highlight") return "";
      const allowed = getFilteredTaxaSet();
      return allowed.has(taxonName) ? "filter-match" : "filter-dimmed";
    }

    function isTreeFilteringActive() {
      const colDef = getActiveFilterColumnDef();
      let colFiltering = false;
      if (colDef && colDef.type === "continuous") {
        colFiltering = (filterState.minVal !== null && filterState.minVal > (colDef.min || 0)) || 
                       (filterState.maxVal !== null && filterState.maxVal < (colDef.max || 100));
      } else if (colDef) {
        const allValsCount = (colDef.values || []).length;
        colFiltering = filterState.selectedCategories.size < allValsCount;
      }
      const msaCovFiltering = !!(typeof msaState !== 'undefined' && msaState.filterCoverage);
      const msaPartFiltering = !!(typeof msaState !== 'undefined' && msaState.activePartition && msaState.activePartition !== "all");
      return colFiltering || msaCovFiltering || msaPartFiltering;
    }

    function getFilteredTaxaSet() {
      const allTaxa = Object.keys(TAXA_METADATA);
      const colDef = getActiveFilterColumnDef();
      const res = new Set();

      // Check metadata column filter
      let isColFiltering = false;
      if (colDef) {
        if (colDef.type === "continuous") {
          isColFiltering = (filterState.minVal !== null && filterState.minVal > (colDef.min || 0)) || 
                           (filterState.maxVal !== null && filterState.maxVal < (colDef.max || 100));
        } else {
          const allValsCount = (colDef.values || []).length;
          isColFiltering = filterState.selectedCategories.size < allValsCount;
        }
      }

      // Check MSA coverage filter & partition filter
      const align = (typeof getActiveAlignment === 'function') ? getActiveAlignment() : null;
      const minCov = (typeof msaState !== 'undefined' && msaState.coverageThreshold) ? msaState.coverageThreshold : 0.70;
      const isCovFiltering = !!(typeof msaState !== 'undefined' && msaState.filterCoverage && align);
      const isPartFiltering = !!(typeof msaState !== 'undefined' && msaState.activePartition && msaState.activePartition !== "all" && align);

      allTaxa.forEach(name => {
        const m = TAXA_METADATA[name];

        // 1. Metadata column filter check
        if (isColFiltering && colDef) {
          if (colDef.type === "continuous") {
            if (!m || m[colDef.key] === undefined || m[colDef.key] === null) return;
            const v = parseFloat(m[colDef.key]);
            if (isNaN(v)) return;
            const minB = filterState.minVal !== null ? filterState.minVal : (colDef.min || 0);
            const maxB = filterState.maxVal !== null ? filterState.maxVal : (colDef.max || 100);
            if (v < minB || v > maxB) return;
          } else {
            const val = (m && m[colDef.key] !== undefined && m[colDef.key] !== null) ? String(m[colDef.key]) : "Unclassified";
            if (!filterState.selectedCategories.has(val)) return;
          }
        }

        // 2. Alignment partition membership check
        if (isPartFiltering && align) {
          const mode = (typeof msaState !== 'undefined' && msaState.mode) ? msaState.mode : "aa";
          const seqs = align[mode] || {};
          if (!seqs[name]) return;
        }

        // 3. Alignment coverage threshold check
        if (isCovFiltering && align) {
          const mode = (typeof msaState !== 'undefined' && msaState.mode) ? msaState.mode : "aa";
          const cov = getSequenceCoverage(name, align, mode);
          if (cov < minCov) return;
        }

        res.add(name);
      });

      return res;
    }

    function pruneSubtree(node, allowedSet) {
      if (!node) return null;
      if (!node.children || node.children.length === 0) {
        if (allowedSet.has(node.name)) {
          return {
            id: node.id,
            name: node.name,
            length: typeof node.length === 'number' ? node.length : 0.001,
            support: node.support,
            ufboot: node.ufboot,
            alrt: node.alrt,
            supportRaw: node.supportRaw || "",
            children: [],
            _collapsed: false
          };
        }
        return null;
      }

      const keptChildren = [];
      for (let i = 0; i < node.children.length; i++) {
        const pc = pruneSubtree(node.children[i], allowedSet);
        if (pc !== null) keptChildren.push(pc);
      }

      if (keptChildren.length === 0) return null;
      if (keptChildren.length === 1) {
        const single = keptChildren[0];
        single.length = (typeof single.length === 'number' ? single.length : 0.001) + (typeof node.length === 'number' ? node.length : 0.0);
        return single;
      }

      return {
        id: node.id,
        name: node.name || "",
        length: typeof node.length === 'number' ? node.length : 0.0,
        support: node.support,
        ufboot: node.ufboot,
        alrt: node.alrt,
        supportRaw: node.supportRaw || "",
        children: keptChildren,
        _collapsed: node._collapsed || false
      };
    }

    function getActiveFilterColumnDef() {
      const cols = getActiveColumns();
      if (!cols || cols.length === 0) return null;
      return cols.find(c => c.key === filterState.column) || cols[0];
    }

    function setFilterMode(m) {
      filterState.mode = m;
      const btnPrune = document.getElementById("btnFilterModePrune");
      const btnHl = document.getElementById("btnFilterModeHighlight");
      const desc = document.getElementById("filterModeDesc");

      if (m === "prune") {
        btnPrune.className = "py-1 px-2 rounded font-semibold text-center bg-emerald-500 text-white text-[10.5px] transition shadow-sm cursor-pointer";
        btnHl.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10.5px] transition cursor-pointer";
        if (desc) desc.textContent = "Prune contracts the tree to only matching leaves, recalculating branch geometry and evolutionary distances.";
      } else {
        btnHl.className = "py-1 px-2 rounded font-semibold text-center bg-emerald-500 text-white text-[10.5px] transition shadow-sm cursor-pointer";
        btnPrune.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-slate-500/20 text-[10.5px] transition cursor-pointer";
        if (desc) desc.textContent = "Highlight & Dim keeps full tree structure visible, highlighting matching taxa while dimming non-matching branches to 12% opacity.";
      }

      applyTaxaFilter();
    }

    function onFilterColumnChanged(colKey) {
      filterState.column = colKey;
      const colDef = getActiveFilterColumnDef();
      if (!colDef) return;

      if (colDef.type === "continuous") {
        filterState.minVal = colDef.min || 0;
        filterState.maxVal = colDef.max || 100;
      } else {
        filterState.selectedCategories = new Set(colDef.values || []);
      }

      updateFilterUI();
      applyTaxaFilter();
    }

    function updateFilterUI() {
      const sel = document.getElementById("filterColumnSelect");
      if (!sel) return;
      sel.innerHTML = "";
      (activeDataset.columns || []).forEach(c => {
        const opt = document.createElement("option");
        opt.value = c.key;
        opt.textContent = `${c.type === "continuous" ? "📈" : "🏷️"} ${c.label}`;
        if (c.key === filterState.column) opt.selected = true;
        sel.appendChild(opt);
      });

      const colDef = getActiveFilterColumnDef();
      const secCat = document.getElementById("filterCategoricalSection");
      const secNum = document.getElementById("filterContinuousSection");
      const presetsBox = document.getElementById("filterPresetsContainer");

      // Populate Quick Presets
      if (presetsBox) {
        presetsBox.innerHTML = "";
        let presets = [];
        if (currentScale === "1193") {
          presets = [
            { label: "Strong Binders (76)", col: "binding_strength", vals: ["Strong"] },
            { label: "All Binders (115)", col: "binding_strength", vals: ["Strong", "Medium"] },
            { label: "Mainly Alpha (564)", col: "structural_class", vals: ["Mainly Alpha"] },
            { label: "High pLDDT (≥80)", col: "plddt", min: 80, max: 96 }
          ];
        } else if (currentScale === "500") {
          presets = [
            { label: "Rhabdoviridae (120)", col: "family", vals: ["Rhabdoviridae"] },
            { label: "Phenuiviridae (48)", col: "family", vals: ["Phenuiviridae"] },
            { label: "High pLDDT (≥75)", col: "plddt", min: 75, max: 95 }
          ];
        } else {
          presets = [
            { label: "All 6 Taxa", col: colDef ? colDef.key : "family", vals: colDef ? (colDef.values || []) : [] }
          ];
        }

        if (window.ALIGNMENTS_DATA && window.ALIGNMENTS_DATA[currentScale]) {
          presets.push({ label: "Core Align (≥70% Cov)", col: "alignment_coverage", min: 70, max: 100 });
        }

        presets.forEach(p => {
          const chip = document.createElement("button");
          chip.className = "text-[9.5px] px-2 py-0.5 rounded-full bg-emerald-500/10 hover:bg-emerald-500/25 text-emerald-300 border border-emerald-500/30 transition cursor-pointer";
          chip.textContent = p.label;
          chip.onclick = () => {
            filterState.column = p.col;
            if (sel) sel.value = p.col;
            if (p.vals) {
              filterState.selectedCategories = new Set(p.vals);
            } else {
              filterState.minVal = p.min;
              filterState.maxVal = p.max;
            }
            updateFilterUI();
            applyTaxaFilter();
          };
          presetsBox.appendChild(chip);
        });
      }

      if (colDef && colDef.type === "continuous") {
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

        if (inMin) { inMin.min = minV; inMin.max = maxV; inMin.value = filterState.minVal; }
        if (inMax) { inMax.min = minV; inMax.max = maxV; inMax.value = filterState.maxVal; }
        if (slider) { slider.min = minV; slider.max = maxV; slider.value = filterState.minVal; }
        if (legMin) legMin.textContent = `Min: ${minV}`;
        if (legMax) legMax.textContent = `Max: ${maxV}`;
        if (badge) badge.textContent = `[${filterState.minVal} – ${filterState.maxVal}]`;
      } else if (colDef) {
        if (secNum) secNum.classList.add("hidden");
        if (secCat) secCat.classList.remove("hidden");

        const catList = document.getElementById("filterCategoryList");
        if (!catList) return;
        catList.innerHTML = "";

        const counts = {};
        Object.keys(TAXA_METADATA).forEach(name => {
          const m = TAXA_METADATA[name];
          const val = (m && m[colDef.key] !== undefined && m[colDef.key] !== null) ? String(m[colDef.key]) : "Unclassified";
          counts[val] = (counts[val] || 0) + 1;
        });

        if (filterState.selectedCategories.size === 0 && (!colDef.values || colDef.values.length > 0)) {
          filterState.selectedCategories = new Set(colDef.values || Object.keys(counts));
        }

        const q = (filterState.categorySearchQuery || "").toLowerCase();
        const vals = colDef.values || Object.keys(counts);

        vals.forEach(val => {
          if (q && !val.toLowerCase().includes(q)) return;
          const isChecked = filterState.selectedCategories.has(val);
          const color = (colDef.colors && colDef.colors[val]) ? colDef.colors[val] : "#94a3b8";
          const count = counts[val] || 0;

          const row = document.createElement("label");
          row.className = "flex items-center justify-between p-1.5 rounded hover:bg-slate-500/10 cursor-pointer text-xs select-none";
          row.innerHTML = `
            <div class="flex items-center space-x-2 truncate pr-1">
              <input type="checkbox" ${isChecked ? "checked" : ""} class="accent-emerald-500 cursor-pointer">
              <span class="w-2.5 h-2.5 rounded-full shrink-0" style="background-color: ${color}"></span>
              <span class="text-[var(--text-main)] truncate text-[11px] font-medium">${val}</span>
            </div>
            <span class="text-[9.5px] font-mono px-1.5 py-0.2 rounded-full bg-slate-500/20 text-[var(--text-muted)] shrink-0">${count}</span>
          `;

          const cb = row.querySelector("input");
          cb.onchange = (e) => {
            if (e.target.checked) {
              filterState.selectedCategories.add(val);
            } else {
              filterState.selectedCategories.delete(val);
            }
            applyTaxaFilter();
          };

          catList.appendChild(row);
        });
      }
    }

    function onFilterCategorySearch(q) {
      filterState.categorySearchQuery = q;
      updateFilterUI();
    }

    function selectAllFilterCategories() {
      const colDef = getActiveFilterColumnDef();
      if (!colDef || colDef.type === "continuous") return;
      filterState.selectedCategories = new Set(colDef.values || []);
      updateFilterUI();
      applyTaxaFilter();
    }

    function deselectAllFilterCategories() {
      filterState.selectedCategories.clear();
      updateFilterUI();
      applyTaxaFilter();
    }

    function invertFilterCategories() {
      const colDef = getActiveFilterColumnDef();
      if (!colDef || colDef.type === "continuous") return;
      const allVals = colDef.values || [];
      const inverted = new Set();
      allVals.forEach(v => {
        if (!filterState.selectedCategories.has(v)) inverted.add(v);
      });
      filterState.selectedCategories = inverted;
      updateFilterUI();
      applyTaxaFilter();
    }

    function onFilterRangeInputChanged() {
      const inMin = document.getElementById("filterInputMin");
      const inMax = document.getElementById("filterInputMax");
      if (inMin && inMax) {
        filterState.minVal = parseFloat(inMin.value);
        filterState.maxVal = parseFloat(inMax.value);
        const badge = document.getElementById("filterRangeBadge");
        if (badge) badge.textContent = `[${filterState.minVal} – ${filterState.maxVal}]`;
        applyTaxaFilter();
      }
    }

    function onFilterRangeSliderMoved(val) {
      filterState.minVal = parseFloat(val);
      const inMin = document.getElementById("filterInputMin");
      if (inMin) inMin.value = filterState.minVal;
      const badge = document.getElementById("filterRangeBadge");
      if (badge) badge.textContent = `[${filterState.minVal} – ${filterState.maxVal}]`;
      applyTaxaFilter();
    }

    function applyTaxaFilter() {
      const allTaxa = Object.keys(TAXA_METADATA);
      const colDef = getActiveFilterColumnDef();
      const filtered = getFilteredTaxaSet();

      const total = allTaxa.length;
      const matchingCount = filtered.size;

      const isFiltering = isTreeFilteringActive();
      filterState.isActive = isFiltering;

      const pct = total > 0 ? ((matchingCount / total) * 100).toFixed(1) : 100;

      const tabBadge = document.getElementById("tabFilterBadge");
      const statusBadge = document.getElementById("filterStatusBadge");
      const banner = document.getElementById("activeFilterBanner");
      const bannerText = document.getElementById("filterBannerText");

      if (tabBadge) {
        tabBadge.textContent = isFiltering ? `${matchingCount}/${total}` : "All";
        tabBadge.className = isFiltering 
          ? "ml-1 text-[9px] font-mono px-1.5 py-0.5 rounded-full bg-amber-500/20 text-amber-300 font-bold border border-amber-500/40"
          : "ml-1 text-[9px] font-mono px-1.5 py-0.5 rounded-full bg-emerald-500/20 text-emerald-400 font-bold";
      }

      if (statusBadge) {
        statusBadge.textContent = isFiltering ? `Filtered: ${matchingCount} / ${total} (${pct}%)` : `Showing All (${total})`;
        statusBadge.className = isFiltering
          ? "text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-amber-500/10 text-amber-300 border border-amber-500/30"
          : "text-[9.5px] font-mono px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/30";
      }

      if (banner && bannerText) {
        if (isFiltering) {
          banner.classList.remove("hidden");
          banner.classList.add("flex");
          if (settings.scopedClade) {
            banner.classList.add("top-11");
            banner.classList.remove("top-3");
          } else {
            banner.classList.remove("top-11");
            banner.classList.add("top-3");
          }
          const modeLabel = filterState.mode === "prune" ? "Pruned Subtree" : "Highlight";
          let labelDesc = colDef ? colDef.label : "Filter";
          if (typeof msaState !== 'undefined' && msaState.filterCoverage) {
            labelDesc = `≥${Math.round((msaState.coverageThreshold || 0.70) * 100)}% Coverage`;
          }
          if (typeof msaState !== 'undefined' && msaState.activePartition && msaState.activePartition !== "all") {
            const align = (typeof getActiveAlignment === 'function') ? getActiveAlignment() : null;
            const pName = (align && align.name) ? align.name : msaState.activePartition;
            labelDesc += ` (${pName})`;
          }
          bannerText.textContent = `🎯 ${modeLabel}: ${matchingCount} / ${total} taxa (${labelDesc})`;
        } else {
          banner.classList.add("hidden");
          banner.classList.remove("flex");
        }
      }

      applyCurrentRooting();
      updateCongruenceUI();
      if (typeof renderMsa === 'function') {
        if (typeof msaState !== 'undefined') msaState._minimapCacheKey = null;
        renderMsa();
      }
    }

    function clearTaxaFilter() {
      const colDef = getActiveFilterColumnDef();
      if (colDef && colDef.type === "continuous") {
        filterState.minVal = colDef.min || 0;
        filterState.maxVal = colDef.max || 100;
      } else if (colDef) {
        filterState.selectedCategories = new Set(colDef.values || []);
      }
      filterState.isActive = false;
      filterState.categorySearchQuery = "";
      if (typeof msaState !== 'undefined') {
        msaState.filterCoverage = false;
        msaState.activePartition = "all";
        const btn = document.getElementById("btnFilterCoverage");
        if (btn) {
          btn.className = "badge-amber px-1.5 py-0.5 rounded text-[9.5px] font-semibold cursor-pointer transition hover:opacity-90";
          btn.textContent = "Filter";
          btn.title = "Toggle Alignment Coverage Filter on Tree & MSA";
        }
        const partSel = document.getElementById("selectAlignmentPartition");
        if (partSel) partSel.value = "all";
      }
      updateFilterUI();
      applyTaxaFilter();
      if (typeof renderMsa === 'function') {
        if (typeof msaState !== 'undefined') msaState._minimapCacheKey = null;
        renderMsa();
      }
      showCladeToast("Filter reset: displaying full tree dataset.");
    }

    // =========================================================================
    // PHYLOGENETIC CONGRUENCE & MODAL ENGINE
    // =========================================================================
    function cloneTree(node) {
      if (!node) return null;
      const copy = {
        id: node.id,
        name: node.name,
        length: typeof node.length === 'number' ? node.length : 0.0,
        support: node.support,
        ufboot: node.ufboot,
        alrt: node.alrt,
        supportRaw: node.supportRaw || "",
        _collapsed: !!node._collapsed
      };
      if (node.children && node.children.length > 0) {
        copy.children = node.children.map(cloneTree);
      } else {
        copy.children = [];
      }
      return copy;
    }

    function getTreeLeavesInOrder(node) {
      if (!node) return [];
      if (!node.children || node.children.length === 0) return node.name ? [node] : [];
      let res = [];
      node.children.forEach(c => {
        res = res.concat(getTreeLeavesInOrder(c));
      });
      return res;
    }

    function getActiveTanglegramTrees() {
      const compMode = settings.tangleCompare || "3di_vs_aa";
      let leftSrc = rawRoot3Di;
      let rightSrc = rawRootAA;
      let leftTitle = "3Di Structural Tree (FoldMason + Q.3Di.LLM)";
      let rightTitle = "Amino Acid Tree (IQ-TREE LG+G4)";
      let leftColor = "var(--accent)";
      let rightColor = "#a855f7";

      if (compMode === "3di_vs_esm2_cosine" || compMode === "3di_vs_esm2") {
        leftSrc = rawRoot3Di;
        rightSrc = rawTrees["esm2_cosine"] || rawRoot3Di;
        leftTitle = "3Di Structural Tree (FoldMason + Q.3Di.LLM)";
        rightTitle = "ESM-2 PLM (Cosine Distance UPGMA)";
        rightColor = "#10b981";
      } else if (compMode === "3di_vs_esm2_euclidean") {
        leftSrc = rawRoot3Di;
        rightSrc = rawTrees["esm2_euclidean"] || rawRoot3Di;
        leftTitle = "3Di Structural Tree (FoldMason + Q.3Di.LLM)";
        rightTitle = "ESM-2 PLM (Euclidean Distance UPGMA)";
        rightColor = "#06b6d4";
      } else if (compMode === "3di_vs_esm2_l1") {
        leftSrc = rawRoot3Di;
        rightSrc = rawTrees["esm2_l1"] || rawRoot3Di;
        leftTitle = "3Di Structural Tree (FoldMason + Q.3Di.LLM)";
        rightTitle = "ESM-2 PLM (Manhattan / L1 UPGMA)";
        rightColor = "#f59e0b";
      } else if (compMode === "aa_vs_esm2_cosine" || compMode === "aa_vs_esm2") {
        leftSrc = rawRootAA;
        rightSrc = rawTrees["esm2_cosine"] || rawRoot3Di;
        leftTitle = "Amino Acid Tree (IQ-TREE LG+G4)";
        rightTitle = "ESM-2 PLM (Cosine Distance UPGMA)";
        leftColor = "#a855f7";
        rightColor = "#10b981";
      } else if (compMode === "aa_vs_esm2_euclidean") {
        leftSrc = rawRootAA;
        rightSrc = rawTrees["esm2_euclidean"] || rawRoot3Di;
        leftTitle = "Amino Acid Tree (IQ-TREE LG+G4)";
        rightTitle = "ESM-2 PLM (Euclidean Distance UPGMA)";
        leftColor = "#a855f7";
        rightColor = "#06b6d4";
      } else if (compMode === "aa_vs_esm2_l1") {
        leftSrc = rawRootAA;
        rightSrc = rawTrees["esm2_l1"] || rawRoot3Di;
        leftTitle = "Amino Acid Tree (IQ-TREE LG+G4)";
        rightTitle = "ESM-2 PLM (Manhattan / L1 UPGMA)";
        leftColor = "#a855f7";
        rightColor = "#f59e0b";
      } else if (compMode === "esm2_cosine_vs_euclidean") {
        leftSrc = rawTrees["esm2_cosine"] || rawRoot3Di;
        rightSrc = rawTrees["esm2_euclidean"] || rawRoot3Di;
        leftTitle = "ESM-2 PLM (Cosine Distance)";
        rightTitle = "ESM-2 PLM (Euclidean Distance)";
        leftColor = "#10b981";
        rightColor = "#06b6d4";
      } else if (compMode === "esm2_cosine_vs_l1") {
        leftSrc = rawTrees["esm2_cosine"] || rawRoot3Di;
        rightSrc = rawTrees["esm2_l1"] || rawRoot3Di;
        leftTitle = "ESM-2 PLM (Cosine Distance)";
        rightTitle = "ESM-2 PLM (Manhattan / L1 Distance)";
        leftColor = "#10b981";
        rightColor = "#f59e0b";
      }

      let tRootLeft = leftSrc ? cloneTree(leftSrc) : null;
      let tRootRight = rightSrc ? cloneTree(rightSrc) : null;

      if (filterState.isActive && filterState.mode === "prune") {
        const allowed = getFilteredTaxaSet();
        if (allowed && allowed.size > 0 && tRootLeft && tRootRight) {
          const pL = pruneSubtree(tRootLeft, allowed);
          const pR = pruneSubtree(tRootRight, allowed);
          if (pL && pR) {
            tRootLeft = pL;
            tRootRight = pR;
          }
        }
      }

      if (settings.scopedClade && settings.scopedClade.taxa && settings.scopedClade.taxa.size > 0 && tRootLeft && tRootRight) {
        const pL = pruneSubtree(tRootLeft, settings.scopedClade.taxa);
        const pR = pruneSubtree(tRootRight, settings.scopedClade.taxa);
        if (pL && pR) {
          tRootLeft = pL;
          tRootRight = pR;
        }
      }

      // Ensure both trees contain identical common taxa
      if (tRootLeft && tRootRight) {
        const leavesL = getAllLeaves(tRootLeft).map(l => l.name).filter(Boolean);
        const leavesR = getAllLeaves(tRootRight).map(l => l.name).filter(Boolean);
        const setL = new Set(leavesL);
        const setR = new Set(leavesR);
        if (setL.size !== setR.size || leavesL.some(n => !setR.has(n))) {
          const common = new Set(leavesL.filter(n => setR.has(n)));
          if (common.size > 0) {
            tRootLeft = pruneSubtree(tRootLeft, common);
            tRootRight = pruneSubtree(tRootRight, common);
          }
        }
      }

      return {
        leftTree: tRootLeft,
        rightTree: tRootRight,
        leftTitle: leftTitle,
        rightTitle: rightTitle,
        leftColor: leftColor,
        rightColor: rightColor,
        compMode: compMode
      };
    }

    // Dynamic Robinson-Foulds calculation
    function computeRF(t1, t2) {
      const leaves1 = getAllLeaves(t1).map(l => l.name).filter(Boolean);
      const leaves2 = getAllLeaves(t2).map(l => l.name).filter(Boolean);
      const set2 = new Set(leaves2);
      const common = leaves1.filter(n => set2.has(n));
      const N = common.length;
      if (N < 3) {
        return { shared: 0, maxRf: 0, rfDist: 0, pct: 100, totalTaxa: N };
      }

      const taxonToIdx = {};
      common.forEach((name, idx) => { taxonToIdx[name] = idx; });

      function extractSplits(root) {
        const splits = new Set();
        function postOrder(node) {
          let leafIndices = [];
          if (!node.children || node.children.length === 0) {
            if (node.name && taxonToIdx[node.name] !== undefined) {
              leafIndices.push(taxonToIdx[node.name]);
            }
            return leafIndices;
          }
          node.children.forEach(c => {
            const cl = postOrder(c);
            for (let i = 0; i < cl.length; i++) leafIndices.push(cl[i]);
          });
          const k = leafIndices.length;
          if (k > 1 && k < N) {
            let bitArr;
            if (leafIndices.includes(0)) {
              const subSet = new Set(leafIndices);
              bitArr = [];
              for (let i = 0; i < N; i++) {
                if (!subSet.has(i)) bitArr.push(i);
              }
            } else {
              bitArr = leafIndices.slice().sort((a, b) => a - b);
            }
            splits.add(bitArr.join(","));
          }
          return leafIndices;
        }
        postOrder(root);
        return splits;
      }

      const s1 = extractSplits(t1);
      const s2 = extractSplits(t2);

      let shared = 0;
      s1.forEach(sp => {
        if (s2.has(sp)) shared++;
      });

      const totalPossible = Math.max(1, Math.max(s1.size, s2.size, N - 3));
      const maxRf = (s1.size + s2.size) > 0 ? (s1.size + s2.size) : 2 * Math.max(1, N - 3);
      const rfDist = (s1.size + s2.size) - 2 * shared;
      const pct = Math.max(0, Math.min(100, (shared / totalPossible) * 100));

      return { shared, maxRf, rfDist, pct, totalTaxa: N };
    }

    // Fast Cophenetic Pearson r
    function computeCopheneticR(t1, t2) {
      const leaves1 = getAllLeaves(t1).map(l => l.name).filter(Boolean);
      const leaves2 = getAllLeaves(t2).map(l => l.name).filter(Boolean);
      const set2 = new Set(leaves2);
      const common = leaves1.filter(n => set2.has(n));
      const N = common.length;
      if (N < 3) return 1.0;

      const taxonToIdx = {};
      common.forEach((n, i) => { taxonToIdx[n] = i; });

      function getPairwiseDistances(root) {
        const dist = new Float32Array(N * N);
        const leafDepths = new Float32Array(N);

        function assignDepthsRec(node, currentD) {
          const d = currentD + (typeof node.length === 'number' && !isNaN(node.length) ? Math.max(0.0001, node.length) : 0.01);
          if (!node.children || node.children.length === 0) {
            if (node.name && taxonToIdx[node.name] !== undefined) {
              leafDepths[taxonToIdx[node.name]] = d;
            }
            return [taxonToIdx[node.name]];
          }
          const childLeaves = [];
          node.children.forEach(c => {
            const sub = assignDepthsRec(c, d);
            if (sub && sub.length > 0) childLeaves.push(sub.filter(idx => idx !== undefined));
          });
          for (let a = 0; a < childLeaves.length; a++) {
            for (let b = a + 1; b < childLeaves.length; b++) {
              const grpA = childLeaves[a];
              const grpB = childLeaves[b];
              for (let i = 0; i < grpA.length; i++) {
                const u = grpA[i];
                for (let j = 0; j < grpB.length; j++) {
                  const v = grpB[j];
                  dist[u * N + v] = d;
                  dist[v * N + u] = d;
                }
              }
            }
          }
          const allSub = [];
          childLeaves.forEach(g => {
            for (let i = 0; i < g.length; i++) allSub.push(g[i]);
          });
          return allSub;
        }

        assignDepthsRec(root, 0);

        for (let i = 0; i < N; i++) {
          for (let j = i + 1; j < N; j++) {
            const lcaDepth = dist[i * N + j];
            const patristic = (leafDepths[i] + leafDepths[j]) - 2 * lcaDepth;
            dist[i * N + j] = Math.max(0.0001, patristic);
          }
        }
        return dist;
      }

      const d1 = getPairwiseDistances(t1);
      const d2 = getPairwiseDistances(t2);

      let nPairs = 0;
      let sumX = 0, sumY = 0, sumX2 = 0, sumY2 = 0, sumXY = 0;
      for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
          const x = d1[i * N + j];
          const y = d2[i * N + j];
          sumX += x;
          sumY += y;
          sumX2 += x * x;
          sumY2 += y * y;
          sumXY += x * y;
          nPairs++;
        }
      }

      if (nPairs < 2) return 1.0;
      const num = nPairs * sumXY - sumX * sumY;
      const den = Math.sqrt(Math.max(0, (nPairs * sumX2 - sumX * sumX) * (nPairs * sumY2 - sumY * sumY)));
      if (den === 0) return 0.0;
      return Math.max(-1.0, Math.min(1.0, num / den));
    }

    // Discordance & Taxa Shifts
    function computeDiscordanceAndTaxaShifts(tLeft, tRight, metadata) {
      const leavesL = getTreeLeavesInOrder(tLeft).map(l => l.name).filter(Boolean);
      const leavesR = getTreeLeavesInOrder(tRight).map(l => l.name).filter(Boolean);
      const setR = new Set(leavesR);
      const common = leavesL.filter(n => setR.has(n));
      const N = common.length;
      if (N < 2) {
        return { nativeCrossings: 0, nativePct: 0, untangledCrossings: 0, untangledPct: 0, discordantTaxa: [] };
      }

      const rankL = {};
      leavesL.forEach((n, idx) => { rankL[n] = idx; });
      const rankR = {};
      leavesR.forEach((n, idx) => { rankR[n] = idx; });

      const shifts = common.map(id => {
        const r1 = rankL[id];
        const r2 = rankR[id];
        const shift = Math.abs(r1 - r2);
        const m = (metadata && metadata[id]) || {};
        const cat = m.family || m.lineage || m.structural_class || m.binding_strength || m.genus || "Taxon";
        return { id, category: cat, r1: r1 + 1, r2: r2 + 1, shift };
      });
      shifts.sort((a, b) => b.shift - a.shift);

      const maxCross = (N * (N - 1)) / 2;
      let nativeCross = 0;
      const yArr = common.map(id => rankR[id]);
      for (let i = 0; i < yArr.length; i++) {
        for (let j = i + 1; j < yArr.length; j++) {
          if (yArr[i] > yArr[j]) nativeCross++;
        }
      }
      const nativePct = maxCross > 0 ? (nativeCross / maxCross) * 100 : 0;

      // Untangled rotation
      const tRightClone = cloneTree(tRight);
      const mapL = rankL;
      function getAvgL(node) {
        const lf = getTreeLeavesInOrder(node);
        const idxs = lf.map(l => mapL[l.name] !== undefined ? mapL[l.name] : 0);
        return idxs.reduce((a, b) => a + b, 0) / (idxs.length || 1);
      }
      function rotateSubtrees(node) {
        if (!node.children || node.children.length <= 1) return;
        node.children.forEach(rotateSubtrees);
        node.children.sort((a, b) => getAvgL(a) - getAvgL(b));
      }
      rotateSubtrees(tRightClone);
      const untangledLeavesR = getTreeLeavesInOrder(tRightClone).map(l => l.name).filter(n => rankL[n] !== undefined);
      const untangledRankR = {};
      untangledLeavesR.forEach((n, idx) => { untangledRankR[n] = idx; });
      const uArr = common.map(id => untangledRankR[id]);
      let untangledCross = 0;
      for (let i = 0; i < uArr.length; i++) {
        for (let j = i + 1; j < uArr.length; j++) {
          if (uArr[i] > uArr[j]) untangledCross++;
        }
      }
      const untangledPct = maxCross > 0 ? (untangledCross / maxCross) * 100 : 0;

      return {
        nativeCrossings: nativeCross,
        nativePct: nativePct,
        untangledCrossings: untangledCross,
        untangledPct: untangledPct,
        discordantTaxa: shifts.slice(0, 15)
      };
    }

    const congruenceCache = new Map();

    function computeActiveCongruence() {
      const treesInfo = getActiveTanglegramTrees();
      const { leftTree, rightTree, leftTitle, rightTitle, compMode } = treesInfo;

      if (!leftTree || !rightTree) {
        return {
          taxa_count: 0,
          shared_splits: 0,
          max_rf: 0,
          rf_distance: 0,
          congruence_pct: 0,
          cophenetic_r: 0,
          discordance_native: 0,
          discordance_untangled: 0,
          discordant_taxa: [],
          interpretation: "No overlapping taxa found for the current tree pair and filter.",
          leftTitle: leftTitle || "Left Tree",
          rightTitle: rightTitle || "Right Tree",
          filterDesc: ""
        };
      }

      // Compute descriptive filter/partition label
      let filterDesc = "";
      if (typeof msaState !== 'undefined' && msaState.filterCoverage) {
        filterDesc += `≥${Math.round((msaState.coverageThreshold || 0.70) * 100)}% Coverage`;
      }
      if (typeof msaState !== 'undefined' && msaState.activePartition && msaState.activePartition !== "all") {
        const align = (typeof getActiveAlignment === 'function') ? getActiveAlignment() : null;
        const pName = (align && align.name) ? align.name : msaState.activePartition;
        filterDesc += (filterDesc ? ", " : "") + `Partition: ${pName}`;
      }
      if (filterState.isActive && filterState.mode === "prune") {
        const colDef = getActiveFilterColumnDef();
        if (colDef) {
          filterDesc += (filterDesc ? ", " : "") + `${colDef.label} filtered`;
        }
      }
      if (settings.scopedClade && settings.scopedClade.name) {
        filterDesc += (filterDesc ? ", " : "") + `Clade: ${settings.scopedClade.name}`;
      }

      const leaves = getAllLeaves(leftTree);
      const taxaCount = leaves.length;
      const cacheKey = `${currentScale || "active"}:${compMode}:${taxaCount}:${filterDesc}`;
      if (congruenceCache.has(cacheKey)) {
        return congruenceCache.get(cacheKey);
      }

      const dsCong = (activeDataset && activeDataset.congruence && activeDataset.congruence[compMode]) || null;
      const isFiltered = (filterState.isActive && filterState.mode === "prune") || 
                         (typeof msaState !== 'undefined' && (msaState.filterCoverage || (msaState.activePartition && msaState.activePartition !== "all"))) ||
                         (settings.scopedClade && settings.scopedClade.taxa && settings.scopedClade.taxa.size > 0);

      const rf = computeRF(leftTree, rightTree);
      const copheneticR = computeCopheneticR(leftTree, rightTree);
      const disc = computeDiscordanceAndTaxaShifts(leftTree, rightTree, TAXA_METADATA);

      let interp = "";
      if (!isFiltered && dsCong && dsCong.interpretation) {
        interp = dsCong.interpretation;
      } else if (isFiltered && currentScale === "1193" && filterState.column === "binding_strength" && filterState.selectedCategories.has("Strong") && filterState.selectedCategories.size === 1) {
        interp = "Sub-cohort analysis of the 76 Strong Binders reveals enhanced cophylogenetic concordance: cophenetic distance correlation jumps from r = +0.317 in the full library up to r = +0.578 among strong binders. This indicates that tight henipavirus receptor engagement imposes rigid structural constraints that co-evolve with sequence signatures.";
      } else {
        interp = `Congruence analysis across ${taxaCount} taxa reveals ${rf.pct.toFixed(1)}% Robinson-Foulds topological congruence (${rf.shared} shared clades, RF distance ${rf.rfDist}) and cophenetic distance correlation r = ${copheneticR >= 0 ? "+" : ""}${copheneticR.toFixed(3)}. Planar crossing discordance is ${disc.nativePct.toFixed(1)}% (reduced to ${disc.untangledPct.toFixed(1)}% upon clade-preserving untangling).`;
      }

      const result = {
        taxa_count: taxaCount,
        shared_splits: rf.shared,
        max_rf: rf.maxRf,
        rf_distance: rf.rfDist,
        congruence_pct: rf.pct,
        cophenetic_r: copheneticR,
        discordance_native: disc.nativePct,
        discordance_untangled: disc.untangledPct,
        discordant_taxa: disc.discordantTaxa,
        interpretation: interp,
        leftTitle,
        rightTitle,
        filterDesc
      };

      congruenceCache.set(cacheKey, result);
      return result;
    }

    function updateCongruenceUI() {
      const activeCong = computeActiveCongruence();

      const rfPct = activeCong.congruence_pct;
      const sharedSplits = activeCong.shared_splits;
      const rfDist = activeCong.rf_distance;
      const copheneticR = activeCong.cophenetic_r;
      const discordance = settings.tangleMode === "min_crossings" 
        ? activeCong.discordance_untangled 
        : (settings.tangleMode === "aligned" ? 0 : activeCong.discordance_native);

      const elRf = document.getElementById("tangleRfCongruence");
      const elRfBar = document.getElementById("tangleRfProgressBar");
      const elShared = document.getElementById("tangleRfShared");
      const elDist = document.getElementById("tangleRfDist");
      const elR = document.getElementById("tangleCopheneticR");
      const elDisc = document.getElementById("tangleDiscordance");

      if (elRf) elRf.textContent = `${rfPct.toFixed(1)}% (${sharedSplits} clades)`;
      if (elRfBar) elRfBar.style.width = `${Math.min(100, Math.max(5, rfPct))}%`;
      if (elShared) elShared.textContent = `${sharedSplits} shared clades`;
      if (elDist) elDist.textContent = `RF dist: ${rfDist}`;
      if (elR) elR.textContent = `r = ${copheneticR >= 0 ? "+" : ""}${copheneticR.toFixed(3)}`;
      if (elDisc) {
        const modeNote = settings.tangleMode === "min_crossings" ? " (untangled)" : (settings.tangleMode === "aligned" ? " (aligned)" : "");
        elDisc.textContent = `${discordance.toFixed(1)}% discordance${modeNote}`;
      }
    }

    function openCongruenceModal() {
      const modal = document.getElementById("congruenceModal");
      if (!modal) return;
      modal.classList.remove("hidden");

      const activeCong = computeActiveCongruence();

      const rfPct = activeCong.congruence_pct;
      const sharedSplits = activeCong.shared_splits;
      const rfDist = activeCong.rf_distance;
      const maxRf = activeCong.max_rf;
      const copheneticR = activeCong.cophenetic_r;
      const discordance = activeCong.discordance_native;
      const untangledDisc = activeCong.discordance_untangled;
      const interp = activeCong.interpretation;

      const mTitle = document.getElementById("modalCongruenceTitle");
      const mSub = document.getElementById("modalCongruenceSubtitle");
      const mRf = document.getElementById("modalRfScore");
      const mRfSub = document.getElementById("modalRfSub");
      const mCoph = document.getElementById("modalCopheneticScore");
      const mCophSub = document.getElementById("modalCopheneticSub");
      const mDisc = document.getElementById("modalDiscordanceScore");
      const mDiscSub = document.getElementById("modalDiscordanceSub");
      const mInterp = document.getElementById("modalInterpretationText");
      const mTable = document.getElementById("modalDiscordantTableBody");

      if (mTitle) {
        mTitle.textContent = `${activeCong.leftTitle} vs ${activeCong.rightTitle}`;
      }
      if (mSub) {
        const filterStr = activeCong.filterDesc ? ` • Filter: ${activeCong.filterDesc}` : "";
        mSub.textContent = `Comparing phylogenetic congruence across ${activeCong.taxa_count} taxa${filterStr}`;
      }
      if (mRf) mRf.textContent = `${rfPct.toFixed(1)}%`;
      if (mRfSub) mRfSub.textContent = `${sharedSplits} shared clades out of ${Math.round(maxRf/2)} internal bipartitions (RF distance: ${rfDist})`;
      if (mCoph) mCoph.textContent = `r = ${copheneticR >= 0 ? "+" : ""}${copheneticR.toFixed(3)}`;
      if (mCophSub) mCophSub.textContent = `Pearson correlation on pairwise patristic branch substitution distances`;
      if (mDisc) mDisc.textContent = `${discordance.toFixed(1)}%`;
      if (mDiscSub) mDiscSub.textContent = `Native planar crossing discordance (reduced to ${untangledDisc.toFixed(1)}% upon untangling)`;
      if (mInterp) mInterp.textContent = interp;

      if (mTable) {
        mTable.innerHTML = "";
        const discordant = activeCong.discordant_taxa || [];
        discordant.forEach((item) => {
          const row = document.createElement("tr");
          row.className = "hover:bg-slate-500/10 cursor-pointer transition";
          row.onclick = () => {
            closeCongruenceModal();
            setLayout("tanglegram");
            selectTaxon(item.id);
            highlightTangleTaxon(item.id);
          };
          row.innerHTML = `
            <td class="py-2 px-3 font-mono font-bold text-sky-400 truncate max-w-[170px]">${item.id}</td>
            <td class="py-2 px-2 text-[var(--text-muted)] truncate max-w-[140px]">${item.category}</td>
            <td class="py-2 px-2 text-center font-mono text-emerald-400">#${item.r1}</td>
            <td class="py-2 px-2 text-center font-mono text-purple-400">#${item.r2}</td>
            <td class="py-2 px-3 text-right font-mono font-bold text-amber-400">&plusmn;${item.shift} ranks</td>
          `;
          mTable.appendChild(row);
        });
      }
    }

    function closeCongruenceModal() {
      const modal = document.getElementById("congruenceModal");
      if (modal) modal.classList.add("hidden");
    }

    
    // =========================================================================
    // MULTIPLE SEQUENCE ALIGNMENT (MSA) VIEWER ENGINE
    // =========================================================================
    const msaState = {
      mode: "aa", // "aa" or "3di"
      colorScheme: "clustal", // "clustal", "zappo", "hydro", "identity", "foldstate"
      sortOrder: "tree", // "tree", "selected", "alpha"
      syncWithTree: true, // Synced with tree selection by default
      showMinimap: true, // Docked 2D whole-alignment radar overview
      scrollX: 0, // Column offset
      scrollY: 0, // Row offset
      cellWidth: 16,
      cellHeight: 20,
      heightMode: "standard", // "compact", "standard", "tall"
      isMinimized: false,
      hoveredCol: null,
      hoveredRow: null,
      stripGapCols: true, // Automatically strip all-gap columns in filtered / scoped sets
      activeKeptCols: null, // Cache of currently visible column indices
      _minimapCacheKey: null,
      _minimapOffscreen: null,
      coverageThreshold: 0.70, // Default 70% alignment coverage threshold
      filterCoverage: false,   // Whether to filter out sequences below coverage threshold
      activePartition: "all"   // "all" or specific partition ID (e.g. "cluster_1_cov70")
    };

    // High-DPI (Retina) Canvas Initialization Helper
    function initHighDpiCanvas(canvas) {
      if (!canvas) return null;
      const parent = canvas.parentElement;
      if (!parent) return null;
      const rect = parent.getBoundingClientRect();
      const dpr = Math.max(1, window.devicePixelRatio || 1);
      const w = Math.max(1, Math.floor(rect.width));
      const h = Math.max(1, Math.floor(rect.height));

      const targetW = Math.round(w * dpr);
      const targetH = Math.round(h * dpr);

      if (canvas.width !== targetW || canvas.height !== targetH) {
        canvas.width = targetW;
        canvas.height = targetH;
        canvas.style.width = `${w}px`;
        canvas.style.height = `${h}px`;
      }

      const ctx = canvas.getContext("2d");
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.scale(dpr, dpr);
      ctx.clearRect(0, 0, w, h);
      ctx.imageSmoothingEnabled = true;
      ctx.imageSmoothingQuality = "high";
      return { ctx, width: w, height: h, dpr };
    }

    // High-Contrast Text Color Calculator
    function getContrastTextColor(hexBg) {
      if (!hexBg || hexBg.length < 7) return "#ffffff";
      const r = parseInt(hexBg.slice(1, 3), 16) || 0;
      const g = parseInt(hexBg.slice(3, 5), 16) || 0;
      const b = parseInt(hexBg.slice(5, 7), 16) || 0;
      const lum = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
      return lum > 0.58 ? "#0f172a" : "#ffffff";
    }

    function toggleMsaSync() {
      msaState.syncWithTree = !msaState.syncWithTree;
      const btn = document.getElementById("btnMsaSync");
      const icon = document.getElementById("msaSyncIcon");
      const label = document.getElementById("msaSyncLabel");
      if (msaState.syncWithTree) {
        if (btn) btn.className = "badge-sky px-2 py-0.5 rounded text-[10px] font-semibold transition flex items-center space-x-1 cursor-pointer hover:opacity-90";
        if (icon) icon.textContent = "🔗";
        if (label) label.textContent = "Synced";
        if (settings.selectedTaxon) {
          const taxa = getMsaTaxaList();
          const idx = taxa.indexOf(settings.selectedTaxon);
          if (idx >= 0) {
            msaState.scrollY = Math.max(0, idx - 2);
            renderMsa();
          }
        }
      } else {
        if (btn) btn.className = "px-2 py-0.5 rounded text-[10px] font-semibold bg-slate-800 hover:bg-slate-700 text-slate-300 border border-[var(--border-color)] transition flex items-center space-x-1 cursor-pointer";
        if (icon) icon.textContent = "🔓";
        if (label) label.textContent = "Decoupled";
      }
    }

    // Standard Amino Acid Color Matrices
    const AA_CLUSTAL_COLORS = {
      'A': '#2563eb', 'I': '#2563eb', 'L': '#2563eb', 'M': '#2563eb', 'F': '#2563eb', 'W': '#2563eb', 'V': '#2563eb',
      'K': '#dc2626', 'R': '#dc2626',
      'D': '#db2777', 'E': '#db2777',
      'S': '#059669', 'T': '#059669', 'N': '#059669', 'Q': '#059669',
      'C': '#e11d48',
      'G': '#ea580c',
      'P': '#ca8a04',
      'H': '#0891b2', 'Y': '#0891b2',
      '-': '#1e293b'
    };

    const AA_ZAPPO_COLORS = {
      'I': '#fb7185', 'L': '#fb7185', 'V': '#fb7185', 'A': '#fb7185', 'M': '#fb7185', // Aliphatic/Hydrophobic (Rose)
      'F': '#f59e0b', 'W': '#f59e0b', 'Y': '#f59e0b', // Aromatic (Amber)
      'K': '#3b82f6', 'R': '#3b82f6', 'H': '#0284c7', // Basic/Positive (Blue)
      'D': '#ef4444', 'E': '#dc2626', // Acidic/Negative (Red)
      'S': '#10b981', 'T': '#059669', 'N': '#14b8a6', 'Q': '#0d9488', // Polar/Conformation (Emerald/Teal)
      'P': '#eab308', // Proline (Gold)
      'G': '#f97316', // Glycine (Orange)
      'C': '#ec4899', // Cysteine (Pink)
      '-': '#94a3b8'
    };

    const AA_HYDRO_COLORS = {
      'I': '#1e3a8a', 'V': '#1d4ed8', 'L': '#2563eb', 'F': '#3b82f6', 'C': '#60a5fa',
      'M': '#93c5fd', 'A': '#bfdbfe', 'W': '#dbeafe', 'G': '#64748b', 'T': '#e2e8f0',
      'S': '#f1f5f9', 'Y': '#f8fafc', 'P': '#fed7aa', 'H': '#fdba74', 'N': '#fb923c',
      'D': '#f97316', 'Q': '#ea580c', 'E': '#c2410c', 'K': '#9a3412', 'R': '#7c2d12',
      '-': '#1e293b'
    };

    // FoldMason 3Di Structural Alphabet Colors (20 tertiary conformation states)
    const STRUCT_3DI_COLORS = {
      'a': '#0284c7', 'b': '#0369a1', 'c': '#0ea5e9', 'd': '#38bdf8', // Alpha-helix core
      'e': '#d97706', 'f': '#b45309', 'g': '#f59e0b', 'h': '#fbbf24', 'i': '#fde68a', // Beta-sheet strands
      'k': '#16a34a', 'l': '#15803d', 'm': '#22c55e', 'n': '#4ade80', // Turns & sharp bends
      'p': '#9333ea', 'q': '#7e22ce', 'r': '#a855f7', 's': '#c084fc', 't': '#e9d5ff', 'v': '#6b21a8', // Loops & coils
      '-': '#1e293b'
    };

    function getActiveAlignment() {
      if (!window.ALIGNMENTS_DATA || !window.ALIGNMENTS_DATA[currentScale]) return null;
      const base = window.ALIGNMENTS_DATA[currentScale];
      if (msaState.activePartition && msaState.activePartition !== "all" && base.partitions && base.partitions[msaState.activePartition]) {
        return base.partitions[msaState.activePartition];
      }
      return base;
    }

    function getSequenceCoverage(taxon, align, mode) {
      if (!align) align = getActiveAlignment();
      if (!align) return 0.0;
      if (align.coverage && typeof align.coverage[taxon] === "number") {
        return align.coverage[taxon];
      }
      const seqs = align[mode || msaState.mode] || {};
      const seq = seqs[taxon];
      if (!seq || seq.length === 0) return 0.0;
      let nonGaps = 0;
      for (let i = 0; i < seq.length; i++) {
        if (seq[i] !== "-") nonGaps++;
      }
      return nonGaps / seq.length;
    }

    function setMsaCoverageThreshold(val) {
      const num = parseFloat(val);
      if (!isNaN(num)) {
        msaState.coverageThreshold = Math.max(0.01, Math.min(1.0, num / 100));
        const btn = document.getElementById("btnFilterCoverage");
        if (btn && msaState.filterCoverage) {
          btn.textContent = `≥${Math.round(msaState.coverageThreshold * 100)}% ON`;
        }
        if (msaState.filterCoverage) {
          applyTaxaFilter();
        } else {
          renderMsa();
        }
      }
    }

    function toggleFilterCoverage() {
      msaState.filterCoverage = !msaState.filterCoverage;
      const btn = document.getElementById("btnFilterCoverage");
      if (btn) {
        if (msaState.filterCoverage) {
          btn.className = "px-1.5 py-0.5 rounded text-[9.5px] font-bold bg-amber-500 text-white shadow-sm cursor-pointer";
          btn.textContent = `≥${Math.round(msaState.coverageThreshold * 100)}% ON`;
          btn.title = `Alignment Coverage Filter Active (≥${Math.round(msaState.coverageThreshold * 100)}%): Tree and MSA filtered`;
        } else {
          btn.className = "badge-amber px-1.5 py-0.5 rounded text-[9.5px] font-semibold cursor-pointer";
          btn.textContent = "Filter";
          btn.title = "Toggle Alignment Coverage Filter on Tree & MSA";
        }
      }
      applyTaxaFilter();
      if (typeof showToastNotification === 'function') {
        showToastNotification(msaState.filterCoverage
          ? `🎯 Tree & MSA Filter: Showing taxa with ≥${Math.round(msaState.coverageThreshold * 100)}% alignment coverage`
          : `Showing all taxa (Coverage filter disabled)`);
      }
    }

    function setMsaPartition(val) {
      msaState.activePartition = val || "all";
      msaState.scrollX = 0;
      msaState.scrollY = 0;
      msaState._minimapCacheKey = null;
      applyTaxaFilter();
      if (typeof showToastNotification === 'function' && val && val !== "all") {
        const align = (typeof getActiveAlignment === 'function') ? getActiveAlignment() : null;
        const pName = (align && align.name) ? align.name : val;
        showToastNotification(`🔀 Tree & MSA synced to alignment partition: <strong>${pName}</strong>`);
      }
    }

    function getMsaTaxaList() {
      const align = getActiveAlignment();
      if (!align) return [];

      const seqs = align[msaState.mode] || {};
      let taxa = Object.keys(seqs);

      // Filter by active metadata filter
      if (filterState.isActive) {
        const allowed = getFilteredTaxaSet();
        taxa = taxa.filter(t => allowed.has(t));
      }

      // Filter by active scoped clade
      if (settings.scopedClade && settings.scopedClade.taxa) {
        taxa = taxa.filter(t => settings.scopedClade.taxa.has(t));
      }

      // Filter by alignment coverage threshold
      if (msaState.filterCoverage) {
        const minCov = msaState.coverageThreshold || 0.70;
        taxa = taxa.filter(t => getSequenceCoverage(t, align, msaState.mode) >= minCov);
      }

      if (msaState.sortOrder === "tree") {
        // Follow top-to-bottom vertical leaf order in tree
        const leaves = getAllLeaves(activeTreeRoot);
        const order = {};
        leaves.forEach((l, idx) => { order[l.name] = idx; });
        taxa.sort((a, b) => {
          const idxA = order[a] !== undefined ? order[a] : 99999;
          const idxB = order[b] !== undefined ? order[b] : 99999;
          return idxA - idxB;
        });
      } else if (msaState.sortOrder === "selected") {
        taxa.sort((a, b) => {
          if (a === settings.selectedTaxon) return -1;
          if (b === settings.selectedTaxon) return 1;
          return a.localeCompare(b);
        });
      } else {
        taxa.sort((a, b) => a.localeCompare(b));
      }

      return taxa;
    }

    function toggleMsaStripGaps() {
      msaState.stripGapCols = !msaState.stripGapCols;
      const btn = document.getElementById("btnMsaStripGaps");
      const lbl = document.getElementById("msaStripGapsLabel");
      if (btn && lbl) {
        if (msaState.stripGapCols) {
          btn.className = "badge-emerald px-2 py-1 rounded text-[10px] font-semibold transition flex items-center space-x-1 cursor-pointer whitespace-nowrap hover:opacity-90";
          lbl.textContent = "Strip Gaps: ON";
        } else {
          btn.className = "px-2 py-1 rounded text-[10px] font-semibold bg-slate-800 hover:bg-slate-700 text-slate-400 border border-[var(--border-color)] transition flex items-center space-x-1 cursor-pointer whitespace-nowrap";
          lbl.textContent = "Strip Gaps: OFF";
        }
      }
      msaState._minimapCacheKey = null;
      renderMsa();
    }

    function getMsaKeptColumns(taxa, align) {
      if (!align) return [];
      const origLen = align.length || 533;
      if (!msaState.stripGapCols || !taxa || taxa.length === 0) {
        const all = new Array(origLen);
        for (let i = 0; i < origLen; i++) all[i] = i;
        return all;
      }

      const seqs = align[msaState.mode] || {};
      const kept = [];
      for (let c = 0; c < origLen; c++) {
        let hasRes = false;
        for (let i = 0; i < taxa.length; i++) {
          const s = seqs[taxa[i]];
          if (s && c < s.length) {
            const ch = s[c];
            if (ch !== '-' && ch !== '.' && ch !== ' ' && ch !== '?') {
              hasRes = true;
              break;
            }
          }
        }
        if (hasRes) kept.push(c);
      }
      return kept;
    }

    function computeColumnConsensus(taxaList, keptCols) {
      const align = getActiveAlignment();
      if (!align || taxaList.length === 0) return { consensus: "", conservation: [] };

      const seqs = align[msaState.mode] || {};
      let consensus = "";
      const conservation = [];

      for (let idx = 0; idx < keptCols.length; idx++) {
        const c = keptCols[idx];
        const counts = {};
        let totalNonGap = 0;

        for (let i = 0; i < taxaList.length; i++) {
          const seq = seqs[taxaList[i]];
          if (seq && c < seq.length) {
            const ch = seq[c];
            if (ch !== '-' && ch !== '.' && ch !== ' ' && ch !== '?') {
              counts[ch] = (counts[ch] || 0) + 1;
              totalNonGap++;
            }
          }
        }

        let maxChar = "-";
        let maxCount = 0;
        for (let ch in counts) {
          if (counts[ch] > maxCount) {
            maxCount = counts[ch];
            maxChar = ch;
          }
        }

        consensus += maxChar;
        const score = taxaList.length > 0 ? (maxCount / taxaList.length) : 0;
        conservation.push(score);
      }

      return { consensus, conservation };
    }

    function setMsaMode(mode) {
      msaState.mode = mode;
      const btnAA = document.getElementById("btnMsaModeAA");
      const btn3Di = document.getElementById("btnMsaMode3Di");
      const selColor = document.getElementById("msaColorSelect");

      const prevScheme = msaState.colorScheme;
      if (mode === "aa") {
        if (btnAA) btnAA.className = "px-2 py-0.5 rounded font-bold bg-sky-500 text-white transition shadow-sm cursor-pointer";
        if (btn3Di) btn3Di.className = "px-2 py-0.5 rounded font-medium text-[var(--text-muted)] hover:text-white transition cursor-pointer";
        if (selColor) {
          const valid = ["clustal", "custom", "identity", "zappo", "hydro", "taxon"];
          const targetScheme = valid.includes(prevScheme) ? prevScheme : "clustal";
          selColor.innerHTML = `
            <option value="clustal"${targetScheme === "clustal" ? " selected" : ""}>🎨 ClustalX (Classic)</option>
            <option value="custom"${targetScheme === "custom" ? " selected" : ""}>🌿 Custom Palette (Chemistry)</option>
            <option value="identity"${targetScheme === "identity" ? " selected" : ""}>🎯 Conservation Identity</option>
            <option value="zappo"${targetScheme === "zappo" ? " selected" : ""}>🌈 Zappo (Physicochemical)</option>
            <option value="hydro"${targetScheme === "hydro" ? " selected" : ""}>💧 Hydrophobicity</option>
            <option value="taxon"${targetScheme === "taxon" ? " selected" : ""}>🏷️ Taxon Metadata Color</option>
          `;
          msaState.colorScheme = targetScheme;
        }
      } else {
        if (btn3Di) btn3Di.className = "px-2 py-0.5 rounded font-bold bg-purple-500 text-white transition shadow-sm cursor-pointer";
        if (btnAA) btnAA.className = "px-2 py-0.5 rounded font-medium text-[var(--text-muted)] hover:text-white transition cursor-pointer";
        if (selColor) {
          const valid = ["foldstate", "custom", "identity", "taxon"];
          const targetScheme = valid.includes(prevScheme) ? prevScheme : "foldstate";
          selColor.innerHTML = `
            <option value="foldstate"${targetScheme === "foldstate" ? " selected" : ""}>🧊 3Di Secondary Structure</option>
            <option value="custom"${targetScheme === "custom" ? " selected" : ""}>🌿 Custom Palette (3Di Geometry)</option>
            <option value="identity"${targetScheme === "identity" ? " selected" : ""}>🎯 Conservation Identity</option>
            <option value="taxon"${targetScheme === "taxon" ? " selected" : ""}>🏷️ Taxon Metadata Color</option>
          `;
          msaState.colorScheme = targetScheme;
        }
      }

      msaState._minimapCacheKey = null;
      renderMsa();
    }

    function setMsaColorScheme(cs) {
      msaState.colorScheme = cs;
      msaState._minimapCacheKey = null;
      renderMsa();
    }

    function setMsaSort(so) {
      msaState.sortOrder = so;
      renderMsa();
    }

    function jumpMsaToPosition(col) {
      const c = parseInt(col, 10);
      const align = getActiveAlignment();
      const maxCol = (msaState.activeKeptCols && msaState.activeKeptCols.length) ? msaState.activeKeptCols.length : (align ? align.length : 533);
      if (!isNaN(c) && c >= 1 && c <= maxCol) {
        msaState.scrollX = Math.max(0, c - 1);
        renderMsa();
      }
    }

    function zoomMsa(delta) {
      msaState.cellWidth = Math.max(10, Math.min(32, msaState.cellWidth + delta));
      msaState.cellHeight = Math.round(msaState.cellWidth * 1.25);
      const badge = document.getElementById("msaZoomBadge");
      if (badge) badge.textContent = `${msaState.cellWidth}px`;
      renderMsa();
    }

    function toggleMsaHeight() {
      const drawer = document.getElementById("alignmentDrawer");
      const icon = document.getElementById("msaHeightIcon");
      if (!drawer) return;

      if (msaState.heightMode === "standard") {
        msaState.heightMode = "tall";
        drawer.style.height = "380px";
        if (icon) icon.textContent = "↕ 380px";
      } else if (msaState.heightMode === "tall") {
        msaState.heightMode = "compact";
        drawer.style.height = "160px";
        if (icon) icon.textContent = "↕ 160px";
      } else {
        msaState.heightMode = "standard";
        drawer.style.height = "240px";
        if (icon) icon.textContent = "↕ 240px";
      }

      renderMsa();
      setTimeout(() => updateViewportTransform(false), 150);
    }

    function toggleMsaMinimap() {
      msaState.showMinimap = !msaState.showMinimap;
      const panel = document.getElementById("msaMinimapContainer");
      const btn = document.getElementById("btnMsaMinimapToggle");
      if (panel) {
        if (msaState.showMinimap) {
          panel.classList.remove("hidden");
          panel.classList.add("flex");
        } else {
          panel.classList.add("hidden");
          panel.classList.remove("flex");
        }
      }
      if (btn) {
        if (msaState.showMinimap) {
          btn.className = "badge-sky px-2 py-1 rounded text-[10px] font-semibold transition flex items-center space-x-1 cursor-pointer whitespace-nowrap hover:opacity-90";
        } else {
          btn.className = "px-2 py-1 rounded text-[10px] font-semibold bg-slate-800 hover:bg-slate-700 text-slate-300 border border-[var(--border-color)] transition flex items-center space-x-1 cursor-pointer whitespace-nowrap";
        }
      }
      renderMsa();
    }

    function toggleMsaDrawer() {
      const drawer = document.getElementById("alignmentDrawer");
      const icon = document.getElementById("msaDrawerIcon");
      const label = document.getElementById("msaDrawerLabel");
      if (!drawer) return;

      if (!msaState.isMinimized) {
        msaState.prevHeight = drawer.style.height || "240px";
        drawer.style.height = "40px";
        msaState.isMinimized = true;
        if (icon) icon.textContent = "▲";
        if (label) label.textContent = "Exp";
      } else {
        drawer.style.height = msaState.prevHeight;
        msaState.isMinimized = false;
        if (icon) icon.textContent = "▼";
        if (label) label.textContent = "Min";
        renderMsa();
      }

      setTimeout(() => updateViewportTransform(false), 150);
    }

    function onMsaWheel(e) {
      e.preventDefault();
      e.stopPropagation(); // Decoupled: prevent wheel from moving tree
      const align = getActiveAlignment();
      if (!align) return;

      const taxa = getMsaTaxaList();
      const maxCols = (msaState.activeKeptCols && msaState.activeKeptCols.length) ? msaState.activeKeptCols.length : (align.length || 533);
      const maxRows = taxa.length;

      // Handle scrolling on left taxa list specifically
      if (e.currentTarget && e.currentTarget.id === "msaTaxaList") {
        const rowStep = (e.deltaY || e.deltaX) / (msaState.cellHeight * 0.85);
        msaState.scrollY = Math.max(0, Math.min(maxRows - 2, msaState.scrollY + rowStep));
        renderMsa();
        return;
      }

      // Proportional smooth scrolling on matrix
      const dx = e.shiftKey ? e.deltaY : e.deltaX;
      const dy = e.shiftKey ? 0 : e.deltaY;

      if (Math.abs(dx) > 0.4) {
        const colStep = dx / (msaState.cellWidth * 0.85);
        msaState.scrollX = Math.max(0, Math.min(maxCols - 5, msaState.scrollX + colStep));
      }

      if (Math.abs(dy) > 0.4) {
        const rowStep = dy / (msaState.cellHeight * 0.85);
        msaState.scrollY = Math.max(0, Math.min(maxRows - 2, msaState.scrollY + rowStep));
      }

      const posInput = document.getElementById("msaPosInput");
      if (posInput) posInput.value = Math.floor(msaState.scrollX) + 1;

      renderMsa();
    }

    function renderMsa() {
      if (msaState.isMinimized) return;
      const align = getActiveAlignment();
      if (!align) return;

      const taxa = getMsaTaxaList();
      const origAlignLen = align.length || 533;
      const keptCols = getMsaKeptColumns(taxa, align);
      msaState.activeKeptCols = keptCols;
      const alignLen = keptCols.length;
      const strippedCount = origAlignLen - alignLen;
      const seqs = align[msaState.mode] || {};

      if (msaState.scrollX >= alignLen) {
        msaState.scrollX = Math.max(0, alignLen - 5);
      }

      // Update Summary Header Badge
      const badge = document.getElementById("msaSummaryBadge");
      const visCount = document.getElementById("msaVisibleCount");
      const maxColLbl = document.getElementById("msaMaxColLabel");
      const posInput = document.getElementById("msaPosInput");

      // Update Partition Selector Options if partitions exist
      const partSelect = document.getElementById("selectAlignmentPartition");
      if (partSelect) {
        const baseAlign = window.ALIGNMENTS_DATA ? window.ALIGNMENTS_DATA[currentScale] : null;
        const parts = baseAlign && baseAlign.partitions ? baseAlign.partitions : {};
        const partKeys = Object.keys(parts);
        const currentVal = msaState.activePartition || "all";
        let html = `<option value="all"${currentVal === "all" ? " selected" : ""}>📦 Full Alignment (${baseAlign ? baseAlign.taxa_count : 0})</option>`;
        partKeys.forEach(k => {
          const p = parts[k];
          html += `<option value="${k}"${currentVal === k ? " selected" : ""}>🔀 ${p.name || k}</option>`;
        });
        if (partSelect.getAttribute("data-scale") !== currentScale || partSelect.options.length !== (partKeys.length + 1)) {
          partSelect.innerHTML = html;
          partSelect.setAttribute("data-scale", currentScale);
        }
      }

      if (badge) {
        const stripNote = strippedCount > 0 ? ` &bull; <span class="text-emerald-300 font-semibold">${strippedCount} gap-only cols removed</span>` : "";
        const covFilterNote = msaState.filterCoverage ? ` &bull; <span class="text-amber-400 font-semibold">≥${Math.round(msaState.coverageThreshold * 100)}% cov</span>` : "";
        badge.innerHTML = `${taxa.length.toLocaleString()} seqs &bull; ${alignLen} cols${stripNote}${covFilterNote} (${msaState.mode.toUpperCase()})`;
      }
      if (visCount) visCount.textContent = taxa.length.toLocaleString();
      if (maxColLbl) maxColLbl.textContent = `/ ${alignLen}`;
      if (posInput) posInput.max = alignLen;

      // Compute Column Consensus & Conservation over kept columns
      const { consensus, conservation } = computeColumnConsensus(taxa, keptCols);

      // Average Conservation Badge
      const avgScore = conservation.length > 0 ? (conservation.reduce((a, b) => a + b, 0) / conservation.length * 100).toFixed(1) : 0;
      const consBadge = document.getElementById("msaConsensusScore");
      if (consBadge) consBadge.textContent = `${avgScore}% Avg`;

      // 1. RENDER LEFT TAXA COLUMN
      const taxaListEl = document.getElementById("msaTaxaList");
      const isDark = isDarkTheme();
      if (taxaListEl) {
        taxaListEl.innerHTML = "";
        const visibleRowCount = Math.ceil(taxaListEl.clientHeight / msaState.cellHeight) + 1;
        const startR = Math.floor(msaState.scrollY);
        const endR = Math.min(taxa.length, startR + visibleRowCount);

        for (let i = startR; i < endR; i++) {
          const tName = taxa[i];
          const m = TAXA_METADATA[tName] || {};
          const isSelected = (tName === settings.selectedTaxon);

          // Category color dot: honors custom palette, colorColumn, and categorical metadata
          const dotColor = getNodeColor({ name: tName });

          const rowEl = document.createElement("div");
          rowEl.className = `flex items-center justify-between px-2.5 text-[11px] cursor-pointer truncate select-none border-b border-[var(--border-color)] transition-all ${
            isSelected 
              ? (isDark ? "bg-sky-500/25 text-sky-300 font-bold border-l-2 border-sky-400" : "bg-sky-500/20 text-sky-800 font-bold border-l-2 border-sky-500") 
              : "hover:bg-slate-500/10 text-[var(--text-main)]"
          }`;
          rowEl.style.height = `${msaState.cellHeight}px`;
          rowEl.style.lineHeight = `${msaState.cellHeight}px`;
          rowEl.title = `${tName} - Click to select in tree`;

          rowEl.setAttribute("data-taxon", tName);
          rowEl.innerHTML = `
            <div class="flex items-center space-x-1.5 truncate">
              <span class="w-2 h-2 rounded-full shrink-0 shadow-sm" style="background-color: ${dotColor}"></span>
              <span class="truncate font-mono font-medium" title="${tName}">${getLeafLabelText(tName)}</span>
            </div>
            <div class="flex items-center space-x-1 shrink-0">
              <span class="text-[8px] font-mono px-1 rounded ${Math.round(getSequenceCoverage(tName, align, msaState.mode) * 100) >= Math.round(msaState.coverageThreshold * 100) ? "text-emerald-400 bg-emerald-500/15 border border-emerald-500/30" : "text-amber-400 bg-amber-500/15 border border-amber-500/30"}" title="Alignment Coverage: ${Math.round(getSequenceCoverage(tName, align, msaState.mode) * 100)}%">${Math.round(getSequenceCoverage(tName, align, msaState.mode) * 100)}%</span>
              <span class="text-[8.5px] font-mono text-[var(--text-muted)]">#${i + 1}</span>
            </div>
          `;

          // DECOUPLED SELECTION: Click selects taxon without moving/panning the tree canvas
          rowEl.onclick = (e) => {
            e.stopPropagation();
            selectTaxon(tName, false); // false = do not move tree
            highlightTangleTaxon(tName);
            renderMsa();
          };
          rowEl.onmouseenter = () => highlightTangleTaxon(tName);
          rowEl.onmouseleave = () => clearTangleHighlight();

          taxaListEl.appendChild(rowEl);
        }
      }

      // 2. RENDER RULER CANVAS (Retina High-DPI Crisp Vector Rendering)
      const rulerCanvas = document.getElementById("msaRulerCanvas");
      const rSetup = initHighDpiCanvas(rulerCanvas);
      if (rSetup) {
        const { ctx: rctx, width: rWidth } = rSetup;
        const cW = msaState.cellWidth;
        const startC = Math.floor(msaState.scrollX);
        const colCount = Math.ceil(rWidth / cW) + 1;
        const endC = Math.min(alignLen, startC + colCount);

        // Theme-aware ruler fill and text
        rctx.fillStyle = isDark ? "rgba(15, 23, 42, 0.4)" : "rgba(241, 245, 249, 0.95)";
        rctx.fillRect(0, 0, rWidth, 24);

        rctx.fillStyle = isDark ? "#94a3b8" : "#0f172a";
        rctx.font = "600 10px ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Monaco, Consolas, monospace";
        rctx.textBaseline = "middle";

        for (let c = startC; c < endC; c++) {
          const origCol = keptCols[c];
          const x = (c - msaState.scrollX) * cW;
          const posNum = origCol + 1;

          if (posNum % 10 === 0 || posNum === 1) {
            rctx.strokeStyle = isDark ? "#475569" : "#64748b";
            rctx.lineWidth = 1;
            rctx.beginPath();
            const tickX = Math.floor(x) + 0.5;
            rctx.moveTo(tickX, 13);
            rctx.lineTo(tickX, 24);
            rctx.stroke();
            rctx.fillText(String(posNum), Math.floor(x) + 3, 7);
          } else if (posNum % 5 === 0) {
            rctx.strokeStyle = isDark ? "#334155" : "#94a3b8";
            rctx.lineWidth = 1;
            rctx.beginPath();
            const tickX = Math.floor(x) + 0.5;
            rctx.moveTo(tickX, 17);
            rctx.lineTo(tickX, 24);
            rctx.stroke();
          }
        }
      }

      // 3. RENDER RESIDUE MATRIX CANVAS (Theme-Aware, Publication-Grade Clean Rendering)
      const matrixCanvas = document.getElementById("msaMatrixCanvas");
      const mSetup = initHighDpiCanvas(matrixCanvas);
      if (mSetup) {
        const { ctx: mctx, width: mWidth, height: mHeight } = mSetup;
        const cW = msaState.cellWidth;
        const cH = msaState.cellHeight;

        // Clean Canvas Background (No harsh black margins)
        const canvasBg = isDark ? "#090d16" : "#ffffff";
        mctx.fillStyle = canvasBg;
        mctx.fillRect(0, 0, mWidth, mHeight);

        const startC = Math.floor(msaState.scrollX);
        const colCount = Math.ceil(mWidth / cW) + 1;
        const endC = Math.min(alignLen, startC + colCount);

        const startR = Math.floor(msaState.scrollY);
        const rowCount = Math.ceil(mHeight / cH) + 1;
        const endR = Math.min(taxa.length, startR + rowCount);

        const fontSize = Math.max(9, Math.min(14, cW - 3));
        mctx.textAlign = "center";
        mctx.textBaseline = "middle";

        const customPal = (customPaletteState.isActive && customPaletteState.palette.length > 0)
          ? customPaletteState.palette
          : ["#B9554E", "#627B08", "#267567", "#294719", "#72A183"];

        for (let r = startR; r < endR; r++) {
          const tName = taxa[r];
          const seq = seqs[tName] || "";
          const y = (r - msaState.scrollY) * cH;
          const isSelected = (tName === settings.selectedTaxon);
          const ry = Math.floor(y);
          const rh = Math.max(1, Math.floor(cH - 1));

          for (let c = startC; c < endC; c++) {
            const origCol = keptCols[c];
            const x = (c - msaState.scrollX) * cW;
            const ch = (origCol < seq.length) ? seq[origCol] : "-";
            const rx = Math.floor(x);
            const rw = Math.max(1, Math.floor(cW - 1));

            if (ch === "-" || ch === ".") {
              // SUBTLE, CLEAN GAP RENDERING (NO distracting black boxes or thick cages)
              if (cW >= 7) {
                mctx.font = `600 ${Math.max(9, Math.min(13, cW - 2))}px ui-monospace, SFMono-Regular, monospace`;
                mctx.fillStyle = isDark ? "#334155" : "#94a3b8";
                mctx.fillText("–", Math.floor(x + cW / 2), Math.floor(y + cH / 2));
              }
              continue;
            }

            // Calculate Residue Color Scheme
            let bg = isDark ? "#1e293b" : "#f1f5f9";
            if (msaState.colorScheme === "taxon") {
              bg = getNodeColor({ name: tName });
            } else if (msaState.colorScheme === "custom") {
              if (msaState.mode === "3di") {
                const chL = ch.toLowerCase();
                if ("abcd".includes(chL)) bg = customPal[0];
                else if ("efghi".includes(chL)) bg = customPal[1 % customPal.length];
                else if ("klmn".includes(chL)) bg = customPal[2 % customPal.length];
                else if ("pqrstuv".includes(chL)) bg = customPal[3 % customPal.length];
                else bg = customPal[4 % customPal.length];
              } else {
                const chU = ch.toUpperCase();
                if ("AVLIMF W".includes(chU)) bg = customPal[0];
                else if ("KRH".includes(chU)) bg = customPal[1 % customPal.length];
                else if ("DE".includes(chU)) bg = customPal[2 % customPal.length];
                else if ("STNQ".includes(chU)) bg = customPal[3 % customPal.length];
                else bg = customPal[4 % customPal.length];
              }
            } else if (msaState.colorScheme === "identity") {
              const isCons = (ch === consensus[c]);
              const sc = conservation[c] || 0;
              if (isCons && sc >= 0.8) bg = isDark ? "#10b981" : "#059669";
              else if (isCons || sc >= 0.5) bg = isDark ? "#38bdf8" : "#0284c7";
              else bg = isDark ? "#1e293b" : "#e2e8f0";
            } else if (msaState.mode === "3di") {
              bg = STRUCT_3DI_COLORS[ch.toLowerCase()] || "#9333ea";
            } else {
              if (msaState.colorScheme === "zappo") {
                bg = AA_ZAPPO_COLORS[ch.toUpperCase()] || (isDark ? "#1e293b" : "#f1f5f9");
              } else if (msaState.colorScheme === "hydro") {
                bg = AA_HYDRO_COLORS[ch.toUpperCase()] || (isDark ? "#1e293b" : "#f1f5f9");
              } else {
                bg = AA_CLUSTAL_COLORS[ch.toUpperCase()] || (isDark ? "#1e293b" : "#f1f5f9");
              }
            }

            mctx.fillStyle = bg;
            if (cW >= 14 && typeof mctx.roundRect === "function") {
              mctx.beginPath();
              mctx.roundRect(rx, ry, rw, rh, 2.5);
              mctx.fill();
            } else {
              mctx.fillRect(rx, ry, rw, rh);
            }

            // High-Contrast Residue Character Rendering
            if (cW >= 8) {
              mctx.font = `700 ${fontSize}px ui-monospace, SFMono-Regular, "SF Mono", Menlo, Monaco, Consolas, monospace`;
              mctx.fillStyle = getContrastTextColor(bg);
              mctx.fillText(ch, Math.floor(x + cW / 2), Math.floor(y + cH / 2));
            }
          }

          if (isSelected) {
            mctx.strokeStyle = isDark ? "#38bdf8" : "#0284c7";
            mctx.lineWidth = 2;
            mctx.strokeRect(0.5, ry + 0.5, mWidth - 1, rh);
          }
        }
      }

      // 4. RENDER CONSENSUS & CONSERVATION BAR CANVAS (Theme-Aware Contrast)
      const consCanvas = document.getElementById("msaConsensusCanvas");
      const cSetup = initHighDpiCanvas(consCanvas);
      if (cSetup) {
        const { ctx: cctx, width: cWidth } = cSetup;
        const cW = msaState.cellWidth;
        const startC = Math.floor(msaState.scrollX);
        const colCount = Math.ceil(cWidth / cW) + 1;
        const endC = Math.min(alignLen, startC + colCount);

        cctx.fillStyle = isDark ? "rgba(15, 23, 42, 0.4)" : "rgba(241, 245, 249, 0.95)";
        cctx.fillRect(0, 0, cWidth, 28);

        cctx.textAlign = "center";
        cctx.textBaseline = "middle";

        for (let c = startC; c < endC; c++) {
          const x = (c - msaState.scrollX) * cW;
          const ch = consensus[c] || "-";
          const score = conservation[c] || 0;

          // Conservation Bar (bottom half)
          const barH = Math.round(score * 15);
          if (isDark) {
            cctx.fillStyle = score >= 0.8 ? "#10b981" : (score >= 0.5 ? "#38bdf8" : "#f59e0b");
          } else {
            cctx.fillStyle = score >= 0.8 ? "#059669" : (score >= 0.5 ? "#0284c7" : "#d97706");
          }
          cctx.fillRect(Math.floor(x), 28 - barH, Math.max(1, Math.floor(cW - 1)), barH);

          // Consensus Character (top half with 7:1+ contrast)
          if (cW >= 8 && ch !== "-" && ch !== ".") {
            cctx.font = `700 ${Math.max(9, Math.min(13, cW - 3))}px ui-monospace, SFMono-Regular, "SF Mono", Menlo, Monaco, Consolas, monospace`;
            cctx.fillStyle = isDark ? "#f8fafc" : "#0f172a";
            cctx.fillText(ch, Math.floor(x + cW / 2), 7);
          }
        }
      }

      // 5. RENDER DOCKED ALIGNMENT MINIMAP / RADAR
      renderMsaMinimap(taxa, alignLen, consensus, conservation, keptCols);
    }

    // ALIGNMENT MINIMAP 2D OVERVIEW RENDERER
    function renderMsaMinimap(taxa, alignLen, consensus, conservation, keptCols) {
      if (!msaState.showMinimap) return;
      const canvas = document.getElementById("msaMinimapCanvas");
      const vpBox = document.getElementById("msaMinimapViewport");
      const infoBadge = document.getElementById("msaMinimapInfo");
      if (!canvas || !vpBox) return;

      const miniSetup = initHighDpiCanvas(canvas);
      if (!miniSetup) return;
      const { ctx: miniCtx, width: mw, height: mh } = miniSetup;
      if (mw <= 0 || mh <= 0 || !taxa || taxa.length === 0 || alignLen <= 0) return;

      const stripKey = msaState.stripGapCols ? "strip" : "raw";
      const isDark = isDarkTheme();
      const palSig = (msaState.colorScheme === "custom" && customPaletteState.isActive)
        ? customPaletteState.palette.join("")
        : msaState.colorScheme;
      const themeSig = isDark ? "dark" : "light";
      const cacheKey = `${currentScale}_${msaState.mode}_${msaState.sortOrder}_${taxa.length}_${alignLen}_${stripKey}_${themeSig}_${palSig}_${mw}_${mh}`;

      if (msaState._minimapCacheKey !== cacheKey || !msaState._minimapOffscreen) {
        const offscreen = document.createElement("canvas");
        offscreen.width = mw;
        offscreen.height = mh;
        const octx = offscreen.getContext("2d");
        const imgData = octx.createImageData(mw, mh);
        const data = imgData.data;

        const align = getActiveAlignment();
        const seqs = align ? (align[msaState.mode] || {}) : {};
        const customPal = (customPaletteState.isActive && customPaletteState.palette.length > 0)
          ? customPaletteState.palette
          : ["#B9554E", "#627B08", "#267567", "#294719", "#72A183"];

        for (let py = 0; py < mh; py++) {
          const r = Math.min(taxa.length - 1, Math.floor((py / mh) * taxa.length));
          const tName = taxa[r];
          const seq = seqs[tName] || "";

          for (let px = 0; px < mw; px++) {
            const c = Math.min(alignLen - 1, Math.floor((px / mw) * alignLen));
            const origCol = (keptCols && c < keptCols.length) ? keptCols[c] : c;
            const idx = (py * mw + px) * 4;

            if (origCol >= seq.length || seq[origCol] === "-" || seq[origCol] === ".") {
              if (isDark) {
                data[idx] = 15;     // Gap dark: deep slate #0f172a
                data[idx + 1] = 23;
                data[idx + 2] = 42;
                data[idx + 3] = 255;
              } else {
                data[idx] = 241;    // Gap light: soft paper slate #f1f5f9
                data[idx + 1] = 245;
                data[idx + 2] = 249;
                data[idx + 3] = 255;
              }
            } else {
              const ch = seq[origCol];
              const isCons = (ch === consensus[c]);
              const sc = conservation[c] || 0;

              if (msaState.colorScheme === "custom") {
                let hex = customPal[0];
                if (msaState.mode === "3di") {
                  const chL = ch.toLowerCase();
                  if ("abcd".includes(chL)) hex = customPal[0];
                  else if ("efghi".includes(chL)) hex = customPal[1 % customPal.length];
                  else if ("klmn".includes(chL)) hex = customPal[2 % customPal.length];
                  else if ("pqrstuv".includes(chL)) hex = customPal[3 % customPal.length];
                  else hex = customPal[4 % customPal.length];
                } else {
                  const chU = ch.toUpperCase();
                  if ("AVLIMF W".includes(chU)) hex = customPal[0];
                  else if ("KRH".includes(chU)) hex = customPal[1 % customPal.length];
                  else if ("DE".includes(chU)) hex = customPal[2 % customPal.length];
                  else if ("STNQ".includes(chU)) hex = customPal[3 % customPal.length];
                  else hex = customPal[4 % customPal.length];
                }
                data[idx] = parseInt(hex.slice(1, 3), 16) || 56;
                data[idx + 1] = parseInt(hex.slice(3, 5), 16) || 189;
                data[idx + 2] = parseInt(hex.slice(5, 7), 16) || 248;
                data[idx + 3] = 255;
              } else if (isCons && sc >= 0.8) {
                // Core conservation >= 80%
                data[idx] = isDark ? 16 : 5;
                data[idx + 1] = isDark ? 185 : 150;
                data[idx + 2] = isDark ? 129 : 105;
                data[idx + 3] = 255;
              } else if (isCons || sc >= 0.5) {
                // Moderate conservation >= 50%
                data[idx] = isDark ? 56 : 2;
                data[idx + 1] = isDark ? 189 : 132;
                data[idx + 2] = isDark ? 248 : 199;
                data[idx + 3] = 255;
              } else {
                // Variable / background residue
                data[idx] = isDark ? 71 : 148;
                data[idx + 1] = isDark ? 85 : 163;
                data[idx + 2] = isDark ? 105 : 184;
                data[idx + 3] = 255;
              }
            }
          }
        }

        octx.putImageData(imgData, 0, 0);
        msaState._minimapOffscreen = offscreen;
        msaState._minimapCacheKey = cacheKey;
      }

      miniCtx.drawImage(msaState._minimapOffscreen, 0, 0, mw, mh);

      // Synchronize Draggable Viewport Box
      const matrixCanvas = document.getElementById("msaMatrixCanvas");
      const matW = matrixCanvas ? matrixCanvas.clientWidth : 600;
      const matH = matrixCanvas ? matrixCanvas.clientHeight : 160;

      const visCols = Math.max(1, matW / msaState.cellWidth);
      const visRows = Math.max(1, matH / msaState.cellHeight);

      const vx = (msaState.scrollX / alignLen) * mw;
      const vy = (msaState.scrollY / taxa.length) * mh;
      const vw = Math.max(6, Math.min(mw, (visCols / alignLen) * mw));
      const vh = Math.max(6, Math.min(mh, (visRows / taxa.length) * mh));

      const boundedX = Math.max(0, Math.min(mw - vw, vx));
      const boundedY = Math.max(0, Math.min(mh - vh, vy));

      vpBox.style.left = `${Math.round(boundedX)}px`;
      vpBox.style.top = `${Math.round(boundedY)}px`;
      vpBox.style.width = `${Math.round(vw)}px`;
      vpBox.style.height = `${Math.round(vh)}px`;
      vpBox.style.display = "block";

      if (infoBadge) {
        const startC = Math.floor(msaState.scrollX) + 1;
        const endC = Math.min(alignLen, Math.floor(msaState.scrollX + visCols));
        infoBadge.textContent = `${startC}–${endC}`;
      }
    }

        // CARTESIAN LAYOUT (Phylogram & Cladogram with Triangular Clades)
    function renderCartesianTree(g, visibleLeaves, maxDepth, defs) {
      const spacing = settings.verticalSpacing;
      const xSpan = (visibleLeaves.length > 50) ? 650 : 480;

      visibleLeaves.forEach((leaf, idx) => {
        leaf.y = 50 + idx * spacing;
      });

      function computeInternalCoords(node) {
        if (!node.children || node.children.length === 0 || node._collapsed) return;
        node.children.forEach(computeInternalCoords);
        const validChildY = node.children.map(c => c.y).filter(y => typeof y === 'number' && !isNaN(y));
        if (validChildY.length > 0) {
          node.y = validChildY.reduce((acc, y) => acc + y, 0) / validChildY.length;
        } else {
          node.y = 50.0;
        }
      }
      computeInternalCoords(activeTreeRoot);

      function drawBranches(node, currentX) {
        node.x = currentX;

        if (node._collapsed) {
          drawCollapsedTriangle(g, node, maxDepth, xSpan, spacing);
          drawNodePoint(g, node);
          return;
        }

        if (!node.children || node.children.length === 0) return;

        const validChildY = node.children.map(c => c.y).filter(y => typeof y === 'number' && !isNaN(y));
        const minY = validChildY.length > 0 ? Math.min(...validChildY) : node.y;
        const maxY = validChildY.length > 0 ? Math.max(...validChildY) : node.y;

        const vLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
        vLine.setAttribute("x1", node.x);
        vLine.setAttribute("y1", minY);
        vLine.setAttribute("x2", node.x);
        vLine.setAttribute("y2", maxY);
        vLine.setAttribute("class", "branch-path");
        applyVerticalBranchGradient(vLine, node, minY, maxY, defs);
        vLine.addEventListener("mouseenter", (e) => showCladeTooltip(e, node));
        vLine.addEventListener("mouseleave", () => scheduleHideTooltip(250));
        g.appendChild(vLine);

        node.children.forEach(child => {
          const bLen = settings.branchLengths ? (typeof child.length === 'number' ? child.length : 1.0) : 1.0;
          const childX = node.x + (bLen / maxDepth) * xSpan;
          child.x = childX;

          const hLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          hLine.setAttribute("x1", node.x);
          hLine.setAttribute("y1", child.y);
          hLine.setAttribute("x2", childX);
          hLine.setAttribute("y2", child.y);
          hLine.setAttribute("class", "branch-path");
          const bCol = getBranchStrokeColor(child);
          if (bCol) hLine.style.stroke = bCol;

          hLine.addEventListener("mouseenter", (e) => showBranchTooltip(e, child));
          hLine.addEventListener("mouseleave", () => scheduleHideTooltip(250));
          g.appendChild(hLine);

          drawBranches(child, childX);
        });

        drawNodePoint(g, node);
      }

      drawBranches(activeTreeRoot, 50);

      // Render Leaf Labels
      const validX = visibleLeaves.map(l => l.x).filter(x => typeof x === 'number' && !isNaN(x));
      const maxX = validX.length > 0 ? Math.max(...validX) : 200;

      visibleLeaves.forEach(leaf => {
        if (leaf._collapsed) return;

        drawNodePoint(g, leaf);

        const nodeR = typeof settings.nodeRadius === 'number' ? settings.nodeRadius : 3;
        const lx = settings.alignLabels ? maxX + 14 : leaf.x + nodeR + 4;
        const ly = leaf.y;

        if (settings.alignLabels && lx > leaf.x + nodeR + 6) {
          const guideLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          guideLine.setAttribute("x1", leaf.x + nodeR + 3);
          guideLine.setAttribute("y1", leaf.y);
          guideLine.setAttribute("x2", lx - 4);
          guideLine.setAttribute("y2", leaf.y);
          guideLine.setAttribute("stroke", "var(--border-color)");
          guideLine.setAttribute("stroke-dasharray", "2,3");
          guideLine.setAttribute("stroke-width", "0.75");
          guideLine.setAttribute("opacity", "0.45");
          g.appendChild(guideLine);
        }

        const meta = TAXA_METADATA[leaf.name];
        const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt.setAttribute("x", lx);
        txt.setAttribute("y", ly);
        txt.setAttribute("dominant-baseline", "central");
        txt.setAttribute("class", `tip-label ${leaf.name === settings.selectedTaxon ? "selected" : ""} ${getLeafFilterClass(leaf.name)}`);
        txt.setAttribute("data-taxon", leaf.name);
        txt.textContent = getLeafLabelText(leaf.name);
        txt.onclick = (e) => { e.stopPropagation(); selectTaxon(leaf.name); };
        txt.onmouseenter = (e) => showNodeTooltip(e, leaf);
        txt.onmouseleave = () => scheduleHideTooltip(250);
        g.appendChild(txt);

        if (meta && settings.showMeta && visibleLeaves.length <= 65) {
          const colDef = getActiveColorColumnDef();
          const subText = colDef && meta[colDef.key] !== undefined ? `${meta[colDef.key]}` : "";
          if (subText) {
            const sub = document.createElementNS("http://www.w3.org/2000/svg", "text");
            sub.setAttribute("x", lx + 140);
            sub.setAttribute("y", ly);
            sub.setAttribute("dominant-baseline", "central");
            sub.setAttribute("class", "tip-label-sub");
            sub.textContent = subText;
            sub.onclick = (e) => { e.stopPropagation(); selectTaxon(leaf.name); };
            g.appendChild(sub);
          }
        }
      });
    }

    // TRIANGULAR CLADE RENDERING
    function drawCollapsedTriangle(g, node, maxDepth, xSpan, spacing) {
      const leaves = getAllLeaves(node);
      const info = getCladeInfo(node);
      const cladeColor = info.cladeColor || "#38bdf8";

      const validChildDepths = leaves.map(l => settings.branchLengths ? l.depth : l.cladoDepth);
      const minChildDepth = Math.min(...validChildDepths, node.depth || 0);
      const maxChildDepth = Math.max(...validChildDepths, node.depth || 0);

      const deltaDepth = Math.max(maxChildDepth - (node.depth || 0), 0.08);
      const triWidth = Math.max((deltaDepth / maxDepth) * xSpan, 65);

      const triHeight = Math.min(Math.max(leaves.length * spacing * 0.55, 18), 120);

      const apexX = node.x;
      const apexY = node.y;
      const baseX = apexX + triWidth;
      const topY = apexY - triHeight / 2;
      const botY = apexY + triHeight / 2;

      const pathData = `M ${apexX} ${apexY} L ${baseX} ${topY} L ${baseX} ${botY} Z`;

      const wedge = document.createElementNS("http://www.w3.org/2000/svg", "path");
      wedge.setAttribute("d", pathData);
      wedge.setAttribute("class", "clade-wedge");
      wedge.style.fill = cladeColor;
      wedge.style.fillOpacity = "0.45";
      wedge.style.stroke = cladeColor;
      wedge.style.strokeWidth = "1.5px";

      wedge.addEventListener("mouseenter", (e) => showCollapsedTooltip(e, node, info));
      wedge.addEventListener("mouseleave", () => scheduleHideTooltip(250));
      wedge.addEventListener("click", (e) => {
        e.stopPropagation();
        toggleCladeCollapse(node);
      });
      g.appendChild(wedge);

      // Clade Title Label
      const titleTxt = document.createElementNS("http://www.w3.org/2000/svg", "text");
      titleTxt.setAttribute("x", baseX + 10);
      titleTxt.setAttribute("y", apexY - 2);
      titleTxt.setAttribute("class", "clade-summary-label");
      titleTxt.style.fill = cladeColor;
      titleTxt.textContent = `${info.dominantCategory} [${leaves.length} taxa]`;
      titleTxt.onclick = (e) => { e.stopPropagation(); toggleCladeCollapse(node); };
      g.appendChild(titleTxt);

      // Sub-label with homogeneity
      const subTxt = document.createElementNS("http://www.w3.org/2000/svg", "text");
      subTxt.setAttribute("x", baseX + 10);
      subTxt.setAttribute("y", apexY + 12);
      subTxt.setAttribute("class", "clade-summary-sub");
      subTxt.textContent = `${info.homogeneity.toFixed(0)}% ${info.colLabel} &bull; Click to Expand`;
      subTxt.onclick = (e) => { e.stopPropagation(); toggleCladeCollapse(node); };
      g.appendChild(subTxt);
    }

    // RADAR MINIMAP & REFINED NODE POINT RENDERING
    function drawNodePoint(g, node) {
      const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      circle.setAttribute("cx", node.x);
      circle.setAttribute("cy", node.y);

      const isCollapsible = node.children && node.children.length > 0 && node !== activeTreeRoot;
      const radius = node._collapsed ? settings.nodeRadius * 1.3 : settings.nodeRadius;
      circle.setAttribute("r", radius);
      circle.setAttribute("class", `node-dot ${getLeafFilterClass(node.name)}`);
      circle.style.fill = getNodeColor(node);

      if (node._collapsed) {
        circle.style.stroke = isDarkTheme() ? "#ffffff" : "#0f172a";
        circle.style.strokeWidth = "1.0px";
      } else if (isCollapsible) {
        circle.style.stroke = "var(--accent)";
        circle.style.strokeWidth = "0.8px";
      } else {
        circle.style.stroke = isDarkTheme() ? "rgba(255, 255, 255, 0.4)" : "rgba(15, 23, 42, 0.4)";
        circle.style.strokeWidth = "0.65px";
      }

      circle.addEventListener("mouseenter", (e) => {
        if (node._collapsed) {
          showCollapsedTooltip(e, node, getCladeInfo(node));
        } else if (node.name) {
          highlightTangleConnector(node.name);
          showNodeTooltip(e, node);
        } else if (isCollapsible) {
          showCladeTooltip(e, node);
        }
      });
      circle.addEventListener("mouseleave", () => {
        if (node.name) clearTangleHighlight();
        scheduleHideTooltip(250);
      });
      circle.addEventListener("click", (e) => {
        e.stopPropagation();
        if (node._collapsed || isCollapsible) {
          toggleCladeCollapse(node);
        } else if (node.name) {
          selectTaxon(node.name);
        }
      });

      g.appendChild(circle);
    }

    // RADIAL TREE LAYOUT (iTOL / FigTree Enhanced Readability Engine)
    function renderRadialTree(g, visibleLeaves, maxDepth) {
      const totalLeaves = visibleLeaves.length;
      if (totalLeaves === 0) return;

      const centerX = 450;
      const centerY = 450;

      const arcDeg = settings.radialArc !== undefined ? settings.radialArc : 360;
      const rotDeg = settings.treeRotation !== undefined ? settings.treeRotation : 0;
      const startAngle = (rotDeg * Math.PI) / 180;
      const totalArcRad = (arcDeg / 360) * (2 * Math.PI);
      const angleStep = totalArcRad / (totalLeaves > 1 ? (arcDeg === 360 ? totalLeaves : totalLeaves - 1) : 1);

      const radiusScale = typeof settings.radialRadiusScale === 'number' ? settings.radialRadiusScale : 1.0;
      const baseMaxRadius = Math.min(380, 50 + totalLeaves * 2.2);
      const maxRadius = baseMaxRadius * radiusScale;

      visibleLeaves.forEach((leaf, idx) => {
        leaf.angle = startAngle + idx * angleStep;
      });

      function computeInternalAngles(node) {
        if (!node.children || node.children.length === 0 || node._collapsed) return;
        node.children.forEach(computeInternalAngles);
        const childAngles = node.children.map(c => c.angle).filter(a => typeof a === 'number');
        if (childAngles.length > 0) {
          node.angle = childAngles.reduce((acc, a) => acc + a, 0) / childAngles.length;
        } else {
          node.angle = 0;
        }
      }
      computeInternalAngles(activeTreeRoot);

      // 1. Clade Sector Halos (iTOL style background wedges)
      if (settings.cladeSectors && settings.colorColumn !== "solid" && visibleLeaves.length > 1) {
        const colDef = getActiveColorColumnDef();
        if (colDef && colDef.type === "categorical") {
          let currCat = null;
          let segStartAngle = null;
          let segEndAngle = null;
          const sectors = [];

          visibleLeaves.forEach((leaf, idx) => {
            const meta = TAXA_METADATA[leaf.name] || {};
            const cat = meta[colDef.key] || "Unclassified";
            if (currCat === null) {
              currCat = cat;
              segStartAngle = leaf.angle - angleStep / 2;
              segEndAngle = leaf.angle + angleStep / 2;
            } else if (cat === currCat) {
              segEndAngle = leaf.angle + angleStep / 2;
            } else {
              sectors.push({ cat: currCat, start: segStartAngle, end: segEndAngle });
              currCat = cat;
              segStartAngle = leaf.angle - angleStep / 2;
              segEndAngle = leaf.angle + angleStep / 2;
            }
            if (idx === visibleLeaves.length - 1) {
              sectors.push({ cat: currCat, start: segStartAngle, end: segEndAngle });
            }
          });

          const wedgesGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
          wedgesGroup.setAttribute("class", "clade-sectors");
          sectors.forEach(sec => {
            const secColor = (colDef.colors && colDef.colors[sec.cat]) || "#94a3b8";
            const rOuter = maxRadius + 18;
            const rInner = 20;
            const sa = sec.start;
            const ea = sec.end;
            const diffA = ea - sa;
            if (diffA <= 0.001) return;
            const largeArc = diffA > Math.PI ? 1 : 0;

            const p1x = centerX + rInner * Math.cos(sa);
            const p1y = centerY + rInner * Math.sin(sa);
            const p2x = centerX + rOuter * Math.cos(sa);
            const p2y = centerY + rOuter * Math.sin(sa);
            const p3x = centerX + rOuter * Math.cos(ea);
            const p3y = centerY + rOuter * Math.sin(ea);
            const p4x = centerX + rInner * Math.cos(ea);
            const p4y = centerY + rInner * Math.sin(ea);

            const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
            path.setAttribute("d", `M ${p1x} ${p1y} L ${p2x} ${p2y} A ${rOuter} ${rOuter} 0 ${largeArc} 1 ${p3x} ${p3y} L ${p4x} ${p4y} A ${rInner} ${rInner} 0 ${largeArc} 0 ${p1x} ${p1y} Z`);
            path.setAttribute("fill", secColor);
            path.setAttribute("fill-opacity", isDarkTheme() ? "0.09" : "0.07");
            path.setAttribute("stroke", secColor);
            path.setAttribute("stroke-opacity", "0.25");
            path.setAttribute("stroke-width", "0.75px");
            wedgesGroup.appendChild(path);
          });
          g.appendChild(wedgesGroup);
        }
      }

      // 2. Concentric Distance Scale Rings
      if (settings.concentricRings) {
        const ringsGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
        ringsGroup.setAttribute("class", "concentric-rings");
        const ringFractions = [0.25, 0.5, 0.75, 1.0];
        ringFractions.forEach(frac => {
          const r = 20 + (maxRadius - 20) * frac;
          const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
          circle.setAttribute("cx", centerX);
          circle.setAttribute("cy", centerY);
          circle.setAttribute("r", r);
          circle.setAttribute("fill", "none");
          circle.setAttribute("stroke", "var(--border-color)");
          circle.setAttribute("stroke-width", "0.6px");
          circle.setAttribute("stroke-dasharray", "3,3");
          circle.setAttribute("opacity", isDarkTheme() ? "0.35" : "0.45");
          ringsGroup.appendChild(circle);

          if (settings.branchLengths && maxDepth > 0) {
            const distLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
            distLabel.setAttribute("x", centerX + 4);
            distLabel.setAttribute("y", centerY - r - 3);
            distLabel.setAttribute("fill", "var(--text-muted)");
            distLabel.setAttribute("font-size", "8.5px");
            distLabel.setAttribute("font-family", "ui-monospace, monospace");
            distLabel.setAttribute("opacity", "0.75");
            distLabel.textContent = (maxDepth * frac).toFixed(3);
            ringsGroup.appendChild(distLabel);
          }
        });
        g.appendChild(ringsGroup);
      }

      // 3. Draw Radial Branches
      function drawRadialBranches(node, currRadius) {
        node.radius = currRadius;
        node.x = centerX + currRadius * Math.cos(node.angle);
        node.y = centerY + currRadius * Math.sin(node.angle);

        if (node._collapsed) {
          drawNodePoint(g, node);
          return;
        }

        if (!node.children || node.children.length === 0) return;

        node.children.forEach(child => {
          const bLen = settings.branchLengths ? (typeof child.length === 'number' ? child.length : 1.0) : 1.0;
          const childRadius = node.radius + (bLen / Math.max(0.0001, maxDepth)) * (maxRadius - 20);
          child.radius = childRadius;
          child.x = centerX + childRadius * Math.cos(child.angle);
          child.y = centerY + childRadius * Math.sin(child.angle);

          const arcPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
          const midX = centerX + node.radius * Math.cos(child.angle);
          const midY = centerY + node.radius * Math.sin(child.angle);
          const d = `M ${node.x} ${node.y} A ${node.radius} ${node.radius} 0 0 ${child.angle > node.angle ? 1 : 0} ${midX} ${midY} L ${child.x} ${child.y}`;
          arcPath.setAttribute("d", d);
          arcPath.setAttribute("class", "branch-path");
          const bCol = getBranchStrokeColor(child);
          if (bCol) arcPath.style.stroke = bCol;

          arcPath.addEventListener("mouseenter", (e) => showBranchTooltip(e, child));
          arcPath.addEventListener("mouseleave", () => scheduleHideTooltip(250));
          g.appendChild(arcPath);

          drawRadialBranches(child, childRadius);
        });

        drawNodePoint(g, node);
      }

      drawRadialBranches(activeTreeRoot, 20);

      // 4. Draw Tip Labels & Dotted Circular Alignment Guidelines
      visibleLeaves.forEach(leaf => {
        if (leaf._collapsed) return;
        drawNodePoint(g, leaf);

        const nodeR = typeof settings.nodeRadius === 'number' ? settings.nodeRadius : 3;
        const rLabel = settings.alignLabels ? (maxRadius + nodeR + 10) : (leaf.radius + nodeR + 4);
        const lx = centerX + rLabel * Math.cos(leaf.angle);
        const ly = centerY + rLabel * Math.sin(leaf.angle);

        // Circular Alignment Dotted Hairline
        if (settings.alignLabels && rLabel > leaf.radius + nodeR + 2) {
          const guideLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          guideLine.setAttribute("x1", leaf.x);
          guideLine.setAttribute("y1", leaf.y);
          guideLine.setAttribute("x2", centerX + (rLabel - 4) * Math.cos(leaf.angle));
          guideLine.setAttribute("y2", centerY + (rLabel - 4) * Math.sin(leaf.angle));
          guideLine.setAttribute("stroke", "var(--border-color)");
          guideLine.setAttribute("stroke-dasharray", "2,3");
          guideLine.setAttribute("stroke-width", "0.75");
          guideLine.setAttribute("opacity", "0.5");
          g.appendChild(guideLine);
        }

        const normAngle = ((leaf.angle % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
        const deg = (normAngle * 180) / Math.PI;
        const isFlipped = deg > 90 && deg < 270;
        const rotDeg = isFlipped ? deg + 180 : deg;

        const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt.setAttribute("x", lx);
        txt.setAttribute("y", ly);
        txt.setAttribute("dominant-baseline", "central");

        if (settings.labelOrientation === "horizontal") {
          const isRight = Math.cos(leaf.angle) >= 0;
          txt.setAttribute("text-anchor", isRight ? "start" : "end");
        } else {
          txt.setAttribute("transform", `rotate(${rotDeg}, ${lx}, ${ly})`);
          txt.setAttribute("text-anchor", isFlipped ? "end" : "start");
        }

        txt.setAttribute("class", `tip-label tip-label-radial ${leaf.name === settings.selectedTaxon ? "selected" : ""} ${getLeafFilterClass(leaf.name)}`);
        txt.setAttribute("data-taxon", leaf.name);
        txt.textContent = getLeafLabelText(leaf.name);
        txt.onclick = (e) => { e.stopPropagation(); selectTaxon(leaf.name); };
        txt.onmouseenter = (e) => showNodeTooltip(e, leaf);
        txt.onmouseleave = () => scheduleHideTooltip(250);
        g.appendChild(txt);
      });
    }

    // UNROOTED EQUAL-ANGLE & EQUAL-DAYLIGHT STAR TREE (Large-Cohort Disentangled Engine)
    function renderUnrootedTree(g, visibleLeaves, maxDepth) {
      const centerX = 450;
      const centerY = 450;
      const unrootedScale = typeof settings.unrootedScale === 'number' ? settings.unrootedScale : 1.0;
      const maxSpan = 380 * unrootedScale;

      function countLeaves(node) {
        if (!node.children || node.children.length === 0 || node._collapsed) return 1;
        return node.children.reduce((acc, c) => acc + countLeaves(c), 0);
      }

      // Concentric Scale Rings for Unrooted View
      if (settings.concentricRings) {
        const ringsGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
        ringsGroup.setAttribute("class", "unrooted-concentric-rings");
        [0.25, 0.5, 0.75, 1.0].forEach(frac => {
          const r = maxSpan * frac * 0.85;
          const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
          circle.setAttribute("cx", centerX);
          circle.setAttribute("cy", centerY);
          circle.setAttribute("r", r);
          circle.setAttribute("fill", "none");
          circle.setAttribute("stroke", "var(--border-color)");
          circle.setAttribute("stroke-width", "0.6px");
          circle.setAttribute("stroke-dasharray", "3,3");
          circle.setAttribute("opacity", isDarkTheme() ? "0.3" : "0.4");
          ringsGroup.appendChild(circle);
        });
        g.appendChild(ringsGroup);
      }

      function layoutEqualAngle(node, startA, endA, curX, curY) {
        node.x = curX;
        node.y = curY;

        if (node._collapsed) {
          drawNodePoint(g, node);
          return;
        }
        if (!node.children || node.children.length === 0) return;

        const span = endA - startA;
        let currA = startA;

        // Equal-Daylight vs Equal-Angle weight calculation
        const useDaylight = Boolean(settings.unrootedDaylight);
        const getChildWeight = (child) => {
          const lc = countLeaves(child);
          if (!useDaylight) return lc;
          // Sublinear power damping (N^0.72 + 0.35) distributes daylight evenly to prevent hairballing
          return Math.pow(lc, 0.72) + 0.35;
        };

        const totalWeight = node.children.reduce((acc, c) => acc + getChildWeight(c), 0);

        node.children.forEach(child => {
          const cWeight = getChildWeight(child);
          const childSpan = (cWeight / Math.max(0.0001, totalWeight)) * span;
          const childA = currA + childSpan / 2;

          const bLen = settings.branchLengths ? (typeof child.length === 'number' && !isNaN(child.length) ? child.length : 1.0) : 1.0;
          const safeMaxDepth = Math.max(0.0001, maxDepth);

          // Branch Length Spreading Transformations:
          let r;
          const lenMode = settings.unrootedLengthMode || "sqrt";
          if (lenMode === "sqrt") {
            // Square-root root-to-tip expansion: expands short internal branches to untangle the core
            const frac = Math.sqrt(Math.max(0.00001, bLen)) / Math.sqrt(safeMaxDepth);
            r = Math.max(12, Math.min(maxSpan * 0.85, frac * maxSpan * 0.75));
          } else if (lenMode === "equal") {
            // Equal topological steps (Cladogram-style star): complete elimination of center crowding
            r = Math.max(16, maxSpan / (Math.max(3, safeMaxDepth) * 1.35));
          } else {
            // Strict linear evolutionary distances
            r = Math.max(8, Math.min(maxSpan * 0.75, (bLen / safeMaxDepth) * maxSpan));
          }

          const nextX = curX + r * Math.cos(childA);
          const nextY = curY + r * Math.sin(childA);

          const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
          line.setAttribute("x1", curX);
          line.setAttribute("y1", curY);
          line.setAttribute("x2", nextX);
          line.setAttribute("y2", nextY);
          line.setAttribute("class", "branch-path");
          const bCol = getBranchStrokeColor(child);
          if (bCol) line.style.stroke = bCol;

          line.addEventListener("mouseenter", (e) => showBranchTooltip(e, child));
          line.addEventListener("mouseleave", () => scheduleHideTooltip(250));
          g.appendChild(line);

          layoutEqualAngle(child, currA, currA + childSpan, nextX, nextY);
          currA += childSpan;
        });

        if (node.children && node.children.length > 0 && !node._collapsed) {
          drawNodePoint(g, node);
        }
      }

      const rotRad = ((settings.treeRotation || 0) * Math.PI) / 180;
      layoutEqualAngle(activeTreeRoot, rotRad, rotRad + 2 * Math.PI, centerX, centerY);

      // Smart label collision decluttering for high-tip counts
      const labelFilter = settings.unrootedLabelFilter || "smart";
      let lastRenderedAngle = -999;
      const minAngularGap = (visibleLeaves.length > 150) ? (2 * Math.PI / 80) : (visibleLeaves.length > 60 ? (2 * Math.PI / 65) : 0);

      visibleLeaves.forEach((leaf, idx) => {
        drawNodePoint(g, leaf);

        const dx = leaf.x - centerX;
        const dy = leaf.y - centerY;
        const angle = Math.atan2(dy, dx);
        const normAngle = ((angle % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);

        // Check if label should be rendered or decluttered
        const isSelected = leaf.name === settings.selectedTaxon;
        const searchInput = document.getElementById("searchInput");
        const hasSearch = searchInput && searchInput.value.trim().length > 0;
        const matchesSearch = hasSearch && getLeafFilterClass(leaf.name) === "";

        let shouldRenderText = true;
        if (labelFilter === "none") {
          shouldRenderText = isSelected || matchesSearch;
        } else if (labelFilter === "smart" && minAngularGap > 0) {
          if (!isSelected && !matchesSearch) {
            if (lastRenderedAngle !== -999 && Math.abs(normAngle - lastRenderedAngle) < minAngularGap) {
              shouldRenderText = false;
            }
          }
        }

        if (!shouldRenderText) return;
        lastRenderedAngle = normAngle;

        const nodeR = typeof settings.nodeRadius === 'number' ? settings.nodeRadius : 3;
        // Radial staggering (+14px on alternating tips) prevents horizontal text collision
        const staggerOffset = (settings.staggerLabels && (idx % 2 === 1)) ? 14 : 0;
        const offset = nodeR + 4 + staggerOffset;
        const lx = leaf.x + offset * Math.cos(angle);
        const ly = leaf.y + offset * Math.sin(angle);

        const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt.setAttribute("x", lx);
        txt.setAttribute("y", ly);
        txt.setAttribute("dominant-baseline", "central");

        if (settings.labelOrientation === "radial") {
          const deg = (normAngle * 180) / Math.PI;
          const isFlipped = deg > 90 && deg < 270;
          const rotDeg = isFlipped ? deg + 180 : deg;
          txt.setAttribute("transform", `rotate(${rotDeg}, ${lx}, ${ly})`);
          txt.setAttribute("text-anchor", isFlipped ? "end" : "start");
        } else {
          txt.setAttribute("text-anchor", dx < 0 ? "end" : "start");
        }

        txt.setAttribute("class", `tip-label tip-label-unrooted ${isSelected ? "selected" : ""} ${getLeafFilterClass(leaf.name)}`);
        txt.setAttribute("data-taxon", leaf.name);
        txt.textContent = getLeafLabelText(leaf.name);
        txt.onclick = (e) => { e.stopPropagation(); selectTaxon(leaf.name); };
        txt.onmouseenter = (e) => showNodeTooltip(e, leaf);
        txt.onmouseleave = () => scheduleHideTooltip(250);
        g.appendChild(txt);
      });
    }

    // DUAL TANGLEGRAM LAYOUT (With True Topology, Min-Crossing Untangle, and Aligned Modes)
    function renderTanglegram(g, defs) {
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

      const treesInfo = getActiveTanglegramTrees();
      const compMode = treesInfo.compMode;
      let tRoot3Di = treesInfo.leftTree;
      let tRootAA = treesInfo.rightTree;
      const leftTitle = treesInfo.leftTitle;
      const rightTitle = treesInfo.rightTitle;
      const leftColor = treesInfo.leftColor;
      const rightColor = treesInfo.rightColor;

      if (!tRoot3Di || !tRootAA) {
        const txt = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt.setAttribute("x", 200);
        txt.setAttribute("y", 150);
        txt.setAttribute("fill", "var(--text-muted)");
        txt.textContent = "No overlapping taxa found for the current active filter/partition.";
        g.appendChild(txt);
        return;
      }

      assignDepths(tRoot3Di, 0);
      assignDepths(tRootAA, 0);

      // 1. Natural DFS post-order leaf traversal for 3Di
      function getTreeLeavesInOrder(node) {
        if (!node.children || node.children.length === 0) return [node];
        let res = [];
        node.children.forEach(c => {
          res = res.concat(getTreeLeavesInOrder(c));
        });
        return res;
      }

      const leaves3Di = getTreeLeavesInOrder(tRoot3Di);
      leaves3Di.forEach((leaf, idx) => {
        leaf.y = 50 + idx * spacing;
      });

      const map3DiIndex = {};
      leaves3Di.forEach((l, idx) => {
        map3DiIndex[l.name] = idx;
      });

      // 2. Order AA Tree Leaves based on settings.tangleMode
      let leavesAA = [];

      if (settings.tangleMode === "aligned") {
        // Aligned / Parallel: order right-hand leaves to match left-hand leaves
        leavesAA = [...getAllLeaves(tRootAA)].sort((a, b) => {
          const idxA = map3DiIndex[a.name] !== undefined ? map3DiIndex[a.name] : 9999;
          const idxB = map3DiIndex[b.name] !== undefined ? map3DiIndex[b.name] : 9999;
          return idxA - idxB;
        });
        leavesAA.forEach((leaf, idx) => {
          leaf.y = 50 + idx * spacing;
        });
      } else if (settings.tangleMode === "min_crossings") {
        // Optimal Subtree Rotations (Barycenter Heuristic): rotate children around internal nodes
        // to minimize crossings while strictly preserving phylogenetic clades
        function getAvg3Di(node) {
          const lf = getTreeLeavesInOrder(node);
          const idxs = lf.map(l => map3DiIndex[l.name] !== undefined ? map3DiIndex[l.name] : 0);
          return idxs.reduce((a, b) => a + b, 0) / (idxs.length || 1);
        }
        function rotateSubtrees(node) {
          if (!node.children || node.children.length <= 1) return;
          node.children.forEach(rotateSubtrees);
          node.children.sort((a, b) => getAvg3Di(a) - getAvg3Di(b));
        }
        rotateSubtrees(tRootAA);
        leavesAA = getTreeLeavesInOrder(tRootAA);
        leavesAA.forEach((leaf, idx) => {
          leaf.y = 50 + idx * spacing;
        });
      } else {
        // "true_topology": Independent DFS post-order traversal of the AA tree (Native IQ-TREE order)
        // Reveals true topological discordance and crossing tangles!
        leavesAA = getTreeLeavesInOrder(tRootAA);
        leavesAA.forEach((leaf, idx) => {
          leaf.y = 50 + idx * spacing;
        });
      }

      // Compute internal Y coordinates for both trees
      function computeInternal(node) {
        if (!node.children || node.children.length === 0) return;
        node.children.forEach(computeInternal);
        const childY = node.children.map(c => c.y).filter(y => typeof y === 'number' && !isNaN(y));
        if (childY.length > 0) {
          node.y = childY.reduce((a, b) => a + b, 0) / childY.length;
        }
      }
      computeInternal(tRoot3Di);
      computeInternal(tRootAA);

      const maxDepth3Di = Math.max(...leaves3Di.map(l => settings.branchLengths ? l.depth : l.cladoDepth)) || 1.0;
      const maxDepthAA = Math.max(...leavesAA.map(l => settings.branchLengths ? l.depth : l.cladoDepth)) || 1.0;

      // Draw Left Tree (3Di Structural, branching rightwards towards center)
      function drawLeft(node, currentX) {
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
        applyVerticalBranchGradient(vLine, node, minY, maxY, defs);
        g.appendChild(vLine);

        node.children.forEach(child => {
          const bLen = settings.branchLengths ? (typeof child.length === 'number' ? child.length : 1.0) : 1.0;
          const childX = node.x + (bLen / maxDepth3Di) * leftTreeSpan;
          child.x = childX;

          const hLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          hLine.setAttribute("x1", node.x);
          hLine.setAttribute("y1", child.y);
          hLine.setAttribute("x2", childX);
          hLine.setAttribute("y2", child.y);
          hLine.setAttribute("class", "branch-path");
          const bCol = getBranchStrokeColor(child);
          if (bCol) hLine.style.stroke = bCol;
          hLine.addEventListener("mouseenter", (e) => showBranchTooltip(e, child));
          hLine.addEventListener("mouseleave", () => scheduleHideTooltip(250));
          g.appendChild(hLine);

          drawLeft(child, childX);
        });
        drawNodePoint(g, node);
      }

      // Draw Right Tree (AA Sequence, branching leftwards towards center)
      function drawRight(node, currentX) {
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
        applyVerticalBranchGradient(vLine, node, minY, maxY, defs);
        g.appendChild(vLine);

        node.children.forEach(child => {
          const bLen = settings.branchLengths ? (typeof child.length === 'number' ? child.length : 1.0) : 1.0;
          const childX = node.x - (bLen / maxDepthAA) * rightTreeSpan;
          child.x = childX;

          const hLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
          hLine.setAttribute("x1", node.x);
          hLine.setAttribute("y1", child.y);
          hLine.setAttribute("x2", childX);
          hLine.setAttribute("y2", child.y);
          hLine.setAttribute("class", "branch-path");
          const bCol = getBranchStrokeColor(child);
          if (bCol) hLine.style.stroke = bCol;
          hLine.addEventListener("mouseenter", (e) => showBranchTooltip(e, child));
          hLine.addEventListener("mouseleave", () => scheduleHideTooltip(250));
          g.appendChild(hLine);

          drawRight(child, childX);
        });
        drawNodePoint(g, node);
      }

      drawLeft(tRoot3Di, leftTreeRootX);
      drawRight(tRootAA, rightTreeRootX);

      // Section Header Badges
      const leftHeader = document.createElementNS("http://www.w3.org/2000/svg", "text");
      leftHeader.setAttribute("x", leftTreeRootX);
      leftHeader.setAttribute("y", 25);
      leftHeader.setAttribute("fill", leftColor);
      leftHeader.setAttribute("font-weight", "bold");
      leftHeader.setAttribute("font-size", "12px");
      leftHeader.textContent = leftTitle;
      g.appendChild(leftHeader);

      const rightHeader = document.createElementNS("http://www.w3.org/2000/svg", "text");
      rightHeader.setAttribute("x", rightTreeRootX);
      rightHeader.setAttribute("y", 25);
      rightHeader.setAttribute("text-anchor", "end");
      rightHeader.setAttribute("fill", rightColor);
      rightHeader.setAttribute("font-weight", "bold");
      rightHeader.setAttribute("font-size", "12px");
      rightHeader.textContent = rightTitle;
      g.appendChild(rightHeader);

      // Count Inversions / Crossings
      const mapAA = {};
      leavesAA.forEach(l => { mapAA[l.name] = l; });

      const arrCross = [];
      leaves3Di.forEach(l1 => {
        const l2 = mapAA[l1.name];
        if (l2 && typeof l2.y === 'number') arrCross.push(l2.y);
      });
      let crossings = 0;
      for (let i = 0; i < arrCross.length; i++) {
        for (let j = i + 1; j < arrCross.length; j++) {
          if (arrCross[i] > arrCross[j]) crossings++;
        }
      }
      const badgeCross = document.getElementById("tangleCrossingBadge");
      if (badgeCross) {
        badgeCross.textContent = `Crossings: ${crossings.toLocaleString()}`;
      }

      // Draw Tanglegram Connecting Curves and Labels
      leaves3Di.forEach(l1 => {
        drawNodePoint(g, l1);
        const l2 = mapAA[l1.name];
        if (!l2) return;
        drawNodePoint(g, l2);

        const cleanId = l1.name.replace(/[^a-zA-Z0-9]/g, "_");

        // Hairline extension from 3Di leaf to left label
        if (leftLabelX - 5 > l1.x) {
          const extL = document.createElementNS("http://www.w3.org/2000/svg", "line");
          extL.setAttribute("x1", l1.x);
          extL.setAttribute("y1", l1.y);
          extL.setAttribute("x2", leftLabelX - 5);
          extL.setAttribute("y2", l1.y);
          extL.setAttribute("stroke", "var(--grid-line)");
          extL.setAttribute("stroke-width", "0.6px");
          extL.setAttribute("stroke-dasharray", "2,2");
          g.appendChild(extL);
        }

        // Left Label
        const txt1 = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt1.setAttribute("id", "tangle_label_left_" + cleanId);
        txt1.setAttribute("x", leftLabelX);
        txt1.setAttribute("y", l1.y + 3.5);
        txt1.setAttribute("class", `tip-label ${l1.name === settings.selectedTaxon ? "selected" : ""}`);
        txt1.setAttribute("data-taxon", l1.name);
        txt1.textContent = getLeafLabelText(l1.name);
        txt1.onclick = () => selectTaxon(l1.name);
        txt1.onmouseenter = () => highlightTangleTaxon(l1.name);
        txt1.onmouseleave = () => clearTangleHighlight();
        g.appendChild(txt1);

        // Hairline extension from AA leaf to right label
        if (l2.x > rightLabelX + 5) {
          const extR = document.createElementNS("http://www.w3.org/2000/svg", "line");
          extR.setAttribute("x1", rightLabelX + 5);
          extR.setAttribute("y1", l2.y);
          extR.setAttribute("x2", l2.x);
          extR.setAttribute("y2", l2.y);
          extR.setAttribute("stroke", "var(--grid-line)");
          extR.setAttribute("stroke-width", "0.6px");
          extR.setAttribute("stroke-dasharray", "2,2");
          g.appendChild(extR);
        }

        // Right Label
        const txt2 = document.createElementNS("http://www.w3.org/2000/svg", "text");
        txt2.setAttribute("id", "tangle_label_right_" + cleanId);
        txt2.setAttribute("x", rightLabelX);
        txt2.setAttribute("y", l2.y + 3.5);
        txt2.setAttribute("text-anchor", "end");
        txt2.setAttribute("class", `tip-label ${l1.name === settings.selectedTaxon ? "selected" : ""}`);
        txt2.setAttribute("data-taxon", l2.name);
        txt2.textContent = getLeafLabelText(l2.name);
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
        path.setAttribute("d", `M ${x1} ${y1} C ${cp1x} ${y1}, ${cp2x} ${y2}, ${x2} ${y2}`);
        path.setAttribute("class", `tangle-connector tangle_conn_${cleanId} ${getLeafFilterClass(l1.name)}`);
        path.style.stroke = getNodeColor(l1);

        path.onmouseenter = () => highlightTangleTaxon(l1.name);
        path.onmouseleave = () => clearTangleHighlight();
        path.onclick = () => selectTaxon(l1.name);
        g.appendChild(path);
      });

      // Record tanglegram geometry and active roots for the radar minimap
      tangleState.activeLeftRoot = tRoot3Di;
      tangleState.activeRightRoot = tRootAA;
      tangleState.leftLeaves = leaves3Di;
      tangleState.rightLeaves = leavesAA;
      tangleState.leftTreeRootX = leftTreeRootX;
      tangleState.rightTreeRootX = rightTreeRootX;
      tangleState.leftConnectorX = leftConnectorX;
      tangleState.rightConnectorX = rightConnectorX;
      tangleState.leftColor = leftColor;
      tangleState.rightColor = rightColor;
    }

    // 11. TOOLTIPS & ACTIONS
    const tooltip = document.getElementById("treeTooltip");
    let hideTimer = null;

    function scheduleHideTooltip(delay = 250) {
      clearTimeout(hideTimer);
      hideTimer = setTimeout(() => {
        tooltip.classList.add("hidden");
      }, delay);
    }

    function hideTooltipImmediate() {
      clearTimeout(hideTimer);
      tooltip.classList.add("hidden");
    }

    tooltip.addEventListener("mouseenter", () => {
      clearTimeout(hideTimer);
    });
    tooltip.addEventListener("mouseleave", () => {
      scheduleHideTooltip(180);
    });

    function showNodeTooltip(e, leaf) {
      clearTimeout(hideTimer);
      pendingHoverNode = leaf;

      const meta = TAXA_METADATA[leaf.name] || {};
      const colDef = getActiveColorColumnDef();

      document.getElementById("tooltipId").textContent = leaf.name;
      const primaryVal = (colDef && meta[colDef.key] !== undefined) ? String(meta[colDef.key]) : "-";
      document.getElementById("tooltipBadge").textContent = primaryVal;

      const metaList = document.getElementById("tooltipMetaList");
      metaList.innerHTML = "";

      Object.keys(meta).forEach(k => {
        if (k === "color") return;
        const cleanK = k.replace(/_/g, " ").replace(/\\b\\w/g, l => l.toUpperCase());
        const div = document.createElement("div");
        div.className = "flex justify-between items-center text-[10px]";
        div.innerHTML = `<span>${cleanK}:</span><strong class="text-[var(--text-main)] font-mono ml-2">${meta[k]}</strong>`;
        metaList.appendChild(div);
      });

      document.getElementById("btnTooltipReroot").textContent = "⚓ Reroot at Leaf Edge";
      document.getElementById("tooltipViewerWrapper").classList.remove("hidden");

      positionTooltip(e);
      tooltip.classList.remove("hidden");
      tooltipViewer.loadStructure(leaf.name);
    }

    function showCladeTooltip(e, node) {
      clearTimeout(hideTimer);
      pendingHoverNode = node;

      const info = getCladeInfo(node);
      document.getElementById("tooltipId").textContent = `Clade (${info.count} Taxa)`;
      document.getElementById("tooltipBadge").textContent = `${info.homogeneity.toFixed(0)}% ${info.dominantCategory}`;

      const metaList = document.getElementById("tooltipMetaList");
      metaList.innerHTML = `
        <div class="flex justify-between"><span>Leaves:</span><strong class="text-[var(--text-main)] font-mono">${info.count} taxa</strong></div>
        <div class="flex justify-between"><span>Dominant ${info.colLabel}:</span><strong class="text-[var(--text-main)] font-mono">${info.dominantCategory}</strong></div>
        <div class="flex justify-between"><span>Purity:</span><strong class="text-emerald-400 font-mono">${info.homogeneity.toFixed(1)}%</strong></div>
        <div class="flex justify-between"><span>Internal Depth:</span><strong class="text-sky-400 font-mono">${node.cladoDepth || 0}</strong></div>
      `;

      document.getElementById("btnTooltipReroot").textContent = "⚓ Reroot at Clade Stem";
      document.getElementById("tooltipViewerWrapper").classList.add("hidden");

      positionTooltip(e);
      tooltip.classList.remove("hidden");
    }

    function showBranchTooltip(e, childNode) {
      clearTimeout(hideTimer);
      pendingHoverNode = childNode;

      const bLen = typeof childNode.length === 'number' ? childNode.length.toFixed(5) : "-";
      document.getElementById("tooltipId").textContent = childNode.name ? `Branch: ${childNode.name}` : "Phylogenetic Branch";
      document.getElementById("tooltipBadge").textContent = `${bLen} subs/site`;

      const metaList = document.getElementById("tooltipMetaList");
      let supportHtml = "";
      if (childNode.ufboot !== undefined && childNode.ufboot !== null && !isNaN(childNode.ufboot)) {
        const ufVal = childNode.ufboot <= 1.0 && childNode.ufboot > 0.0 ? (childNode.ufboot * 100).toFixed(1) : childNode.ufboot;
        supportHtml += `<div class="flex justify-between"><span>UFboot Support:</span><strong class="text-sky-400 font-mono">${ufVal}%</strong></div>`;
      }
      if (childNode.alrt !== undefined && childNode.alrt !== null && !isNaN(childNode.alrt)) {
        const alrtVal = childNode.alrt <= 1.0 && childNode.alrt > 0.0 ? (childNode.alrt * 100).toFixed(1) : childNode.alrt;
        supportHtml += `<div class="flex justify-between"><span>SH-aLRT Support:</span><strong class="text-indigo-400 font-mono">${alrtVal}%</strong></div>`;
      }
      if (!supportHtml) {
        const supVal = childNode.support !== null && !isNaN(childNode.support) ? (childNode.support <= 1.0 && childNode.support > 0.0 ? (childNode.support * 100).toFixed(1) + "%" : childNode.support) : "N/A";
        supportHtml = `<div class="flex justify-between"><span>Support:</span><strong class="text-sky-400 font-mono">${supVal}</strong></div>`;
      }

      metaList.innerHTML = `
        <div class="flex justify-between"><span>Branch Length:</span><strong class="text-emerald-400 font-mono">${bLen}</strong></div>
        ${supportHtml}
      `;

      document.getElementById("btnTooltipReroot").textContent = "⚓ Reroot at this Branch";
      document.getElementById("tooltipViewerWrapper").classList.add("hidden");

      positionTooltip(e);
      tooltip.classList.remove("hidden");
    }

    function showCollapsedTooltip(e, node, info) {
      clearTimeout(hideTimer);
      pendingHoverNode = node;

      document.getElementById("tooltipId").textContent = `Collapsed Clade: ${info.dominantCategory}`;
      document.getElementById("tooltipBadge").textContent = `${info.count} taxa`;

      const metaList = document.getElementById("tooltipMetaList");
      metaList.innerHTML = `
        <div class="flex justify-between"><span>Collapsed Leaves:</span><strong class="text-[var(--text-main)] font-mono">${info.count} taxa</strong></div>
        <div class="flex justify-between"><span>Majority Category:</span><strong class="text-sky-400 font-mono">${info.dominantCategory}</strong></div>
        <div class="flex justify-between"><span>Homogeneity:</span><strong class="text-emerald-400 font-mono">${info.homogeneity.toFixed(1)}%</strong></div>
        <p class="text-[9.5px] text-amber-400 pt-1">Click wedge to expand back into full sub-tree.</p>
      `;

      document.getElementById("btnTooltipReroot").textContent = "⚓ Reroot at Clade Stem";
      document.getElementById("tooltipViewerWrapper").classList.add("hidden");

      positionTooltip(e);
      tooltip.classList.remove("hidden");
    }

    function positionTooltip(e) {
      const pad = 16;
      let x = e.clientX + pad;
      let y = e.clientY + pad;

      const maxW = window.innerWidth - 300;
      const maxH = window.innerHeight - 260;

      if (x > maxW) x = e.clientX - 290;
      if (y > maxH) y = e.clientY - 250;

      tooltip.style.left = Math.max(10, x) + "px";
      tooltip.style.top = Math.max(10, y) + "px";
    }

    // PINNED STRUCTURE SELECTION
    function selectTaxon(taxName, moveTree = false) {
      settings.selectedTaxon = taxName;
      // Only scroll MSA if sync is explicitly enabled by user
      if (msaState.syncWithTree && typeof getMsaTaxaList === 'function') {
        const msaTaxa = getMsaTaxaList();
        const idx = msaTaxa.indexOf(taxName);
        if (idx >= 0) {
          msaState.scrollY = Math.max(0, idx - 2);
        }
      }
      if (typeof renderMsa === 'function') {
        renderMsa();
      }
      document.querySelectorAll(".tip-label").forEach(el => {
        const taxon = el.getAttribute("data-taxon") || el.textContent;
        if (taxon === taxName) {
          el.classList.add("selected");
        } else {
          el.classList.remove("selected");
        }
      });

      const card = document.getElementById("selectedCard");
      if (!card) return;
      card.classList.remove("hidden");
      document.getElementById("cardId").textContent = taxName;

      const meta = TAXA_METADATA[taxName] || {};
      const colDef = getActiveColorColumnDef();
      const primaryVal = (colDef && meta[colDef.key] !== undefined) ? String(meta[colDef.key]) : "-";
      const primaryColor = (colDef && colDef.colors && colDef.colors[primaryVal]) ? colDef.colors[primaryVal] : "#38bdf8";

      const badgeCategory = document.getElementById("cardBadgeCategory");
      if (badgeCategory) {
        badgeCategory.textContent = primaryVal;
        badgeCategory.style.backgroundColor = primaryColor + "22";
        badgeCategory.style.borderColor = primaryColor + "66";
        badgeCategory.style.color = primaryColor;
      }

      const cardPlddt = document.getElementById("cardPlddt");
      if (cardPlddt) cardPlddt.textContent = meta.plddt !== undefined ? meta.plddt : "-";

      const cardLen = document.getElementById("cardLen");
      if (cardLen) cardLen.textContent = meta.length !== undefined ? (meta.length + " aa") : "-";

      const attrTable = document.getElementById("cardAttributesTable");
      if (attrTable) {
        attrTable.innerHTML = "";
        const skipKeys = new Set(["id", "color", colDef ? colDef.key : ""]);
        const displayKeys = Object.keys(meta).filter(k => !skipKeys.has(k)).slice(0, 8);

        displayKeys.forEach(k => {
          const row = document.createElement("div");
          row.className = "flex justify-between items-center text-[9.5px] border-b border-slate-800/40 py-0.5";
          const cleanK = k.replace(/_/g, " ").replace(/\b\w/g, l => l.toUpperCase());
          row.innerHTML = `
            <span class="text-[var(--text-muted)] truncate max-w-[105px]">${cleanK}:</span>
            <span class="font-medium text-[var(--text-main)] truncate max-w-[155px] text-right" title="${meta[k]}">${meta[k]}</span>
          `;
          attrTable.appendChild(row);
        });
      }

      sidebarViewer.loadStructure(taxName);
    }

    function highlightTangleTaxon(taxName) {
      const cleanId = taxName.replace(/[^a-zA-Z0-9]/g, "_");
      const el = document.getElementById("tangle_" + cleanId);
      if (el) el.classList.add("highlighted");
      const leftLbl = document.getElementById("tangle_label_left_" + cleanId);
      if (leftLbl) leftLbl.classList.add("selected");
      const rightLbl = document.getElementById("tangle_label_right_" + cleanId);
      if (rightLbl) rightLbl.classList.add("selected");
    }

    function highlightTangleConnector(taxName) {
      highlightTangleTaxon(taxName);
    }

    function clearTangleHighlight() {
      document.querySelectorAll(".tangle-connector.highlighted").forEach(el => el.classList.remove("highlighted"));
      document.querySelectorAll(".tip-label.selected").forEach(el => {
        const taxon = el.getAttribute("data-taxon") || el.textContent;
        if (taxon !== settings.selectedTaxon) {
          el.classList.remove("selected");
        }
      });
    }

    // 12. NATURAL SMOOTH NAVIGATION & VIEWPORT PANNING
    const container = document.getElementById("treeCanvasWrapper") || document.getElementById("treeContainer");
    let isPanning = false;
    let startPan = { x: 0, y: 0 };

    container.addEventListener("mousedown", (e) => {
      if (e.target.closest("button") || e.target.closest("#minimapBox") || e.target.closest("#treeTooltip") || e.target.closest("#alignmentDrawer")) return;
      isPanning = true;
      startPan = { x: e.clientX - settings.zoom.x, y: e.clientY - settings.zoom.y };
    });

    window.addEventListener("mouseup", () => { isPanning = false; });
    window.addEventListener("mousemove", (e) => {
      if (!isPanning) return;
      settings.zoom.x = e.clientX - startPan.x;
      settings.zoom.y = e.clientY - startPan.y;
      updateViewportTransform(false);
      updateMinimap();
    });

    // NATURAL NAVIGATION: 2-Finger Trackpad / Wheel Pan + Ctrl/Meta Pinch-to-Zoom
    container.addEventListener("wheel", (e) => {
      if (e.target.closest("#alignmentDrawer")) return; // Decoupled from alignment drawer
      e.preventDefault();
      const rect = container.getBoundingClientRect();
      const mouseX = e.clientX - rect.left;
      const mouseY = e.clientY - rect.top;

      if (e.ctrlKey || e.metaKey) {
        // Pinch-to-zoom or Ctrl+wheel zoom
        const factor = Math.exp(-e.deltaY * 0.005);
        const zoomFactor = Math.min(Math.max(factor, 0.85), 1.15);
        const oldK = settings.zoom.k;
        const newK = Math.max(0.04, Math.min(oldK * zoomFactor, 5.0));

        settings.zoom.x = mouseX - (mouseX - settings.zoom.x) * (newK / oldK);
        settings.zoom.y = mouseY - (mouseY - settings.zoom.y) * (newK / oldK);
        settings.zoom.k = newK;
      } else {
        // Smooth 2D panning via trackpad scroll or mouse wheel
        settings.zoom.x -= e.deltaX;
        settings.zoom.y -= e.deltaY;
      }

      updateViewportTransform(false);
      updateMinimap();
    }, { passive: false });

    // ISOLATE ALIGNMENT DRAWER EVENTS & ADD MATRIX DRAG-TO-PAN
    function setupMsaInteractions() {
      const drawer = document.getElementById("alignmentDrawer");
      if (drawer) {
        drawer.addEventListener("wheel", (e) => {
          e.stopPropagation();
        }, { passive: false });
        drawer.addEventListener("mousedown", (e) => {
          e.stopPropagation();
        });
        drawer.addEventListener("pointerdown", (e) => {
          e.stopPropagation();
        });
      }

      const matrixWrap = document.getElementById("msaMatrixWrapper");
      if (matrixWrap) {
        let isMsaDragging = false;
        let msaDragStart = { x: 0, y: 0, scrollX: 0, scrollY: 0 };

        matrixWrap.addEventListener("mousedown", (e) => {
          if (e.button !== 0) return;
          e.stopPropagation();
          isMsaDragging = true;
          msaDragStart = {
            x: e.clientX,
            y: e.clientY,
            scrollX: msaState.scrollX,
            scrollY: msaState.scrollY
          };
          matrixWrap.style.cursor = "grabbing";
        });

        window.addEventListener("mousemove", (e) => {
          if (!isMsaDragging) return;
          const dx = (e.clientX - msaDragStart.x) / msaState.cellWidth;
          const dy = (e.clientY - msaDragStart.y) / msaState.cellHeight;
          const align = getActiveAlignment();
          const maxCols = (msaState.activeKeptCols && msaState.activeKeptCols.length) ? msaState.activeKeptCols.length : (align ? align.length : 533);
          const taxa = getMsaTaxaList();
          const maxRows = taxa.length;

          msaState.scrollX = Math.max(0, Math.min(maxCols - 5, msaDragStart.scrollX - dx));
          msaState.scrollY = Math.max(0, Math.min(maxRows - 2, msaDragStart.scrollY - dy));

          const posInput = document.getElementById("msaPosInput");
          if (posInput) posInput.value = Math.floor(msaState.scrollX) + 1;

          renderMsa();
        });

        window.addEventListener("mouseup", () => {
          if (isMsaDragging) {
            isMsaDragging = false;
            if (matrixWrap) matrixWrap.style.cursor = "crosshair";
          }
        });
      }

      // Interactive Alignment Minimap Drag & Pan Navigation
      const minimapWrap = document.getElementById("msaMinimapWrapper");
      if (minimapWrap) {
        let isMiniDragging = false;
        function handleMinimapNav(e) {
          const rect = minimapWrap.getBoundingClientRect();
          if (rect.width <= 0 || rect.height <= 0) return;
          const x = Math.max(0, Math.min(rect.width, e.clientX - rect.left));
          const y = Math.max(0, Math.min(rect.height, e.clientY - rect.top));
          const align = getActiveAlignment();
          if (!align) return;
          const alignLen = (msaState.activeKeptCols && msaState.activeKeptCols.length) ? msaState.activeKeptCols.length : (align.length || 533);
          const taxa = getMsaTaxaList();
          if (!taxa.length) return;

          const matrixCanvas = document.getElementById("msaMatrixCanvas");
          const visCols = matrixCanvas ? (matrixCanvas.clientWidth / msaState.cellWidth) : 20;
          const visRows = matrixCanvas ? (matrixCanvas.clientHeight / msaState.cellHeight) : 10;

          const targetCol = (x / rect.width) * alignLen;
          const targetRow = (y / rect.height) * taxa.length;

          msaState.scrollX = Math.max(0, Math.min(alignLen - 5, targetCol - visCols / 2));
          msaState.scrollY = Math.max(0, Math.min(taxa.length - 2, targetRow - visRows / 2));

          const posInput = document.getElementById("msaPosInput");
          if (posInput) posInput.value = Math.floor(msaState.scrollX) + 1;

          renderMsa();
        }

        minimapWrap.addEventListener("mousedown", (e) => {
          if (e.button !== 0) return;
          e.stopPropagation();
          e.preventDefault();
          isMiniDragging = true;
          handleMinimapNav(e);
        });

        window.addEventListener("mousemove", (e) => {
          if (!isMiniDragging) return;
          e.stopPropagation();
          handleMinimapNav(e);
        });

        window.addEventListener("mouseup", () => {
          if (isMiniDragging) isMiniDragging = false;
        });
      }
    }

    // DYNAMIC ZOOM-ADAPTIVE LABEL SCALING (RADIAL & UNROOTED)
    function updateLabelScaling() {
      const base = typeof settings.labelSize === 'number' ? settings.labelSize : 10;
      const k = Math.max(0.01, (settings.zoom && settings.zoom.k) || 1);

      if (base === 0) {
        document.documentElement.setAttribute("data-labels-hidden", "true");
        return;
      } else {
        document.documentElement.removeAttribute("data-labels-hidden");
      }

      let radialSize;
      if (settings.zoomAdaptiveLabels) {
        // Shrink font size on zoom in so dense branches in radial & unrooted layouts remain unobscured
        if (k <= 1.0) {
          radialSize = Math.min(20, base / Math.pow(k, 0.55));
        } else {
          // As zoom k increases, shrink font size in screen pixels
          radialSize = Math.max(2.2, base / Math.pow(k, 1.25));
        }
      } else {
        radialSize = base;
      }

      document.documentElement.style.setProperty('--radial-tip-size', `${radialSize.toFixed(2)}px`);
      document.documentElement.style.setProperty('--tree-tip-size', `${base}px`);
    }

    function updateViewportTransform(smooth = false) {
      const vp = document.getElementById("mainViewport");
      if (vp) {
        vp.style.transition = smooth ? "transform 0.22s ease-out" : "none";
        vp.setAttribute("transform", `translate(${settings.zoom.x}, ${settings.zoom.y}) scale(${settings.zoom.k})`);
      }
      updateLabelScaling();
    }

    function zoomStep(factor) {
      const rect = container.getBoundingClientRect();
      const cx = rect.width / 2;
      const cy = rect.height / 2;

      const oldK = settings.zoom.k;
      const newK = Math.max(0.04, Math.min(oldK * factor, 5.0));

      settings.zoom.x = cx - (cx - settings.zoom.x) * (newK / oldK);
      settings.zoom.y = cy - (cy - settings.zoom.y) * (newK / oldK);
      settings.zoom.k = newK;

      updateViewportTransform(true);
      updateMinimap();
    }

    // AUTO-FIT TREE TO SCREEN (LARGE COHORT OPTIMIZED)
    function fitTreeToScreen(smooth = true) {
      if (settings.layout === "umap") {
        const cW = container.clientWidth || 900;
        const cH = container.clientHeight || 700;
        const minX = 60;
        const maxX = 1180;
        const minY = 40;
        const maxY = 780;
        const spanW = maxX - minX;
        const spanH = maxY - minY;

        const pad = 30;
        const scaleX = (cW - pad * 2) / spanW;
        const scaleY = (cH - pad * 2) / spanH;
        const k = Math.max(0.25, Math.min(scaleX, scaleY, 1.4));

        settings.zoom.k = k;
        settings.zoom.x = (cW - spanW * k) / 2 - minX * k;
        settings.zoom.y = (cH - spanH * k) / 2 - minY * k;

        updateViewportTransform(smooth);
        updateMinimap();
        return;
      }

      if (settings.layout === "tanglegram") {
        const cW = container.clientWidth || 900;
        const cH = container.clientHeight || 700;
        const minX = 40;
        const maxX = 1170;
        const treeW = maxX - minX;

        const leftLeaves = tangleState.leftLeaves.length > 0 ? tangleState.leftLeaves : getAllLeaves(rawRoot3Di);
        const leavesCount = (leftLeaves && leftLeaves.length) || (activeDataset && Object.keys(activeDataset.taxa).length) || 6;
        const treeH = Math.max(100, 50 + leavesCount * settings.verticalSpacing + 40);

        const pad = 30;
        const scaleX = (cW - pad * 2) / treeW;
        const scaleY = (cH - pad * 2) / treeH;

        if (leavesCount <= 12) {
          const k = Math.max(0.20, Math.min(scaleX, scaleY, 1.4));
          settings.zoom.k = k;
          settings.zoom.x = (cW - treeW * k) / 2 - minX * k;
          settings.zoom.y = (cH - treeH * k) / 2 - 30 * k;
        } else {
          const k = Math.max(0.25, Math.min(scaleX, 0.95));
          settings.zoom.k = k;
          settings.zoom.x = (cW - treeW * k) / 2 - minX * k;
          settings.zoom.y = 35;
        }

        updateViewportTransform(smooth);
        updateMinimap();
        return;
      }

      const visibleLeaves = getVisibleLeaves(activeTreeRoot);
      if (visibleLeaves.length === 0) return;

      const validX = visibleLeaves.map(l => l.x).filter(x => typeof x === 'number' && !isNaN(x));
      const validY = visibleLeaves.map(l => l.y).filter(y => typeof y === 'number' && !isNaN(y));

      if (validX.length === 0 || validY.length === 0) return;

      let minX, maxX, minY, maxY;
      if (settings.layout === "radial" || settings.layout === "unrooted") {
        const labelPadX = 140;
        const labelPadY = 35;
        minX = Math.min(...validX) - labelPadX;
        maxX = Math.max(...validX) + labelPadX;
        minY = Math.min(...validY) - labelPadY;
        maxY = Math.max(...validY) + labelPadY;
      } else {
        minX = Math.min(...validX, activeTreeRoot.x || 50);
        maxX = Math.max(...validX) + 140;
        minY = Math.min(...validY);
        maxY = Math.max(...validY);
      }

      const treeW = Math.max(maxX - minX, 100);
      const treeH = Math.max(maxY - minY, 100);

      const cW = container.clientWidth || 900;
      const cH = container.clientHeight || 700;

      const pad = 40;
      const scaleX = (cW - pad * 2) / treeW;
      const scaleY = (cH - pad * 2) / treeH;

      if (settings.layout === "radial" || settings.layout === "unrooted") {
        const k = Math.max(0.12, Math.min(scaleX, scaleY, 1.8));
        settings.zoom.k = k;
        settings.zoom.x = (cW - treeW * k) / 2 - minX * k;
        settings.zoom.y = (cH - treeH * k) / 2 - minY * k;
      } else {
        if (visibleLeaves.length <= 80) {
          const k = Math.max(0.25, Math.min(scaleX, scaleY, 1.6));
          settings.zoom.k = k;
          settings.zoom.x = (cW - treeW * k) / 2 - minX * k;
          settings.zoom.y = (cH - treeH * k) / 2 - minY * k;
        } else {
          const k = Math.max(0.40, Math.min(scaleX, 0.85));
          settings.zoom.k = k;
          settings.zoom.x = 40;
          settings.zoom.y = 35;
        }
      }

      updateViewportTransform(smooth);
      updateMinimap();
    }

    function centerOnSelection() {
      if (!settings.selectedTaxon) {
        fitTreeToScreen(true);
        return;
      }
      if (settings.layout === "umap") {
        const target = umapPointCoords[settings.selectedTaxon];
        if (!target || typeof target.x !== 'number' || typeof target.y !== 'number') {
          fitTreeToScreen(true);
          return;
        }
        const W = container.clientWidth || 900;
        const H = container.clientHeight || 700;
        settings.zoom.x = W / 2 - target.x * settings.zoom.k;
        settings.zoom.y = H / 2 - target.y * settings.zoom.k;
        updateViewportTransform(true);
        updateMinimap();
        return;
      }
      const visibleLeaves = getVisibleLeaves(activeTreeRoot);
      const target = visibleLeaves.find(l => l.name === settings.selectedTaxon);
      if (!target || typeof target.x !== 'number' || typeof target.y !== 'number') {
        fitTreeToScreen(true);
        return;
      }

      const W = container.clientWidth || 900;
      const H = container.clientHeight || 700;

      settings.zoom.x = W / 2 - target.x * settings.zoom.k;
      settings.zoom.y = H / 2 - target.y * settings.zoom.k;

      updateViewportTransform(true);
      updateMinimap();
    }

    // RADAR OVERVIEW MINIMAP
    const minimapBox = document.getElementById("minimapBox");
    const minimapCanvas = document.getElementById("minimapCanvas");
    const minimapCtx = minimapCanvas.getContext("2d");
    const minimapViewport = document.getElementById("minimapViewport");
    let minimapBounds = {};

    function updateMinimap() {
      if (!minimapCanvas) return;
      const mw = minimapCanvas.width;
      const mh = minimapCanvas.height;
      minimapCtx.clearRect(0, 0, mw, mh);

      const isDark = isDarkTheme();
      minimapCtx.fillStyle = isDark ? (settings.theme === "obsidian" ? "#030712" : (settings.theme === "forest" ? "#041f16" : "#0b1120")) : (settings.theme === "solarized" ? "#fdf6e3" : (settings.theme === "nordic" ? "#eceff4" : "#ffffff"));
      minimapCtx.fillRect(0, 0, mw, mh);

      const titleEl = document.getElementById("minimapTitle");
      if (titleEl) {
        titleEl.textContent = (settings.layout === "tanglegram") ? "TANGLEGRAM RADAR" : ((settings.layout === "umap") ? "UMAP RADAR" : "RADAR OVERVIEW");
      }

      const pad = 10;
      let minX, maxX, minY, maxY, scaleX, scaleY;

      if (settings.layout === "umap") {
        minX = 100;
        maxX = 1140;
        minY = 50;
        maxY = 740;
        const spanX = Math.max(maxX - minX, 100);
        const spanY = Math.max(maxY - minY, 60);

        scaleX = (mw - pad * 2) / spanX;
        scaleY = (mh - pad * 2) / spanY;
        minimapBounds = { minX, maxX, minY, maxY, scaleX, scaleY, pad };

        minimapCtx.strokeStyle = isDark ? "rgba(168, 85, 247, 0.4)" : "rgba(147, 51, 234, 0.4)";
        minimapCtx.lineWidth = 0.8;
        minimapCtx.strokeRect(pad + (120 - minX) * scaleX, pad + (70 - minY) * scaleY, 1000 * scaleX, 650 * scaleY);

        minimapCtx.fillStyle = isDark ? "#c084fc" : "#9333ea";
        const taxaEntries = Object.entries(umapPointCoords);
        const step = taxaEntries.length > 300 ? Math.ceil(taxaEntries.length / 250) : 1;
        for (let i = 0; i < taxaEntries.length; i += step) {
          const [tax, pt] = taxaEntries[i];
          if (pt && typeof pt.x === 'number' && typeof pt.y === 'number') {
            const nx = pad + (pt.x - minX) * scaleX;
            const ny = pad + (pt.y - minY) * scaleY;
            minimapCtx.fillRect(nx - 1, ny - 1, 2, 2);
          }
        }

      } else if (settings.layout === "tanglegram") {
        const leftRoot = tangleState.activeLeftRoot || rawRoot3Di;
        const rightRoot = tangleState.activeRightRoot || rawRootAA;
        const leftLeaves = (tangleState.leftLeaves && tangleState.leftLeaves.length > 0) ? tangleState.leftLeaves : getAllLeaves(leftRoot);
        const rightLeaves = (tangleState.rightLeaves && tangleState.rightLeaves.length > 0) ? tangleState.rightLeaves : getAllLeaves(rightRoot);

        if (!leftLeaves || leftLeaves.length === 0) return;

        minX = 40;
        maxX = 1170;
        const validY = [...leftLeaves, ...rightLeaves].map(l => l.y).filter(y => typeof y === 'number' && !isNaN(y));
        minY = (validY.length > 0) ? Math.min(...validY) - 10 : 40;
        maxY = (validY.length > 0) ? Math.max(...validY) + 10 : 600;

        const spanX = Math.max(maxX - minX, 100);
        const spanY = Math.max(maxY - minY, 60);

        scaleX = (mw - pad * 2) / spanX;
        scaleY = (mh - pad * 2) / spanY;

        minimapBounds = { minX, maxX, minY, maxY, scaleX, scaleY, pad };

        // 1. Draw Sampled Connecting Curves in Middle Channel
        const mapRight = {};
        rightLeaves.forEach(l => { mapRight[l.name] = l; });
        
        minimapCtx.strokeStyle = isDark ? "rgba(168, 85, 247, 0.3)" : "rgba(147, 51, 234, 0.4)";
        minimapCtx.lineWidth = 0.6;
        
        const step = leftLeaves.length > 150 ? Math.ceil(leftLeaves.length / 100) : 1;
        for (let i = 0; i < leftLeaves.length; i += step) {
          const l1 = leftLeaves[i];
          const l2 = mapRight[l1.name];
          if (l2 && typeof l1.y === 'number' && typeof l2.y === 'number') {
            const x1 = pad + (tangleState.leftConnectorX - minX) * scaleX;
            const y1 = pad + (l1.y - minY) * scaleY;
            const x2 = pad + (tangleState.rightConnectorX - minX) * scaleX;
            const y2 = pad + (l2.y - minY) * scaleY;
            minimapCtx.beginPath();
            minimapCtx.moveTo(x1, y1);
            minimapCtx.bezierCurveTo(x1 + (x2 - x1) * 0.5, y1, x2 - (x2 - x1) * 0.5, y2, x2, y2);
            minimapCtx.stroke();
          }
        }

        // 2. Draw Left Tree
        minimapCtx.strokeStyle = isDark ? "rgba(56, 189, 248, 0.75)" : "rgba(2, 132, 199, 0.8)";
        minimapCtx.lineWidth = 0.85;
        function drawMiniLeft(node) {
          if (!node || typeof node.x !== 'number' || typeof node.y !== 'number') return;
          const nx = pad + (node.x - minX) * scaleX;
          const ny = pad + (node.y - minY) * scaleY;

          if (!node.children || node.children.length === 0) {
            minimapCtx.fillStyle = isDark ? "#38bdf8" : "#0284c7";
            minimapCtx.beginPath();
            minimapCtx.arc(nx, ny, leftLeaves.length <= 12 ? 2.5 : 1.2, 0, 2 * Math.PI);
            minimapCtx.fill();
            return;
          }

          node.children.forEach(c => {
            if (typeof c.x === 'number' && typeof c.y === 'number') {
              const cx = pad + (c.x - minX) * scaleX;
              const cy = pad + (c.y - minY) * scaleY;
              minimapCtx.beginPath();
              minimapCtx.moveTo(nx, ny);
              minimapCtx.lineTo(nx, cy);
              minimapCtx.lineTo(cx, cy);
              minimapCtx.stroke();
              drawMiniLeft(c);
            }
          });
        }
        if (leftRoot) drawMiniLeft(leftRoot);

        // 3. Draw Right Tree
        minimapCtx.strokeStyle = isDark ? "rgba(168, 85, 247, 0.75)" : "rgba(147, 51, 234, 0.8)";
        minimapCtx.lineWidth = 0.85;
        function drawMiniRight(node) {
          if (!node || typeof node.x !== 'number' || typeof node.y !== 'number') return;
          const nx = pad + (node.x - minX) * scaleX;
          const ny = pad + (node.y - minY) * scaleY;

          if (!node.children || node.children.length === 0) {
            minimapCtx.fillStyle = isDark ? "#a855f7" : "#9333ea";
            minimapCtx.beginPath();
            minimapCtx.arc(nx, ny, rightLeaves.length <= 12 ? 2.5 : 1.2, 0, 2 * Math.PI);
            minimapCtx.fill();
            return;
          }

          node.children.forEach(c => {
            if (typeof c.x === 'number' && typeof c.y === 'number') {
              const cx = pad + (c.x - minX) * scaleX;
              const cy = pad + (c.y - minY) * scaleY;
              minimapCtx.beginPath();
              minimapCtx.moveTo(nx, ny);
              minimapCtx.lineTo(nx, cy);
              minimapCtx.lineTo(cx, cy);
              minimapCtx.stroke();
              drawMiniRight(c);
            }
          });
        }
        if (rightRoot) drawMiniRight(rightRoot);

      } else {
        // Normal single-tree minimap (Cartesian, Radial, Unrooted)
        let rootNode = activeTreeRoot;
        let visibleLeaves = getVisibleLeaves(rootNode);
        if (!visibleLeaves || visibleLeaves.length === 0) return;

        const validX = visibleLeaves.map(l => l.x).filter(x => typeof x === 'number' && !isNaN(x));
        const validY = visibleLeaves.map(l => l.y).filter(y => typeof y === 'number' && !isNaN(y));
        if (validX.length === 0 || validY.length === 0) return;

        minX = Math.min(...validX, rootNode.x || 50);
        maxX = Math.max(...validX) + 120;
        minY = Math.min(...validY);
        maxY = Math.max(...validY);

        const spanX = Math.max(maxX - minX, 60);
        const spanY = Math.max(maxY - minY, 60);

        if (settings.layout === "radial" || settings.layout === "unrooted") {
          const uniformScale = Math.min((mw - pad * 2) / spanX, (mh - pad * 2) / spanY);
          scaleX = uniformScale;
          scaleY = uniformScale;
        } else {
          scaleX = (mw - pad * 2) / spanX;
          scaleY = (mh - pad * 2) / spanY;
        }

        minimapBounds = { minX, maxX, minY, maxY, scaleX, scaleY, pad };

        minimapCtx.strokeStyle = isDark ? "rgba(100, 116, 139, 0.45)" : "rgba(148, 163, 184, 0.6)";
        minimapCtx.lineWidth = 1.0;

        function drawMiniBranch(node) {
          if (!node || typeof node.x !== 'number' || typeof node.y !== 'number') return;
          const nx = pad + (node.x - minX) * scaleX;
          const ny = pad + (node.y - minY) * scaleY;

          if (node._collapsed) {
            minimapCtx.fillStyle = getNodeColor(node);
            minimapCtx.beginPath();
            minimapCtx.arc(nx, ny, 3.0, 0, 2 * Math.PI);
            minimapCtx.fill();
            return;
          }

          if (!node.children || node.children.length === 0) {
            minimapCtx.fillStyle = isDark ? "#38bdf8" : "#0284c7";
            minimapCtx.beginPath();
            minimapCtx.arc(nx, ny, visibleLeaves.length <= 10 ? 3.0 : 1.5, 0, 2 * Math.PI);
            minimapCtx.fill();
            return;
          }

          node.children.forEach(c => {
            if (typeof c.x === 'number' && typeof c.y === 'number') {
              const cx = pad + (c.x - minX) * scaleX;
              const cy = pad + (c.y - minY) * scaleY;
              minimapCtx.beginPath();
              minimapCtx.moveTo(nx, ny);
              minimapCtx.lineTo(nx, cy);
              minimapCtx.lineTo(cx, cy);
              minimapCtx.stroke();
              drawMiniBranch(c);
            }
          });
        }
        drawMiniBranch(rootNode);
      }

      // Viewport Rectangle (Shared Across All Layouts)
      const treeWrapper = document.getElementById("treeCanvasWrapper") || document.getElementById("treeContainer");
      const W = treeWrapper ? (treeWrapper.clientWidth || 900) : 900;
      const H = treeWrapper ? (treeWrapper.clientHeight || 700) : 700;

      const viewTreeLeft = (0 - settings.zoom.x) / settings.zoom.k;
      const viewTreeTop = (0 - settings.zoom.y) / settings.zoom.k;
      const viewTreeRight = (W - settings.zoom.x) / settings.zoom.k;
      const viewTreeBottom = (H - settings.zoom.y) / settings.zoom.k;

      const miniVpX = pad + (viewTreeLeft - minX) * scaleX;
      const miniVpY = pad + (viewTreeTop - minY) * scaleY;
      const miniVpW = (viewTreeRight - viewTreeLeft) * scaleX;
      const miniVpH = (viewTreeBottom - viewTreeTop) * scaleY;

      const leafCount = (settings.layout === "tanglegram")
        ? ((tangleState.leftLeaves && tangleState.leftLeaves.length) || 6)
        : (settings.layout === "umap" ? Object.keys(umapPointCoords).length : getVisibleLeaves(activeTreeRoot).length);

      const allInView = (miniVpX <= pad && miniVpY <= pad && (miniVpX + miniVpW) >= (mw - pad) && (miniVpY + miniVpH) >= (mh - pad));

      if (allInView && leafCount <= 12) {
        minimapViewport.style.left = "4px";
        minimapViewport.style.top = "4px";
        minimapViewport.style.width = (mw - 8) + "px";
        minimapViewport.style.height = (mh - 8) + "px";
        minimapViewport.style.borderColor = "rgba(56, 189, 248, 0.4)";
        minimapViewport.style.backgroundColor = "transparent";
      } else {
        const clampedX = Math.max(0, Math.min(miniVpX, mw));
        const clampedY = Math.max(0, Math.min(miniVpY, mh));
        const clampedW = Math.max(10, Math.min(miniVpW, mw - clampedX));
        const clampedH = Math.max(10, Math.min(miniVpH, mh - clampedY));

        minimapViewport.style.left = clampedX + "px";
        minimapViewport.style.top = clampedY + "px";
        minimapViewport.style.width = clampedW + "px";
        minimapViewport.style.height = clampedH + "px";
        minimapViewport.style.borderColor = "#38bdf8";
        minimapViewport.style.backgroundColor = "rgba(56, 189, 248, 0.2)";
      }
    }

    // Minimap Click & Drag Navigation
    let isMinimapNavigating = false;
    function navigateFromMinimap(e) {
      const rect = minimapBox.getBoundingClientRect();
      const clickX = e.clientX - rect.left;
      const clickY = e.clientY - rect.top;

      const b = minimapBounds;
      if (!b || !b.scaleX || !b.scaleY) return;

      const targetTreeX = b.minX + (clickX - b.pad) / b.scaleX;
      const targetTreeY = b.minY + (clickY - b.pad) / b.scaleY;

      const treeWrapper = document.getElementById("treeCanvasWrapper") || document.getElementById("treeContainer");
      const W = treeWrapper ? (treeWrapper.clientWidth || 900) : 900;
      const H = treeWrapper ? (treeWrapper.clientHeight || 700) : 700;

      settings.zoom.x = W / 2 - targetTreeX * settings.zoom.k;
      settings.zoom.y = H / 2 - targetTreeY * settings.zoom.k;

      updateViewportTransform(false);
      updateMinimap();
    }

    minimapBox.addEventListener("mousedown", (e) => {
      if (e.button !== 0) return;
      isMinimapNavigating = true;
      navigateFromMinimap(e);
    });

    window.addEventListener("mousemove", (e) => {
      if (isMinimapNavigating) {
        navigateFromMinimap(e);
      }
    });

    window.addEventListener("mouseup", () => {
      if (isMinimapNavigating) {
        isMinimapNavigating = false;
        updateMinimap();
      }
    });

    // 13. DATASET SWITCHING & CONTROLS INTERACTION
    function setLayout(layoutName) {
      settings.layout = layoutName;
      ["btnRect", "btnClado", "btnRadial", "btnUnrooted", "btnTangle", "btnUmap"].forEach(id => {
        const btn = document.getElementById(id);
        if (btn) btn.className = "py-1.5 rounded font-medium text-center hover:bg-slate-500/20 text-[var(--text-muted)] text-[10.5px] transition";
      });
      const activeBtnMap = {
        "rectangular": "btnRect",
        "cladogram": "btnClado",
        "radial": "btnRadial",
        "unrooted": "btnUnrooted",
        "tanglegram": "btnTangle",
        "umap": "btnUmap"
      };
      if (activeBtnMap[layoutName]) {
        const activeBtn = document.getElementById(activeBtnMap[layoutName]);
        if (activeBtn) {
          activeBtn.className = (layoutName === "umap")
            ? "py-1.5 rounded font-medium text-center bg-purple-600 text-white text-[10.5px] transition shadow-sm cursor-pointer"
            : "py-1.5 rounded font-medium text-center bg-sky-500 text-white text-[10.5px] transition";
        }
      }

      // Sync ESM-2 representation toggle buttons
      const btnEsmTree = document.getElementById("btnEsmViewTree");
      const btnEsmUmap = document.getElementById("btnEsmViewUmap");
      if (btnEsmTree && btnEsmUmap) {
        if (layoutName === "umap") {
          btnEsmTree.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-purple-500/20 text-[10.5px] transition cursor-pointer flex items-center justify-center space-x-1";
          btnEsmUmap.className = "py-1 px-2 rounded font-semibold text-center bg-purple-600 text-white text-[10.5px] transition shadow-sm cursor-pointer flex items-center justify-center space-x-1";
        } else {
          btnEsmTree.className = "py-1 px-2 rounded font-semibold text-center bg-purple-600 text-white text-[10.5px] transition shadow-sm cursor-pointer flex items-center justify-center space-x-1";
          btnEsmUmap.className = "py-1 px-2 rounded font-medium text-center text-[var(--text-muted)] hover:text-white hover:bg-purple-500/20 text-[10.5px] transition cursor-pointer flex items-center justify-center space-x-1";
        }
      }

      // Toggle single-tree dataset selector vs tanglegram alignment options
      const dsSec = document.getElementById("treeDatasetSection");
      const tangleOpt = document.getElementById("tanglegramOptions");
      if (layoutName === "tanglegram") {
        if (dsSec) dsSec.classList.add("hidden");
        if (tangleOpt) tangleOpt.classList.remove("hidden");
      } else {
        if (dsSec) dsSec.classList.remove("hidden");
        if (tangleOpt) tangleOpt.classList.add("hidden");
      }

      updateLayoutSpecificControls();
      renderTree();
      fitTreeToScreen(true);
    }

    function updateLayoutSpecificControls() {
      const mode = settings.layout;
      const rectCard = document.getElementById("rectControlsCard");
      const radialCard = document.getElementById("radialControlsCard");
      const unrootedCard = document.getElementById("unrootedControlsCard");
      const umapCard = document.getElementById("umapControlsCard");

      if (rectCard) {
        if (mode === "rectangular" || mode === "cladogram") rectCard.classList.remove("hidden");
        else rectCard.classList.add("hidden");
      }
      if (radialCard) {
        if (mode === "radial") radialCard.classList.remove("hidden");
        else radialCard.classList.add("hidden");
      }
      if (unrootedCard) {
        if (mode === "unrooted") unrootedCard.classList.remove("hidden");
        else unrootedCard.classList.add("hidden");
      }
      if (umapCard) {
        if (mode === "umap") umapCard.classList.remove("hidden");
        else umapCard.classList.add("hidden");
      }
    }

    function setTangleMode(mode) {
      settings.tangleMode = mode;
      const descMap = {
        "true_topology": "True topology preserves native IQ-TREE branch order on both sides to expose structural vs sequence discordance.",
        "min_crossings": "Rotates internal clades of the sequence tree to minimize crossing tangles while preserving 100% of phylogenetic clades.",
        "aligned": "Orders right-hand leaves directly alongside matching left-hand leaves for parallel 1-to-1 visual comparison."
      };
      const descEl = document.getElementById("tangleModeDesc");
      if (descEl) descEl.textContent = descMap[mode] || "";
      updateCongruenceUI();
      renderTree();
      updateMinimap();
    }

    function switchDataset(dsName) {
      settings.dataset = dsName;
      if (dsName === "3di" && typeof msaState !== "undefined" && msaState.mode !== "3di") {
        setMsaMode("3di");
      } else if (dsName === "aa" && typeof msaState !== "undefined" && msaState.mode !== "aa") {
        setMsaMode("aa");
      }
      let label = "3Di Tree";
      const embedSub = document.getElementById("embedMetricSubSection");
      if (dsName === "aa") {
        label = "AA Tree";
        if (embedSub) embedSub.classList.add("hidden");
      } else if (dsName === "esm2") {
        const metricName = (settings.embedMetric === "euclidean") ? "Euclidean" : (settings.embedMetric === "l1" ? "L1 / Manhattan" : "Cosine");
        label = `ESM-2 (${metricName})`;
        if (embedSub) embedSub.classList.remove("hidden");
      } else {
        if (embedSub) embedSub.classList.add("hidden");
      }
      const b = document.getElementById("badgeTree");
      if (b) b.textContent = label;
      applyCurrentRooting();
    }

    function switchEmbedMetric(metric) {
      settings.embedMetric = metric;
      const descMap = {
        "cosine": "Cosine distance measures angular alignment in the 1280-dim PLM representation space, robust to overall norm shifts.",
        "euclidean": "Euclidean distance measures absolute L2 geometric vector displacement between PLM sequence representations.",
        "l1": "Manhattan / L1 distance measures the sum of absolute coordinate differences across all 1280 PLM embedding dimensions."
      };
      const descEl = document.getElementById("embedMetricDesc");
      if (descEl) descEl.textContent = descMap[metric] || "";
      if (settings.dataset === "esm2") {
        const metricName = (metric === "euclidean") ? "Euclidean" : (metric === "l1" ? "L1 / Manhattan" : "Cosine");
        const b = document.getElementById("badgeTree");
        if (b) b.textContent = `ESM-2 (${metricName})`;
        applyCurrentRooting();
      }
    }

    function setTangleCompare(compareMode) {
      settings.tangleCompare = compareMode;
      updateCongruenceUI(); // Fix bug: call updateCongruenceUI directly!
      renderTree();
      updateMinimap();
    }

    function switchDatasetScale(scale) {
      currentScale = scale;
      const ds = DATASETS[scale];
      activeDataset = ds; // Update activeDataset reference!
      NEWICK_3DI = ds.newick_3di;
      NEWICK_AA = ds.newick_aa;
      NEWICK_ESM2 = ds.newick_esm2 || ds.newick_esm2_cosine;
      NEWICK_ESM2_COSINE = ds.newick_esm2_cosine;
      NEWICK_ESM2_EUCLIDEAN = ds.newick_esm2_euclidean;
      NEWICK_ESM2_L1 = ds.newick_esm2_l1;
      if (typeof syncDatasetAlignmentCoverage === 'function') {
        syncDatasetAlignmentCoverage(scale);
      }
      TAXA_METADATA = ds.taxa;

      parseAllActiveTrees();
      updateModalityOptions();
      updatePaletteMiniStrip();

      settings.verticalSpacing = ds.defaultSpacing;
      settings.nodeRadius = ds.defaultRadius;
      settings.selectedTaxon = null;
      settings.colorColumn = ds.defaultColorCol || (ds.columns && ds.columns[0].key) || "family";
      settings.cladeGroupColumn = ds.defaultCladeCol || (ds.columns && ds.columns.find(c => c.type === 'categorical')?.key) || "family";

      // Reset scopedClade cleanly on scale change
      if (settings.scopedClade) {
        settings.scopedClade = null;
        const banner = document.getElementById("scopedCladeBanner");
        if (banner) {
          banner.classList.add("hidden");
          banner.classList.remove("flex");
        }
        const scopeBadge = document.getElementById("cladeScopeActiveBadge");
        if (scopeBadge) {
          scopeBadge.textContent = "Full Cohort";
          scopeBadge.className = "text-[9px] font-mono px-2 py-0.5 rounded-full bg-slate-700/60 text-slate-400 border border-slate-600/40";
        }
      }
      const activeSil = ds["silhouette_" + (cladePartitionState.source || "esm2")] || ds.silhouette;
      if (activeSil && activeSil.best_k) {
        cladePartitionState.k = activeSil.best_k;
      } else {
        cladePartitionState.k = 3;
      }
      const kSlider = document.getElementById("cladeKSlider");
      if (kSlider) {
        const maxK = activeSil && activeSil.profile ? Math.max(...activeSil.profile.map(p => p.k)) : (scale === "6" ? 5 : 35);
        kSlider.max = maxK;
        kSlider.value = cladePartitionState.k;
      }
      const kValEl = document.getElementById("cladeKVal");
      if (kValEl) kValEl.textContent = `k = ${cladePartitionState.k}`;

      // Reset filterState cleanly to active dataset columns
      filterState.isActive = false;
      const activeCols = ds.columns || [];
      const firstCatCol = activeCols.find(c => c.type === "categorical") || activeCols[0];
      filterState.column = firstCatCol ? firstCatCol.key : "family";
      filterState.selectedCategories = new Set();
      filterState.minVal = null;
      filterState.maxVal = null;
      filterState.categorySearchQuery = "";
      const banner = document.getElementById("activeFilterBanner");
      if (banner) banner.classList.add("hidden");

      document.getElementById("spacingSlider").value = settings.verticalSpacing;
      document.getElementById("spacingVal").textContent = settings.verticalSpacing + "px";
      document.getElementById("radiusSlider").value = settings.nodeRadius;
      document.getElementById("radiusVal").textContent = settings.nodeRadius + "px";

      // Adapt tip label size to cohort scale
      if (scale === "1193") {
        settings.labelSize = 8;
      } else if (scale === "500") {
        settings.labelSize = 9;
      } else {
        settings.labelSize = 10;
      }
      const lblSlider = document.getElementById("labelSizeSlider");
      if (lblSlider) lblSlider.value = settings.labelSize;
      const lblVal = document.getElementById("labelSizeVal");
      if (lblVal) lblVal.textContent = settings.labelSize + "px";
      updateLabelScaling();
      msaState._minimapCacheKey = null;

      if (scale === "1193") {
        document.getElementById("badgeTaxa").textContent = "1,193 ESMFold Designs";
      } else if (scale === "500") {
        document.getElementById("badgeTaxa").textContent = "500 Viral Structures";
      } else if (scale === "100") {
        document.getElementById("badgeTaxa").textContent = "100 RdRp Structures";
      } else {
        document.getElementById("badgeTaxa").textContent = "6 Benchmark Taxa";
      }

      const scaleSel = document.getElementById("scaleSelect");
      if (scaleSel && scaleSel.value !== scale) scaleSel.value = scale;

      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      updateFilterUI();
      applyCurrentRooting();
      updateCladeManagementUI();
      updateCongruenceUI();
      msaState.scrollX = 0;
      msaState.scrollY = 0;
      renderMsa();

      const firstTaxon = Object.keys(TAXA_METADATA)[0];
      if (firstTaxon) selectTaxon(firstTaxon);
    }

    function setTreeRotation(deg) {
      settings.treeRotation = parseFloat(deg) || 0;
      const rBadge = document.getElementById("treeRotationVal");
      if (rBadge) rBadge.textContent = Math.round(settings.treeRotation) + "°";
      const uBadge = document.getElementById("unrootedRotationVal");
      if (uBadge) uBadge.textContent = Math.round(settings.treeRotation) + "°";
      renderTree();
      updateMinimap();
    }

    function setRadialArc(deg) {
      settings.radialArc = parseFloat(deg) || 360;
      const aBadge = document.getElementById("radialArcVal");
      if (aBadge) aBadge.textContent = Math.round(settings.radialArc) + "°";
      renderTree();
      updateMinimap();
    }

    function setRadialRadiusScale(val) {
      settings.radialRadiusScale = parseFloat(val) || 1.0;
      const sBadge = document.getElementById("radialRadiusScaleVal");
      if (sBadge) sBadge.textContent = settings.radialRadiusScale.toFixed(1) + "x";
      renderTree();
      updateMinimap();
    }

    function setUnrootedScale(val) {
      settings.unrootedScale = parseFloat(val) || 1.0;
      const sBadge = document.getElementById("unrootedScaleVal");
      if (sBadge) sBadge.textContent = settings.unrootedScale.toFixed(1) + "x";
      renderTree();
      updateMinimap();
    }

    function setLabelOrientation(mode) {
      settings.labelOrientation = mode;
      const selR = document.getElementById("radialLabelOrientationSelect");
      if (selR) selR.value = mode;
      const selU = document.getElementById("unrootedLabelOrientationSelect");
      if (selU) selU.value = mode;
      renderTree();
    }

    function setBranchWidth(w) {
      settings.branchWidth = parseFloat(w) || 1.4;
      const bBadge = document.getElementById("branchWidthVal");
      if (bBadge) bBadge.textContent = settings.branchWidth.toFixed(1) + "px";
      document.documentElement.style.setProperty("--branch-width", settings.branchWidth + "px");
      renderTree();
    }

    function setSpacing(v) {
      settings.verticalSpacing = parseInt(v);
      document.getElementById("spacingVal").textContent = v + "px";
      renderTree();
    }

    function setNodeRadius(r) {
      settings.nodeRadius = parseFloat(r);
      document.getElementById("radiusVal").textContent = r + "px";
      renderTree();
    }

    function setLabelSize(v) {
      settings.labelSize = parseFloat(v);
      const valEl = document.getElementById("labelSizeVal");
      if (valEl) valEl.textContent = settings.labelSize === 0 ? "Hidden" : `${settings.labelSize}px`;
      updateLabelScaling();
    }

    function toggleBackgroundGrid(forceVal) {
      if (typeof forceVal === "boolean") {
        settings.showGrid = forceVal;
      } else {
        settings.showGrid = !settings.showGrid;
      }
      const chk = document.getElementById("toggleShowGrid");
      if (chk) chk.checked = settings.showGrid;

      const svgEl = document.getElementById("treeSvg");
      if (svgEl) {
        if (settings.showGrid) {
          svgEl.classList.remove("no-grid");
        } else {
          svgEl.classList.add("no-grid");
        }
      }

      const btnIcon = document.getElementById("btnToggleGridIcon");
      const btnLabel = document.getElementById("btnToggleGridLabel");
      const btn = document.getElementById("btnToggleGrid");
      if (btn) {
        if (settings.showGrid) {
          btn.classList.remove("text-rose-400");
          if (btnIcon) btnIcon.textContent = "🌐";
          if (btnLabel) btnLabel.textContent = "Grid";
          btn.title = "Hide background grid lines";
        } else {
          btn.classList.add("text-rose-400");
          if (btnIcon) btnIcon.textContent = "🚫";
          if (btnLabel) btnLabel.textContent = "No Grid";
          btn.title = "Show background grid lines";
        }
      }
    }

    function setUnrootedLengthMode(mode) {
      settings.unrootedLengthMode = mode;
      renderTree();
    }

    function toggleUnrootedDaylight(val) {
      settings.unrootedDaylight = Boolean(val);
      renderTree();
    }

    function setUnrootedLabelFilter(val) {
      settings.unrootedLabelFilter = val;
      renderTree();
    }

    function toggleStaggerLabels(val) {
      settings.staggerLabels = Boolean(val);
      renderTree();
    }

    function toggleSetting(key, val) {
      settings[key] = val;
      if (key === 'zoomAdaptiveLabels') {
        updateLabelScaling();
      }
      renderTree();
    }

    function toggleSupportColoring(val) {
      settings.colorBranchesBySupport = Boolean(val);
      const ctrl = document.getElementById("supportControls");
      const badge = document.getElementById("supportBadge");
      if (ctrl) {
        if (settings.colorBranchesBySupport) {
          ctrl.classList.remove("hidden");
        } else {
          ctrl.classList.add("hidden");
        }
      }
      if (badge) {
        badge.textContent = settings.colorBranchesBySupport ? "Active" : "Off";
        if (settings.colorBranchesBySupport) {
          badge.className = "text-[9px] font-mono px-1.5 py-0.5 rounded bg-emerald-500/20 text-emerald-400 font-semibold";
        } else {
          badge.className = "text-[9px] font-mono px-1.5 py-0.5 rounded bg-slate-700/50 text-[var(--text-muted)]";
        }
      }
      updateLegend();
      renderTree();
    }

    function setSupportMetric(m) {
      settings.supportMetric = m;
      updateLegend();
      renderTree();
    }

    function setSupportPalette(p) {
      settings.supportPalette = p;
      const preview = document.getElementById("supportPalettePreview");
      const pLow = document.getElementById("previewLow");
      const pMid = document.getElementById("previewMid");
      const pHigh = document.getElementById("previewHigh");
      if (preview) {
        if (p === "traffic") {
          if (pLow) pLow.textContent = "<70%";
          if (pMid) pMid.textContent = "70-94%";
          if (pHigh) pHigh.textContent = "≥95%";
          preview.innerHTML = `
            <div class="h-full bg-rose-500" style="width:33.3%"></div>
            <div class="h-full bg-amber-400" style="width:33.3%"></div>
            <div class="h-full bg-emerald-500" style="width:33.4%"></div>
          `;
        } else if (p === "viridis") {
          if (pLow) pLow.textContent = "0%";
          if (pMid) pMid.textContent = "50%";
          if (pHigh) pHigh.textContent = "100%";
          preview.innerHTML = `
            <div class="h-full w-full" style="background: linear-gradient(90deg, rgb(68,1,84), rgb(59,82,139), rgb(33,145,140), rgb(94,201,98), rgb(253,231,37));"></div>
          `;
        } else if (p === "grayscale") {
          if (pLow) pLow.textContent = "Fade";
          if (pMid) pMid.textContent = "Mid";
          if (pHigh) pHigh.textContent = "Full";
          preview.innerHTML = `
            <div class="h-full w-full" style="background: linear-gradient(90deg, rgba(56,189,248,0.18), rgba(56,189,248,1.0));"></div>
          `;
        }
      }
      updateLegend();
      renderTree();
    }

    function handleSearch(query) {
      const q = query.trim().toLowerCase();
      if (!q) {
        document.querySelectorAll(".tip-label").forEach(el => {
          el.style.opacity = "1";
          el.classList.remove("selected");
        });
        return;
      }
      let firstMatch = null;
      document.querySelectorAll(".tip-label").forEach(el => {
        const taxonName = el.getAttribute("data-taxon") || el.textContent;
        const displayText = el.textContent;
        const meta = TAXA_METADATA[taxonName] || {};
        
        let match = taxonName.toLowerCase().includes(q) || displayText.toLowerCase().includes(q);
        if (!match) {
          for (let key in meta) {
            if (String(meta[key]).toLowerCase().includes(q)) {
              match = true;
              break;
            }
          }
        }

        if (match) {
          el.style.opacity = "1";
          el.classList.add("selected");
          if (!firstMatch) firstMatch = taxonName;
        } else {
          el.style.opacity = "0.2";
          el.classList.remove("selected");
        }
      });
      if (firstMatch) selectTaxon(firstMatch);
    }

    function populateOutgroupSelect() {
      const sel = document.getElementById("outgroupSelect");
      if (!sel) return;
      sel.innerHTML = "";
      const names = Object.keys(TAXA_METADATA);
      names.forEach(n => {
        const opt = document.createElement("option");
        opt.value = n;
        const m = TAXA_METADATA[n];
        const colDef = getActiveColorColumnDef();
        const primaryVal = (colDef && m && m[colDef.key]) ? ` (${m[colDef.key]})` : "";
        opt.textContent = `${n}${primaryVal}`;
        sel.appendChild(opt);
      });
    }

    function copyTextWithFallback(text, onSuccess, onFallback) {
      if (navigator.clipboard && typeof navigator.clipboard.writeText === 'function') {
        navigator.clipboard.writeText(text).then(() => {
          if (typeof onSuccess === 'function') onSuccess();
        }).catch(() => {
          fallbackExecCopy(text, onSuccess, onFallback);
        });
      } else {
        fallbackExecCopy(text, onSuccess, onFallback);
      }
    }

    function fallbackExecCopy(text, onSuccess, onFallback) {
      try {
        const ta = document.createElement("textarea");
        ta.value = text;
        ta.style.position = "fixed";
        ta.style.left = "-9999px";
        ta.style.top = "-9999px";
        document.body.appendChild(ta);
        ta.focus();
        ta.select();
        const success = document.execCommand("copy");
        document.body.removeChild(ta);
        if (success && typeof onSuccess === 'function') {
          onSuccess();
          return;
        }
      } catch (e) {}
      if (typeof onFallback === 'function') onFallback();
    }

    function getActiveNewickString() {
      if (settings.scopedClade && settings.scopedClade.taxa && settings.scopedClade.taxa.size > 0) {
        const rawRoot = (settings.dataset === "3di") ? rawRoot3Di : ((settings.dataset === "aa") ? rawRootAA : rawRootESM2);
        if (rawRoot && typeof pruneSubtree === 'function') {
          const pruned = pruneSubtree(rawRoot, settings.scopedClade.taxa);
          if (pruned) return serializeNewick(pruned);
        }
      }
      if (settings.dataset === "3di") return NEWICK_3DI || "";
      if (settings.dataset === "aa") return NEWICK_AA || "";
      if (settings.dataset === "esm2_euclidean") return (typeof NEWICK_ESM2_EUCLIDEAN !== 'undefined' && NEWICK_ESM2_EUCLIDEAN) || (typeof NEWICK_ESM2 !== 'undefined' && NEWICK_ESM2) || "";
      if (settings.dataset === "esm2_l1") return (typeof NEWICK_ESM2_L1 !== 'undefined' && NEWICK_ESM2_L1) || (typeof NEWICK_ESM2 !== 'undefined' && NEWICK_ESM2) || "";
      return (typeof NEWICK_ESM2_COSINE !== 'undefined' && NEWICK_ESM2_COSINE) || (typeof NEWICK_ESM2 !== 'undefined' && NEWICK_ESM2) || "";
    }

    function exportNewick() {
      try {
        const nwk = getActiveNewickString();
        if (!nwk) {
          showToastNotification("⚠️ No Newick tree data available for active modality.");
          return;
        }

        const filename = `tree_${settings.dataset}_cohort_${currentScale}${settings.scopedClade ? '_subclade' : ''}.nwk`;
        const blob = new Blob([nwk], { type: "text/plain;charset=utf-8" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        setTimeout(() => URL.revokeObjectURL(url), 1000);

        copyTextWithFallback(nwk, () => {
          showToastNotification(`📋 Exported <strong>${filename}</strong> & copied Newick to clipboard!`);
        }, () => {
          showToastNotification(`📋 Exported and downloaded <strong>${filename}</strong>!`);
        });
      } catch (err) {
        console.error("Newick export error:", err);
        showToastNotification("⚠️ Error exporting Newick: " + err.message);
      }
    }

    function copyNewick() {
      exportNewick();
    }

    function exportSVG() {
      const s = new XMLSerializer().serializeToString(svg);
      const blob = new Blob([s], { type: "image/svg+xml;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `viral_tree_${settings.dataset}_${settings.layout}.svg`;
      a.click();
      URL.revokeObjectURL(url);
    }

    // INITIALIZATION
    toggleBackgroundGrid(settings.showGrid);
    window.addEventListener("DOMContentLoaded", () => {
      try {
        const savedPal = localStorage.getItem("phylo_custom_palette");
        if (savedPal) {
          const parsed = JSON.parse(savedPal);
          if (parsed && Array.isArray(parsed.palette) && parsed.palette.length >= 2) {
            customPaletteState = Object.assign(customPaletteState, parsed);
          }
        }
      } catch(e) {}
      initHeaderLogo();
      updatePaletteMiniStrip();
      initCustomThemeFromStorage();
      setTheme(settings.theme);
      updateModalityOptions();
      populateMetadataSelectors();
      populateOutgroupSelect();
      updateLegend();
      updateFilterUI();
      applyCurrentRooting();
      updateCladeManagementUI();
      updateCongruenceUI();
      updateLabelScaling();
      updateLayoutSpecificControls();
      updatePipelineCommandPreview();
      setupMsaInteractions();
      renderMsa();

      const firstTaxon = Object.keys(TAXA_METADATA)[0];
      if (firstTaxon) {
        selectTaxon(firstTaxon);
      }
    });


    // =========================================================================
    // PURE JAVASCRIPT ZERO-DEPENDENCY ZIP ARCHIVE GENERATOR (STORE METHOD)
    // =========================================================================
    function createZipArchive(files) {
      const crcTable = new Uint32Array(256);
      for (let i = 0; i < 256; i++) {
        let c = i;
        for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
        crcTable[i] = c >>> 0;
      }
      function calcCrc(buf) {
        let crc = 0xFFFFFFFF;
        for (let i = 0; i < buf.length; i++) crc = (crc >>> 8) ^ crcTable[(crc ^ buf[i]) & 0xFF];
        return (crc ^ 0xFFFFFFFF) >>> 0;
      }

      const encoder = new TextEncoder();
      const fileEntries = [];
      let offset = 0;
      const parts = [];

      for (const f of files) {
        const nameBytes = encoder.encode(f.name);
        const dataBytes = (typeof f.content === 'string') ? encoder.encode(f.content) : f.content;
        const crc = calcCrc(dataBytes);
        const size = dataBytes.length;

        // Local file header (30 bytes + name length)
        const lh = new Uint8Array(30 + nameBytes.length);
        const dv = new DataView(lh.buffer);
        dv.setUint32(0, 0x04034b50, true);
        dv.setUint16(4, 20, true);
        dv.setUint16(6, 0, true);
        dv.setUint16(8, 0, true);
        dv.setUint16(10, 0x4821, true);
        dv.setUint16(12, 0x5821, true);
        dv.setUint32(14, crc, true);
        dv.setUint32(18, size, true);
        dv.setUint32(22, size, true);
        dv.setUint16(26, nameBytes.length, true);
        dv.setUint16(28, 0, true);
        lh.set(nameBytes, 30);

        parts.push(lh);
        parts.push(dataBytes);

        fileEntries.push({ nameBytes, crc, size, offset });
        offset += lh.length + dataBytes.length;
      }

      const cdStart = offset;
      for (const fe of fileEntries) {
        const cd = new Uint8Array(46 + fe.nameBytes.length);
        const dv = new DataView(cd.buffer);
        dv.setUint32(0, 0x02014b50, true);
        dv.setUint16(4, 20, true);
        dv.setUint16(6, 20, true);
        dv.setUint16(8, 0, true);
        dv.setUint16(10, 0, true);
        dv.setUint16(12, 0x4821, true);
        dv.setUint16(14, 0x5821, true);
        dv.setUint32(16, fe.crc, true);
        dv.setUint32(20, fe.size, true);
        dv.setUint32(24, fe.size, true);
        dv.setUint16(28, fe.nameBytes.length, true);
        dv.setUint16(30, 0, true);
        dv.setUint16(32, 0, true);
        dv.setUint16(34, 0, true);
        dv.setUint16(36, 0, true);
        dv.setUint32(38, 0x81a40000, true);
        dv.setUint32(42, fe.offset, true);
        cd.set(fe.nameBytes, 46);
        parts.push(cd);
        offset += cd.length;
      }

      const cdSize = offset - cdStart;
      const eocd = new Uint8Array(22);
      const dvEocd = new DataView(eocd.buffer);
      dvEocd.setUint32(0, 0x06054b50, true);
      dvEocd.setUint16(4, 0, true);
      dvEocd.setUint16(6, 0, true);
      dvEocd.setUint16(8, fileEntries.length, true);
      dvEocd.setUint16(10, fileEntries.length, true);
      dvEocd.setUint32(12, cdSize, true);
      dvEocd.setUint32(16, cdStart, true);
      dvEocd.setUint16(20, 0, true);
      parts.push(eocd);

      return new Blob(parts, { type: 'application/zip' });
    }

    function serializeNewick(node) {
      if (!node) return ";";
      function stringify(n) {
        if (!n.children || n.children.length === 0) {
          let s = n.name || "";
          if (n.length !== undefined && n.length !== null) s += ":" + (typeof n.length === 'number' ? n.length.toFixed(6) : n.length);
          return s;
        }
        const childrenStr = n.children.map(stringify).join(",");
        let s = "(" + childrenStr + ")";
        if (n.name) s += n.name;
        if (n.length !== undefined && n.length !== null) s += ":" + (typeof n.length === 'number' ? n.length.toFixed(6) : n.length);
        return s;
      }
      return stringify(node) + ";";
    }

    function showToastNotification(msg, duration = 3500) {
      let toast = document.getElementById("globalToast");
      if (!toast) {
        toast = document.createElement("div");
        toast.id = "globalToast";
        toast.className = "fixed bottom-6 left-1/2 -translate-x-1/2 z-50 px-4 py-2.5 rounded-xl bg-slate-900/95 border border-sky-400 text-sky-200 text-xs font-semibold shadow-2xl backdrop-blur-md transition-all duration-300 opacity-0 pointer-events-none transform translate-y-2";
        document.body.appendChild(toast);
      }
      toast.innerHTML = msg;
      toast.classList.remove("opacity-0", "translate-y-2", "pointer-events-none");
      toast.classList.add("opacity-100", "translate-y-0");
      clearTimeout(toast._timer);
      toast._timer = setTimeout(() => {
        toast.classList.remove("opacity-100", "translate-y-0");
        toast.classList.add("opacity-0", "translate-y-2", "pointer-events-none");
      }, duration);
    }

    function generateStandaloneSubcladeViewer(label, targetTaxa, subsetStructs, metaObj) {
      const taxaJson = JSON.stringify(targetTaxa);
      const structsJson = JSON.stringify(subsetStructs);
      const metaJson = JSON.stringify(metaObj);
      const openScriptTag = "<" + "script>";
      const closeScriptTag = "<" + "/script>";
      return `<!DOCTYPE html>
<html lang="en" class="dark">
<head>
  <meta charset="UTF-8">
  <title>Subclade Analysis Viewer - ${label}</title>
  <script src="https://cdn.tailwindcss.com">${closeScriptTag}
  <style>
    body { background-color: #0b1120; color: #f8fafc; font-family: ui-sans-serif, system-ui, sans-serif; }
    .custom-scroll::-webkit-scrollbar { width: 6px; height: 6px; }
    .custom-scroll::-webkit-scrollbar-thumb { background: #334155; border-radius: 4px; }
  </style>
</head>
<body class="h-screen flex flex-col overflow-hidden">
  <header class="p-3 bg-slate-900 border-b border-slate-800 flex items-center justify-between shrink-0 shadow-md">
    <div class="flex items-center space-x-3">
      <span class="w-3 h-3 rounded-full bg-sky-400 animate-pulse"></span>
      <h1 class="text-sm font-bold text-white tracking-wide">Subclade Package Viewer: <span class="text-sky-400">${label}</span></h1>
      <span class="text-xs px-2.5 py-0.5 rounded-full bg-sky-500/20 text-sky-300 border border-sky-500/30 font-mono">${targetTaxa.length} Taxa</span>
    </div>
    <div class="text-[11px] text-slate-400">Offline Standalone Bundle</div>
  </header>
  <div class="flex-1 flex overflow-hidden">
    <div class="w-80 border-r border-slate-800 bg-slate-900/60 p-3 flex flex-col space-y-2 shrink-0">
      <input id="subcladeSearch" type="text" placeholder="Search taxa ID / metadata..." class="w-full bg-slate-950 border border-slate-700 rounded-lg px-2.5 py-1.5 text-xs text-white placeholder-slate-500 focus:outline-none focus:border-sky-500">
      <div id="taxaList" class="flex-1 overflow-y-auto space-y-1 custom-scroll pr-1"></div>
    </div>
    <div class="flex-1 flex flex-col p-4 bg-slate-950 overflow-hidden space-y-3">
      <div class="flex items-center justify-between">
        <div>
          <h2 id="activeTaxonTitle" class="text-sm font-bold text-white">Select a taxon to inspect 3D structure</h2>
          <div id="activeTaxonSubtitle" class="text-xs text-slate-400 mt-0.5">Drag to rotate &bull; Scroll to zoom</div>
        </div>
        <div id="activeTaxonBadge" class="text-xs font-mono px-3 py-1 rounded-full bg-slate-800 text-slate-300"></div>
      </div>
      <div class="flex-1 rounded-xl border border-slate-800 bg-slate-900/80 relative overflow-hidden flex items-center justify-center">
        <canvas id="viewer3dCanvas" class="w-full h-full block cursor-grab active:cursor-grabbing"></canvas>
      </div>
    </div>
  </div>
  ${openScriptTag}
    const TAXA = ${taxaJson};
    const STRUCTS = ${structsJson};
    const META = ${metaJson};
    let selectedTaxon = TAXA[0];

    const listEl = document.getElementById("taxaList");
    function renderList(query = "") {
      listEl.innerHTML = "";
      const q = query.toLowerCase();
      TAXA.filter(t => t.toLowerCase().includes(q) || (META[t] && JSON.stringify(META[t]).toLowerCase().includes(q))).forEach(t => {
        const item = document.createElement("div");
        const m = META[t] || {};
        const isSel = t === selectedTaxon;
        item.className = "p-2 rounded-lg text-xs cursor-pointer transition flex items-center justify-between " + (isSel ? "bg-sky-500/20 border border-sky-400/50 text-white font-bold" : "hover:bg-slate-800/80 text-slate-300 border border-transparent");
        item.innerHTML = "<span class='truncate'>" + t + "</span><span class='text-[10px] font-mono text-slate-400'>" + (m.structural_class || m.Family || "") + "</span>";
        item.onclick = () => selectTaxon(t);
        listEl.appendChild(item);
      });
    }

    document.getElementById("subcladeSearch").addEventListener("input", e => renderList(e.target.value));

    const canvas = document.getElementById("viewer3dCanvas");
    const ctx = canvas.getContext("2d");
    let rotX = 0.3, rotY = 0.5, zoom = 1.0;
    let isDragging = false, lastMouseX = 0, lastMouseY = 0;

    canvas.addEventListener("mousedown", e => { isDragging = true; lastMouseX = e.clientX; lastMouseY = e.clientY; });
    window.addEventListener("mouseup", () => { isDragging = false; });
    window.addEventListener("mousemove", e => {
      if (!isDragging) return;
      rotY += (e.clientX - lastMouseX) * 0.01;
      rotX += (e.clientY - lastMouseY) * 0.01;
      lastMouseX = e.clientX; lastMouseY = e.clientY;
      render3D();
    });
    canvas.addEventListener("wheel", e => {
      e.preventDefault();
      zoom *= (e.deltaY < 0 ? 1.08 : 0.92);
      render3D();
    });

    function selectTaxon(tid) {
      selectedTaxon = tid;
      renderList(document.getElementById("subcladeSearch").value);
      document.getElementById("activeTaxonTitle").innerText = tid;
      const m = META[tid] || {};
      document.getElementById("activeTaxonBadge").innerText = (m.structural_class || m.Family || "Taxon");
      render3D();
    }

    function render3D() {
      const rect = canvas.getBoundingClientRect();
      canvas.width = rect.width; canvas.height = rect.height;
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      const str = STRUCTS[selectedTaxon];
      if (!str || !str.atoms || str.atoms.length === 0) {
        ctx.fillStyle = "#64748b"; ctx.font = "12px sans-serif"; ctx.textAlign = "center";
        ctx.fillText("No 3D structure data available for " + selectedTaxon, canvas.width / 2, canvas.height / 2);
        return;
      }
      const cx = canvas.width / 2, cy = canvas.height / 2;
      const scale = (Math.min(canvas.width, canvas.height) / 130) * zoom;
      const cosY = Math.cos(rotY), sinY = Math.sin(rotY);
      const cosX = Math.cos(rotX), sinX = Math.sin(rotX);

      const proj = str.atoms.map(a => {
        let x1 = a.x * cosY + a.z * sinY;
        let z1 = -a.x * sinY + a.z * cosY;
        let y2 = a.y * cosX - z1 * sinX;
        let z2 = a.y * sinX + z1 * cosX;
        return { px: cx + x1 * scale, py: cy + y2 * scale, pz: z2, plddt: a.plddt };
      });

      for (let i = 0; i < proj.length - 1; i++) {
        const p1 = proj[i], p2 = proj[i+1];
        ctx.beginPath();
        ctx.moveTo(p1.px, p1.py); ctx.lineTo(p2.px, p2.py);
        const plddt = p1.plddt || 70;
        ctx.strokeStyle = plddt >= 90 ? "#2563eb" : plddt >= 70 ? "#38bdf8" : plddt >= 50 ? "#facc15" : "#f97316";
        ctx.lineWidth = 2.0; ctx.stroke();
      }
    }

    if (TAXA.length > 0) selectTaxon(TAXA[0]);
    window.addEventListener("resize", render3D);
  ${closeScriptTag}
</body>
</html>`;
    }

    function exportSubcladePackage() {
      try {
        let targetTaxa = null;
        let label = "subclade";

        if (settings.scopedClade && settings.scopedClade.taxa && settings.scopedClade.taxa.size > 0) {
          targetTaxa = Array.from(settings.scopedClade.taxa);
          label = (settings.scopedClade.name || "clade").toLowerCase().replace(/[^a-z0-9]/g, "_");
        } else if (typeof filterState !== 'undefined' && filterState.isActive) {
          const filteredSet = getFilteredTaxaSet();
          targetTaxa = Array.from(filteredSet);
          label = "filtered_" + (filterState.column || "subset").toLowerCase().replace(/[^a-z0-9]/g, "_");
        } else {
          targetTaxa = (typeof getAllLeaves === 'function' && activeTreeRoot) ? getAllLeaves(activeTreeRoot).map(l => l.name) : Object.keys(TAXA_METADATA);
          label = "cohort_" + currentScale;
        }

        if (!targetTaxa || targetTaxa.length === 0) {
          showToastNotification("⚠️ No active taxa available to export.");
          return;
        }

        const taxaSet = new Set(targetTaxa);
        const files = [];

        // A. Alignments (with dynamic gap stripping)
        const alnData = window.ALIGNMENTS_DATA ? window.ALIGNMENTS_DATA[currentScale] : null;
        if (alnData) {
          if (alnData.aa) {
            let faAa = "";
            const NL = String.fromCharCode(10);
            for (const tid of targetTaxa) {
              if (alnData.aa[tid]) faAa += ">" + tid + NL + alnData.aa[tid] + NL;
            }
            if (faAa) files.push({ name: "alignments/" + label + "_aa.fasta", content: faAa });
          }
          if (alnData["3di"]) {
            let fa3di = "";
            const NL = String.fromCharCode(10);
            for (const tid of targetTaxa) {
              if (alnData["3di"][tid]) fa3di += ">" + tid + NL + alnData["3di"][tid] + NL;
            }
            if (fa3di) files.push({ name: "alignments/" + label + "_3di.fasta", content: fa3di });
          }
        }

        // B. Newick Trees
        const esmNwk = (typeof NEWICK_ESM2 !== 'undefined' && NEWICK_ESM2) || 
                       (typeof NEWICK_ESM2_COSINE !== 'undefined' && NEWICK_ESM2_COSINE) || 
                       (typeof NEWICK_ESM2_EUCLIDEAN !== 'undefined' && NEWICK_ESM2_EUCLIDEAN) || 
                       (typeof NEWICK_ESM2_L1 !== 'undefined' && NEWICK_ESM2_L1) || null;
        if (targetTaxa.length >= Object.keys(TAXA_METADATA).length) {
          if (NEWICK_3DI) files.push({ name: "trees/" + label + "_3di.nwk", content: NEWICK_3DI });
          if (NEWICK_AA) files.push({ name: "trees/" + label + "_aa.nwk", content: NEWICK_AA });
          if (esmNwk) files.push({ name: "trees/" + label + "_esm2.nwk", content: esmNwk });
        } else {
          if (rawTrees && rawTrees["3di"]) {
            const p3 = pruneSubtree(rawTrees["3di"], taxaSet);
            if (p3) files.push({ name: "trees/" + label + "_3di.nwk", content: serializeNewick(p3) });
          }
          if (rawTrees && rawTrees["aa"]) {
            const pa = pruneSubtree(rawTrees["aa"], taxaSet);
            if (pa) files.push({ name: "trees/" + label + "_aa.nwk", content: serializeNewick(pa) });
          }
          const esmKey = "esm2_" + (settings.tanglegramEsmMetric || "cosine");
          if (rawTrees && rawTrees[esmKey]) {
            const pe = pruneSubtree(rawTrees[esmKey], taxaSet);
            if (pe) files.push({ name: "trees/" + label + "_" + esmKey + ".nwk", content: serializeNewick(pe) });
          }
        }

        // C. Metadata (TSV + JSON)
        const metaObj = {};
        const NL = String.fromCharCode(10);
        const TAB = String.fromCharCode(9);
        let tsvContent = "taxa_id";
        const firstTaxon = targetTaxa.find(t => TAXA_METADATA && TAXA_METADATA[t]);
        const headers = firstTaxon ? Object.keys(TAXA_METADATA[firstTaxon]) : [];
        if (headers.length > 0) tsvContent += TAB + headers.join(TAB) + NL;
        else tsvContent += NL;

        for (const tid of targetTaxa) {
          const m = (TAXA_METADATA && TAXA_METADATA[tid]) ? TAXA_METADATA[tid] : {};
          metaObj[tid] = m;
          const row = [tid];
          for (const h of headers) {
            const cell = (m[h] !== undefined && m[h] !== null) ? String(m[h]) : "";
            row.push(cell.split(TAB).join(" ").split(NL).join(" "));
          }
          tsvContent += row.join(TAB) + NL;
        }
        files.push({ name: "metadata/" + label + "_metadata.tsv", content: tsvContent });
        files.push({ name: "metadata/" + label + "_metadata.json", content: JSON.stringify(metaObj, null, 2) });

        // D. 3D Structures (C-alpha traces)
        const structDict = window.CA_STRUCTURES || window.CA_500_STRUCTURES || window.CA_100_STRUCTURES || {};
        const subsetStructs = {};
        if (structDict) {
          for (const tid of targetTaxa) {
            if (structDict[tid]) subsetStructs[tid] = structDict[tid];
          }
        }
        files.push({ name: "structures/" + label + "_ca_traces.json", content: JSON.stringify(subsetStructs) });

      // E. Summary & Silhouette Info
      const dateStr = new Date().toISOString();
      const summary = {
        label: label,
        cohort: currentScale,
        taxa_count: targetTaxa.length,
        taxa: targetTaxa,
        exported_at: dateStr,
        silhouette: settings.scopedClade ? settings.scopedClade.score : null
      };
      files.push({ name: "embeddings/" + label + "_summary.json", content: JSON.stringify(summary, null, 2) });

      // F. Reproducibility Scripts & Documentation
      const reproduceSh = `#!/usr/bin/env bash
# ==============================================================================
# REPRODUCIBILITY EXECUTION SCRIPT
# Subclade Analysis: ${label}
# Cohort Source: ${currentScale}
# Taxa Count: ${targetTaxa.length}
# Export Timestamp: ${dateStr}
# ==============================================================================

set -euo pipefail

echo "========================================================================"
echo "Reproducing Viral Structural Phylogenetics Analysis for: ${label}"
echo "Taxa in subset: ${targetTaxa.length}"
echo "========================================================================"

# 1. Verification of Prerequisite Tools
echo "==> [1/3] Checking environment and prerequisite binaries..."
command -v foldmason >/dev/null 2>&1 || { echo "ERROR: foldmason not found in PATH. Run ./install.sh to setup environment."; exit 1; }
command -v iqtree >/dev/null 2>&1 || { echo "ERROR: iqtree not found in PATH. Run ./install.sh to setup environment."; exit 1; }

OUTDIR="reproduced_results/${label}"
mkdir -p "$OUTDIR/phylogeny" "$OUTDIR/embeddings"

# 2. Maximum Likelihood Phylogenetic Inference (3Di Tertiary + AA Primary)
echo "==> [2/3] Inferring dual 3Di and AA phylogenies with IQ-TREE..."
python3 scripts/viral_phylogenetics.py tree \
  --alignment "alignments/${label}_3di.fasta" \
  --alignment-aa "alignments/${label}_aa.fasta" \
  --tree-type both \
  --method iqtree \
  --matrix alphafold \
  --rate-heterogeneity auto \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads AUTO \
  --output-dir "$OUTDIR/phylogeny" \
  --prefix "${label}_tree"

# 3. PLM Protein Language Model Embedding & UPGMA Hierarchical Clustering
echo "==> [3/3] Calculating PLM embeddings and hierarchical distance trees..."
python3 scripts/embed_and_cluster.py \
  --fasta "alignments/${label}_aa.fasta" \
  --model esm2_t33_650M_UR50D \
  --metric cosine \
  --output-dir "$OUTDIR/embeddings" \
  --prefix "${label}"

echo "========================================================================"
echo "✅ Reproduction complete!"
echo "Outputs stored in: $OUTDIR"
echo "========================================================================"
`;
      files.push({ name: "REPRODUCE.sh", content: reproduceSh });

      const fence = "```";
      const tick = "`";
      const reproduceMd = `# Reproducibility Report: ${label}

**Exported:** ${dateStr}  
**Parent Cohort:** ${currentScale} Taxa  
**Subset Size:** ${targetTaxa.length} Taxa  
${settings.scopedClade ? `**Silhouette Score:** S = ${settings.scopedClade.score || "N/A"}\n` : ""}

---

## 1. Quick Reproduction
To re-run the entire structural phylogenetics and PLM clustering pipeline on this exact subset:

${fence}bash
# Ensure Conda environment is active:
conda activate spt

# Execute the bundled reproduction script:
chmod +x REPRODUCE.sh
./REPRODUCE.sh
${fence}

---

## 2. Exact Commands & Parameters

### Step A: Maximum Likelihood Phylogeny Inference
Infers both 3Di structural and amino acid sequence trees:
${fence}bash
python3 scripts/viral_phylogenetics.py tree \\
  --alignment alignments/${label}_3di.fasta \\
  --alignment-aa alignments/${label}_aa.fasta \\
  --tree-type both \\
  --method iqtree \\
  --matrix alphafold \\
  --rate-heterogeneity auto \\
  --bootstrap 1000 \\
  --alrt 1000 \\
  --threads AUTO \\
  --output-dir results/reproduced/${label}/phylogeny \\
  --prefix ${label}_tree
${fence}

- **Substitution Model (3Di)**: Empirical AlphaFold matrix (${tick}Q.3Di.AF${tick}) with automated rate heterogeneity selection (${tick}+G4${tick}, ${tick}+I${tick}, ${tick}+R${tick}).
- **Substitution Model (AA)**: ModelFinder Plus automatic selection.
- **Resampling**: 1,000 Ultrafast Bootstrap (${tick}UFboot${tick}) and 1,000 SH-aLRT replicates.

### Step B: ESM-2 Protein Language Model Clustering
Generates mean-pooled sequence embeddings and constructs UPGMA distance trees:
${fence}bash
python3 scripts/embed_and_cluster.py \\
  --fasta alignments/${label}_aa.fasta \\
  --model esm2_t33_650M_UR50D \\
  --metric cosine \\
  --output-dir results/reproduced/${label}/embeddings \\
  --prefix ${label}
${fence}

---

## 3. Included Dataset Components
- ${tick}alignments/${label}_aa.fasta${tick}: Primary amino acid multiple sequence alignment.
- ${tick}alignments/${label}_3di.fasta${tick}: Tertiary 3Di structural multiple sequence alignment.
- ${tick}trees/${tick}: Pruned Newick trees (${tick}3di${tick}, ${tick}aa${tick}, ${tick}esm2${tick}).
- ${tick}structures/${label}_ca_traces.json${tick}: 3D backbone coordinates and pLDDT scores.
- ${tick}metadata/${label}_metadata.tsv${tick}: Tab-separated metadata table.
- ${tick}${label}_viewer.html${tick}: Standalone zero-dependency interactive 3D viewer.
- ${tick}REPRODUCE.sh${tick}: Executable reproduction script.
`;
      files.push({ name: "REPRODUCIBILITY.md", content: reproduceMd });

      // G. Standalone Subclade Viewer HTML
      const viewerHtml = generateStandaloneSubcladeViewer(label, targetTaxa, subsetStructs, metaObj);
      files.push({ name: label + "_viewer.html", content: viewerHtml });

      // H. Build ZIP and trigger download
      const zipBlob = createZipArchive(files);
      const zipName = label + "_" + targetTaxa.length + "_taxa_package.zip";
      const url = URL.createObjectURL(zipBlob);
      const a = document.createElement("a");
      a.href = url;
      a.download = zipName;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      setTimeout(() => URL.revokeObjectURL(url), 1000);

      showToastNotification("📦 Exported <strong>" + zipName + "</strong> (" + files.length + " files packaged)");
    } catch (err) {
      console.error("ZIP Export Error:", err);
      showToastNotification("⚠️ Error packaging ZIP: " + err.message);
    }
  }

    window.addEventListener("resize", () => {
      renderTree();
      renderMsa();
    });
