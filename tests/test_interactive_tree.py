import json
import unittest
from pathlib import Path

from viral_phylo.web.report import assemble_html


def _cohort(title, has_esm, n_taxa=60, **extra):
    """Build a synthetic cohort large enough that the chunk-window assertions below
    stay inside a single dataset entry, as they do for real cohorts."""
    cohort = {
        "title": title,
        "has_esm": has_esm,
        "taxa": {
            f"taxon_{i:03d}": {"taxon_id": f"taxon_{i:03d}", "family": "Demo"}
            for i in range(n_taxa)
        },
    }
    cohort.update(extra)
    return cohort


class TestInteractiveTree(unittest.TestCase):
    """Assert against the page assembled from the tracked templates in
    viral_phylo/web/template/, NOT against a built interactive_tree.html.

    The built artifact is generated output and is deliberately untracked, so
    reading it here would make the suite depend on whether a pipeline has been
    run. assemble_html() is the same function the builder uses, so these tests
    still cover the real template, CSS, viewer JS and dataset wiring.
    """

    DATASETS_FIXTURE = {
        "1193": _cohort("\U0001f9ec 1,193 Nipah ESMFold Structures", True,
                        esm2_umap={"x": [0.0], "y": [0.0], "taxa": ["taxon_000"]}),
        "500": _cohort("\U0001f9ec 500 Viral Glycoproteins", True),
        "6": _cohort("\U0001f9ec 6 Benchmark Glycoproteins", True,
                     congruence={"3di_vs_esm2_cosine": {"rf_distance": 0}}),
        "100": _cohort("\U0001f9ec 100 RNA-dependent RNA Polymerases", False),
    }

    def setUp(self):
        self.repo_root = Path(__file__).resolve().parent.parent
        self.content = assemble_html(
            json.dumps(self.DATASETS_FIXTURE),
            "v1.0.0 (abc1234)",
            "",
        )

    def test_interactive_tree_exists_and_valid(self):
        content = self.content

        # Assembled page carries the full inlined CSS and viewer JS
        self.assertGreater(len(content), 400_000)

        # The builder substitutes the data script tags into this placeholder;
        # assert the contract exists rather than the substituted result.
        self.assertIn("<!-- DATA_SCRIPTS_PLACEHOLDER -->", content)

        # Verify cohorts embedded
        self.assertIn('"1193":', content)
        self.assertIn('"500":', content)
        self.assertIn('"100":', content)
        self.assertIn('"6":', content)

        # Verify ZIP exporter and reproducibility components
        self.assertIn("function exportSubcladePackage()", content)
        self.assertIn("function createZipArchive(files)", content)
        self.assertIn("REPRODUCE.sh", content)
        self.assertIn("REPRODUCIBILITY.md", content)

    def test_example_pdbs_exist(self):
        """Verify that the 6 benchmark example PDB structures exist in tests/fixtures/structures."""
        structures_dir = self.repo_root / "tests/fixtures/structures"
        self.assertTrue(structures_dir.exists(), "tests/fixtures/structures must exist")
        pdbs = list(structures_dir.glob("*.pdb"))
        self.assertEqual(len(pdbs), 6, f"Expected exactly 6 example PDB structures, found {len(pdbs)}")
        for pdb in pdbs:
            self.assertGreater(pdb.stat().st_size, 10_000, f"{pdb.name} is too small to be a valid PDB")

    def test_uncomputed_esm_options_greyed_out(self):
        """Verify that datasets without ESM embeddings are flagged and uncomputed options are handled."""
        content = self.content
        self.assertIn("function updateModalityOptions()", content)
        self.assertIn("🤖 ESM-2 PLM Tree (Not Run in Pipeline)", content)
        self.assertIn("(Not Run)", content)

        # Check cohort 100 has has_esm: false
        idx_100 = content.find('"100":')
        self.assertNotEqual(idx_100, -1)
        chunk_100 = content[idx_100:idx_100 + 1500]
        self.assertIn('"has_esm": false', chunk_100)
        self.assertNotIn("3di_vs_esm2_cosine", chunk_100)

        # Check cohort 1193 has has_esm: true
        self.assertIn('"1193": {"title": "\\ud83e\\uddec 1,193 Nipah ESMFold Structures", "has_esm": true', content)

        # Check cohort 6 has has_esm: true and full ESM congruence profiles
        idx_6 = content.find('"6":')
        self.assertNotEqual(idx_6, -1)
        chunk_6 = content[idx_6:]
        self.assertIn('"has_esm": true', chunk_6)
        self.assertIn("3di_vs_esm2_cosine", chunk_6)



    def test_themes_and_scoped_banner_layout(self):
        """Verify that multi-theme switching and compact scoped clade positioning are present."""
        content = self.content

        # Verify expanded theme collection exists in CSS
        self.assertIn('[data-theme="dark"]', content)
        self.assertIn('[data-theme="obsidian"]', content)
        self.assertIn('[data-theme="forest"]', content)
        self.assertIn('[data-theme="dracula"]', content)
        self.assertIn('[data-theme="cyberpunk"]', content)
        self.assertIn('[data-theme="steel"]', content)
        self.assertIn('[data-theme="espresso"]', content)
        self.assertIn('[data-theme="light"]', content)
        self.assertIn('[data-theme="solarized"]', content)
        self.assertIn('[data-theme="nordic"]', content)
        self.assertIn('[data-theme="parchment"]', content)
        self.assertIn('[data-theme="mint"]', content)
        self.assertIn('[data-theme="high_contrast"]', content)

        # Verify high-contrast clade badge CSS classes
        self.assertIn(".badge-sky", content)
        self.assertIn(".badge-emerald", content)
        self.assertIn(".badge-rose", content)

        # Verify high-contrast branch strokes for publication white
        self.assertIn("--branch-stroke: #1e293b", content)
        self.assertIn("--tip-label: #0f172a", content)

        # Verify theme switching JS functions
        self.assertIn("THEMES_META", content)
        self.assertIn("function isDarkTheme(", content)
        self.assertIn("function toggleThemeMenu(", content)

        # Verify scoped clade banner and active filter banner are docked at top-3 left-4
        self.assertIn('id="scopedCladeBanner" class="absolute top-3 left-4', content)
        self.assertIn('id="activeFilterBanner" class="absolute top-3 left-4', content)


    def test_palette_creator_studio(self):
        """Verify that the custom Palette Creator Studio, hex string parsing, and colour wheel picker are embedded."""
        content = self.content

        # Verify Palette Modal & UI buttons
        self.assertIn('id="paletteModal"', content)
        self.assertIn('id="headerPaletteBtn"', content)
        self.assertIn('id="paletteMiniStrip"', content)
        self.assertIn('id="paletteHexInput"', content)

        # Verify user requested botanical earth palette preset (#B9554E, #627B08, #267567, #294719, #72A183)
        self.assertIn("#B9554E", content)
        self.assertIn("#627B08", content)
        self.assertIn("#267567", content)
        self.assertIn("#294719", content)
        self.assertIn("#72A183", content)

        # Verify color wheel input and JS handler functions
        self.assertIn('input type="color"', content)
        self.assertIn("function openPaletteModal(", content)
        self.assertIn("function closePaletteModal(", content)
        self.assertIn("function onColorWheelChange(", content)
        self.assertIn("function parseHexList(", content)
        self.assertIn("function interpolatePalette(", content)
        self.assertIn("function applyPaletteStudio(", content)
        self.assertIn("function resetToDefaultPalette(", content)
        self.assertIn("function getCategoryColorFromCustomPalette(", content)


    def test_alignment_viewer_synchronization(self):
        """Verify that the alignment viewer updates synchronously with tree filtering, rooting, scoping, and themes."""
        content = self.content

        # Verify default sync with tree selection is enabled
        self.assertIn("syncWithTree: true", content)

        # Verify theme-aware MSA rendering and gap handling
        self.assertIn("data[idx] = 241;    // Gap light: soft paper slate #f1f5f9", content)
        self.assertIn('msaState.colorScheme === "custom"', content)

        # Verify renderMsa is called on rooting, scoping, filtering, and theme change
        self.assertIn("function setMsaMode(", content)
        self.assertIn("function setMsaColorScheme(", content)
        self.assertIn("function renderMsa(", content)

    def test_custom_theme_studio(self):
        """Verify that the interactive Custom Theme Studio modal, presets, and JS engine are properly integrated."""
        content = self.content

        # Verify [data-theme="custom"] CSS definition
        self.assertIn('[data-theme="custom"]', content)

        # Verify #themeModal and trigger button
        self.assertIn('id="themeModal"', content)
        self.assertIn('openThemeModal()', content)
        self.assertIn('applyCustomThemeStudio()', content)

        # Verify Custom Theme Studio engine functions
        self.assertIn("function openThemeModal(", content)
        self.assertIn("function closeThemeModal(", content)
        self.assertIn("function loadThemePreset(", content)
        self.assertIn("function onThemeColorWheelChange(", content)
        self.assertIn("function onThemeHexChange(", content)
        self.assertIn("function applyCustomThemeStudio(", content)
        self.assertIn("function resetCustomThemeToDefault(", content)
        self.assertIn("function copyThemeJson(", content)
        self.assertIn("function importThemeJson(", content)
        self.assertIn("function initCustomThemeFromStorage(", content)

        # Verify curated presets exist
        self.assertIn("Linear Crisp Light", content)
        self.assertIn("GitHub Minimal Light", content)
        self.assertIn("Editorial Parchment", content)
        self.assertIn("Cafe Latte Cream", content)
        self.assertIn("Obsidian Velvet", content)

    def test_light_mode_aesthetics_and_contrast(self):
        """Verify that light themes have crisp publication surfaces and high-contrast badges."""
        content = self.content

        # Verify .badge-purple class and variables
        self.assertIn(".badge-purple", content)
        self.assertIn("--badge-purple-bg:", content)
        self.assertIn("--badge-purple-text:", content)

        # Verify publication theme variables are crisp white/charcoal
        self.assertIn("--card-bg: #ffffff;", content)
        self.assertIn("--text-main: #0f172a;", content)

    def test_toolkit_branding_and_logo(self):
        """Verify that the header is renamed to The Structural Phylogenetics Toolkit and has C-alpha logo + toolkit emoji."""
        content = self.content

        # Verify page title and header h1
        self.assertIn("<title>The Structural Phylogenetics Toolkit</title>", content)
        self.assertIn("The Structural Phylogenetics Toolkit</h1>", content)

        # Verify C-alpha protein structure logo elements
        self.assertIn('id="headerLogoCanvas"', content)
        self.assertIn('class="MiniHeaderLogoViewer"', content.replace("class MiniHeaderLogoViewer", 'class="MiniHeaderLogoViewer"'))  # or assert class definition
        self.assertIn("class MiniHeaderLogoViewer", content)
        self.assertIn("function initHeaderLogo(", content)
        self.assertIn("initHeaderLogo()", content)

        # Verify C-alpha SVG backbone gradient fallback
        self.assertIn('id="logoBackboneGrad"', content)

        # Verify stylized SPT monogram, toolkit handle, and internal tree & structure tool bays
        self.assertIn(">SPT</span>", content)
        self.assertIn('title="Radial Phylogenetic Tree"', content)
        self.assertIn('title="Live 3D C-alpha Protein Structure"', content)
        self.assertIn("Toolbox Carry Handle", content)
        self.assertNotIn("🧰", content)

        # Verify removal of Viro3D bubble and subtitle, and presence of git-linked version
        self.assertNotIn("Viro3D &amp; Local", content)
        self.assertNotIn("Viro3D & Local", content)
        self.assertNotIn("Dual IQ-TREE 3Di &amp; AA Inference", content)
        self.assertRegex(content, r'v\d+\.\d+\.\d+\s*\([0-9a-fA-F]+\)')

    def test_tip_label_metadata_selection(self):
        """Verify dynamic tip label selector, metadata resolution, and data-taxon attribute decoupling."""
        content = self.content

        # Verify selector DOM elements exist
        self.assertIn('id="tipLabelColumnSelect"', content)
        self.assertIn('id="tipLabelActiveBadge"', content)
        self.assertIn('onchange="setTipLabelColumn(this.value)"', content)

        # Verify JavaScript state and functions exist
        self.assertIn('tipLabelColumn: "taxon_id"', content)
        self.assertIn('function getLeafLabelText(leafName)', content)
        self.assertIn('function setTipLabelColumn(colKey)', content)

        # Verify data-taxon attribute decoupling across layouts
        self.assertIn('txt.setAttribute("data-taxon", leaf.name);', content)
        self.assertIn('txt1.setAttribute("data-taxon", l1.name);', content)
        self.assertIn('txt2.setAttribute("data-taxon", l2.name);', content)
        self.assertIn('rowEl.setAttribute("data-taxon", tName);', content)

        # Verify selection and search robustness
        self.assertIn('const taxon = el.getAttribute("data-taxon") || el.textContent;', content)
        self.assertIn('const taxonName = el.getAttribute("data-taxon") || el.textContent;', content)

    def test_export_newick_and_zip_package(self):
        """Verify Newick export and ZIP archive package export buttons, functions, and error resilience."""
        content = self.content

        # Verify buttons in header
        self.assertIn('onclick="exportNewick()"', content)
        self.assertIn('onclick="exportSubcladePackage()"', content)
        self.assertIn('Export Newick', content)
        self.assertIn('Export ZIP', content)

        # Verify Newick export functions and multi-strategy copy
        self.assertIn('function exportNewick()', content)
        self.assertIn('function copyNewick()', content)
        self.assertIn('function getActiveNewickString()', content)
        self.assertIn('function copyTextWithFallback(', content)

        # Verify ZIP export implementation and no undefined variable references
        self.assertIn('function exportSubcladePackage()', content)
        self.assertIn('function createZipArchive(files)', content)
        self.assertNotIn('currentLeaves', content)
        self.assertNotIn('currentMeta', content)

        # Verify toast notifications wired into export flows
        self.assertIn('function showToastNotification(', content)
        self.assertIn('showToastNotification(msg);', content)

    def test_radial_and_unrooted_readability_controls(self):
        """Verify advanced iTOL & FigTree inspired readability controls for Radial and Unrooted trees."""
        content = self.content

        # Verify layout-specific control cards exist in HTML
        self.assertIn('id="radialControlsCard"', content)
        self.assertIn('id="unrootedControlsCard"', content)
        self.assertIn('id="rectControlsCard"', content)

        # Verify interactive sliders and toggles exist
        self.assertIn('id="treeRotationSlider"', content)
        self.assertIn('id="radialArcSlider"', content)
        self.assertIn('id="radialRadiusScaleSlider"', content)
        self.assertIn('id="unrootedRotationSlider"', content)
        self.assertIn('id="unrootedScaleSlider"', content)
        self.assertIn('id="branchWidthSlider"', content)
        self.assertIn('id="radialLabelOrientationSelect"', content)
        self.assertIn('id="unrootedLabelOrientationSelect"', content)
        self.assertIn('id="toggleRadialAlign"', content)
        self.assertIn('id="toggleConcentricRings"', content)
        self.assertIn('id="toggleCladeSectors"', content)

        # Verify JavaScript state properties exist
        self.assertIn('treeRotation: 0', content)
        self.assertIn('radialArc: 360', content)
        self.assertIn('radialRadiusScale: 1.0', content)
        self.assertIn('unrootedScale: 1.0', content)
        self.assertIn('labelOrientation: "radial"', content)
        self.assertIn('concentricRings: true', content)
        self.assertIn('cladeSectors: true', content)
        self.assertIn('branchWidth: 1.4', content)

        # Verify setter functions exist
        self.assertIn('function setTreeRotation(deg)', content)
        self.assertIn('function setRadialArc(deg)', content)
        self.assertIn('function setRadialRadiusScale(val)', content)
        self.assertIn('function setUnrootedScale(val)', content)
        self.assertIn('function setLabelOrientation(mode)', content)
        self.assertIn('function setBranchWidth(w)', content)
        self.assertIn('function updateLayoutSpecificControls()', content)

        # Verify dynamic rendering features in Radial & Unrooted layout engines
        self.assertIn('"concentric-rings"', content)
        self.assertIn('"clade-sectors"', content)
        self.assertIn('settings.labelOrientation === "horizontal"', content)

    def test_alphafold_database_querying_studio(self):
        """Verify that the Pipeline & Structure Studio supports querying AlphaFold DB (AFDB)."""
        content = self.content

        # Verify studio title and source selector buttons
        self.assertIn("Pipeline &amp; Structure Studio", content)
        self.assertIn('id="btnPipeSourceAlphaFold"', content)
        self.assertIn("setPipelineSource('alphafold')", content)

        # Verify AlphaFold DB input section, format selector, and PAE toggle
        self.assertIn('id="pipeAlphaFoldSection"', content)
        self.assertIn('id="pipeAfdbInput"', content)
        self.assertIn('id="pipeAfdbFormatSelect"', content)
        self.assertIn('id="pipeAfdbPaeToggle"', content)

        # Verify check button and JS client query functions
        self.assertIn('id="btnAfdbCheck"', content)
        self.assertIn('id="afdbCheckStatus"', content)
        self.assertIn("function queryAlphaFoldApi()", content)
        self.assertIn("function setAlphaFoldPreset(", content)

        # Verify AFDB API endpoint integration
        self.assertIn("https://alphafold.ebi.ac.uk/api/prediction/", content)
        self.assertIn("--source alphafold", content)

    def test_bootstrap_support_branch_coloring(self):
        """Verify that branch bootstrap support coloring and IQ-TREE support parsing are present."""
        content = self.content

        # Verify UI controls in Tree Display panel
        self.assertIn('id="toggleSupportColor"', content)
        self.assertIn('id="selectSupportMetric"', content)
        self.assertIn('id="selectSupportPalette"', content)
        self.assertIn('id="supportPalettePreview"', content)
        self.assertIn('id="branchSupportLegend"', content)
        self.assertIn('id="branchSupportLegendBody"', content)

        # Verify default settings
        self.assertIn('colorBranchesBySupport: false', content)
        self.assertIn('supportMetric: "ufboot"', content)
        self.assertIn('supportPalette: "traffic"', content)

        # Verify support functions
        self.assertIn('function getBranchStrokeColor(node)', content)
        self.assertIn('function toggleSupportColoring(val)', content)
        self.assertIn('function setSupportMetric(m)', content)
        self.assertIn('function setSupportPalette(p)', content)

        # Verify composite Newick support parsing (UFboot and SH-aLRT)
        self.assertIn('node.supportRaw', content)
        self.assertIn('node.ufboot', content)
        self.assertIn('node.alrt', content)

        # Verify tooltip breakdown for bootstrap values
        self.assertIn('UFboot Support:', content)
        self.assertIn('SH-aLRT Support:', content)

        # Verify vertical branch gradient between horizontal branches
        self.assertIn('function applyVerticalBranchGradient(vLine, node, minY, maxY, defs)', content)
        self.assertIn('gradientUnits', content)
        self.assertIn('userSpaceOnUse', content)
        self.assertIn('applyVerticalBranchGradient(vLine', content)


    def test_background_grid_toggle_and_unrooted_enhancements(self):
        """Verify background grid toggle (UI + HUD) and high-density unrooted tree features."""
        content = self.content

        # 1. Background grid toggle CSS and elements
        self.assertIn("#treeSvg.no-grid", content)
        self.assertIn('id="toggleShowGrid"', content)
        self.assertIn('id="btnToggleGrid"', content)
        self.assertIn('id="btnToggleGridLabel"', content)
        self.assertIn("function toggleBackgroundGrid(", content)
        self.assertIn("showGrid: true", content)

        # 2. Unrooted tree controls in display tab
        self.assertIn('id="unrootedControlsCard"', content)
        self.assertIn('id="selectUnrootedLenMode"', content)
        self.assertIn('id="toggleUnrootedDaylight"', content)
        self.assertIn('id="selectUnrootedLabelFilter"', content)
        self.assertIn('id="toggleStaggerLabels"', content)

        # 3. Unrooted tree settings and helper functions
        self.assertIn('unrootedLengthMode: "sqrt"', content)
        self.assertIn('unrootedDaylight: true', content)
        self.assertIn('unrootedLabelFilter: "smart"', content)
        self.assertIn('staggerLabels: true', content)
        self.assertIn("function setUnrootedLengthMode(", content)
        self.assertIn("function toggleUnrootedDaylight(", content)
        self.assertIn("function setUnrootedLabelFilter(", content)
        self.assertIn("function toggleStaggerLabels(", content)

        # 4. Equal-daylight fan spreading and smart label decluttering logic
        self.assertIn("layoutEqualAngle", content)
        self.assertIn("Equal-Daylight vs Equal-Angle weight calculation", content)
        self.assertIn("Math.pow(lc, 0.72) + 0.35", content)
        self.assertIn("settings.unrootedLabelFilter", content)
        self.assertIn("settings.staggerLabels", content)


    def test_alignment_coverage_threshold_and_partitions_ui(self):
        """Verify 70% coverage threshold UI controls, JS helpers, and multi-alignment partitions."""
        content = self.content

        # 1. Alignment Drawer UI controls for coverage threshold & partitions
        self.assertIn('id="selectAlignmentPartition"', content)
        self.assertIn('id="msaCoverageThreshold"', content)
        self.assertIn('id="btnFilterCoverage"', content)

        # 2. Run/Fetch Pipeline Generator controls
        self.assertIn('id="pipeCoverageInput"', content)
        self.assertIn('id="pipeMultiAlnToggle"', content)

        # 3. JavaScript state & functions
        self.assertIn("coverageThreshold: 0.70", content)
        self.assertIn("filterCoverage: false", content)
        self.assertIn("activePartition: \"all\"", content)
        self.assertIn("function getSequenceCoverage(", content)
        self.assertIn("function setMsaCoverageThreshold(", content)
        self.assertIn("function toggleFilterCoverage(", content)
        self.assertIn("function setMsaPartition(", content)

        # 4. Command line generator includes coverage options
        self.assertIn("--min-coverage", content)
        self.assertIn("--multi-alignment", content)

        # 5. Tree Synchronization with Coverage Filter & Partitions
        self.assertIn("function isTreeFilteringActive()", content)
        self.assertIn("function syncDatasetAlignmentCoverage(", content)
        self.assertIn("isCovFiltering", content)
        self.assertIn("isPartFiltering", content)
        self.assertIn("Core Align (≥70% Cov)", content)
        self.assertIn("Alignment Coverage (%)", content)


    def test_esm2_umap_scatter_view(self):
        """Verify ESM-2 2D UMAP scatter projection view, controls, and radar synchronization."""
        content = self.content

        # 1. Embedded UMAP projection datasets
        self.assertIn('"esm2_umap":', content)
        self.assertIn("UMAP-1 (ESM-2 Latent Dimension 1)", content)
        self.assertIn("UMAP-2 (ESM-2 Latent Dimension 2)", content)

        # 2. UI view switcher buttons
        self.assertIn('id="btnUmap"', content)
        self.assertIn('id="esmViewModeSection"', content)
        self.assertIn('id="btnEsmViewTree"', content)
        self.assertIn('id="btnEsmViewUmap"', content)
        self.assertIn('id="umapControlsCard"', content)

        # 3. JavaScript handlers & render pipeline
        self.assertIn("function setEsmViewMode(", content)
        self.assertIn("function renderUmapScatter(", content)
        self.assertIn('settings.layout === "umap"', content)
        self.assertIn("UMAP RADAR", content)
        self.assertIn("umapPointCoords", content)


if __name__ == "__main__":
    unittest.main()








class TestAlphabetLengths(unittest.TestCase):
    """Regression: the viewer used one length (from the 3Di alignment) for both
    alphabets. MAFFT aligns AA and 3Di separately, so their lengths differ and AA
    mode was silently truncated or over-read."""

    def test_bundle_records_a_length_per_alphabet(self):
        from viral_phylo.web.alignments import annotate_alphabet_lengths
        data = {"run": {"length": 5, "aa": {"t1": "ACDEFGH"}, "3di": {"t1": "ACDEF"},
                        "partitions": {"p1": {"aa": {"t1": "AC-D"}, "3di": {"t1": "ACD"}}}}}
        annotate_alphabet_lengths(data)
        self.assertEqual(data["run"]["lengths"], {"aa": 7, "3di": 5})
        self.assertEqual(data["run"]["length"], 5, "legacy field must be left unchanged")
        self.assertEqual(data["run"]["partitions"]["p1"]["lengths"], {"aa": 4, "3di": 3})

    def test_viewer_reads_the_length_for_the_active_alphabet(self):
        content = assemble_html(json.dumps(TestInteractiveTree.DATASETS_FIXTURE), "v1.0.0 (abc1234)", "")
        self.assertIn("function alignLength(align)", content)
        self.assertIn("align.lengths[mode]", content)
        self.assertNotIn("origLen = align.length", content)
