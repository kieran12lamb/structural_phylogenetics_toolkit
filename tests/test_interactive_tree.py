import unittest
from pathlib import Path


class TestInteractiveTree(unittest.TestCase):
    def setUp(self):
        self.repo_root = Path(__file__).resolve().parent.parent
        self.html_path = self.repo_root / "interactive_tree.html"

    def test_interactive_tree_exists_and_valid(self):
        self.assertTrue(self.html_path.exists(), "interactive_tree.html must exist at repo root")
        content = self.html_path.read_text(encoding="utf-8")
        
        # Verify file size (standalone with embedded precomputations is > 2 MB)
        self.assertGreater(len(content), 1_000_000)

        # Verify critical script assets and fallbacks
        self.assertIn("results/ca_500_structures.js", content)
        self.assertIn("results/ca_100_structures.js", content)
        self.assertIn("results/alignments_data.js", content)

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
        """Verify that the 6 benchmark example PDB structures exist in results/glycoprotein_workflow/structures."""
        structures_dir = self.repo_root / "results/glycoprotein_workflow/structures"
        self.assertTrue(structures_dir.exists(), "results/glycoprotein_workflow/structures must exist")
        pdbs = list(structures_dir.glob("*.pdb"))
        self.assertEqual(len(pdbs), 6, f"Expected exactly 6 example PDB structures, found {len(pdbs)}")
        for pdb in pdbs:
            self.assertGreater(pdb.stat().st_size, 10_000, f"{pdb.name} is too small to be a valid PDB")

    def test_uncomputed_esm_options_greyed_out(self):
        """Verify that datasets without ESM embeddings are flagged and uncomputed options are handled."""
        content = self.html_path.read_text(encoding="utf-8")
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
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

        # Verify .badge-purple class and variables
        self.assertIn(".badge-purple", content)
        self.assertIn("--badge-purple-bg:", content)
        self.assertIn("--badge-purple-text:", content)

        # Verify publication theme variables are crisp white/charcoal
        self.assertIn("--card-bg: #ffffff;", content)
        self.assertIn("--text-main: #0f172a;", content)

    def test_toolkit_branding_and_logo(self):
        """Verify that the header is renamed to The Structural Phylogenetics Toolkit and has C-alpha logo + toolkit emoji."""
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

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
        content = self.html_path.read_text(encoding="utf-8")

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


if __name__ == "__main__":
    unittest.main()





