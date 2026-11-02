"""Tests for the QC, gating and partition-selection logic in viral_phylo.clustering.

These cover the pure-Python decision making only, so they run in CI without mmseqs
or foldseek installed. Two regressions are pinned deliberately:

  * selection must NOT maximise silhouette when the constraints cannot be met -
    doing so picks the most fragmented partition, because near-singleton clusters
    are trivially well separated;
  * the suggested alternative must prefer more usable clusters over a fraction of
    a percentage point of extra retention.
"""

import unittest

import numpy as np

from viral_phylo.clustering import (
    alignment_gap_fraction,
    clusters_to_labels,
    full_set_alignment_is_usable,
    partition_shape,
    qc_filter_sequences,
    select_partition,
    silhouette_for_partition,
    suggest_alternative,
)


def _candidate(threshold, coverage, sizes, silhouette):
    """Build a scored-candidate dict of the shape sweep_clusterings emits."""
    clusters = {}
    taxon = 0
    for index, size in enumerate(sizes):
        members = [f"t{taxon + i}" for i in range(size)]
        taxon += size
        clusters[f"rep{index}"] = members
    labels = np.array([i for i, size in enumerate(sizes) for _ in range(size)], dtype=int)
    shape = partition_shape(labels, min_cluster_size=10)
    shape.update({
        "tool": "foldseek", "threshold": threshold, "threshold_kind": "evalue",
        "coverage": coverage, "clusters": clusters,
        "silhouette": {"structural": silhouette},
    })
    return shape


def _sized(threshold, coverage, n_usable, retained_fraction, silhouette, total=1134):
    """Build a candidate hitting a target retention and usable-cluster count.

    Selection and suggestion depend only on those two quantities plus silhouette,
    so constructing them directly keeps the fixtures faithful to the real sweep
    without hand-balancing long lists of cluster sizes.
    """
    retained = int(round(total * retained_fraction))
    base, remainder = divmod(retained, n_usable)
    sizes = [base + (1 if i < remainder else 0) for i in range(n_usable)]
    assert min(sizes) >= 10, "usable clusters must be at least min_cluster_size"
    sizes += [1] * (total - retained)
    return _candidate(threshold, coverage, sizes, silhouette)


class TestQualityControl(unittest.TestCase):
    def test_short_fragments_are_excluded_with_a_reason(self):
        seqs = {"tiny": "A" * 9, "short": "A" * 35, "normal": "A" * 120, "long": "A" * 259}
        kept, excluded = qc_filter_sequences(seqs, min_length=50, min_length_frac=0.4)
        self.assertEqual(set(kept), {"normal", "long"})
        self.assertEqual(len(excluded), 2)
        for item in excluded:
            self.assertIn("reason", item)
            self.assertIn("floor", item["reason"])

    def test_relative_floor_scales_with_median(self):
        seqs = {f"s{i}": "A" * 1000 for i in range(10)}
        seqs["half"] = "A" * 300  # below 0.4 x median of 1000
        kept, excluded = qc_filter_sequences(seqs, min_length=50, min_length_frac=0.4)
        self.assertNotIn("half", kept)
        self.assertEqual(excluded[0]["taxon_id"], "half")

    def test_gaps_do_not_count_towards_length(self):
        seqs = {"gappy": "-" * 200 + "A" * 60, "solid": "A" * 60}
        kept, _ = qc_filter_sequences(seqs, min_length=50, min_length_frac=0.4)
        self.assertEqual(set(kept), {"gappy", "solid"})


class TestFullSetGate(unittest.TestCase):
    def test_rejects_a_sparse_over_expanded_alignment(self):
        """The 25 Sep Nipah alignment: 1458 columns for a median length of 120."""
        alignment = {f"t{i}": "A" * 120 + "-" * 1338 for i in range(20)}
        stats = alignment_gap_fraction(alignment)
        self.assertGreater(stats["gap_fraction"], 0.9)
        self.assertGreater(stats["expansion"], 10)
        usable, why = full_set_alignment_is_usable(stats)
        self.assertFalse(usable)
        self.assertIn("gap fraction", why)

    def test_accepts_a_compact_alignment(self):
        alignment = {f"t{i}": "ACDEFGHIKL" * 12 for i in range(20)}
        usable, why = full_set_alignment_is_usable(alignment_gap_fraction(alignment))
        self.assertTrue(usable, why)


class TestPartitionShape(unittest.TestCase):
    def test_counts_sizes_retention_and_singletons(self):
        labels = np.array([0] * 50 + [1] * 30 + [2] * 5 + [3] + [4], dtype=int)
        shape = partition_shape(labels, min_cluster_size=10)
        self.assertEqual(shape["n_clusters"], 5)
        self.assertEqual(shape["clusters_at_or_above_min"], 2)
        self.assertEqual(shape["singletons"], 2)
        self.assertEqual(shape["taxa_retained"], 80)
        self.assertAlmostEqual(shape["retained_fraction"], 80 / 87, places=3)

    def test_unassigned_taxa_are_tracked(self):
        labels = np.array([0] * 10 + [-1] * 3, dtype=int)
        self.assertEqual(partition_shape(labels, 10)["unassigned"], 3)

    def test_clusters_to_labels_marks_missing_taxa(self):
        labels = clusters_to_labels({"r": ["a", "b"]}, ["a", "b", "c"])
        self.assertEqual(labels[2], -1)


class TestSilhouetteScoping(unittest.TestCase):
    def test_small_clusters_are_excluded_from_the_score(self):
        """A partition of mostly singletons must not score well."""
        rng = np.random.default_rng(0)
        points = np.vstack([rng.normal(0, 0.1, (20, 2)), rng.normal(5, 0.1, (20, 2))])
        distance = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
        good = np.array([0] * 20 + [1] * 20)
        self.assertGreater(silhouette_for_partition(good, distance, 10), 0.8)
        # All singletons: nothing survives the min-size filter, so no score.
        singletons = np.arange(40)
        self.assertIsNone(silhouette_for_partition(singletons, distance, 10))


class TestSelection(unittest.TestCase):
    def test_maximises_silhouette_within_the_feasible_set(self):
        candidates = [
            _candidate(0.4, 0.7, [400, 300, 200, 100], 0.20),
            _candidate(0.5, 0.7, [350, 300, 250, 100], 0.45),
        ]
        result = select_partition(candidates, min_cluster_size=10, min_retained=0.70,
                                  primary_space="structural")
        self.assertFalse(result["relaxed_constraints"])
        self.assertEqual(result["selected"]["threshold"], 0.5)

    def test_falls_back_to_retention_not_silhouette_when_infeasible(self):
        """Regression: maximising silhouette here picks the most fragmented partition."""
        fragmented = _candidate(0.9, 0.8, [42, 20] + [1] * 800, 0.71)
        inclusive = _candidate(0.3, 0.5, [200, 150, 100] + [1] * 400, 0.21)
        result = select_partition([fragmented, inclusive], min_cluster_size=10,
                                  min_retained=0.70, primary_space="structural")
        self.assertTrue(result["relaxed_constraints"])
        self.assertEqual(result["selected"]["threshold"], 0.3,
                         "must prefer the partition retaining more taxa, not the higher silhouette")
        self.assertIn("recommendation", result)
        self.assertIn("too divergent", result["recommendation"])


class TestAlternativeSuggestion(unittest.TestCase):
    def setUp(self):
        # Synthetic, with deliberate headroom above the 10-point gain threshold so
        # that this exercises the retention tie band rather than the gain cutoff.
        self.selected = _sized(1.0, 0.5, n_usable=21, retained_fraction=0.75, silhouette=0.131)
        # c=0.5 retains 2 points more, but c=0.7 resolves 3 more usable clusters;
        # within the tie band the cluster count must decide.
        self.wide_c50 = _sized(10.0, 0.5, n_usable=15, retained_fraction=0.95, silhouette=0.115)
        self.wide_c70 = _sized(10.0, 0.7, n_usable=18, retained_fraction=0.93, silhouette=0.072)

    def test_prefers_more_usable_clusters_over_a_hair_more_retention(self):
        """Regression: within the retention tie band, more usable clusters must win."""
        suggestion = suggest_alternative([self.selected, self.wide_c50, self.wide_c70], self.selected)
        self.assertIsNotNone(suggestion)
        self.assertEqual(suggestion["coverage"], 0.7)
        self.assertEqual(suggestion["clusters_at_or_above_min"], 18)
        self.assertIn("--cluster-evalue 10.0", suggestion["how_to_use"])

    def test_no_suggestion_when_nothing_is_meaningfully_more_inclusive(self):
        marginal = _sized(0.5, 0.5, n_usable=21, retained_fraction=0.78, silhouette=0.20)
        self.assertIsNone(suggest_alternative([self.selected, marginal], self.selected))

    def test_suggestion_explains_the_trade_off(self):
        suggestion = suggest_alternative([self.selected, self.wide_c50, self.wide_c70], self.selected)
        self.assertIn("retains", suggestion["rationale"])
        self.assertIn("silhouette", suggestion["rationale"])


if __name__ == "__main__":
    unittest.main()


class TestRecursionSchedule(unittest.TestCase):
    """Regression: the sweep grid jumped 1.0 -> 0.01, shattering an 88-taxon cluster
    into 55 clusters (41 singletons) and keeping only 28% of its taxa in trees."""

    def test_evalue_schedule_starts_at_current_and_steps_gently(self):
        from viral_phylo.clustering import recursion_schedule
        self.assertEqual(recursion_schedule("evalue", 1.0), [1.0, 0.3, 0.1])

    def test_evalue_schedule_never_jumps_by_100x_in_one_step(self):
        from viral_phylo.clustering import recursion_schedule
        schedule = recursion_schedule("evalue", 1.0, max_attempts=6)
        for a, b in zip(schedule, schedule[1:]):
            self.assertLessEqual(a / b, 10, f"step {a} -> {b} is too large")

    def test_tmscore_schedule_is_strictly_stricter(self):
        from viral_phylo.clustering import recursion_schedule
        schedule = recursion_schedule("tmscore", 0.5)
        self.assertTrue(schedule and all(t > 0.5 for t in schedule))


class TestShatterGuard(unittest.TestCase):
    def test_fraction_kept_counts_only_tree_eligible_subclusters(self):
        from viral_phylo.clustering import fraction_kept
        subs = {"a": ["x"] * 40, "b": ["y"] * 3, "c": ["z"]}
        self.assertAlmostEqual(fraction_kept(subs, 44), 40 / 44)

    def test_shattered_split_is_below_the_default_bar(self):
        """The real e=0.01 split of the 88-taxon cluster kept 25 of 88 taxa."""
        from viral_phylo.clustering import RECURSE_MIN_KEPT, fraction_kept
        subs = {f"s{i}": ["t"] for i in range(41)}
        subs.update({"k1": ["t"] * 10, "k2": ["t"] * 8, "k3": ["t"] * 4, "k4": ["t"] * 3})
        subs.update({f"m{i}": ["t"] * 2 for i in range(10)})
        self.assertLess(fraction_kept(subs, 88), RECURSE_MIN_KEPT)


class TestSplitJudgement(unittest.TestCase):
    """Splits are judged on 3Di by default; AA only when asked (--recurse-on both)."""

    # 3Di fine, AA gappy: the typical MAFFT case.
    STATS = {"gap_fraction": 0.40, "expansion": 1.8, "aa": {"gap_fraction": 0.62, "expansion": 2.9}}

    def test_default_judges_3di_only(self):
        from viral_phylo.clustering import cluster_needs_split
        split, _ = cluster_needs_split(self.STATS, max_gap_fraction=0.55, max_expansion=2.5)
        self.assertFalse(split, "a gappy AA alignment must not trigger a split by default")

    def test_both_also_splits_on_aa(self):
        from viral_phylo.clustering import cluster_needs_split
        split, why = cluster_needs_split(self.STATS, recurse_on="both",
                                         max_gap_fraction=0.55, max_expansion=2.5)
        self.assertTrue(split)
        self.assertIn("AA:", why)
        self.assertNotIn("3Di:", why)

    def test_gappy_3di_splits_under_either_setting(self):
        from viral_phylo.clustering import cluster_needs_split
        stats = {"gap_fraction": 0.70, "expansion": 3.0, "aa": {"gap_fraction": 0.30, "expansion": 1.2}}
        for mode in ("3di", "both"):
            self.assertTrue(cluster_needs_split(stats, recurse_on=mode,
                                                max_gap_fraction=0.55, max_expansion=2.5)[0])


class TestSelectionWithoutSilhouette(unittest.TestCase):
    """Regression: with no distance space (a sweep on a fresh run, before embeddings
    exist) every candidate scored the same and the first in the grid was chosen."""

    def test_falls_back_to_retention_when_no_silhouette(self):
        a = _sized(0.001, 0.8, n_usable=15, retained_fraction=0.72, silhouette=None)
        b = _sized(1.0, 0.5, n_usable=21, retained_fraction=0.87, silhouette=None)
        for c in (a, b):
            c["silhouette"] = {}
        result = select_partition([a, b], min_cluster_size=10, min_retained=0.70)
        self.assertEqual(result["selected"]["threshold"], 1.0)


class TestTwoSpaceSelection(unittest.TestCase):
    """The sweep scores partitions in two spaces: ESM-2 when embeddings exist, and a
    structural (1 - TM-score) or identity space computed from the data."""

    def _pair(self):
        a = _sized(1.0, 0.5, n_usable=21, retained_fraction=0.87, silhouette=None)
        b = _sized(10.0, 0.5, n_usable=15, retained_fraction=0.93, silhouette=None)
        return a, b

    def test_selects_on_the_primary_space_and_records_disagreement(self):
        a, b = self._pair()
        a["silhouette"] = {"esm2": 0.30, "structural": 0.10}
        b["silhouette"] = {"esm2": 0.10, "structural": 0.40}
        result = select_partition([a, b], min_cluster_size=10, min_retained=0.70, primary_space="esm2")
        self.assertEqual(result["selected"]["threshold"], 1.0)
        self.assertIsNotNone(result["disagreement"])
        self.assertEqual(result["disagreement"]["space"], "structural")
        self.assertEqual(result["disagreement"]["alternative_choice"]["threshold"], 10.0)

    def test_no_disagreement_when_the_spaces_agree(self):
        a, b = self._pair()
        a["silhouette"] = {"esm2": 0.30, "structural": 0.40}
        b["silhouette"] = {"esm2": 0.10, "structural": 0.10}
        result = select_partition([a, b], min_cluster_size=10, min_retained=0.70, primary_space="esm2")
        self.assertIsNone(result["disagreement"])

    def test_falls_back_to_the_other_space_when_primary_is_missing(self):
        """A fresh run has no ESM-2 embeddings yet; the structural space must decide."""
        a, b = self._pair()
        a["silhouette"] = {"structural": 0.10}
        b["silhouette"] = {"structural": 0.40}
        result = select_partition([a, b], min_cluster_size=10, min_retained=0.70, primary_space="esm2")
        self.assertEqual(result["selected"]["threshold"], 10.0)


def _foldseek():
    from viral_phylo.binaries import tool_is_available
    return tool_is_available("foldseek")


@unittest.skipUnless(_foldseek(), "foldseek required")
class TestStructuralDistanceSpace(unittest.TestCase):
    def test_structural_space_is_a_valid_distance_matrix(self):
        import os
        import tempfile
        from pathlib import Path
        from viral_phylo.clustering import identity_distance_matrix
        fixtures = Path(__file__).resolve().parent / "fixtures" / "structures"
        taxa = sorted(p.stem for p in fixtures.glob("*.pdb"))
        with tempfile.TemporaryDirectory() as d:
            got, D = identity_distance_matrix(str(fixtures), taxa, os.path.join(d, "s"), tool="foldseek",
                                              score="alntmscore", evalue="10", alignment_type=2, threads=2)
        self.assertEqual(got, taxa)
        self.assertEqual(D.shape, (6, 6))
        self.assertTrue(np.allclose(D, D.T), "must be symmetric")
        self.assertTrue(np.all(np.diag(D) == 0))
        self.assertTrue(np.all((D >= 0) & (D <= 1)))
        self.assertTrue(np.any(D < 1), "some pairs must have a real TM-score hit")


class TestPLMClustering(unittest.TestCase):
    """Tests for PLM embedding-space clustering (k-means, agglomerative, ward, spectral)
    and automatic optimal-k selection by silhouette score."""

    def setUp(self):
        rng = np.random.default_rng(42)
        # 3 well-separated clusters in 16D embedding space (10 taxa each = 30 total)
        c1 = rng.normal(loc=[5.0, 0.0, 0.0] + [0.0] * 13, scale=0.05, size=(10, 16))
        c2 = rng.normal(loc=[0.0, 5.0, 0.0] + [0.0] * 13, scale=0.05, size=(10, 16))
        c3 = rng.normal(loc=[0.0, 0.0, 5.0] + [0.0] * 13, scale=0.05, size=(10, 16))
        self.embeddings = np.vstack([c1, c2, c3])
        self.taxa = [f"g1_{i}" for i in range(10)] + [f"g2_{i}" for i in range(10)] + [f"g3_{i}" for i in range(10)]

        norms = np.linalg.norm(self.embeddings, axis=1, keepdims=True)
        unit = self.embeddings / np.where(norms == 0, 1e-12, norms)
        self.distance = 1.0 - np.clip(unit @ unit.T, -1.0, 1.0)
        np.fill_diagonal(self.distance, 0.0)

    def test_plm_cluster_recovers_exact_groups_across_algorithms(self):
        from viral_phylo.clustering import plm_cluster
        expected_sets = {
            frozenset(f"g1_{i}" for i in range(10)),
            frozenset(f"g2_{i}" for i in range(10)),
            frozenset(f"g3_{i}" for i in range(10)),
        }
        for alg in ("kmeans", "agglomerative", "ward", "spectral"):
            clusters = plm_cluster(self.taxa, self.embeddings, k=3, algorithm=alg, metric="cosine")
            self.assertEqual(len(clusters), 3, f"{alg} should produce 3 clusters")
            for rep, members in clusters.items():
                self.assertIn(rep, members, "representative medoid must be a member of its cluster")
                self.assertEqual(members[0], rep, "representative medoid must be first in member list")
            got_sets = {frozenset(m) for m in clusters.values()}
            self.assertEqual(got_sets, expected_sets, f"{alg} did not recover the 3 ground-truth clusters")

    def test_sweep_plm_clusterings_picks_optimal_k_by_silhouette(self):
        from viral_phylo.clustering import sweep_plm_clusterings
        result = sweep_plm_clusterings(
            self.taxa,
            self.embeddings,
            self.distance,
            algorithm="kmeans",
            metric="cosine",
            k_min=2,
            k_max=7,
            min_cluster_size=5,
            min_retained=0.70,
            verbose=False,
        )
        selected = result["selected"]
        self.assertIsNotNone(selected)
        self.assertEqual(selected["k"], 3, "should auto-select k=3 as optimal silhouette score")
        self.assertEqual(selected["n_clusters"], 3)
        self.assertGreater(selected["silhouette_all"], 0.90)
        self.assertIn("optimal PLM (kmeans) silhouette score", result["reason"])

    def test_partition_structures_by_similarity_plm_mode_auto_and_fixed_k(self):
        import json
        import os
        import tempfile
        from unittest.mock import patch
        from viral_phylo.clustering import partition_structures_by_similarity

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_dir = os.path.join(tmpdir, "structures")
            os.makedirs(pdb_dir, exist_ok=True)
            fake_seqs = {}
            for t in self.taxa:
                with open(os.path.join(pdb_dir, f"{t}.pdb"), "w", encoding="utf-8") as fh:
                    fh.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n")
                fake_seqs[t] = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 residues

            npz_path = os.path.join(tmpdir, "embeddings.npz")
            np.savez_compressed(npz_path, taxa=np.array(self.taxa), embeddings=self.embeddings)

            def _fake_extract(structure_files, work_dir, **kwargs):
                fa = os.path.join(work_dir, "unaligned_aa.fa")
                with open(fa, "w", encoding="utf-8") as fh:
                    for t, s in fake_seqs.items():
                        fh.write(f">{t}\n{s}\n")
                return fa, fa

            def _fake_align(struct_dir, out_dir, **kwargs):
                os.makedirs(out_dir, exist_ok=True)
                members = [
                    os.path.splitext(f)[0]
                    for f in os.listdir(struct_dir)
                    if f.endswith(".pdb")
                ]
                for fname in ("foldmason.fasta_3di.fa", "foldmason.fasta_aa.fa"):
                    with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as fh:
                        for m in sorted(members):
                            fh.write(f">{m}\n{fake_seqs[m]}\n")

            with patch("viral_phylo.alignment.extract_unaligned_sequences", side_effect=_fake_extract), \
                 patch("viral_phylo.alignment.align_structures", side_effect=_fake_align):
                # 1. Default auto-k (optimal silhouette score -> k=3)
                out_auto = os.path.join(tmpdir, "aln_auto")
                summary_auto = partition_structures_by_similarity(
                    pdb_dir,
                    out_auto,
                    mode="plm",
                    embeddings_npz=npz_path,
                    plm_algorithm="kmeans",
                    plm_k=None,
                    plm_k_min=2,
                    plm_k_max=6,
                    min_cluster_size=5,
                )
                self.assertEqual(summary_auto["mode"], "plm")
                self.assertEqual(summary_auto["operating_point"]["threshold_kind"], "k")
                self.assertEqual(summary_auto["operating_point"]["threshold"], 3)
                self.assertEqual(summary_auto["operating_point"]["plm_k_selection"], "silhouette_auto")
                self.assertEqual(summary_auto["num_clusters"], 3)
                self.assertTrue(os.path.isdir(os.path.join(out_auto, "cluster_1_k3", "structures")))
                self.assertTrue(os.path.isfile(os.path.join(out_auto, "multi_alignment_summary.json")))

                # 2. Fixed k=4 provided beforehand
                out_fixed = os.path.join(tmpdir, "aln_fixed")
                summary_fixed = partition_structures_by_similarity(
                    pdb_dir,
                    out_fixed,
                    mode="plm",
                    embeddings_npz=npz_path,
                    plm_algorithm="kmeans",
                    plm_k=4,
                    min_cluster_size=5,
                )
                self.assertEqual(summary_fixed["operating_point"]["threshold"], 4)
                self.assertEqual(summary_fixed["operating_point"]["plm_k_selection"], "fixed")
                self.assertEqual(summary_fixed["num_clusters"], 4)
                self.assertTrue(os.path.isdir(os.path.join(out_fixed, "cluster_1_k4", "structures")))

    def test_partition_structures_by_similarity_both_mmseqs_and_plm(self):
        import os
        import tempfile
        from unittest.mock import patch
        from viral_phylo.clustering import partition_structures_by_similarity, resolve_cluster_modes

        self.assertEqual(resolve_cluster_modes("both"), ["sequence", "plm"])
        self.assertEqual(resolve_cluster_modes("sequence_plm"), ["sequence", "plm"])
        self.assertEqual(resolve_cluster_modes("structural_plm"), ["structural", "plm"])
        self.assertEqual(resolve_cluster_modes("all"), ["structural", "sequence", "plm"])

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_dir = os.path.join(tmpdir, "structures")
            os.makedirs(pdb_dir, exist_ok=True)
            fake_seqs = {}
            for t in self.taxa:
                with open(os.path.join(pdb_dir, f"{t}.pdb"), "w", encoding="utf-8") as fh:
                    fh.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n")
                fake_seqs[t] = "ACDEFGHIKLMNPQRSTVWY" * 5

            npz_path = os.path.join(tmpdir, "embeddings.npz")
            np.savez_compressed(npz_path, taxa=np.array(self.taxa), embeddings=self.embeddings)

            def _fake_extract(structure_files, work_dir, **kwargs):
                fa = os.path.join(work_dir, "unaligned_aa.fa")
                with open(fa, "w", encoding="utf-8") as fh:
                    for t, s in fake_seqs.items():
                        fh.write(f">{t}\n{s}\n")
                return fa, fa

            def _fake_align(struct_dir, out_dir, **kwargs):
                os.makedirs(out_dir, exist_ok=True)
                members = [os.path.splitext(f)[0] for f in os.listdir(struct_dir) if f.endswith(".pdb")]
                for fname in ("foldmason.fasta_3di.fa", "foldmason.fasta_aa.fa"):
                    with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as fh:
                        for m in sorted(members):
                            fh.write(f">{m}\n{fake_seqs[m]}\n")

            def _fake_mmseqs(fasta_path, workdir, min_seq_id, coverage, **kwargs):
                return {
                    "g1_0": [f"g1_{i}" for i in range(10)],
                    "g2_0": [f"g2_{i}" for i in range(10)],
                    "g3_0": [f"g3_{i}" for i in range(10)],
                }

            with patch("viral_phylo.alignment.extract_unaligned_sequences", side_effect=_fake_extract), \
                 patch("viral_phylo.alignment.align_structures", side_effect=_fake_align), \
                 patch("viral_phylo.clustering.mmseqs_cluster", side_effect=_fake_mmseqs):
                out_both = os.path.join(tmpdir, "aln_both")
                summary_both = partition_structures_by_similarity(
                    pdb_dir,
                    out_both,
                    mode="sequence_plm",
                    embeddings_npz=npz_path,
                    plm_algorithm="kmeans",
                    plm_k=None,
                    plm_k_min=2,
                    plm_k_max=5,
                    min_cluster_size=5,
                )
                self.assertEqual(summary_both["modes_run"], ["sequence", "plm"])
                self.assertEqual(summary_both["num_clusters"], 6)
                self.assertTrue(os.path.isdir(os.path.join(out_both, "cluster_mmseqs_1_id0.5", "structures")))
                self.assertTrue(os.path.isdir(os.path.join(out_both, "cluster_plm_1_k3", "structures")))
                self.assertTrue(os.path.isfile(os.path.join(out_both, "multi_alignment_summary_sequence.json")))
                self.assertTrue(os.path.isfile(os.path.join(out_both, "multi_alignment_summary_plm.json")))
                comp = summary_both["comparison"]
                self.assertIn("sequence", comp["modes"])
                self.assertIn("plm", comp["modes"])
                self.assertAlmostEqual(
                    comp["pairwise_agreement"]["sequence_vs_plm"]["adjusted_rand_index"], 1.0
                )
