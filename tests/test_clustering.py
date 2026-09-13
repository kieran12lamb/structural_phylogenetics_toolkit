import unittest
import numpy as np

def cosine_distance_matrix(vectors):
    """Compute pairwise cosine distances."""
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1e-12, norms)
    normed = vectors / norms
    similarity = np.dot(normed, normed.T)
    dist = 1.0 - similarity
    return np.clip(dist, 0.0, 2.0)

def euclidean_distance_matrix(vectors):
    """Compute pairwise Euclidean distances."""
    diff = vectors[:, np.newaxis, :] - vectors[np.newaxis, :, :]
    return np.sqrt(np.sum(diff ** 2, axis=-1))

def l1_distance_matrix(vectors):
    """Compute pairwise L1 / Manhattan distances."""
    diff = vectors[:, np.newaxis, :] - vectors[np.newaxis, :, :]
    return np.sum(np.abs(diff), axis=-1)


class TestClustering(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.taxa = ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"]
        vecs = np.random.randn(4, 16)
        self.vectors = vecs / np.linalg.norm(vecs, axis=1, keepdims=True)

    def test_distance_metrics(self):
        # Cosine
        d_cos = cosine_distance_matrix(self.vectors)
        self.assertEqual(d_cos.shape, (4, 4))
        np.testing.assert_allclose(np.diag(d_cos), 0.0, atol=1e-6)
        np.testing.assert_allclose(d_cos, d_cos.T, atol=1e-6)

        # Euclidean
        d_euc = euclidean_distance_matrix(self.vectors)
        self.assertEqual(d_euc.shape, (4, 4))
        np.testing.assert_allclose(np.diag(d_euc), 0.0, atol=1e-6)
        np.testing.assert_allclose(d_euc, d_euc.T, atol=1e-6)

        # L1
        d_l1 = l1_distance_matrix(self.vectors)
        self.assertEqual(d_l1.shape, (4, 4))
        np.testing.assert_allclose(np.diag(d_l1), 0.0, atol=1e-6)
        np.testing.assert_allclose(d_l1, d_l1.T, atol=1e-6)

    def test_scipy_clustering_if_available(self):
        """If scipy is installed, test embed_and_cluster hierarchical tree and silhouette."""
        try:
            from scripts.embed_and_cluster import build_hierarchical_tree, compute_silhouette_profile
        except ImportError:
            self.skipTest("SciPy or embed_and_cluster dependencies not available")

        newick, Z, condensed = build_hierarchical_tree(self.taxa, self.vectors, clustering="upgma", metric="cosine")
        
        self.assertIsInstance(newick, str)
        self.assertTrue(newick.endswith(";"))
        for t in self.taxa:
            self.assertIn(t, newick)

        profile_data = compute_silhouette_profile(Z, condensed, self.taxa, max_k=3)
        self.assertIn("profile", profile_data)
        self.assertIn("best_k", profile_data)
        self.assertGreater(len(profile_data["profile"]), 0)
        
        for item in profile_data["profile"]:
            self.assertIn("k", item)
            self.assertIn("score", item)
            self.assertTrue(-1.0 <= item["score"] <= 1.0)


if __name__ == "__main__":
    unittest.main()
