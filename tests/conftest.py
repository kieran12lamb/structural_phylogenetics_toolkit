import pytest
import numpy as np
from pathlib import Path

@pytest.fixture
def sample_aa_records():
    return {
        "Taxon_A": "MKVLLLLALLSTLTLL",
        "Taxon_B": "MKVLLLLALLST-TLL",
        "Taxon_C": "MK-LLLLALLST-TLL",
        "Taxon_D": "MK-LLLLALLST-TLL"
    }

@pytest.fixture
def sample_3di_records():
    return {
        "Taxon_A": "DPCDVPDDPDPCDVPD",
        "Taxon_B": "DPCDVPDDPDP--VPD",
        "Taxon_C": "DP-DVPDDPDP--VPD",
        "Taxon_D": "DP-DVPDDPDP--VPD"
    }

@pytest.fixture
def sample_gapped_alignment():
    return {
        "Taxon_1": "M-K--L-V",
        "Taxon_2": "A-K--I-V",
        "Taxon_3": "V-R--L-I"
    }

@pytest.fixture
def sample_embeddings():
    # 4 synthetic 16-dimensional embedding vectors
    np.random.seed(42)
    taxa = ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"]
    vectors = np.random.randn(4, 16)
    # Normalize vectors
    vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    return taxa, vectors

@pytest.fixture
def sample_metadata():
    return {
        "Taxon_A": {"Family": "Rhabdoviridae", "pLDDT": 85.5, "length": 450},
        "Taxon_B": {"Family": "Rhabdoviridae", "pLDDT": 78.2, "length": 448},
        "Taxon_C": {"Family": "Flaviviridae", "pLDDT": 91.0, "length": 500},
        "Taxon_D": {"Family": "Flaviviridae", "pLDDT": 62.4, "length": 510}
    }

@pytest.fixture
def repo_root():
    return Path(__file__).resolve().parent.parent
