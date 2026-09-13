#!/usr/bin/env python3
"""Protein Language Model (ESM-2 / ESM-C) Embedding & Hierarchical Tree Construction.

This module embeds protein amino acid sequences using high-capacity Protein
Language Models (PLMs) such as ESM-2 650M (facebook/esm2_t33_650M_UR50D) or
ESM-C 600M (biohub/ESMC-600M), computes pairwise semantic distance matrices,
and constructs hierarchical clustering trees (UPGMA / Neighbor-Joining) in Newick format.

Compatible with IQ-TREE outputs and the Viral Structural Phylogenetics Interactive Suite.
"""

import argparse
import io
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np


def find_torch_python() -> Optional[str]:
    """Find a Python interpreter to use for ML/embeddings.

    Priority:
    1. The install script environment 'spt' (e.g. envs/spt/bin/python)
    2. The user's active Conda environment ($CONDA_PREFIX/bin/python)
    3. The currently executing Python interpreter (sys.executable)
    4. Fallback candidate environments
    """
    def is_valid_python(path: str) -> bool:
        return bool(path and os.path.isfile(path) and os.access(path, os.X_OK))

    # 1. Primary default: The install script environment 'spt'
    spt_candidates = []
    if os.environ.get("CONDA_DEFAULT_ENV") == "spt" and os.environ.get("CONDA_PREFIX"):
        spt_candidates.append(os.path.join(os.environ["CONDA_PREFIX"], "bin", "python"))

    conda_base_hints = [
        os.environ.get("CONDA_PREFIX", ""),
        os.path.expanduser("~/miniconda3"),
        os.path.expanduser("~/anaconda3"),
        os.path.expanduser("~/miniforge3"),
        os.path.expanduser("~/micromamba"),
        "/opt/conda",
        "/opt/homebrew/Caskroom/miniconda/base"
    ]
    for hint in conda_base_hints:
        if not hint:
            continue
        if "/envs/" in hint:
            base_dir = hint.split("/envs/")[0]
            spt_candidates.append(os.path.join(base_dir, "envs", "spt", "bin", "python"))
        else:
            spt_candidates.append(os.path.join(hint, "envs", "spt", "bin", "python"))

    for cand in spt_candidates:
        if is_valid_python(cand):
            return cand

    # 2. Secondary default: Assume the user's active Conda environment
    if os.environ.get("CONDA_PREFIX"):
        active_py = os.path.join(os.environ["CONDA_PREFIX"], "bin", "python")
        if is_valid_python(active_py):
            return active_py

    # 3. Third priority: Currently executing Python interpreter
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
        import scipy  # noqa: F401
        return sys.executable
    except ImportError:
        pass

    # 4. Fallback candidate environments
    fallbacks = [
        os.path.expanduser("~/miniconda3/envs/nipah/bin/python"),
        os.path.expanduser("~/miniconda3/bin/python"),
        "/usr/local/bin/python3",
        sys.executable,
    ]
    for cand in fallbacks:
        if is_valid_python(cand):
            return cand

    return sys.executable


def extract_sequences_from_fasta(fasta_path: str) -> List[Tuple[str, str]]:
    """Read FASTA and return list of (taxon_id, clean_amino_acid_sequence).

    Gaps ('-') and whitespace are automatically removed to present the native
    unaligned biological context to the protein language model.
    """
    records = []
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    clean_seq = "".join(seq_parts).replace("-", "").replace(".", "").upper()
                    if clean_seq:
                        records.append((header, clean_seq))
                header = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            clean_seq = "".join(seq_parts).replace("-", "").replace(".", "").upper()
            if clean_seq:
                records.append((header, clean_seq))

    return records


def extract_sequences_from_pdb_dir(pdb_dir: str) -> List[Tuple[str, str]]:
    """Extract amino acid sequences from PDB structure files in a directory."""
    import glob
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'MSE': 'M'
    }

    files = []
    for ext in ("*.pdb", "*.cif", "*.mmcif"):
        files.extend(glob.glob(os.path.join(pdb_dir, ext)))
    files = sorted(set(files))

    records = []
    for p in files:
        basename = os.path.basename(p)
        taxon_id = os.path.splitext(basename)[0]
        seq_res = []
        seen_res = set()
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("ATOM  "):
                    atom_name = line[12:16].strip()
                    if atom_name == "CA":
                        res_name = line[17:20].strip()
                        chain_id = line[21:22].strip()
                        res_seq = line[22:27].strip()
                        key = (chain_id, res_seq)
                        if key not in seen_res:
                            seen_res.add(key)
                            seq_res.append(three_to_one.get(res_name, 'X'))
        clean_seq = "".join(seq_res)
        if clean_seq:
            records.append((taxon_id, clean_seq))

    return records


def generate_esm_embeddings(
    records: List[Tuple[str, str]],
    model_name_or_key: str = "esm2",
    batch_size: int = 1,
    max_length: int = 1024,
    device: Optional[str] = None
) -> Tuple[List[str], np.ndarray]:
    """Generate mean-pooled residue embeddings for each sequence using ESM-2 or ESM-C.

    Args:
        records: List of (taxon_id, sequence)
        model_name_or_key: 'esm2', 'esmc', or HuggingFace model ID
        batch_size: Inference batch size
        max_length: Maximum sequence length to truncate
        device: 'cuda', 'mps', or 'cpu' (auto-selected if None)

    Returns:
        (taxa_ids, embeddings_array of shape [N, hidden_dim])
    """
    import torch
    from transformers import AutoTokenizer, AutoModel

    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"

    print(f"[PLM] Computing embeddings using PyTorch on device: {device.upper()}")

    # Model resolution
    key_norm = model_name_or_key.lower().strip()
    if key_norm in ("esm2", "esm-2", "esm2_650m", "esm2_t33"):
        hf_model_id = "facebook/esm2_t33_650M_UR50D"
    elif key_norm in ("esmc", "esm-c", "esmc_600m", "esmc-600m"):
        hf_model_id = "biohub/ESMC-600M"
    else:
        hf_model_id = model_name_or_key

    print(f"[PLM] Loading model checkpoint: '{hf_model_id}'...")
    tokenizer = AutoTokenizer.from_pretrained(hf_model_id, trust_remote_code=True)
    try:
        model = AutoModel.from_pretrained(hf_model_id, trust_remote_code=True)
    except Exception as e:
        if "esmc" in key_norm:
            # Fallback for ESM-C: try evolutionaryscale package
            try:
                from esm.models.esmc import ESMC
                from esm.sdk.api import ESMProtein
                print("[PLM] Loading ESM-C via EvolutionaryScale esm SDK...")
                client = ESMC.from_pretrained("esmc_600m")
                names = []
                embs = []
                for idx, (tid, seq) in enumerate(records):
                    names.append(tid)
                    prot = ESMProtein(sequence=seq[:max_length])
                    p_tensor = client.encode(prot)
                    out = client.logits(p_tensor)
                    rep = out.logits.squeeze(0).mean(dim=0).cpu().numpy()
                    embs.append(rep)
                return names, np.array(embs)
            except Exception as esm_err:
                print(f"[PLM Warning] Could not load ESMC via esm package ({esm_err}). Falling back to ESM-2 650M.")
                hf_model_id = "facebook/esm2_t33_650M_UR50D"
                tokenizer = AutoTokenizer.from_pretrained(hf_model_id)
                model = AutoModel.from_pretrained(hf_model_id)
        else:
            raise e

    model = model.to(device)
    model.eval()

    names = []
    embeddings = []

    total = len(records)
    print(f"[PLM] Embedding {total} sequences with {hf_model_id}...")

    with torch.no_grad():
        for i in range(0, total, batch_size):
            batch = records[i : i + batch_size]
            batch_names = [item[0] for item in batch]
            batch_seqs = [item[1] for item in batch]

            if batch_size == 1:
                # Fast single-sequence path without padding to avoid MPS attention mask bugs
                inputs = tokenizer(
                    batch_seqs[0],
                    return_tensors="pt",
                    truncation=True,
                    max_length=max_length
                )
                inputs = {k: v.to(device) for k, v in inputs.items()}
                outputs = model(**inputs)
                hidden_states = outputs.last_hidden_state  # shape: [1, L, D]
                # Exclude [CLS] and [EOS] tokens
                if hidden_states.shape[1] > 2:
                    pooled = hidden_states[:, 1:-1, :].mean(dim=1).cpu().numpy()
                else:
                    pooled = hidden_states.mean(dim=1).cpu().numpy()
                names.append(batch_names[0])
                embeddings.append(pooled[0])
            else:
                inputs = tokenizer(
                    batch_seqs,
                    return_tensors="pt",
                    padding=True,
                    truncation=True,
                    max_length=max_length
                )
                inputs = {k: v.to(device) for k, v in inputs.items()}
                outputs = model(**inputs)
                hidden_states = outputs.last_hidden_state

                attention_mask = inputs.get("attention_mask")
                if attention_mask is not None:
                    mask = attention_mask.clone()
                    for b_idx in range(len(batch_seqs)):
                        seq_len = int(attention_mask[b_idx].sum().item())
                        if seq_len > 2:
                            mask[b_idx, 0] = 0
                            mask[b_idx, seq_len - 1] = 0
                        elif seq_len > 0:
                            mask[b_idx, :] = 1
                    mask = mask.unsqueeze(-1)
                    sum_embeddings = (hidden_states * mask).sum(dim=1)
                    sum_mask = mask.sum(dim=1).clamp(min=1e-9)
                    pooled = (sum_embeddings / sum_mask).cpu().numpy()
                else:
                    if hidden_states.shape[1] > 2:
                        pooled = hidden_states[:, 1:-1, :].mean(dim=1).cpu().numpy()
                    else:
                        pooled = hidden_states.mean(dim=1).cpu().numpy()

                for name, emb in zip(batch_names, pooled):
                    names.append(name)
                    embeddings.append(emb)

            curr_done = min(i + len(batch), total)
            if curr_done % 25 == 0 or curr_done == total:
                print(f"  Progress: {curr_done}/{total} sequences embedded", flush=True)

    emb_arr = np.array(embeddings, dtype=np.float32)
    print(f"[PLM] Generated embedding matrix: shape {emb_arr.shape} ({emb_arr.nbytes / 1024 / 1024:.2f} MB)\n")
    return names, emb_arr


def build_hierarchical_tree(
    taxa_names: List[str],
    embeddings: np.ndarray,
    clustering: str = "upgma",
    metric: str = "cosine"
) -> str:
    """Construct a phylogenetic/hierarchical clustering tree from embeddings in Newick format."""
    from scipy.spatial.distance import pdist, squareform

    metric_lower = metric.lower().strip()
    if metric_lower == "cosine":
        norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
        norms = np.where(norms == 0, 1e-12, norms)
        norm_embs = embeddings / norms
        condensed_dist = pdist(norm_embs, metric="cosine")
        condensed_dist = np.clip(condensed_dist, 0.0, 2.0)
    elif metric_lower in ("l1", "manhattan", "cityblock"):
        condensed_dist = pdist(embeddings, metric="cityblock")
    elif metric_lower == "euclidean":
        condensed_dist = pdist(embeddings, metric="euclidean")
    else:
        condensed_dist = pdist(embeddings, metric=metric)

    dist_matrix = squareform(condensed_dist)
    n = len(taxa_names)

    clustering_lower = clustering.lower().strip()

    # For Neighbor-Joining (NJ), use BioPython if available
    if clustering_lower == "nj":
        try:
            from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
            from Bio import Phylo

            matrix_lower = []
            for i in range(n):
                matrix_lower.append([float(dist_matrix[i, j]) for j in range(i + 1)])

            dm = DistanceMatrix(names=taxa_names, matrix=matrix_lower)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)

            out_buf = io.StringIO()
            Phylo.write(tree, out_buf, "newick")
            return out_buf.getvalue().strip()
        except Exception as e:
            print(f"[PLM Warning] NJ tree construction fallback to UPGMA: {e}")

    # Fast C-accelerated SciPy linkage for UPGMA / Average / Complete / Ward
    from scipy.cluster.hierarchy import linkage, to_tree

    method_map = {
        "upgma": "average",
        "average": "average",
        "nj": "average",
        "complete": "complete",
        "single": "single",
        "ward": "ward"
    }
    scipy_method = method_map.get(clustering_lower, "average")
    Z = linkage(condensed_dist, method=scipy_method)
    root_node, _ = to_tree(Z, rd=True)

    def node_to_newick(node):
        if node.is_leaf():
            return f"{taxa_names[node.id]}:{node.dist:.5f}"
        left_str = node_to_newick(node.left)
        right_str = node_to_newick(node.right)
        branch_len = max(0.0, node.dist - max(node.left.dist, node.right.dist))
        return f"({left_str},{right_str}):{branch_len:.5f}"

    newick_str = f"({node_to_newick(root_node.left)},{node_to_newick(root_node.right)}):0.00000;"
    return newick_str, Z, condensed_dist


def compute_silhouette_profile(
    Z: np.ndarray,
    condensed_dist: np.ndarray,
    taxa_names: List[str],
    max_k: int = 20
) -> Dict:
    """Compute silhouette scores for hierarchical cuts from k=2 up to max_k.

    Identifies 'reasonable peaks' (local maxima in the silhouette profile)
    indicating optimal, cohesive biological sub-clades.
    """
    from scipy.cluster.hierarchy import fcluster
    from scipy.spatial.distance import squareform

    dist_mat = squareform(condensed_dist)
    n = len(taxa_names)
    max_eval_k = min(max_k, n - 1)
    if max_eval_k < 2:
        return {"profile": [], "best_k": 2, "best_score": 0.0, "peaks": []}

    profile = []
    scores = []

    for k in range(2, max_eval_k + 1):
        clusts = fcluster(Z, t=k, criterion="maxclust")
        unique_c = np.unique(clusts)
        if len(unique_c) < 2:
            continue

        s_scores = np.zeros(n)
        for i in range(n):
            ci = clusts[i]
            same_mask = (clusts == ci)
            same_count = np.sum(same_mask)
            if same_count <= 1:
                s_scores[i] = 0.0
                continue

            a_i = np.sum(dist_mat[i, same_mask]) / (same_count - 1)
            b_i = min(np.mean(dist_mat[i, clusts == oc]) for oc in unique_c if oc != ci)
            denom = max(a_i, b_i)
            s_scores[i] = (b_i - a_i) / denom if denom > 0 else 0.0

        mean_s = float(np.mean(s_scores))
        scores.append((k, mean_s, clusts))

    if not scores:
        return {"profile": [], "best_k": 2, "best_score": 0.0, "peaks": []}

    peaks = []
    for idx in range(len(scores)):
        k_val, s_val, clusts = scores[idx]
        is_left = (idx == 0) or (s_val >= scores[idx - 1][1])
        is_right = (idx == len(scores) - 1) or (s_val >= scores[idx + 1][1])
        is_peak = is_left and is_right and (s_val > 0.35 or idx == 0)

        sizes = [int(np.sum(clusts == c)) for c in np.unique(clusts)]
        entry = {
            "k": int(k_val),
            "score": round(float(s_val), 4),
            "is_peak": bool(is_peak),
            "cluster_sizes": sorted(sizes, reverse=True)
        }
        profile.append(entry)
        if is_peak:
            peaks.append({"k": int(k_val), "score": round(float(s_val), 4)})

    peaks = sorted(peaks, key=lambda p: p["score"], reverse=True)
    best = max(scores, key=lambda x: x[1])

    return {
        "n_taxa": int(n),
        "best_k": int(best[0]),
        "best_score": round(float(best[1]), 4),
        "peaks": peaks,
        "profile": profile
    }


def run_embedding_pipeline(
    input_source: str,
    output_dir: str = "plm_embeddings",
    model: str = "esm2",
    clustering: str = "upgma",
    metric: str = "cosine",
    prefix: str = "viral_plm",
    batch_size: int = 1,
    max_length: int = 1024,
    compute_silhouette: bool = True
) -> Tuple[str, str]:
    """Execute complete embedding, tree construction, and silhouette analysis workflow."""
    import json
    os.makedirs(output_dir, exist_ok=True)

    if os.path.isfile(input_source):
        records = extract_sequences_from_fasta(input_source)
    elif os.path.isdir(input_source):
        candidate_fa = os.path.join(input_source, "foldmason.fasta_aa.fa")
        if os.path.isfile(candidate_fa):
            records = extract_sequences_from_fasta(candidate_fa)
        else:
            records = extract_sequences_from_pdb_dir(input_source)
    else:
        raise FileNotFoundError(f"Input source '{input_source}' not found.")

    if len(records) < 2:
        raise ValueError(f"Need at least 2 sequences to build hierarchical tree, found {len(records)}.")

    print(f"[PLM] Successfully loaded {len(records)} protein sequences from '{input_source}'.")

    taxa_names, embeddings = generate_esm_embeddings(
        records=records,
        model_name_or_key=model,
        batch_size=batch_size,
        max_length=max_length
    )

    tree_res = build_hierarchical_tree(
        taxa_names=taxa_names,
        embeddings=embeddings,
        clustering=clustering,
        metric=metric
    )
    if isinstance(tree_res, tuple):
        newick, Z, condensed_dist = tree_res
    else:
        newick = tree_res
        Z, condensed_dist = None, None

    # Compute Silhouette Profile across increasing cut heights
    sil_info = None
    if compute_silhouette and Z is not None and condensed_dist is not None:
        try:
            sil_info = compute_silhouette_profile(Z, condensed_dist, taxa_names, max_k=min(25, len(taxa_names)-1))
            sil_path = os.path.join(output_dir, f"{prefix}_silhouette_{model}.json")
            with open(sil_path, "w", encoding="utf-8") as f:
                json.dump(sil_info, f, indent=2)
            print(f"[PLM] Silhouette score profile saved to: '{sil_path}' (Optimal k={sil_info['best_k']} with S={sil_info['best_score']})")
        except Exception as e:
            print(f"[PLM Warning] Could not compute silhouette profile: {e}")

    npz_path = os.path.join(output_dir, f"{prefix}_embeddings_{model}.npz")
    np.savez_compressed(
        npz_path,
        taxa=np.array(taxa_names),
        embeddings=embeddings,
        model=model,
        clustering=clustering,
        metric=metric,
        silhouette=json.dumps(sil_info) if sil_info else ""
    )
    print(f"[PLM] Embeddings saved to: '{npz_path}'")

    treefile_path = os.path.join(output_dir, f"{prefix}_tree_{model}.treefile")
    with open(treefile_path, "w", encoding="utf-8") as f:
        f.write(newick + "\n")

    print(f"[PLM] Hierarchical clustering tree ({clustering.upper()}, {metric} distance) saved to: '{treefile_path}'.\n")
    return treefile_path, npz_path


def main():
    parser = argparse.ArgumentParser(
        description="Embed protein sequences with ESM-2/ESM-C PLM and construct hierarchical clustering trees in Newick format."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input FASTA file (e.g. foldmason.fasta_aa.fa) or directory of PDB structures."
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="plm_results",
        help="Directory to store resulting treefile and embeddings (default: plm_results)."
    )
    parser.add_argument(
        "-m", "--model",
        default="esm2",
        choices=["esm2", "esmc", "facebook/esm2_t33_650M_UR50D", "biohub/ESMC-600M"],
        help="Protein Language Model to use (default: esm2 [650M])."
    )
    parser.add_argument(
        "-c", "--clustering",
        default="upgma",
        choices=["upgma", "nj", "average", "complete", "single", "ward"],
        help="Hierarchical clustering algorithm (default: upgma)."
    )
    parser.add_argument(
        "--metric",
        default="cosine",
        choices=["cosine", "euclidean", "l1", "cityblock", "manhattan"],
        help="Pairwise distance metric (default: cosine, choices: cosine, euclidean, l1, manhattan)."
    )
    parser.add_argument(
        "-p", "--prefix",
        default="viral_plm",
        help="Output file prefix (default: viral_plm)."
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1,
        help="Inference batch size (default: 1)."
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=1024,
        help="Max sequence length for transformer truncation (default: 1024)."
    )

    args = parser.parse_args()

    # Check if we need to dispatch to another python interpreter
    try:
        import torch
        import transformers
        import scipy
    except ImportError:
        target_py = find_torch_python()
        if target_py and target_py != sys.executable:
            print(f"[Environment] Dispatching embedding task to Python with PyTorch: {target_py}")
            import subprocess
            cmd = [target_py, __file__] + sys.argv[1:]
            res = subprocess.run(cmd)
            sys.exit(res.returncode)
        else:
            print("ERROR: PyTorch and Transformers are required for PLM embedding extraction.")
            print("Please run within a Python environment having 'torch', 'transformers', and 'scipy'.")
            sys.exit(1)

    run_embedding_pipeline(
        input_source=args.input,
        output_dir=args.output_dir,
        model=args.model,
        clustering=args.clustering,
        metric=args.metric,
        prefix=args.prefix,
        batch_size=args.batch_size,
        max_length=args.max_length
    )


if __name__ == "__main__":
    main()
