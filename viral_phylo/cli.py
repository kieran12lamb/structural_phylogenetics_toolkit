"""Command Line Interface for the Viral Structural Phylogenetics Toolkit."""

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

from viral_phylo.binaries import (
    find_iqtree_bin,
    find_foldmason_bin,
    find_mafft_bin,
    find_torch_python,
)
from viral_phylo.matrices import ensure_matrix_file, ensure_3di_matrix
from viral_phylo.fetch import fetch_structures, fetch_alphafold_structures
from viral_phylo.alignment import (
    align_structures,
    filter_alignment_by_coverage,
    partition_structures_by_coverage,
    compute_alignment_coverage,
)
from viral_phylo.tree import build_tree

def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Viral Structural Phylogenetics & PLM Embedding Suite.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Subcommand to run")

    # Subcommand: fetch
    p_fetch = subparsers.add_parser("fetch", help="Download viral PDB structures from Viro3D or AlphaFold Database.")
    p_fetch.add_argument("--source", choices=["viro3d", "alphafold", "afdb"], default="viro3d", help="Structure repository: 'viro3d' or 'alphafold' (default: viro3d)")
    p_fetch.add_argument("-q", "--qualifier", default="glycoprotein", help="Protein search term (default: glycoprotein)")
    p_fetch.add_argument("-u", "--uniprot", default=None, help="UniProt ID(s) for AlphaFold DB, comma-separated (e.g. 'P00520,P04637')")
    p_fetch.add_argument("--uniprot-file", default=None, help="Path to text file with one UniProt ID per line")
    p_fetch.add_argument("-m", "--max-sequences", type=int, default=6, help="Max structures to download (default: 6)")
    p_fetch.add_argument("--format", choices=["pdb", "cif"], default="pdb", help="File format for AlphaFold DB (default: pdb)")
    p_fetch.add_argument("--download-pae", action="store_true", help="Also download Predicted Aligned Error (PAE) JSON from AlphaFold DB")
    p_fetch.add_argument("-o", "--output-dir", default="viro_3d_structures", help="Output directory for PDBs")

    # Subcommand: align
    p_align = subparsers.add_parser("align", help="Align PDB structures with FoldMason or MAFFT to generate 3Di alignments.")
    p_align.add_argument("-i", "--folder", default="viro_3d_structures", help="Folder containing PDB files (local or fetched)")
    p_align.add_argument("-o", "--output-dir", default="foldmason_alignments", help="Output directory for MSA")
    p_align.add_argument("--aligner", choices=["foldmason", "mafft"], default="foldmason", help="Multiple sequence alignment engine (default: foldmason)")
    p_align.add_argument("--foldmason-bin", default=None, help="Path to foldmason binary")
    p_align.add_argument("--mafft-bin", default=None, help="Path to mafft binary (if using --aligner mafft)")
    p_align.add_argument("--mafft-matrix", default="matrices/mat3di.out", help="3Di substitution matrix for MAFFT (default: matrices/mat3di.out)")
    p_align.add_argument("--min-coverage", type=float, default=0.70, help="Minimum alignment/sequence coverage threshold (0.0 to 1.0, default: 0.70)")
    p_align.add_argument("--multi-alignment", action="store_true", help="Partition dataset into multiple high-coverage sub-alignments (>= min-coverage)")
    p_align.add_argument("--filter-coverage", action="store_true", help="Filter alignment to only retain sequences meeting >= min-coverage")

    # Subcommand: tree
    p_tree = subparsers.add_parser("tree", help="Build structural and/or amino acid phylogeny using IQ-TREE or FoldMason.")
    p_tree.add_argument("-a", "--alignment", default="foldmason_alignments/foldmason.fasta_3di.fa", help="3Di alignment file or directory")
    p_tree.add_argument("--alignment-aa", default=None, help="Amino acid alignment file (default: auto-detected)")
    p_tree.add_argument("--tree-type", choices=["3di", "aa", "both", "tanglegram"], default="both", help="Phylogeny type: '3di' (structural), 'aa' (sequence), or 'both'/'tanglegram' (default: both)")
    p_tree.add_argument("-m", "--method", choices=["iqtree", "foldmason"], default="iqtree", help="Tree method (default: iqtree)")
    p_tree.add_argument("--matrix", choices=["alphafold", "af", "esmfold", "llm", "both", "auto"], default="both", help="3Di substitution matrix")
    p_tree.add_argument("--rate-heterogeneity", default="auto", help="Rate heterogeneity: '+G4', '+I+G4', '+R', or 'auto' (default: auto)")
    p_tree.add_argument("--criterion", choices=["BIC", "AIC", "AICc"], default="BIC", help="Model selection criterion (default: BIC)")
    p_tree.add_argument("-b", "--bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (default: 1000)")
    p_tree.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates (default: 1000)")
    p_tree.add_argument("-t", "--threads", default="AUTO", help="CPU threads (default: AUTO)")
    p_tree.add_argument("--fast", action="store_true", help="Enable fast heuristic search mode (auto-enabled for large datasets)")
    p_tree.add_argument("-o", "--output-dir", default="phylogeny_results", help="Output directory for trees")
    p_tree.add_argument("-p", "--prefix", default="viral_tree", help="Tree output filename prefix")
    p_tree.add_argument("--iqtree-bin", default=None, help="Path to iqtree binary")

    # Subcommand: embed (PLM Embeddings & Hierarchical Clustering)
    p_embed = subparsers.add_parser("embed", help="Extract PLM embeddings (ESM-2 / ESM-C) and construct hierarchical clustering trees.")
    p_embed.add_argument("-i", "--input", "--fasta", dest="input", required=True, help="Path to input FASTA file (e.g. foldmason.fasta_aa.fa) or directory of PDB structures.")
    p_embed.add_argument("-o", "--output-dir", default="plm_results", help="Directory to store resulting treefile and embeddings (default: plm_results).")
    p_embed.add_argument("-m", "--model", default="esm2", choices=["esm2", "esmc", "facebook/esm2_t33_650M_UR50D", "biohub/ESMC-600M"], help="Protein language model (default: esm2 [650M]).")
    p_embed.add_argument("-c", "--clustering", default="upgma", choices=["upgma", "nj", "average", "complete", "single", "ward"], help="Hierarchical clustering method (default: upgma).")
    p_embed.add_argument("--metric", default="cosine", choices=["cosine", "euclidean", "l1", "cityblock", "manhattan"], help="Pairwise distance metric (default: cosine).")
    p_embed.add_argument("-p", "--prefix", default="viral_plm", help="Output file prefix (default: viral_plm).")
    p_embed.add_argument("--batch-size", type=int, default=1, help="Inference batch size (default: 1).")
    p_embed.add_argument("--max-length", type=int, default=1024, help="Max sequence length for truncation (default: 1024).")

    # Subcommand: pipeline (all-in-one)
    p_pipe = subparsers.add_parser("pipeline", help="Run full pipeline: fetch/local -> align -> tree (+ optional PLM embed).")
    p_pipe.add_argument("-i", "--input-folder", default=None, help="Local directory of PDB structures (skips Viro3D fetch if specified)")
    p_pipe.add_argument("--source", choices=["viro3d", "alphafold", "afdb"], default="viro3d", help="Structure repository: 'viro3d' or 'alphafold' (default: viro3d)")
    p_pipe.add_argument("-q", "--qualifier", default="glycoprotein", help="Protein search term if fetching (default: glycoprotein)")
    p_pipe.add_argument("-u", "--uniprot", default=None, help="UniProt ID(s) for AlphaFold DB (comma-separated)")
    p_pipe.add_argument("--uniprot-file", default=None, help="Path to text file with one UniProt ID per line")
    p_pipe.add_argument("-c", "--count", type=int, default=6, help="Number of structures to fetch (default: 6)")
    p_pipe.add_argument("--format", choices=["pdb", "cif"], default="pdb", help="Format for AlphaFold DB: pdb or cif (default: pdb)")
    p_pipe.add_argument("--download-pae", action="store_true", help="Download PAE error matrix from AlphaFold DB")
    p_pipe.add_argument("--aligner", choices=["foldmason", "mafft"], default="foldmason", help="Multiple sequence alignment engine (default: foldmason)")
    p_pipe.add_argument("--mafft-bin", default=None, help="Path to mafft binary (if using --aligner mafft)")
    p_pipe.add_argument("--mafft-matrix", default="matrices/mat3di.out", help="3Di substitution matrix for MAFFT (default: matrices/mat3di.out)")
    p_pipe.add_argument("--min-coverage", type=float, default=0.70, help="Minimum alignment/sequence coverage threshold (0.0 to 1.0, default: 0.70)")
    p_pipe.add_argument("--multi-alignment", action="store_true", help="Partition dataset into multiple high-coverage sub-alignments (>= min-coverage)")
    p_pipe.add_argument("--filter-coverage", action="store_true", help="Filter alignment to only retain sequences meeting >= min-coverage")
    p_pipe.add_argument("--tree-type", choices=["3di", "aa", "both", "tanglegram"], default="both", help="Phylogeny type (default: both)")
    p_pipe.add_argument("-m", "--method", choices=["iqtree", "foldmason"], default="iqtree", help="Tree method (default: iqtree)")
    p_pipe.add_argument("--matrix", choices=["alphafold", "af", "esmfold", "llm", "both", "auto"], default="both", help="3Di matrix (default: both)")
    p_pipe.add_argument("--rate-heterogeneity", default="auto", help="Rate heterogeneity model (default: auto)")
    p_pipe.add_argument("-b", "--bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (default: 1000)")
    p_pipe.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates (default: 1000)")
    p_pipe.add_argument("-t", "--threads", default="AUTO", help="CPU threads (default: AUTO)")
    p_pipe.add_argument("-meta", "--metadata", default=None, help="Path to metadata file (.xlsx, .csv, .tsv, .json). If omitted, inferred from structures or Viro3D API.")
    p_pipe.add_argument("--fast", action="store_true", help="Enable fast search mode")
    p_pipe.add_argument("--embed", action="store_true", help="Extract PLM embeddings and build hierarchical clustering tree (ESM-2 / ESM-C)")
    p_pipe.add_argument("--embed-model", default="esm2", choices=["esm2", "esmc"], help="PLM model for embeddings (default: esm2)")
    p_pipe.add_argument("--embed-clustering", default="upgma", choices=["upgma", "nj", "average", "complete"], help="Clustering method for PLM tree (default: upgma)")
    p_pipe.add_argument("--embed-metric", default="cosine", choices=["cosine", "euclidean", "l1", "cityblock", "manhattan"], help="Pairwise distance metric (default: cosine)")
    p_pipe.add_argument("-o", "--output-dir", default="glycoprotein_workflow", help="Parent output directory")

    return parser


def main():
    parser = build_cli_parser()
    args = parser.parse_args()

    if args.subcommand == "fetch":
        source = getattr(args, "source", "viro3d")
        if source in ("alphafold", "afdb") or getattr(args, "uniprot", None) or getattr(args, "uniprot_file", None):
            u_ids = []
            if getattr(args, "uniprot", None):
                u_ids.extend([x.strip() for x in args.uniprot.split(",") if x.strip()])
            if getattr(args, "uniprot_file", None) and os.path.isfile(args.uniprot_file):
                with open(args.uniprot_file, "r") as uf:
                    u_ids.extend([line.strip() for line in uf if line.strip() and not line.startswith("#")])
            fetch_alphafold_structures(
                uniprot_ids=u_ids if u_ids else None,
                query=args.qualifier if not u_ids else None,
                max_sequences=args.max_sequences,
                output_dir=args.output_dir,
                file_format=getattr(args, "format", "pdb"),
                download_pae=getattr(args, "download_pae", False),
            )
        else:
            fetch_structures(args.qualifier, args.max_sequences, args.output_dir)
    elif args.subcommand == "align":
        align_structures(
            args.folder,
            args.output_dir,
            foldmason_bin=args.foldmason_bin,
            aligner=args.aligner,
            mafft_bin=args.mafft_bin,
            mafft_matrix=args.mafft_matrix,
            min_coverage=getattr(args, "min_coverage", 0.70),
            multi_alignment=getattr(args, "multi_alignment", False),
            filter_coverage=getattr(args, "filter_coverage", False),
        )
    elif args.subcommand == "tree":
        aln_file = args.alignment
        if os.path.isdir(aln_file):
            for cand in ["mafft.fasta_3di.fa", "foldmason.fasta_3di.fa", "foldmason_3di.fa"]:
                if os.path.isfile(os.path.join(aln_file, cand)):
                    aln_file = os.path.join(aln_file, cand)
                    break
            else:
                aln_file = os.path.join(aln_file, "foldmason.fasta_3di.fa")
        build_tree(
            alignment_file=aln_file,
            method=args.method,
            matrix=args.matrix,
            rate_het=args.rate_heterogeneity,
            criterion=args.criterion,
            bootstrap=args.bootstrap,
            alrt=args.alrt,
            threads=args.threads,
            output_dir=args.output_dir,
            prefix=args.prefix,
            iqtree_bin=args.iqtree_bin,
            tree_type=args.tree_type,
            alignment_aa=args.alignment_aa,
            fast=args.fast,
        )
    elif args.subcommand == "embed":
        embed_script = os.path.join(os.path.dirname(__file__), "embed_and_cluster.py")
        if not os.path.isfile(embed_script):
            embed_script = "scripts/embed_and_cluster.py"
        target_py = find_torch_python()
        cmd = [
            target_py, embed_script,
            "-i", args.input,
            "-o", args.output_dir,
            "-m", args.model,
            "-c", args.clustering,
            "--metric", args.metric,
            "-p", args.prefix,
            "--batch-size", str(args.batch_size),
            "--max-length", str(args.max_length),
        ]
        print(f"[PLM] Running embedding extraction using interpreter: {target_py}")
        sub_env = dict(os.environ)
        sub_env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        subprocess.run(cmd, check=True, env=sub_env)
    elif args.subcommand == "pipeline":
        out_base = args.output_dir
        aln_dir = os.path.join(out_base, "alignment")
        phy_dir = os.path.join(out_base, "phylogeny")

        # 1. Fetch or Local Input
        source = getattr(args, "source", "viro3d")
        if args.input_folder:
            if not os.path.isdir(args.input_folder):
                raise FileNotFoundError(f"Specified local structures directory '{args.input_folder}' does not exist.")
            pdb_dir = args.input_folder
            print(f"[Pipeline] Using local structure folder: '{pdb_dir}'")
        elif source in ("alphafold", "afdb") or getattr(args, "uniprot", None) or getattr(args, "uniprot_file", None):
            pdb_dir = os.path.join(out_base, "structures")
            u_ids = []
            if getattr(args, "uniprot", None):
                u_ids.extend([x.strip() for x in args.uniprot.split(",") if x.strip()])
            if getattr(args, "uniprot_file", None) and os.path.isfile(args.uniprot_file):
                with open(args.uniprot_file, "r") as uf:
                    u_ids.extend([line.strip() for line in uf if line.strip() and not line.startswith("#")])
            fetch_alphafold_structures(
                uniprot_ids=u_ids if u_ids else None,
                query=args.qualifier if not u_ids else None,
                max_sequences=args.count,
                output_dir=pdb_dir,
                file_format=getattr(args, "format", "pdb"),
                download_pae=getattr(args, "download_pae", False),
            )
        else:
            pdb_dir = os.path.join(out_base, "structures")
            fetch_structures(args.qualifier, args.count, pdb_dir)

        # 2. Align (with optional Coverage Threshold & Multi-Alignment Partitioning)
        min_cov = getattr(args, "min_coverage", 0.70)
        multi_aln = getattr(args, "multi_alignment", False)
        filter_cov = getattr(args, "filter_coverage", False)

        aln_3di = align_structures(
            pdb_dir,
            aln_dir,
            aligner=args.aligner,
            mafft_bin=getattr(args, "mafft_bin", None),
            mafft_matrix=getattr(args, "mafft_matrix", "matrices/mat3di.out"),
            min_coverage=min_cov,
            multi_alignment=multi_aln,
            filter_coverage=filter_cov,
        )

        multi_summary = None
        if multi_aln:
            summary_f = os.path.join(aln_dir, "multi_alignment_summary.json")
            if os.path.isfile(summary_f):
                with open(summary_f, "r", encoding="utf-8") as f:
                    multi_summary = json.load(f)
                # Copy summary to parent output directory
                shutil.copy2(summary_f, os.path.join(out_base, "multi_alignment_summary.json"))
        # 3. Tree
        if args.input_folder:
            prefix = "local_tree"
        elif source in ("alphafold", "afdb"):
            prefix = "alphafold_tree"
        else:
            prefix = f"{args.qualifier}_tree"

        build_tree(
            alignment_file=aln_3di,
            method=args.method,
            matrix=args.matrix,
            rate_het=args.rate_heterogeneity,
            criterion="BIC",
            bootstrap=args.bootstrap,
            alrt=args.alrt,
            threads=args.threads,
            output_dir=phy_dir,
            prefix=prefix,
            tree_type=args.tree_type,
            fast=args.fast,
        )

        # 3b. Build trees for coverage-partitioned clusters
        if multi_summary and "clusters" in multi_summary:
            print(f"\n[Pipeline] Building phylogenetic trees for {len(multi_summary['clusters'])} coverage-partitioned cluster(s)...")
            for c in multi_summary["clusters"]:
                if c.get("taxa_count", 0) >= 3 and c.get("aln_3di") and os.path.isfile(c["aln_3di"]):
                    c_phy_dir = os.path.join(phy_dir, c["name"])
                    print(f"[Pipeline] Inferring phylogeny for {c['name']} ({c['taxa_count']} taxa)...")
                    try:
                        build_tree(
                            alignment_file=c["aln_3di"],
                            method=args.method,
                            matrix=args.matrix,
                            rate_het=args.rate_heterogeneity,
                            criterion="BIC",
                            bootstrap=args.bootstrap,
                            alrt=args.alrt,
                            threads=args.threads,
                            output_dir=c_phy_dir,
                            prefix=f"{prefix}_{c['name']}",
                            tree_type=args.tree_type,
                            fast=args.fast,
                        )
                    except Exception as e:
                        print(f"  [!] Notice: Tree inference for {c['name']} returned: {e}")

        # 4. Universal Metadata Processing
        if parse_metadata:
            print(f"\n[Pipeline] Analyzing and indexing metadata...")
            meta_src = getattr(args, "metadata", None)
            if not meta_src and not args.input_folder:
                cand_meta = os.path.join(pdb_dir, "taxa_metadata.json")
                if os.path.isfile(cand_meta):
                    meta_src = cand_meta
            meta_res = parse_metadata(meta_src, structure_dir=pdb_dir)
            meta_out_path = os.path.join(out_base, "taxa_metadata.json")
            with open(meta_out_path, "w", encoding="utf-8") as f:
                json.dump(meta_res, f, indent=2)
            print(f"[Pipeline] Saved standardized metadata schema to: '{meta_out_path}' ({len(meta_res['taxa'])} taxa, {len(meta_res['columns'])} columns).")

        # 5. Optional PLM Embedding & Hierarchical Clustering Tree
        if getattr(args, "embed", False):
            print(f"\n[Pipeline] Extracting PLM embeddings ({args.embed_model}) & constructing hierarchical tree ({args.embed_clustering})...")
            embed_script = os.path.join(os.path.dirname(__file__), "embed_and_cluster.py")
            if not os.path.isfile(embed_script):
                embed_script = "scripts/embed_and_cluster.py"
            aa_aln = os.path.join(aln_dir, "foldmason.fasta_aa.fa")
            for cand in ["mafft.fasta_aa.fa", "foldmason.fasta_aa.fa", "foldmason_aa.fa"]:
                if os.path.isfile(os.path.join(aln_dir, cand)):
                    aa_aln = os.path.join(aln_dir, cand)
                    break
            embed_src = aa_aln if os.path.isfile(aa_aln) else pdb_dir
            target_py = find_torch_python()
            cmd = [
                target_py, embed_script,
                "-i", embed_src,
                "-o", phy_dir,
                "-m", args.embed_model,
                "-c", args.embed_clustering,
                "--metric", args.embed_metric,
                "-p", prefix,
                "--batch-size", "1",
            ]
            sub_env = dict(os.environ)
            sub_env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
            subprocess.run(cmd, check=True, env=sub_env)

        # 6. Generate Standalone Interactive HTML Visualizations & Datasets
        print("\n[Pipeline] Compiling interactive HTML suite and supporting assets...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        build_tree_script = os.path.join(script_dir, "build_dynamic_interactive_tree.py")
        build_aln_script = os.path.join(script_dir, "build_alignments_data.py")

        try:
            ca_dict = {}
            pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")) + glob.glob(os.path.join(pdb_dir, "*.cif")))
            for p in pdb_files:
                taxon_id = os.path.splitext(os.path.basename(p))[0]
                ca = []
                with open(p, "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if line.startswith("ATOM") and line[12:16].strip() == "CA":
                            try:
                                resname = line[17:20].strip()
                                resnum = int(line[22:26].strip())
                                x = round(float(line[30:38].strip()), 1)
                                y = round(float(line[38:46].strip()), 1)
                                z = round(float(line[46:54].strip()), 1)
                                b_factor = round(float(line[60:66].strip()), 1)
                                ca.append([x, y, z, b_factor, resnum, resname])
                            except Exception:
                                continue
                if ca:
                    ca_dict[taxon_id] = ca
            if ca_dict:
                ca_js = "window.CA_STRUCTURES = Object.assign(window.CA_STRUCTURES || {}, " + json.dumps(ca_dict) + ");\n"
                with open(os.path.join(out_base, "ca_structures.js"), "w", encoding="utf-8") as f:
                    f.write(ca_js)

            # Also ensure all cohort structure packs exist in out_base so all datasets load cleanly
            for sf in ["ca_500_structures.js", "ca_1193_structures.js", "ca_100_structures.js", "alignments_data.js"]:
                results_sf = os.path.join("results", sf)
                dst_sf = os.path.join(out_base, sf)
                if os.path.isfile(results_sf) and not os.path.isfile(dst_sf):
                    try:
                        shutil.copy2(results_sf, dst_sf)
                    except Exception:
                        pass
        except Exception as e:
            print(f"[Pipeline] Notice: C-alpha extraction encountered: {e}")

        try:
            target_py = sys.executable
            if os.path.isfile(build_aln_script):
                subprocess.run([target_py, build_aln_script], check=False)
            if os.path.isfile(build_tree_script):
                subprocess.run([target_py, build_tree_script], check=False)

            wf_html = os.path.join(out_base, "interactive_tree.html")
            root_html = os.path.abspath("interactive_tree.html")
            print("\n[Pipeline] ========================================================")
            print("[Pipeline] ✅ Full Pipeline Execution & Visualization Complete!")
            if os.path.isfile(wf_html):
                print(f"[Pipeline] 📊 Standalone Visualizer: file://{os.path.abspath(wf_html)}")
            if os.path.isfile(root_html):
                print(f"[Pipeline] 🌐 Dashboard Visualizer:  file://{root_html}")
            print("[Pipeline] ========================================================\n")
        except Exception as e:
            print(f"[Pipeline] Notice: Automated HTML compilation returned: {e}")




if __name__ == "__main__":
    main()
