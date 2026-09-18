"""Structure and metadata fetching from AlphaFold Database (EBI) and Viro3D."""

import json
import os
import urllib.parse
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Dict, Any

BASE_URL = "https://viro3d.cvr.gla.ac.uk/api"


def fetch_structures(qualifier: str, max_sequences: int, output_dir: str, workers: int = 16) -> List[str]:
    """Download viral structures from Viro3D with concurrent worker threads."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n[Viro3D] Querying API for '{qualifier}' (target: {max_sequences} structures)...")
    
    url = f"{BASE_URL}/proteins/protein_name/"
    query = urllib.parse.urlencode({"qualifier": qualifier, "page_size": max_sequences, "page_num": 1})
    req = urllib.request.Request(f"{url}?{query}", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode())

    items = data.get("protein_structures", [])[:max_sequences]
    if not items:
        print(f"No structures found for '{qualifier}'.")
        return []

    # Save metadata JSON for downstream visualizers and tree annotations
    meta_path = os.path.join(output_dir, "taxa_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(items, f, indent=2)
    print(f"[Viro3D] Saved metadata index for {len(items)} structures to '{meta_path}'.")

    downloaded = []
    print(f"[Viro3D] Downloading {len(items)} PDB structures concurrently ({workers} workers)...")

    def _dl_one(item):
        rid = item["record_id"]
        pdb_url = f"{BASE_URL}/pdb/CF-{rid}_relaxed.pdb"
        out_file = os.path.join(output_dir, f"{rid}.pdb")
        if os.path.isfile(out_file) and os.path.getsize(out_file) > 500:
            return rid, out_file, True
        try:
            p_req = urllib.request.Request(pdb_url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(p_req, timeout=30) as r, open(out_file, "wb") as f:
                f.write(r.read())
            return rid, out_file, True
        except Exception as e:
            return rid, str(e), False

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_dl_one, it) for it in items]
        completed = 0
        for fut in as_completed(futures):
            rid, path_or_err, success = fut.result()
            completed += 1
            if success:
                downloaded.append(path_or_err)
            if completed % 50 == 0 or completed == len(items):
                print(f"  Progress: {completed}/{len(items)} downloaded ({len(downloaded)} successful)")

    print(f"[Viro3D] Download complete: {len(downloaded)}/{len(items)} structure(s) saved to '{output_dir}'.\n")
    return downloaded


def query_uniprot_for_accessions(query: str, max_sequences: int = 10) -> List[str]:
    """Query UniProt REST API to resolve keywords/gene names to UniProt accession IDs."""
    print(f"[AlphaFold DB] Resolving search query '{query}' via UniProt API...")
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = urllib.parse.urlencode({
        "query": query,
        "size": max_sequences,
        "fields": "accession,id,gene_names,organism_name,length"
    })
    req = urllib.request.Request(f"{url}?{params}", headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
            results = data.get("results", [])
            accessions = [r.get("primaryAccession") for r in results if r.get("primaryAccession")]
            print(f"[AlphaFold DB] Found {len(accessions)} UniProt accession(s) for '{query}': {', '.join(accessions)}")
            return accessions
    except Exception as e:
        print(f"[!] UniProt search query failed: {e}")
        return []


def fetch_alphafold_structures(
    uniprot_ids: Optional[List[str]] = None,
    query: Optional[str] = None,
    max_sequences: int = 6,
    output_dir: str = "afdb_structures",
    file_format: str = "pdb",
    download_pae: bool = False,
    workers: int = 8,
) -> List[str]:
    """Download AlphaFold structures (.pdb/.cif) and optional PAE matrices from the AlphaFold Database (EBI)."""
    os.makedirs(output_dir, exist_ok=True)

    target_ids = []
    if uniprot_ids:
        for item in uniprot_ids:
            if isinstance(item, str):
                for sub in item.split(","):
                    clean = sub.strip().upper()
                    if clean and clean not in target_ids:
                        target_ids.append(clean)

    if not target_ids and query:
        target_ids = query_uniprot_for_accessions(query, max_sequences=max_sequences)

    if not target_ids:
        print("[!] No UniProt IDs provided or resolved from query. Please specify -u/--uniprot or a valid query term.")
        return []

    target_ids = target_ids[:max_sequences]
    print(f"\n[AlphaFold DB] Querying AFDB API for {len(target_ids)} UniProt ID(s): {', '.join(target_ids)}...")

    metadata_list = []
    download_tasks = []

    for uid in target_ids:
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
        req = urllib.request.Request(api_url, headers={"User-Agent": "Mozilla/5.0"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                entries = json.loads(resp.read().decode())
                if not entries:
                    print(f"  [!] No AlphaFold prediction found for '{uid}'")
                    continue

                # Pick canonical or longest entry
                entry = next((e for e in entries if e.get("uniprotAccession") == uid), None)
                if entry is None:
                    entry = max(entries, key=lambda e: e.get("sequenceEnd", 0))

                entry_acc = entry.get("uniprotAccession", uid)
                plddt = entry.get("globalMetricValue")

                plddt_cat = "Unknown"
                if plddt is not None:
                    if plddt >= 90:
                        plddt_cat = "Very high (>90)"
                    elif plddt >= 70:
                        plddt_cat = "Confident (70-90)"
                    elif plddt >= 50:
                        plddt_cat = "Low (50-70)"
                    else:
                        plddt_cat = "Disordered (<50)"

                prefer_pdb = file_format.lower() == "pdb"
                struct_url = entry.get("pdbUrl") if prefer_pdb else entry.get("cifUrl")
                if not struct_url:
                    struct_url = entry.get("pdbUrl") or entry.get("cifUrl")

                ext = ".pdb" if (struct_url and struct_url.endswith(".pdb")) else ".cif"
                out_struct_file = os.path.join(output_dir, f"{entry_acc}{ext}")

                pae_url = entry.get("paeDocUrl") if download_pae else None
                out_pae_file = os.path.join(output_dir, f"AF-{entry_acc}-F1-predicted_aligned_error.json") if pae_url else None

                meta_entry = {
                    "record_id": entry_acc,
                    "taxon_id": entry_acc,
                    "uniprot_accession": entry_acc,
                    "uniprot_id": entry.get("uniprotId", ""),
                    "gene": entry.get("gene", ""),
                    "organism": entry.get("organismScientificName", "Unknown Organism"),
                    "description": entry.get("uniprotDescription", ""),
                    "mean_plddt": plddt,
                    "plddt_category": plddt_cat,
                    "sequence_length": entry.get("sequenceEnd", 0),
                    "structure_url": struct_url,
                    "pae_url": pae_url,
                    "structure_file": f"{entry_acc}{ext}"
                }
                metadata_list.append(meta_entry)

                if struct_url:
                    download_tasks.append((entry_acc, struct_url, out_struct_file))
                if pae_url and out_pae_file:
                    download_tasks.append((f"{entry_acc}_pae", pae_url, out_pae_file))

        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"  [!] UniProt ID '{uid}' not found in AlphaFold Database (404).")
            else:
                print(f"  [!] HTTP Error {e.code} querying AFDB for '{uid}': {e}")
        except Exception as e:
            print(f"  [!] Error querying AFDB for '{uid}': {e}")

    if not metadata_list:
        print("[!] No AlphaFold entries could be retrieved.")
        return []

    # Save taxa_metadata.json for tree annotation and interactive viewing
    meta_path = os.path.join(output_dir, "taxa_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata_list, f, indent=2)
    print(f"[AlphaFold DB] Saved metadata index for {len(metadata_list)} structures to '{meta_path}'.")

    downloaded = []
    print(f"[AlphaFold DB] Downloading {len(download_tasks)} files concurrently ({workers} workers)...")

    def _dl_file(task):
        name, url, out_path = task
        if os.path.isfile(out_path) and os.path.getsize(out_path) > 300:
            return name, out_path, True
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            with open(out_path, "wb") as f:
                f.write(data)
            return name, out_path, True
        except Exception as e:
            if os.path.isfile(out_path) and os.path.getsize(out_path) == 0:
                try:
                    os.remove(out_path)
                except OSError:
                    pass
            return name, str(e), False

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_dl_file, t) for t in download_tasks]
        completed = 0
        for fut in as_completed(futures):
            name, path_or_err, success = fut.result()
            completed += 1
            if success:
                if not name.endswith("_pae"):
                    downloaded.append(path_or_err)
            else:
                print(f"  [!] Failed to download {name}: {path_or_err}")
            if completed % 5 == 0 or completed == len(download_tasks):
                print(f"  Progress: {completed}/{len(download_tasks)} downloaded ({len(downloaded)} structures)")

    print(f"[AlphaFold DB] Download complete: {len(downloaded)} structure(s) saved to '{output_dir}'.\n")
    return downloaded
