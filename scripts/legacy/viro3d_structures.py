import argparse
import os
import requests

BASE_URL = "https://viro3d.cvr.gla.ac.uk/api"


def main():
    parser = argparse.ArgumentParser(
        description="Search and download protein structures from Viro3D by qualifier."
    )
    parser.add_argument(
        "qualifier",
        nargs="?",
        default=None,
        help="Protein qualifier to search for (default: glycoprotein)",
    )
    parser.add_argument(
        "-q",
        "--qualifier",
        dest="qualifier_flag",
        default=None,
        help="Protein qualifier to search for (optional flag format)",
    )
    parser.add_argument(
        "-m",
        "--max-sequences",
        "--max-downloads",
        dest="max_sequences",
        type=int,
        default=1,
        help="Maximum number of structures/sequences to download (default: 1)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="viro_3d_structures",
        help="Output directory for downloaded structures (default: viro_3d_structures)",
    )
    args = parser.parse_args()
    qualifier = args.qualifier_flag or args.qualifier or "glycoprotein"
    max_sequences = args.max_sequences
    output_dir = args.output_dir

    if max_sequences <= 0:
        print("Maximum number of sequences must be greater than 0.")
        return

    # 1. Search for matching proteins
    query_url = f"{BASE_URL}/proteins/protein_name/"
    response = requests.get(query_url, params={"qualifier": qualifier, "page_size": max_sequences})
    response.raise_for_status()
    data = response.json()

    protein_structures = data.get("protein_structures", [])
    if not protein_structures:
        print(f"No protein structures found for qualifier: {qualifier}")
        return

    os.makedirs(output_dir, exist_ok=True)
    items_to_download = protein_structures[:max_sequences]
    total_found = len(items_to_download)
    print(f"Found {len(protein_structures)} structure(s) (downloading up to {total_found}) for qualifier '{qualifier}':\n")

    for idx, item in enumerate(items_to_download, 1):
        record_id = item["record_id"]
        product = item.get("product", qualifier)
        plddt = item.get("colabfold_json_pLDDT")
        print(f"[{idx}/{total_found}] Found: {product} ({record_id})")
        print(f"ColabFold pLDDT: {plddt}")

        # 2. Construct structure download URL
        pdb_url = f"{BASE_URL}/pdb/CF-{record_id}_relaxed.pdb"
        print(f"Downloading structure from: {pdb_url}")

        # 3. Download the PDB file
        try:
            pdb_res = requests.get(pdb_url)
            pdb_res.raise_for_status()
            output_filename = os.path.join(output_dir, f"{record_id}.pdb")
            with open(output_filename, "w", encoding="utf-8") as f:
                f.write(pdb_res.text)
            print(f"Saved structure to {output_filename} ({len(pdb_res.text)} bytes)\n")
        except requests.RequestException as e:
            print(f"Failed to download structure for {record_id}: {e}\n")


if __name__ == "__main__":
    main()