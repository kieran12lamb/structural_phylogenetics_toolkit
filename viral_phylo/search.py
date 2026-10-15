"""Build a pipeline cohort from a query protein's structural neighbours.

Submits a structure to the hosted Foldseek server, pools the hits across every
database searched, filters them on hit quality, and materialises the survivors as
a structure directory the rest of the pipeline can consume unchanged.

Two design points are worth stating because they are not obvious:

*Filter before capping.* ``--max-hits`` is a ceiling, never a target. Every
database returns rows down to its own E-value floor, so ranking the pooled set and
taking the top N pads a thin result with the random tail: a SARS-CoV-2 ORF8 query
returns 157 genuine hits from BFVD plus 21 hits from CATH50 and PDB100 at E-values
of 3 to 5. Those 21 are noise, and a target of 500 would have pulled them in. What
passes the gate is used; what does not is dropped however short that leaves the
cohort.

*Confidence is not one number.* A predicted structure stores pLDDT in the B-factor
column (0-100, higher is better); an experimental one stores atomic displacement in
the same column (Angstrom squared, higher is worse), and the two ranges overlap, so
reading that column across a pooled cohort silently inverts half the values. Each
taxon therefore carries its confidence *kind* alongside the value, and per-residue
LDDT against the query - which means the same thing whatever the provenance - is
recorded as the scale that is safe to compare.
"""

import glob
import json
import os
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Tuple

SEARCH_HOST = "https://search.foldseek.com"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry"
BFVD_BASE = "https://bfvd.steineggerlab.workers.dev/2025_03_v3"
USER_AGENT = "viral-phylo/1.0 (structural phylogenetics)"

# Databases whose entries are predictions, so their B-factor column is pLDDT.
PREDICTION_DBS = {
    "bfvd_v3", "BFVD_folddisco", "afdb50", "afdb50_folddisco", "afdb-swissprot",
    "afdb-proteome", "afdb-proteome_folddisco", "esm30_folddisco", "mgnify_esm30",
    "cath50", "bfmd",
}
DEFAULT_DATABASES = ["bfvd_v3"]


def _get(url: str, timeout: int = 60, headers: Optional[Dict[str, str]] = None) -> bytes:
    hdrs = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, headers=hdrs)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _get_json(url: str, timeout: int = 60) -> Any:
    return json.loads(_get(url, timeout=timeout).decode("utf-8"))


# ---------------------------------------------------------------------------
# Query resolution
# ---------------------------------------------------------------------------

def resolve_query(spec: str, work_dir: str) -> Dict[str, Any]:
    """Resolve a query to a structure file.

    ``spec`` may be a path to a local structure, a UniProt accession, or an NCBI
    protein accession. The hosted search rejects a bare sequence - submitting FASTA
    returns ``status: ERROR`` - so a sequence-only query has no path here yet and
    is reported as such rather than failing obscurely later.
    """
    os.makedirs(work_dir, exist_ok=True)

    if os.path.isfile(spec):
        return {"kind": "file", "path": spec, "id": os.path.splitext(os.path.basename(spec))[0]}

    acc = spec.strip()
    if re.fullmatch(r"[A-Z]{2}_[0-9]+\.[0-9]+", acc) or re.fullmatch(r"[A-Z]{3}[0-9]{5}\.[0-9]+", acc):
        uni = _refseq_to_uniprot(acc)
        if not uni:
            raise ValueError(
                f"'{acc}' looks like an NCBI accession but no UniProt mapping was found. "
                "Pass a UniProt accession or a structure file instead.")
        print(f"[Search] {acc} -> UniProt {uni}")
        acc = uni

    path, meta = _fetch_afdb_structure(acc, work_dir)
    return {"kind": "afdb", "path": path, "id": acc, **meta}


def _refseq_to_uniprot(accession: str) -> Optional[str]:
    """Map an NCBI protein accession to UniProt via the UniProt ID-mapping search."""
    q = urllib.parse.urlencode({
        "query": f'xref:refseq-{accession} OR "{accession}"',
        "fields": "accession,length",
        "format": "json",
        "size": "5",
    })
    try:
        data = _get_json(f"https://rest.uniprot.org/uniprotkb/search?{q}", timeout=45)
    except (urllib.error.URLError, ValueError, TimeoutError):
        return None
    results = data.get("results") or []
    return results[0].get("primaryAccession") if results else None


def _fetch_afdb_structure(uniprot: str, work_dir: str) -> Tuple[str, Dict[str, Any]]:
    entries = _get_json(f"{AFDB_API}/{uniprot}", timeout=45)
    if not entries:
        raise ValueError(f"No AlphaFold prediction for '{uniprot}'.")
    entry = next((e for e in entries if e.get("uniprotAccession") == uniprot), entries[0])
    cif_url = entry.get("cifUrl")
    if not cif_url:
        raise ValueError(f"AlphaFold entry for '{uniprot}' has no cifUrl.")
    path = os.path.join(work_dir, f"{uniprot}.cif")
    with open(path, "wb") as fh:
        fh.write(_get(cif_url, timeout=120))
    return path, {
        "plddt": entry.get("globalMetricValue"),
        "description": entry.get("uniprotDescription"),
        "length": entry.get("uniprotEnd"),
    }


# ---------------------------------------------------------------------------
# Hosted Foldseek search
# ---------------------------------------------------------------------------

def list_databases() -> List[Dict[str, Any]]:
    return (_get_json(f"{SEARCH_HOST}/api/databases", timeout=45) or {}).get("databases", [])


def submit_search(structure_path: str, databases: List[str], mode: str = "3diaa") -> str:
    """POST a structure and return the ticket id.

    ``q`` has to be a multipart file part; sending it as a plain form field is
    rejected with ``400 http: no such file``.
    """
    boundary = "----viralphylo" + str(int(time.time() * 1000))
    with open(structure_path, "rb") as fh:
        payload = fh.read()

    parts: List[bytes] = []

    def field(name: str, value: str) -> None:
        parts.append(
            f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n'.encode())

    parts.append(
        f'--{boundary}\r\nContent-Disposition: form-data; name="q"; '
        f'filename="{os.path.basename(structure_path)}"\r\n'
        f"Content-Type: application/octet-stream\r\n\r\n".encode())
    parts.append(payload)
    parts.append(b"\r\n")
    field("mode", mode)
    for db in databases:
        field("database[]", db)
    parts.append(f"--{boundary}--\r\n".encode())

    body = b"".join(parts)
    req = urllib.request.Request(
        f"{SEARCH_HOST}/api/ticket", data=body,
        headers={"User-Agent": USER_AGENT,
                 "Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        out = json.loads(resp.read().decode())
    if not out.get("id"):
        raise RuntimeError(f"Foldseek server returned no ticket: {out}")
    return out["id"]


def wait_for_ticket(ticket: str, timeout_s: int = 1800, interval_s: int = 8) -> str:
    """Poll until the job leaves the queue. Returns the terminal status."""
    deadline = time.time() + timeout_s
    last = None
    while time.time() < deadline:
        status = (_get_json(f"{SEARCH_HOST}/api/ticket/{ticket}", timeout=45) or {}).get("status")
        if status != last:
            print(f"[Search] ticket {ticket[:10]}... {status}")
            last = status
        if status in ("COMPLETE", "ERROR"):
            return status
        time.sleep(interval_s)
    raise TimeoutError(f"Foldseek ticket {ticket} did not finish within {timeout_s}s (last status {last}).")


def fetch_results(ticket: str) -> List[Dict[str, Any]]:
    """Return one block per database searched, each with its alignments."""
    data = _get_json(f"{SEARCH_HOST}/api/result/{ticket}/0", timeout=300) or {}
    blocks = []
    for block in data.get("results", []):
        alns = block.get("alignments") or []
        if alns and isinstance(alns[0], list):      # nested one level on some builds
            alns = alns[0]
        blocks.append({"db": block.get("db"), "alignments": alns})
    return blocks


# ---------------------------------------------------------------------------
# Pooling and filtering
# ---------------------------------------------------------------------------

def hit_accession(target: str, db: str) -> str:
    """Reduce a database's target string to a stable identifier.

    The three schemes in play look nothing like each other:
      bfvd_v3  'A0A7M2QE25 ORF8 protein'
      cath50   'af_A0A2R8QEK5_306_422_2.60.120.290'
      pdb100   '2hey-assembly1.cif.gz_F Crystal structure of ...'
    """
    head = str(target).split()[0]
    if db and db.startswith("pdb"):
        m = re.match(r"([0-9a-zA-Z]{4})[-_].*?(?:_([A-Za-z0-9]+))?$", head)
        if m:
            return (m.group(1) + ("_" + m.group(2) if m.group(2) else "")).upper()
    if head.startswith("AF-"):
        m = re.match(r"AF-([A-Z0-9]+)-F\d+", head)
        if m:
            return m.group(1)
    if head.startswith("af_"):
        # CATH domains name an AFDB parent plus the residue window they cover,
        # e.g. af_A0A2R8QEK5_306_422_2.60.120.290. Keep the window: the domain is
        # a slice of the parent, so dropping it would collide distinct domains of
        # the same protein onto one identifier.
        m = re.match(r"af_([A-Z0-9]+)_(\d+)_(\d+)_", head)
        if m:
            return f"{m.group(1)}_{m.group(2)}-{m.group(3)}"
    return head


def cath_parent(target: str) -> Optional[Tuple[str, int, int]]:
    """Parent AFDB accession and residue window for a CATH domain id, if it is one."""
    m = re.match(r"af_([A-Z0-9]+)_(\d+)_(\d+)_", str(target).split()[0])
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def pool_and_filter(
    blocks: List[Dict[str, Any]],
    max_hits: int = 500,
    max_evalue: float = 1e-3,
    min_alnlen: int = 0,
    min_coverage: float = 0.0,
    include_self: bool = True,
) -> Dict[str, Any]:
    """Merge every database's hits, drop what fails the gate, cap at ``max_hits``."""
    pooled: List[Dict[str, Any]] = []
    per_db: Dict[str, Dict[str, int]] = {}

    for block in blocks:
        db = block.get("db")
        stats = per_db.setdefault(db, {"returned": 0, "passed": 0})
        for a in block.get("alignments") or []:
            stats["returned"] += 1
            ev = a.get("eval")
            aln = a.get("alnLength") or 0
            qlen = a.get("qLen") or 0
            cov = (aln / qlen) if qlen else 0.0
            if ev is None or ev > max_evalue:
                continue
            if aln < min_alnlen or cov < min_coverage:
                continue
            stats["passed"] += 1
            pooled.append({
                "accession": hit_accession(a.get("target", ""), db),
                "target": a.get("target"),
                "db": db,
                "evalue": ev,
                "prob": a.get("prob"),
                "seq_id": a.get("seqId"),
                "aln_length": aln,
                "coverage": round(cov, 4),
                "db_len": a.get("dbLen"),
                "tax_id": a.get("taxId"),
                "tax_name": a.get("taxName"),
                "t_seq": a.get("tSeq"),
                "t_ca": a.get("tCa"),
                "predicted": db in PREDICTION_DBS,
            })

    # Same entry can surface from several databases; keep the strongest.
    best: Dict[str, Dict[str, Any]] = {}
    for h in pooled:
        prev = best.get(h["accession"])
        if prev is None or h["evalue"] < prev["evalue"]:
            best[h["accession"]] = h
    hits = sorted(best.values(), key=lambda h: h["evalue"])

    if not include_self:
        hits = [h for h in hits if (h.get("seq_id") or 0) < 99.9]

    capped = hits[:max_hits] if max_hits else hits
    return {
        "hits": capped,
        "pooled": len(pooled),
        "unique": len(hits),
        "selected": len(capped),
        "requested": max_hits,
        "short_of_request": bool(max_hits and len(capped) < max_hits),
        "per_db": per_db,
        "filters": {"max_evalue": max_evalue, "min_alnlen": min_alnlen,
                    "min_coverage": min_coverage, "include_self": include_self},
    }


# ---------------------------------------------------------------------------
# Structure sourcing
#
# tCa gives a complete Calpha trace for every hit, so geometry needs no downloads
# at all. Confidence does: pLDDT lives only in the source file. Each database is
# therefore fetched from its own distribution where an adapter exists, and falls
# back to tCa (geometry only, confidence unavailable) where one does not.
# ---------------------------------------------------------------------------

class BFVDIndex:
    """Random access into BFVD's 220 GB tar via its published sorted index.

    BFVD exposes no per-entry URL, but it serves byte ranges and ships a sorted
    ``accession -> offset -> length`` index. Downloading the whole index is 162 MB;
    binary-searching it over range requests costs about eight small reads instead,
    which is cheaper for the handful of accessions a cohort needs.
    """

    def __init__(self, cache_dir: Optional[str] = None):
        self.url = f"{BFVD_BASE}/bfvd_indexed.tar.index"
        self.tar = f"{BFVD_BASE}/bfvd_indexed.tar"
        self.size = None
        self.cache_path = os.path.join(cache_dir, "bfvd_indexed.tar.index") if cache_dir else None
        self._local = self.cache_path if (self.cache_path and os.path.isfile(self.cache_path)) else None
        self._mem: Dict[str, Tuple[int, int]] = {}

    def _index_size(self) -> int:
        if self.size is None:
            req = urllib.request.Request(self.url, method="HEAD", headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as resp:
                self.size = int(resp.headers.get("Content-Length") or 0)
        return self.size

    def _slice(self, start: int, end: int) -> str:
        data = _get(self.url, timeout=90, headers={"Range": f"bytes={start}-{end}"})
        return data.decode("utf-8", errors="replace")

    def lookup(self, accession: str) -> Optional[Tuple[int, int]]:
        if accession in self._mem:
            return self._mem[accession]
        if self._local:
            with open(self._local, "r") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f and f[0] == accession:
                        self._mem[accession] = (int(f[1]), int(f[2]))
                        return self._mem[accession]
            return None

        lo, hi = 0, self._index_size()
        for _ in range(40):
            if hi - lo <= 400_000:
                break
            mid = (lo + hi) // 2
            chunk = self._slice(mid, mid + 400)
            lines = chunk.split("\n")
            key = lines[1].split("\t")[0] if len(lines) > 1 and lines[1] else None
            if not key:
                break
            if key < accession:
                lo = mid
            else:
                hi = mid
        window = self._slice(lo, min(hi + 4000, self._index_size()))
        for line in window.split("\n"):
            f = line.split("\t")
            if len(f) == 3 and f[0] == accession:
                self._mem[accession] = (int(f[1]), int(f[2]))
                return self._mem[accession]
        return None

    def fetch(self, accession: str, out_path: str) -> Optional[str]:
        import gzip
        loc = self.lookup(accession)
        if not loc:
            return None
        offset, length = loc
        blob = _get(self.tar, timeout=180,
                    headers={"Range": f"bytes={offset}-{offset + length - 1}"})
        try:
            text = gzip.decompress(blob)
        except OSError:
            return None
        with open(out_path, "wb") as fh:
            fh.write(text)
        return out_path


def write_ca_pdb(path: str, t_ca: str, t_seq: str, bfactor: float = 0.0) -> Optional[str]:
    """Write a Calpha-only PDB from the coordinates the search already returned.

    Used where no source adapter exists. The B-factor is written as a constant
    sentinel rather than a plausible-looking number, because the API supplies no
    per-residue confidence and inventing one is how the two incompatible scales get
    confused in the first place.
    """
    nums = [x for x in str(t_ca).split(",") if x.strip()]
    if len(nums) < 3:
        return None
    coords = [tuple(float(v) for v in nums[i:i + 3]) for i in range(0, len(nums) - 2, 3)]
    seq = str(t_seq or "")
    three = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
             "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
             "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
             "Y": "TYR", "V": "VAL"}
    lines = []
    for i, (x, y, z) in enumerate(coords, start=1):
        aa = three.get(seq[i - 1].upper(), "GLY") if i <= len(seq) else "GLY"
        lines.append(
            f"ATOM  {i:>5}  CA  {aa} A{i:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00{bfactor:>6.2f}           C")
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def pdb_confidence(path: str) -> Dict[str, Any]:
    """Per-residue B-factor values from a structure file, with no interpretation."""
    vals: List[float] = []
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
                    try:
                        vals.append(float(line[60:66]))
                    except ValueError:
                        pass
    except OSError:
        return {"values": [], "median": None}
    med = sorted(vals)[len(vals) // 2] if vals else None
    return {"values": vals, "median": med}


def rcsb_metadata(pdb_id: str) -> Dict[str, Any]:
    """Method and resolution for an experimental entry. NMR entries have no resolution."""
    try:
        d = _get_json(f"{RCSB_ENTRY}/{pdb_id[:4].upper()}", timeout=45)
    except (urllib.error.URLError, ValueError, TimeoutError):
        return {}
    res = ((d.get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0]
    method = ((d.get("exptl") or [{}])[0]).get("method")
    return {"method": method, "resolution": res}


def materialise_hits(
    hits: List[Dict[str, Any]],
    out_dir: str,
    cache_dir: Optional[str] = None,
    workers: int = 12,
) -> Dict[str, Any]:
    """Write one structure per hit, recording where its confidence came from.

    Source files are preferred because they carry per-residue confidence; tCa is
    the fallback and yields geometry only. Every taxon records ``confidence_kind``
    so nothing downstream has to guess which scale a B-factor column is on.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    os.makedirs(out_dir, exist_ok=True)
    bfvd = BFVDIndex(cache_dir=cache_dir)
    records: Dict[str, Dict[str, Any]] = {}

    def one(hit: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
        acc = hit["accession"]
        db = hit.get("db") or ""
        dest = os.path.join(out_dir, f"{acc}.pdb")
        rec = {
            "accession": acc, "db": db, "evalue": hit.get("evalue"),
            "seq_id": hit.get("seq_id"), "aln_length": hit.get("aln_length"),
            "coverage": hit.get("coverage"), "tax_id": hit.get("tax_id"),
            "tax_name": hit.get("tax_name"), "db_len": hit.get("db_len"),
            "structure_kind": "predicted" if hit.get("predicted") else "experimental",
            "confidence_kind": "none", "source": None, "resolution": None,
        }
        try:
            if db in ("bfvd_v3", "BFVD_folddisco") and bfvd.fetch(acc, dest):
                rec.update(source="bfvd", confidence_kind="plddt")
            elif db.startswith("afdb") and _try_afdb(acc, dest):
                rec.update(source="afdb", confidence_kind="plddt")
            else:
                if not write_ca_pdb(dest, hit.get("t_ca") or "", hit.get("t_seq") or ""):
                    return acc, {**rec, "error": "no coordinates"}
                rec.update(source="tCa", confidence_kind="none")
                if db.startswith("pdb"):
                    rec.update(rcsb_metadata(acc))
            conf = pdb_confidence(dest)
            rec["confidence"] = conf["values"] if rec["confidence_kind"] != "none" else []
            rec["confidence_median"] = conf["median"] if rec["confidence_kind"] != "none" else None
            rec["n_residues"] = len(conf["values"])
            rec["path"] = dest
        except Exception as exc:                       # one bad hit must not sink the cohort
            rec["error"] = f"{type(exc).__name__}: {exc}"
        return acc, rec

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(one, h) for h in hits]
        for i, fut in enumerate(as_completed(futures), 1):
            acc, rec = fut.result()
            records[acc] = rec
            if i % 25 == 0 or i == len(hits):
                print(f"[Search] materialised {i}/{len(hits)}")

    ok = [r for r in records.values() if not r.get("error")]
    by_source: Dict[str, int] = {}
    for r in ok:
        by_source[r.get("source") or "?"] = by_source.get(r.get("source") or "?", 0) + 1
    return {
        "dir": out_dir,
        "records": records,
        "written": len(ok),
        "failed": [r for r in records.values() if r.get("error")],
        "by_source": by_source,
        "with_confidence": sum(1 for r in ok if r.get("confidence_kind") != "none"),
    }


def _try_afdb(accession: str, dest: str) -> bool:
    try:
        entries = _get_json(f"{AFDB_API}/{accession}", timeout=45)
    except (urllib.error.URLError, ValueError, TimeoutError):
        return False
    if not entries:
        return False
    url = entries[0].get("pdbUrl") or entries[0].get("cifUrl")
    if not url:
        return False
    try:
        with open(dest, "wb") as fh:
            fh.write(_get(url, timeout=120))
        return True
    except (urllib.error.URLError, OSError, TimeoutError):
        return False


def per_residue_lddt(query_structure: str, cohort_dir: str, work_dir: str,
                     threads: int = 8) -> Dict[str, Dict[str, Any]]:
    """Per-residue LDDT of every cohort member against the query.

    This is the one confidence-like quantity that means the same thing for a
    prediction and a crystal structure, so it is what the viewer can safely colour
    by across a mixed cohort. It is similarity to the query, not model quality:
    a flawless structure of a divergent homolog scores low here and is still right.
    """
    from viral_phylo.binaries import find_foldseek_bin
    os.makedirs(work_dir, exist_ok=True)
    out = os.path.join(work_dir, "lddt.m8")
    tmp = os.path.join(work_dir, "tmp_lddt")
    cmd = [find_foldseek_bin(), "easy-search", query_structure, cohort_dir, out, tmp,
           "--exhaustive-search", "1", "--alignment-type", "1", "--exact-tmscore", "1",
           "-e", "inf",
           # lddtfull covers the ALIGNED span, not the whole chain, so the target
           # start/end are needed to place those values on residues. Without them a
           # 76-value array silently looks like a 121-residue track.
           "--format-output", "query,target,alntmscore,lddt,lddtfull,tstart,tend,qstart,qend,tlen",
           "--threads", str(threads)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    scores: Dict[str, Dict[str, Any]] = {}
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            tgt = os.path.splitext(os.path.basename(f[1]))[0]
            try:
                vals = [float(v) for v in f[4].split(",") if v.strip()]
                t_start, t_end = int(f[5]), int(f[6])
                db_len = int(f[9]) if len(f) > 9 and f[9].strip() else None
                # Expand onto the full chain so the array indexes by residue, with
                # None outside the aligned span rather than a misleading 0.
                track: List[Optional[float]] = [None] * (db_len or (t_end))
                for i, v in enumerate(vals):
                    pos = t_start - 1 + i
                    if 0 <= pos < len(track):
                        track[pos] = v
                scores[tgt] = {
                    "tm_to_query": float(f[2]),
                    "lddt_to_query": float(f[3]),
                    "lddt_aligned": vals,
                    "lddt_per_residue": track,
                    "aligned_span": [t_start, t_end],
                    "query_span": [int(f[7]), int(f[8])] if len(f) > 8 else None,
                }
            except (ValueError, IndexError):
                continue
    return scores


# ---------------------------------------------------------------------------
# Taxonomy and host enrichment
# ---------------------------------------------------------------------------

UNIPROT_TAX = "https://rest.uniprot.org/taxonomy"

# Ranks worth promoting to their own metadata column. Viral taxonomy is sparse -
# most entries fill only a few of these - so each is emitted only if something in
# the cohort actually carries it.
TAX_RANKS = ("superkingdom", "kingdom", "phylum", "class", "order",
             "suborder", "family", "subfamily", "genus", "subgenus", "species")


def fetch_taxonomy(tax_ids: List[int], batch: int = 50) -> Dict[int, Dict[str, Any]]:
    """Lineage and virus hosts per taxon id, from UniProt's taxonomy endpoint.

    Keyed on taxon rather than protein on purpose. BFVD ships structures for
    UniProt entries that UniProt has since deleted - A0A7M2QE25 returns
    ``entryType: Inactive, DELETED: Not part of a reference proteome`` - so a
    per-protein lookup silently loses most of a cohort. The taxon id the search
    already returns is stable, and a cohort has far fewer taxa than members:
    157 structures span 20 taxa, 397 span 119.
    """
    out: Dict[int, Dict[str, Any]] = {}
    ids = sorted({int(t) for t in tax_ids if t})
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        q = urllib.parse.urlencode({
            "query": " OR ".join(f"tax_id:{t}" for t in chunk),
            "fields": "id,scientific_name,common_name,rank,hosts,lineage",
            "format": "json",
            "size": str(len(chunk)),
        })
        try:
            data = _get_json(f"{UNIPROT_TAX}/search?{q}", timeout=90)
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            print(f"[Taxonomy] [!] Notice: lookup failed for a batch ({exc}); "
                  "those taxa keep the name the search returned.")
            continue
        for r in data.get("results") or []:
            tid = r.get("taxonId")
            if not tid:
                continue
            lineage = r.get("lineage") or []
            ranks = {}
            for node in lineage:
                rk = node.get("rank")
                if rk and rk in TAX_RANKS and rk not in ranks:
                    ranks[rk] = node.get("scientificName")
            hosts = r.get("hosts") or []
            out[int(tid)] = {
                "scientific_name": r.get("scientificName"),
                "common_name": r.get("commonName"),
                "rank": r.get("rank"),
                "lineage": [n.get("scientificName") for n in lineage if n.get("scientificName")],
                "lineage_str": " > ".join(
                    n.get("scientificName", "") for n in reversed(lineage) if n.get("scientificName")),
                "ranks": ranks,
                "hosts": [{"name": h.get("scientificName"), "common": h.get("commonName"),
                           "tax_id": h.get("taxonId")} for h in hosts],
                "host_names": [h.get("scientificName") for h in hosts if h.get("scientificName")],
                "host_common": [h.get("commonName") for h in hosts if h.get("commonName")],
            }
        print(f"[Taxonomy] resolved {len(out)}/{len(ids)} taxa")
    return out


def enrich_records(records: Dict[str, Dict[str, Any]],
                   taxonomy: Dict[int, Dict[str, Any]]) -> None:
    """Fold taxonomy and host annotation into per-structure records, in place."""
    for rec in records.values():
        tid = rec.get("tax_id")
        tax = taxonomy.get(int(tid)) if tid else None
        if not tax:
            continue
        rec["lineage"] = tax["lineage_str"]
        rec["tax_rank"] = tax.get("rank")
        for rk, val in (tax.get("ranks") or {}).items():
            rec[f"tax_{rk}"] = val
        hosts = tax.get("host_names") or []
        # A virus can have several documented hosts - SARS-CoV lists both masked
        # palm civet and human - so keep the full list and a single primary for
        # colouring, rather than silently dropping the rest.
        rec["hosts"] = "; ".join(hosts) if hosts else None
        rec["host_primary"] = hosts[0] if hosts else None
        rec["host_count"] = len(hosts)
        common = tax.get("host_common") or []
        rec["host_common"] = "; ".join(c for c in common if c) or None
