#!/usr/bin/env python3
import csv, os, re, tarfile, urllib.request, urllib.error

BASE = "https://predictioncenter.org"
CASP = "casp15"
OUT = "casp15_original_inputs"
os.makedirs(OUT, exist_ok=True)

TARGETLIST_CSV = f"{BASE}/{CASP}/targetlist.cgi?type=csv"
LIGANDS_TARBALL = f"{BASE}/download_area/CASP15/targets/casp15.targets.ligands.ALL_09.18.2025.tar.gz"

def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()

def save(path, data: bytes):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(data)

print("==> Downloading CASP15 target list CSV…")
csv_bytes = fetch(TARGETLIST_CSV)
csv_text = csv_bytes.decode("utf-8", errors="replace")

# Write the raw CSV too (for provenance)
save(os.path.join(OUT, "casp15_targetlist.csv"), csv_bytes)

# Parse CSV
targets = []
reader = csv.DictReader(csv_text.splitlines(), delimiter=';')
for row in reader:
    # The CSV columns include "Target" and "Type" among others.
    tid = row.get("Target", "").strip()  # e.g., T1124, H1135, R1136, etc.
    ttype = row.get("Type", "").strip()
    if not tid:
        continue
    targets.append({"id": tid, "type": ttype})

print(f"==> Found {len(targets)} targets in CSV.")

# 1) FASTA for every target: plain-text sequence endpoint
#    Pattern: https://predictioncenter.org/casp15/target.cgi?target=TXXXX&view=sequence
#    Example verified: T1152 sequence plain-text.  (We’ll save each as {id}.fasta)
for t in targets:
    tid = t["id"]
    url = f"{BASE}/{CASP}/target.cgi?target={tid}&view=sequence"
    try:
        data = fetch(url)
        # Sanity check: expect a FASTA header starting with '>'
        if not data.startswith(b">"):
            raise RuntimeError(f"Unexpected content for {tid}")
        save(os.path.join(OUT, "fasta", f"{tid}.fasta"), data)
        print(f"[FASTA] {tid} ✓")
    except Exception as e:
        print(f"[FASTA] {tid} ✗  ({e})")

# 2) SMILES for ligand targets: use the official bulk tarball
lig_dir = os.path.join(OUT, "ligands_tarball")
os.makedirs(lig_dir, exist_ok=True)

print("==> Downloading CASP15 ligands tarball (SMILES)…")
tb_bytes = fetch(LIGANDS_TARBALL)
tb_path = os.path.join(lig_dir, os.path.basename(LIGANDS_TARBALL))
save(tb_path, tb_bytes)

print("==> Extracting ligands tarball…")
with tarfile.open(tb_path, "r:gz") as tar:
    # Extract safely
    def is_within_directory(directory, target):
        abs_directory = os.path.abspath(directory)
        abs_target = os.path.abspath(target)
        return os.path.commonpath([abs_directory]) == os.path.commonpath([abs_directory, abs_target])
    for m in tar.getmembers():
        # optional: only extract .smi/.sdf/.txt that look like SMILES/ligand lists
        if any(m.name.lower().endswith(ext) for ext in (".smi", ".sdf", ".txt", ".tsv", ".csv")):
            out_path = os.path.join(lig_dir, m.name)
            if not is_within_directory(lig_dir, out_path):
                continue
            tar.extract(m, path=lig_dir)
print("==> Done.\n")
print(f"Outputs:\n  - FASTA: {os.path.join(OUT, 'fasta')}/*.fasta\n  - SMILES (from tarball): {lig_dir}/")

