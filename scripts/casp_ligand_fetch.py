#!/usr/bin/env python3
import argparse, os, re, tarfile, io, csv
from pathlib import Path
from urllib.parse import urljoin
import requests

BASE = "https://predictioncenter.org/"
S = requests.Session()
S.headers.update({"User-Agent": "CASP-downloader/1.0 (+for research use)"})


CASP16_PHARMA_DIR = "download_area/CASP16/targets/pharma_ligands/"
CASP16_RESULTS_LIG_DIR = "download_area/CASP16/results/ligands/"
CASP15_PRED_LIG_DIR = "download_area/CASP15/predictions/ligands/"
CASP15_TARGETLIST_LIG = "casp15/targetlist.cgi?view=ligand"

CASP16_SUPERTARGETS = ["L1000","L2000","L3000","L4000"]  # extend if CASP adds more pharma sets


def http_index_list(url):
    """Return list of hrefs from CASP Apache index pages."""
    r = S.get(url, timeout=60)
    r.raise_for_status()
    # Extremely simple: find href="..." entries
    hrefs = re.findall(r'href="([^"]+)"', r.text)
    # Filter out parent links and anchors
    hrefs = [h for h in hrefs if not h.startswith("?") and h not in ("../","./") and not h.startswith("#")]
    return sorted(set(hrefs))


def download(url, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with S.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1<<20):
                if chunk:
                    f.write(chunk)
    return dest


def fetch_casp16_pharma_smiles(outdir: Path):
    """Download L1000/L2000/L3000/L4000 SMILES tarballs and combine their TSVs."""
    index_url = urljoin(BASE, CASP16_PHARMA_DIR)
    entries = http_index_list(index_url)
    tars = [e for e in entries if e.upper().endswith(".SMILES.TAR.GZ")]
    if not tars:
        print("[CASP16] No SMILES tarballs found at", index_url)
        return None

    smiles_dir = outdir / "casp16" / "pharma_smiles"
    smiles_dir.mkdir(parents=True, exist_ok=True)

    combined_rows = []
    combined_header = None

    for tarname in tars:
        tar_url = urljoin(index_url, tarname)
        local_tar = smiles_dir / tarname
        if not local_tar.exists():
            print("[CASP16] downloading", tarname)
            download(tar_url, local_tar)
        else:
            print("[CASP16] exists", tarname)

        # read tar in-memory and pull out .tsv files
        with tarfile.open(local_tar, "r:gz") as tf:
            for m in tf.getmembers():
                if m.isfile() and m.name.lower().endswith(".tsv"):
                    f = tf.extractfile(m)
                    tsv_bytes = f.read()
                    text = tsv_bytes.decode("utf-8", errors="replace")
                    reader = csv.reader(io.StringIO(text), delimiter="\t")
                    header = next(reader)
                    if combined_header is None:
                        combined_header = ["supertarget","source_file"] + header
                    for row in reader:
                        combined_rows.append([tarname.split(".")[0], m.name] + row)

    if combined_header:
        out_tsv = smiles_dir / "CASP16_pharma_SMILES_combined.tsv"
        with open(out_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(combined_header)
            w.writerows(combined_rows)
        print("[CASP16] wrote", out_tsv)


def fetch_casp16_sequences(outdir: Path, supertargets=CASP16_SUPERTARGETS):
    """Grab FASTA for pharma supertargets via the plain-text endpoint."""
    seqdir = outdir / "casp16" / "fastas"
    seqdir.mkdir(parents=True, exist_ok=True)
    for tid in supertargets:
        url = urljoin(BASE, f"casp16/target.cgi?target={tid}&view=sequence")
        dest = seqdir / f"{tid}.fasta"
        if dest.exists():
            print("[CASP16] exists", dest.name)
            continue
        print("[CASP16] sequence", tid)
        r = S.get(url, timeout=60)
        r.raise_for_status()
        text = r.text.strip()
        # Ensure FASTA format has line breaks after header
        if "\n" not in text and "|" in text:
            hdr, seq = text.split("|", 1)
            text = hdr.strip() + "\n" + re.sub(r"\s+", "", seq)
        dest.write_text(text + "\n")


def fetch_casp16_results(outdir: Path):
    """Download pose summary CSVs (optional)."""
    index_url = urljoin(BASE, CASP16_RESULTS_LIG_DIR)
    entries = http_index_list(index_url)
    keep = [e for e in entries if e.endswith(".csv")]
    if not keep:
        print("[CASP16] No results CSVs found at", index_url)
        return
    resdir = outdir / "casp16" / "results"
    resdir.mkdir(parents=True, exist_ok=True)
    for fn in keep:
        url = urljoin(index_url, fn)
        dest = resdir / fn
        if not dest.exists():
            print("[CASP16] results", fn)
            download(url, dest)
        else:
            print("[CASP16] exists", fn)


def fetch_casp15_ligand_fastas(outdir: Path):
    """Iterate the CASP15 ligand target list and fetch FASTA per target via plain-text endpoint."""
    # Pull the ligand view and regex target IDs like T####, H####, R####
    url = urljoin(BASE, CASP15_TARGETLIST_LIG)
    print("[CASP15] loading ligand target list …")
    r = S.get(url, timeout=90)
    r.raise_for_status()
    # Find hrefs like >T1152< etc.
    tids = sorted(set(re.findall(r">([THR]\d{4,5}v?\d*)<", r.text)))
    if not tids:
        # Fallback: also parse explicit "Target:" results page if needed
        print("[CASP15] Warning: no targets parsed from ligand list; page layout may have changed.")
    seqdir = outdir / "casp15" / "fastas"
    seqdir.mkdir(parents=True, exist_ok=True)
    for tid in tids:
        endpoint = urljoin(BASE, f"casp15/target.cgi?target={tid}&view=sequence")
        dest = seqdir / f"{tid}.fasta"
        if dest.exists():
            print("[CASP15] exists", dest.name)
            continue
        try:
            print("[CASP15] sequence", tid)
            rr = S.get(endpoint, timeout=60)
            rr.raise_for_status()
            text = rr.text.strip()
            if "\n" not in text and "|" in text:
                hdr, seq = text.split("|", 1)
                text = hdr.strip() + "\n" + re.sub(r"\s+", "", seq)
            dest.write_text(text + "\n")
        except Exception as e:
            print(f"[CASP15] failed {tid}: {e}")


def fetch_casp15_predictions(outdir: Path):
    """Bulk download all CASP15 ligand prediction tarballs."""
    index_url = urljoin(BASE, CASP15_PRED_LIG_DIR)
    entries = http_index_list(index_url)
    tars = [e for e in entries if e.lower().endswith(".tar.gz")]
    pred_dir = outdir / "casp15" / "predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)
    for t in tars:
        url = urljoin(index_url, t)
        dest = pred_dir / t
        if not dest.exists():
            print("[CASP15] downloading", t)
            download(url, dest)
        else:
            print("[CASP15] exists", t)


def main():
    ap = argparse.ArgumentParser(description="Fetch CASP ligand-docking data (CASP16/CASP15).")
    ap.add_argument("--out", default="./casp_data", help="Output directory")
    ap.add_argument("--skip-casp16-smiles", action="store_true", help="Skip CASP16 pharma SMILES tarballs")
    ap.add_argument("--skip-casp16-fastas", action="store_true", help="Skip CASP16 receptor FASTAs")
    ap.add_argument("--with-casp16-results", action="store_true", help="Also download CASP16 ligand result CSVs")
    ap.add_argument("--skip-casp15-fastas", action="store_true", help="Skip CASP15 ligand FASTAs")
    ap.add_argument("--with-casp15-predictions", action="store_true", help="Also download CASP15 ligand prediction tarballs")

    args = ap.parse_args()
    outdir = Path(args.out)

    if not args.skip_casp16_smiles:
        fetch_casp16_pharma_smiles(outdir)
    if not args.skip_casp16_fastas:
        fetch_casp16_sequences(outdir)
    if args.with_casp16_results:
        fetch_casp16_results(outdir)

    if not args.skip_casp15_fastas:
        fetch_casp15_ligand_fastas(outdir)
    if args.with_casp15_predictions:
        fetch_casp15_predictions(outdir)

    print("\nDone. Tree is rooted at:", outdir.resolve())


if __name__ == "__main__":
    main()

