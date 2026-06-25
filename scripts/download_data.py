#!/usr/bin/env python3
"""Download RNA-seq FASTQ data from Alibaba Cloud OSS.

Reads an Excel sample sheet to extract sample names, lists files in OSS,
matches R1/R2 pairs by sample name, and downloads them to the local
RawData/ directory for consumption by run_upstream.sh / run_pipeline.sh.

Usage:
    export OSS_ACCESS_KEY_ID="your-key-id"
    export OSS_ACCESS_KEY_SECRET="your-key-secret"
    python3 scripts/download_data.py \
        --excel ref/dataset/sample_info.xlsx \
        --project-id FQ260601941

    # Preview without downloading:
    python3 scripts/download_data.py \
        --excel ref/dataset/sample_info.xlsx \
        --project-id FQ260601941 \
        --dry-run
"""

import argparse
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

# --- Constants ---
OSSUTIL = "/home/zhuzp/bulkseq/ossutil"
DEFAULT_OSS_PATH = "oss://jmoss2026-5221/20260620_E260620003_U5228_FQ260601941"
DEFAULT_ENDPOINT = "http://oss-cn-hangzhou.aliyuncs.com"
DEFAULT_BASE_DIR = "/home/zhuzp/bulkseq"

# R1/R2 filename patterns to try when matching OSS files to sample names.
# The first matching pattern wins. Downloaded files are always renamed to
# {sample}_1.fq.gz / {sample}_2.fq.gz regardless of the matched pattern.
PAIRED_PATTERNS = [
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_R1.fq.gz",    "_R2.fq.gz"),
    ("_1.fq.gz",     "_2.fq.gz"),
    ("_1.fastq.gz",  "_2.fastq.gz"),
    ("_1.clean.fq.gz", "_2.clean.fq.gz"),
    ("_L001_R1_001.fastq.gz", "_L001_R2_001.fastq.gz"),
]


class DownloadError(Exception):
    pass


def banner(msg):
    print(f"\n{'='*60}")
    print(f"  {msg}")
    print(f"{'='*60}")


def _safe_print(s):
    """Print, with fallback for terminal encoding issues (e.g. garbled CJK)."""
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("utf-8", errors="replace").decode("utf-8", errors="replace"))


def warn(msg):
    _safe_print(f"  [WARN] {msg}")


def err(msg):
    _safe_print(f"  [ERROR] {msg}")


def info(msg):
    _safe_print(f"  [INFO] {msg}")


def run_cmd(cmd, timeout=600, check=True):
    """Run a shell command and return subprocess.CompletedProcess."""
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout
        )
        if check and result.returncode != 0:
            stderr = result.stderr.strip()
            if "InvalidAccessKeyId" in stderr:
                raise DownloadError("Authentication failed: invalid AccessKey ID or Secret")
            if "NoSuchBucket" in stderr:
                raise DownloadError(f"Bucket not found: verify --oss-path")
            if "AccessDenied" in stderr:
                raise DownloadError(f"Access denied: check credentials and OSS path")
            if "No such file" in stderr:
                raise DownloadError(f"Remote path not found: check --oss-path")
            print(f"  [ossutil stderr] {stderr[:300]}")
            raise DownloadError(f"ossutil exited with code {result.returncode}")
        return result
    except subprocess.TimeoutExpired:
        raise DownloadError(f"ossutil timed out after {timeout}s")


# ===================================================================
# Phase 0: Preflight
# ===================================================================

def preflight_checks(excel_path, output_dir):
    """Verify ossutil and inputs exist before starting.

    Returns resolved Excel path (glob expanded if needed).
    """
    banner("Preflight Checks")

    # ossutil
    if not os.path.exists(OSSUTIL):
        raise DownloadError(f"ossutil not found at {OSSUTIL}")
    if not os.access(OSSUTIL, os.X_OK):
        info("ossutil lacks execute permission, applying chmod +x")
        os.chmod(OSSUTIL, 0o755)

    # Excel (support glob patterns and directories)
    resolved_excel = excel_path
    if os.path.isdir(excel_path):
        matches = glob.glob(os.path.join(excel_path, "*.xlsx"))
        if not matches:
            matches = glob.glob(os.path.join(excel_path, "*.xls"))
        if len(matches) == 1:
            resolved_excel = matches[0]
        elif not matches:
            raise DownloadError(f"No .xlsx/.xls files found in directory: {excel_path}")
        else:
            raise DownloadError(f"Multiple .xlsx files in directory: {matches}")
    elif not os.path.exists(excel_path):
        matches = glob.glob(excel_path)
        if len(matches) == 1:
            resolved_excel = matches[0]
        elif not matches:
            raise DownloadError(f"Excel file not found: {excel_path}")
        else:
            raise DownloadError(f"Multiple files match pattern: {excel_path}")
    info(f"Excel file: {resolved_excel}")

    # Output dir
    raw_dir = os.path.join(output_dir, "RawData")
    os.makedirs(raw_dir, exist_ok=True)
    info(f"Output directory: {raw_dir}")

    # Disk space (simple check: we need at least some space)
    usage = shutil.disk_usage(output_dir)
    free_gb = usage.free / (1024 ** 3)
    info(f"Free disk space: {free_gb:.1f} GB")
    if free_gb < 1:
        raise DownloadError(f"Insufficient disk space ({free_gb:.1f} GB free)")

    return resolved_excel


# ===================================================================
# Phase 1: Argument parsing
# ===================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Download RNA-seq data from Alibaba Cloud OSS"
    )
    p.add_argument("--excel", required=True, metavar="PATH",
                   help="Path to Excel sample info sheet")
    p.add_argument("--project-id", required=True, metavar="ID",
                   help="Project identifier (also used as download subdirectory name)")
    p.add_argument("--oss-path", default=DEFAULT_OSS_PATH, metavar="URL",
                   help=f"OSS bucket path (default: {DEFAULT_OSS_PATH})")
    p.add_argument("--endpoint", default=DEFAULT_ENDPOINT, metavar="URL",
                   help=f"OSS endpoint (default: {DEFAULT_ENDPOINT})")
    p.add_argument("--access-key-id", default=None, metavar="KEY",
                   help="OSS AccessKey ID (env: OSS_ACCESS_KEY_ID)")
    p.add_argument("--access-key-secret", default=None, metavar="SEC",
                   help="OSS AccessKey Secret (env: OSS_ACCESS_KEY_SECRET)")
    p.add_argument("--output-dir", default=None, metavar="DIR",
                   help=f"Download destination (default: {DEFAULT_BASE_DIR}/PROJECT_ID)")
    p.add_argument("--jobs", type=int, default=4, metavar="N",
                   help="Parallel download jobs (default: 4)")
    p.add_argument("--dry-run", action="store_true",
                   help="List matched files only, do not download")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing files without prompting")
    return p.parse_args()


# ===================================================================
# Phase 2: Credential resolution
# ===================================================================

def resolve_credentials(args):
    """Resolve OSS credentials: CLI > env > ossutil default config."""
    access_key_id = args.access_key_id or os.environ.get("OSS_ACCESS_KEY_ID")
    access_key_secret = args.access_key_secret or os.environ.get("OSS_ACCESS_KEY_SECRET")

    if not access_key_id or not access_key_secret:
        raise DownloadError(
            "OSS credentials not found.\n"
            "  Set environment variables:\n"
            "    export OSS_ACCESS_KEY_ID=\"your-key-id\"\n"
            "    export OSS_ACCESS_KEY_SECRET=\"your-key-secret\"\n"
            "  Or pass via CLI:\n"
            "    --access-key-id KEY --access-key-secret SEC"
        )
    return access_key_id, access_key_secret


# ===================================================================
# Phase 3: Excel parsing
# ===================================================================

def parse_excel(excel_path):
    """Extract sample names from the Excel sample info sheet.

    Searches for a column containing '样品名称' or 'EP管标记' across
    all sheets, auto-detecting the header row. Returns (sample_names,
    group_map) where group_map is {sample: group} or None if no group
    column found.
    """
    import pandas as pd

    info(f"Reading: {excel_path}")
    xls = pd.ExcelFile(excel_path)
    info(f"Sheets: {xls.sheet_names}")

    samples = []
    group_map = {}

    for sheet in xls.sheet_names:
        # Try to auto-detect the header row (search rows 0-6)
        df = None
        name_col = None
        for header_row in range(7):
            df = pd.read_excel(xls, sheet_name=sheet, header=header_row)
            if df.empty:
                continue
            for col in df.columns:
                col_str = str(col).strip()
                if "样品名称" in col_str or "EP管标记" in col_str or "EP管" in col_str:
                    name_col = col
                    break
            if name_col is not None:
                info(f"Found header at row {header_row} in sheet '{sheet}'")
                break

        if name_col is None or df is None:
            continue

        # Extract sample names (skip rows where the value IS the header itself)
        vals = df[name_col].dropna().astype(str).str.strip()
        vals = [v for v in vals if v and v.lower() not in ("nan", "null", "")
                and "样品名称" not in v and "EP" not in v]

        if vals:
            samples = vals
            info(f"Found {len(samples)} samples in sheet '{sheet}', column '{name_col}'")

            # Try to find group column
            for col in df.columns:
                col_str = str(col).strip()
                if col_str in ("组名", "group", "Group", "condition", "Condition",
                               "分组", "treatment", "Treatment"):
                    # Map sample -> group for rows where sample name exists
                    for _, row in df.iterrows():
                        s = str(row[name_col]).strip() if pd.notna(row[name_col]) else ""
                        g = str(row[col]).strip() if pd.notna(row[col]) else ""
                        if s and s.lower() not in ("nan", "null", "") and g:
                            group_map[s] = g
                    if group_map:
                        info(f"Found group column '{col}': {set(group_map.values())}")
                    break
            break  # Found the right sheet

    if not samples:
        raise DownloadError(
            "Could not find sample names in Excel.\n"
            "Expected a column containing '样品名称' or 'EP管标记'."
        )

    # Validate
    if len(samples) != len(set(samples)):
        seen = set()
        dups = [s for s in samples if s in seen or seen.add(s)]
        warn(f"Duplicate sample names found: {dups}")

    return samples, (group_map if group_map else None)


def parse_project_info(excel_path):
    """Extract project metadata (species, project name, etc.) from the Excel.

    Returns a dict with keys: species (human/mouse/rat), read_depth, sample_type.
    """
    import pandas as pd

    xls = pd.ExcelFile(excel_path)
    proj = {"species": "human", "read_depth": "6G", "sample_type": ""}

    if "项目信息" not in xls.sheet_names:
        return proj

    df = pd.read_excel(xls, sheet_name="项目信息", header=None)
    if df.empty:
        return proj

    # Flatten the form: scan columns 0-3 for key-value pairs
    key_value = {}
    num_cols = len(df.columns)
    for _, row in df.iterrows():
        for col_idx in range(0, num_cols - 1, 2):
            key = str(row.iloc[col_idx]).strip() if pd.notna(row.iloc[col_idx]) else ""
            val = str(row.iloc[col_idx + 1]).strip() if pd.notna(row.iloc[col_idx + 1]) else ""
            if key and key.lower() not in ("nan", "null"):
                key_value[key] = val

    # Species mapping
    for key, val in key_value.items():
        key_clean = key.replace("*", "").replace(":", "").strip()
        if "物种" in key_clean:
            species_cn = val.strip()
            if "人" in species_cn and "大" not in species_cn:
                proj["species"] = "human"
            elif "大鼠" in species_cn or "rat" in species_cn.lower():
                proj["species"] = "rat"
            elif "小鼠" in species_cn or "鼠" in species_cn or "mouse" in species_cn.lower():
                proj["species"] = "mouse"
            info(f"Species: {species_cn} → {proj['species']}")
        if "深度" in key_clean:
            proj["read_depth"] = val
        if "样本类型" in key_clean or "样品类型" in key_clean:
            proj["sample_type"] = val

    return proj


def generate_config_yaml(output_dir, project_id, species, mode="fastq"):
    """Generate config.yaml by copying the template and filling in auto-detected values.

    Reads config/config.yaml from the pipeline install as the template base.
    """
    banner("Generating config.yaml")

    config_path = os.path.join(output_dir, "config.yaml")
    if os.path.exists(config_path):
        warn(f"config.yaml already exists, skipping: {config_path}")
        return

    # Find the template
    script_dir = os.path.dirname(os.path.abspath(__file__))
    template = os.path.join(script_dir, "..", "config", "config.yaml")

    if not os.path.exists(template):
        warn(f"Template not found: {template}, generating minimal config")
        content = _generate_minimal_config(output_dir, project_id, species, mode)
    else:
        with open(template) as f:
            content = f.read()
        # Fill in auto-detected values
        content = content.replace('project_dir: ""',
                                   f'project_dir: "{output_dir}"')
        content = content.replace('project_id: "Project-XXX"',
                                   f'project_id: "{project_id}"')
        content = content.replace('species: "human"',
                                   f'species: "{species}"')
        content = content.replace('mode: "counts"',
                                   f'mode: "{mode}"')
        content = content.replace('  root: ""',
                                   '  root: "/home/zhuzp/reference"')

    with open(config_path, "w") as f:
        f.write(content)

    info(f"Generated: {config_path}")
    info(f"  species: {species}, mode: {mode}")


def _generate_minimal_config(output_dir, project_id, species, mode):
    """Fallback config if template not found."""
    return f"""# Bulkseq auto-generated config
project_dir: "{output_dir}"
project_id: "{project_id}"
species: "{species}"
mode: "{mode}"

reference:
  root: "/home/zhuzp/reference"
  index: ""
  gtf: ""

upstream:
  trim_front: 15
  skip_qc: false
  multimapping: false

deg:
  logFC_cutoff: 0.5
  p_cutoff: 0.05

comparison_mode: "auto"
comparisons: []

exclude_samples: []

batch:
  correct: false

enrichment:
  go_kegg: true
  gsea: true
  gsea_target_genes: []

checkpoint:
  after_qc: true
  after_alignment: true
"""


# ===================================================================
# Phase 4: OSS file listing and matching
# ===================================================================

def list_oss_files(oss_base, endpoint, access_key_id, access_key_secret):
    """List all objects under the OSS path. Returns list of object keys (full paths)."""
    banner("Listing OSS Files")

    # ossutil ls lists all objects under prefix recursively by default
    cmd = [
        OSSUTIL, "ls", f"{oss_base}/",
        "-s",
        "-e", endpoint,
        "-i", access_key_id,
        "-k", access_key_secret,
    ]
    info(f"Running: ossutil ls {oss_base}/")

    result = run_cmd(cmd, timeout=60)

    files = []
    for line in result.stdout.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        # Skip non-object lines (bucket, prefix headers, summary)
        if line.startswith("oss://"):
            files.append(line)

    # Also parse lines like "oss://bucket/prefix/  (size)" from short format
    # Short format (-s) output: "oss://bucket/key" or "oss://bucket/key   size"
    parsed = []
    for f in files:
        # Take the first token (the path) split by whitespace
        parts = f.split()
        key = parts[0] if parts else f
        if key.endswith("/"):
            continue  # skip directory markers
        parsed.append(key)

    if not parsed:
        raise DownloadError(
            f"No files found at {oss_base}/\n"
            "Check the OSS path and credentials."
        )

    info(f"Found {len(parsed)} objects in OSS")
    return parsed


def match_samples(samples, oss_files):
    """Match sample names to R1/R2 file pairs in the OSS listing.

    Returns (matched, partial, missing) where:
    - matched: {sample: (r1_oss_path, r2_oss_path)}
    - partial: {sample: [paths]} (only one of the pair found)
    - missing: [sample] (no files found)
    """
    banner("Matching Samples to OSS Files")

    matched = {}
    partial = {}
    missing = []

    for sample in samples:
        # Find all OSS files containing this sample name
        candidates = [f for f in oss_files if sample in os.path.basename(f)]

        if not candidates:
            missing.append(sample)
            continue

        # Try each paired pattern
        found = False
        for r1_suffix, r2_suffix in PAIRED_PATTERNS:
            r1_target = f"{sample}{r1_suffix}"
            r2_target = f"{sample}{r2_suffix}"
            r1_match = None
            r2_match = None
            for c in candidates:
                basename = os.path.basename(c)
                # Match: sample name anywhere in filename, suffix at end
                if basename.endswith(r1_target) and sample in basename:
                    r1_match = c
                elif basename.endswith(r2_target) and sample in basename:
                    r2_match = c
            if r1_match and r2_match:
                matched[sample] = (r1_match, r2_match)
                found = True
                break

        if not found:
            # Partial match: found files but not a complete pair
            partial[sample] = candidates

    # Print results table
    print(f"\n  {'Sample':<20s} {'Status':<12s} {'R1':<50s}")
    print(f"  {'-'*20} {'-'*12} {'-'*50}")
    for s in samples:
        if s in matched:
            r1, r2 = matched[s]
            print(f"  {s:<20s} {'MATCHED':<12s} {os.path.basename(r1):<50s}")
        elif s in partial:
            print(f"  {s:<20s} {'PARTIAL':<12s} {os.path.basename(partial[s][0]):<50s}")
        else:
            print(f"  {s:<20s} {'MISSING':<12s} -")

    print(f"\n  Summary: {len(matched)} matched, {len(partial)} partial, {len(missing)} missing")

    if missing:
        warn(f"Missing samples (no files found in OSS): {missing}")
    if partial:
        warn(f"Partial matches (only one of pair found): {list(partial.keys())}")

    if not matched:
        raise DownloadError(
            "No sample-to-file matches found.\n"
            "Try --dry-run first and check the OSS file naming with -s flag."
        )

    return matched, partial, missing


# ===================================================================
# Phase 4b: MD5 checksum
# ===================================================================

def download_md5_file(oss_base, oss_files, endpoint, access_key_id, access_key_secret):
    """Find and download the MD5 checksum file from OSS.

    Returns the local path to the downloaded MD5 file, or None if not found.
    """
    banner("MD5 Checksum File")

    # Look for MD5*.txt in the OSS listing
    md5_oss = None
    for f in oss_files:
        basename = os.path.basename(f)
        if basename.startswith("MD5") and basename.endswith(".txt"):
            md5_oss = f
            break

    if not md5_oss:
        warn("No MD5 checksum file found in OSS")
        return None

    local_md5 = os.path.join("/tmp", os.path.basename(md5_oss))
    info(f"Downloading: {os.path.basename(md5_oss)}")

    cmd = [
        OSSUTIL, "cp", md5_oss, local_md5,
        "-e", endpoint,
        "-i", access_key_id,
        "-k", access_key_secret,
        "-f",
    ]
    try:
        run_cmd(cmd, timeout=30)
        info(f"MD5 file saved to {local_md5}")
        return local_md5
    except DownloadError as e:
        warn(f"Failed to download MD5 file: {e}")
        return None


def parse_md5_file(md5_path):
    """Parse an MD5 checksum file into {oss_relative_path: md5_hash} dict."""
    md5_map = {}
    with open(md5_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                md5_hash = parts[0]
                # The path in the MD5 file is relative: "RawData/SAMPLE/filename.fq.gz"
                rel_path = parts[1]
                md5_map[rel_path] = md5_hash
    info(f"Parsed {len(md5_map)} MD5 entries")
    return md5_map


def compute_md5(filepath):
    """Compute MD5 hash of a local file. Returns hex digest."""
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_file_md5(local_file, oss_file, md5_map, oss_base):
    """Verify a downloaded file's MD5 against the OSS checksum.

    Returns (ok: bool, expected: str, actual: str).
    """
    if not md5_map:
        return True, "", ""

    # Derive the relative path as it appears in the MD5 file
    # OSS file: oss://bucket/prefix/RawData/SAMPLE/filename.fq.gz
    # MD5 entry:        RawData/SAMPLE/filename.fq.gz
    prefix = oss_base.rstrip("/")
    rel_path = oss_file
    if rel_path.startswith(prefix):
        rel_path = rel_path[len(prefix):].lstrip("/")

    expected = md5_map.get(rel_path)
    if not expected:
        # Try matching by basename only (some MD5 files use different prefixes)
        basename = os.path.basename(oss_file)
        for k, v in md5_map.items():
            if k.endswith(basename):
                expected = v
                break

    if not expected:
        return True, "", ""  # Can't verify, skip

    actual = compute_md5(local_file)
    ok = (expected.lower() == actual.lower())
    return ok, expected, actual


# ===================================================================
# Phase 5: Download
# ===================================================================

def download_files(matched, output_dir, endpoint, access_key_id, access_key_secret,
                   jobs, force, oss_base, md5_map=None):
    """Download matched OSS files to local RawData/ directory.

    Files are renamed to {sample}_1.fq.gz / {sample}_2.fq.gz.
    If md5_map is provided, verifies each file's MD5 after download.
    """
    banner("Downloading Files")

    raw_dir = os.path.join(output_dir, "RawData")
    total = len(matched)
    success = 0
    failed = []
    total_bytes = 0
    md5_results = {}  # {sample: (r1_ok, r2_ok)}

    for i, (sample, (r1_oss, r2_oss)) in enumerate(sorted(matched.items()), 1):
        local_r1 = os.path.join(raw_dir, f"{sample}_1.fq.gz")
        local_r2 = os.path.join(raw_dir, f"{sample}_2.fq.gz")
        r1_ok, r2_ok = True, True
        sample_failed = False

        for oss_file, local_file, label in [(r1_oss, local_r1, "R1"),
                                             (r2_oss, local_r2, "R2")]:
            if os.path.exists(local_file) and not force:
                size = os.path.getsize(local_file)
                info(f"[{i}/{total}] {sample} {label}: SKIP (exists, {size/1e9:.1f} GB)")
                total_bytes += size
                # Verify MD5 for existing files too
                if md5_map:
                    ok, _, _ = verify_file_md5(local_file, oss_file, md5_map, oss_base)
                    if not ok:
                        err(f"[{i}/{total}] {sample} {label}: MD5 MISMATCH (existing file)")
                        if label == "R1":
                            r1_ok = False
                        else:
                            r2_ok = False
                continue

            info(f"[{i}/{total}] {sample} {label}: downloading...")
            t0 = time.time()

            cmd = [
                OSSUTIL, "cp", oss_file, local_file,
                "-e", endpoint,
                "-i", access_key_id,
                "-k", access_key_secret,
                "-j", str(jobs),
                "--retry-times", "10",
            ]
            if force:
                cmd.append("-f")

            try:
                run_cmd(cmd, timeout=3600)
                elapsed = time.time() - t0
                size = os.path.getsize(local_file)
                speed = size / elapsed / 1e6 if elapsed > 0 else 0
                total_bytes += size
                info(f"[{i}/{total}] {sample} {label}: OK ({size/1e9:.2f} GB, "
                     f"{elapsed:.0f}s, {speed:.1f} MB/s)")
            except DownloadError as e:
                err(f"[{i}/{total}] {sample} {label}: FAILED - {e}")
                failed.append((sample, label, oss_file))
                sample_failed = True
                break

            # MD5 verification
            if md5_map:
                ok, expected, actual = verify_file_md5(local_file, oss_file, md5_map, oss_base)
                if not ok:
                    err(f"[{i}/{total}] {sample} {label}: MD5 MISMATCH!")
                    err(f"  Expected: {expected}")
                    err(f"  Actual:   {actual}")
                    if label == "R1":
                        r1_ok = False
                    else:
                        r2_ok = False
                else:
                    info(f"[{i}/{total}] {sample} {label}: MD5 OK")

        if not sample_failed:
            success += 1
        md5_results[sample] = (r1_ok, r2_ok)

    return success, failed, total_bytes, md5_results


# ===================================================================
# Phase 6: Verification and summary
# ===================================================================

def verify_and_summarize(samples, matched, output_dir, start_time, success, failed,
                         total_bytes, partial, missing, md5_results=None):
    """Verify downloads and print summary, including MD5 checksums."""
    banner("Download Summary")

    raw_dir = os.path.join(output_dir, "RawData")
    elapsed = time.time() - start_time

    # File existence check
    verified = 0
    for sample in matched:
        r1 = os.path.join(raw_dir, f"{sample}_1.fq.gz")
        r2 = os.path.join(raw_dir, f"{sample}_2.fq.gz")
        r1_ok = os.path.exists(r1) and os.path.getsize(r1) > 0
        r2_ok = os.path.exists(r2) and os.path.getsize(r2) > 0
        if r1_ok and r2_ok:
            verified += 1
        else:
            warn(f"{sample}: incomplete download (R1={'OK' if r1_ok else 'MISS'}, "
                 f"R2={'OK' if r2_ok else 'MISS'})")

    # MD5 summary
    md5_ok = 0
    md5_fail = 0
    md5_skip = 0
    if md5_results:
        for sample, (r1_ok, r2_ok) in md5_results.items():
            if r1_ok and r2_ok:
                md5_ok += 1
            else:
                md5_fail += 1
        md5_skip = len(matched) - md5_ok - md5_fail
    else:
        md5_skip = len(matched)

    print(f"""
  Project:     {os.path.basename(output_dir)}
  Target:      {raw_dir}/
  Samples:     {len(samples)} total, {len(matched)} matched, {len(partial)} partial, {len(missing)} missing
  Downloaded:  {success} samples ({verified} file-exists verified)
  Failed:      {len(failed)}
  Total size:  {total_bytes / 1e9:.2f} GB
  Duration:    {elapsed / 60:.1f} min""")

    if md5_results:
        print(f"  MD5:         {md5_ok} OK, {md5_fail} MISMATCH, {md5_skip} skipped\n")
    else:
        print(f"  MD5:         skipped (no checksum file)\n")

    if md5_fail > 0:
        warn("MD5 mismatches detected! The affected files may be corrupted.")
        for sample, (r1_ok, r2_ok) in md5_results.items():
            if not r1_ok or not r2_ok:
                parts = []
                if not r1_ok:
                    parts.append("R1")
                if not r2_ok:
                    parts.append("R2")
                print(f"    - {sample}: {', '.join(parts)}")

    if failed:
        warn("Failed downloads:")
        for s, label, path in failed:
            print(f"    - {s} {label}: {path}")

    print("  Next steps:")
    print(f"    Run: bash run_pipeline.sh {output_dir} --mode fastq")


# ===================================================================
# Phase 7: Optional sample_sheet.csv generation
# ===================================================================

def generate_sample_sheet(output_dir, samples, group_map):
    """Generate sample_sheet.csv from sample names and optional group mapping."""
    banner("Generating sample_sheet.csv")

    csv_path = os.path.join(output_dir, "sample_sheet.csv")
    if os.path.exists(csv_path):
        warn(f"sample_sheet.csv already exists, skipping: {csv_path}")
        return

    with open(csv_path, "w") as f:
        if group_map:
            f.write("sample,group\n")
            for s in samples:
                g = group_map.get(s, "")
                f.write(f"{s},{g}\n")
        else:
            f.write("sample,group\n")
            for s in samples:
                f.write(f"{s},\n")

    info(f"Generated: {csv_path}")
    if not group_map:
        warn("No group column found in Excel - please fill in group names manually")


# ===================================================================
# Main
# ===================================================================

def main():
    start_time = time.time()

    args = parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(DEFAULT_BASE_DIR, args.project_id)

    info(f"Project ID: {args.project_id}")
    info(f"OSS path:   {args.oss_path}")
    info(f"Output dir: {args.output_dir}")

    # Phase 0
    resolved_excel = preflight_checks(args.excel, args.output_dir)

    # Phase 2
    access_key_id, access_key_secret = resolve_credentials(args)

    # Phase 3
    samples, group_map = parse_excel(resolved_excel)
    banner("Sample List")
    for i, s in enumerate(samples, 1):
        g = f" ({group_map.get(s, '?' )}) " if group_map else ""
        print(f"  [{i}] {s}{g}")

    # Phase 3b: Extract project metadata (species, etc.)
    project_info = parse_project_info(resolved_excel)
    species = project_info["species"]

    # Phase 4
    oss_files = list_oss_files(
        args.oss_path, args.endpoint, access_key_id, access_key_secret
    )
    matched, partial, missing = match_samples(samples, oss_files)

    # Generate sample_sheet.csv and config.yaml (always, even in dry-run)
    generate_sample_sheet(args.output_dir, samples, group_map)
    generate_config_yaml(args.output_dir, args.project_id, species, mode="fastq")

    if args.dry_run:
        banner("Dry Run Complete")
        print("\n  No files downloaded. Remove --dry-run to proceed.")
        return 0

    # Phase 4b: Download MD5 checksum file
    md5_path = download_md5_file(
        args.oss_path, oss_files, args.endpoint, access_key_id, access_key_secret
    )
    md5_map = parse_md5_file(md5_path) if md5_path else None

    # Phase 5
    md5_results = {}
    if matched:
        success, failed, total_bytes, md5_results = download_files(
            matched, args.output_dir, args.endpoint,
            access_key_id, access_key_secret, args.jobs, args.force,
            args.oss_path, md5_map
        )
    else:
        success, failed, total_bytes = 0, [], 0

    # Phase 6
    verify_and_summarize(
        samples, matched, args.output_dir, start_time,
        success, failed, total_bytes, partial, missing, md5_results
    )

    return 0 if not failed else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except DownloadError as e:
        err(str(e))
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n  Interrupted by user.")
        sys.exit(130)
