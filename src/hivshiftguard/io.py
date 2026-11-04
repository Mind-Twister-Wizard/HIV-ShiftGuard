import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

SOURCES = {
    "PI": "https://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
    "NRTI": "https://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt",
    "NNRTI": "https://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt",
    "INI": "https://hivdb.stanford.edu/download/GenoPhenoDatasets/INI_DataSet.txt",
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, allow_nan=False)
    os.replace(tmp, path)


def validate_header(path):
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        first = f.readline().rstrip("\n\r").split("\t")
        second = f.readline()
    if "SeqID" not in first or "P1" not in first or "CompMutList" not in first or len(first) < 20 or not second:
        raise ValueError(f"Invalid genotype-phenotype file: {path}. Check upstream page and filenames.")
    return len(first)


def ensure_data(root, config):
    folder = Path(root) / "data" / "raw"
    folder.mkdir(parents=True, exist_ok=True)
    records = []
    for gene, url in SOURCES.items():
        dest = folder / (gene + "_DataSet.txt")
        if dest.exists():
            try:
                validate_header(dest)
                print(f"[DATA OK] {dest.name} (reusing, sha256 verified in manifest)", flush=True)
            except ValueError:
                invalid = dest.with_suffix(".invalid")
                dest.replace(invalid)
        if not dest.exists():
            print(f"[DOWNLOAD] {url}", flush=True)
            err = None
            for trial in range(config.get("max_download_attempts", 4)):
                tmp = dest.with_suffix(".download")
                try:
                    with requests.get(url, headers={"User-Agent": "HIVShiftGuardAcademicResearch/0.9 (+https://hivdb.stanford.edu/)"},
                                      timeout=config.get("download_timeout_seconds", 60), stream=True) as response:
                        response.raise_for_status()
                        with tmp.open("wb") as out:
                            for block in response.iter_content(chunk_size=1 << 16):
                                if block:
                                    out.write(block)
                    validate_header(tmp)
                    os.replace(tmp, dest)
                    err = None
                    break
                except Exception as exc:
                    err = exc
                    tmp.unlink(missing_ok=True)
                    if trial + 1 < config.get("max_download_attempts", 4):
                        time.sleep(1 + trial * 2)
            if err is not None:
                raise RuntimeError(f"Download failed for {gene}: {err}. Manually download {url} to {dest} and rerun.")
        columns = validate_header(dest)
        records.append({"dataset": gene, "url": url, "filename": dest.name,
                        "sha256": sha256_file(dest), "bytes": dest.stat().st_size,
                        "columns": columns, "checked_utc": datetime.now(timezone.utc).isoformat()})
    return records
