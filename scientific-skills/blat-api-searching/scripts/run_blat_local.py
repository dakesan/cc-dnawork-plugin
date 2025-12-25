#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
import json
from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(add_completion=False)

BLAT_DATA_DIR = Path("~/.local/share/blat").expanduser()
HG38_2BIT_URL = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit"
CHM13_2BIT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit"


def _require_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found in PATH: {tool}")


def _read_single_fasta_sequence(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"FASTA not found: {path}")
    header_count = 0
    sequence_lines = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header_count += 1
                if header_count > 1:
                    raise ValueError("Multiple sequences detected in FASTA.")
                continue
            sequence_lines.append(line)
    if not sequence_lines:
        raise ValueError("No sequence data found in FASTA.")
    return "".join(sequence_lines)


def _write_temp_fasta(sequence: str) -> Path:
    tmp = tempfile.NamedTemporaryFile(prefix="blat_query_", suffix=".fa", delete=False)
    tmp_path = Path(tmp.name)
    with tmp:
        tmp.write(b">query\n")
        tmp.write(sequence.encode("ascii", errors="strict"))
        tmp.write(b"\n")
    return tmp_path


def _parse_psl(path: Path) -> list[dict]:
    records: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return records
    data_lines = lines
    if lines[0].startswith("psLayout") or lines[0].startswith("match"):
        data_lines = lines[5:]
    for line in data_lines:
        parts = line.split("\t")
        if len(parts) < 21:
            continue
        records.append(
            {
                "match": int(parts[0]),
                "mismatch": int(parts[1]),
                "rep_match": int(parts[2]),
                "n_count": int(parts[3]),
                "q_num_insert": int(parts[4]),
                "q_base_insert": int(parts[5]),
                "t_num_insert": int(parts[6]),
                "t_base_insert": int(parts[7]),
                "strand": parts[8],
                "q_name": parts[9],
                "q_size": int(parts[10]),
                "q_start": int(parts[11]),
                "q_end": int(parts[12]),
                "t_name": parts[13],
                "t_size": int(parts[14]),
                "t_start": int(parts[15]),
                "t_end": int(parts[16]),
                "block_count": int(parts[17]),
                "block_sizes": parts[18].strip(","),
                "q_starts": parts[19].strip(","),
                "t_starts": parts[20].strip(","),
            }
        )
    return records


def _ensure_hg38_2bit() -> Path:
    target = BLAT_DATA_DIR / "hg38.2bit"
    if target.exists():
        return target
    BLAT_DATA_DIR.mkdir(parents=True, exist_ok=True)
    _download_2bit(HG38_2BIT_URL, target)
    return target


def _download_2bit(url: str, dest: Path) -> None:
    _require_tool("curl")
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "-L", "-o", str(dest), url]
    subprocess.run(cmd, check=True)


def _ensure_chm13_2bit() -> Path:
    target = BLAT_DATA_DIR / "CHM13.2bit"
    if target.exists():
        return target
    BLAT_DATA_DIR.mkdir(parents=True, exist_ok=True)
    _download_2bit(CHM13_2BIT_URL, target)
    return target


def _resolve_reference(reference: str) -> Path:
    if reference == "hg38":
        return _ensure_hg38_2bit()
    if reference == "CHM13":
        return _ensure_chm13_2bit()
    raise typer.BadParameter("reference must be one of: hg38, CHM13")


@app.command()
def run(
    reference: str = typer.Option(
        ...,
        "--reference",
        help="Reference genome assembly (hg38 or CHM13).",
    ),
    fasta: Optional[Path] = typer.Option(
        None,
        "--fasta",
        help="Path to FASTA file (single sequence only).",
    ),
    sequence: Optional[str] = typer.Option(
        None,
        "--sequence",
        help="Raw DNA sequence string.",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Output path for JSON results.",
    ),
) -> None:
    """
    Run local BLAT using a 2bit reference.
    """
    _require_tool("blat")
    ref_path = _resolve_reference(reference)
    if (fasta is None and sequence is None) or (fasta is not None and sequence is not None):
        raise typer.BadParameter("Provide exactly one of --fasta or --sequence.")

    if fasta is not None:
        _read_single_fasta_sequence(fasta)
        query_fasta = fasta
    else:
        query_seq = (sequence or "").strip()
        if not query_seq:
            raise typer.BadParameter("sequence must be non-empty.")
        query_fasta = _write_temp_fasta(query_seq)

    psl_path = Path(tempfile.mkstemp(prefix="blat_output_", suffix=".psl")[1])

    cmd = ["blat", str(ref_path), str(query_fasta), str(psl_path)]
    subprocess.run(cmd, check=True)

    records = _parse_psl(psl_path)
    json_text = json.dumps(records, indent=2, ensure_ascii=True)
    if output is not None:
        output.write_text(json_text, encoding="utf-8")
    typer.echo(json_text)


if __name__ == "__main__":
    app()
