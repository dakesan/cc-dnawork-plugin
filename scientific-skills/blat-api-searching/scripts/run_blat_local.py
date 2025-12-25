#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
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
        help="Output path for PSL results.",
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
        query_seq = _read_single_fasta_sequence(fasta)
        query_fasta = fasta
    else:
        query_seq = (sequence or "").strip()
        if not query_seq:
            raise typer.BadParameter("sequence must be non-empty.")
        query_fasta = _write_temp_fasta(query_seq)

    if output is None:
        output = Path(tempfile.mkstemp(prefix="blat_output_", suffix=".psl")[1])

    cmd = ["blat", str(ref_path), str(query_fasta), str(output)]
    subprocess.run(cmd, check=True)

    if output and output.exists():
        typer.echo(output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    app()
