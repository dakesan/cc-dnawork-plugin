#!/usr/bin/env python3
import time
from pathlib import Path
from typing import Optional

import requests
import typer

app = typer.Typer(add_completion=False)

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
MIN_POLL_INTERVAL_SEC = 3


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


def _submit_query(
    query: str,
    program: str,
    database: str,
    megablast: bool,
    expect: Optional[float],
    hitlist_size: Optional[int],
    timeout_sec: int,
) -> str:
    payload = {
        "CMD": "Put",
        "QUERY": query,
        "PROGRAM": program,
        "DATABASE": database,
    }
    if megablast:
        payload["MEGABLAST"] = "on"
    if expect is not None:
        payload["EXPECT"] = str(expect)
    if hitlist_size is not None:
        payload["HITLIST_SIZE"] = str(hitlist_size)
    response = requests.post(BLAST_URL, data=payload, timeout=timeout_sec)
    response.raise_for_status()
    rid = None
    for line in response.text.splitlines():
        if line.startswith("RID ="):
            rid = line.split("=", 1)[1].strip()
            break
    if not rid:
        raise RuntimeError("RID not found in BLAST submission response.")
    return rid


def _poll_until_ready(
    rid: str,
    poll_interval_sec: int,
    max_wait_sec: int,
    timeout_sec: int,
) -> None:
    elapsed = 0
    while elapsed <= max_wait_sec:
        params = {"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}
        response = requests.get(BLAST_URL, params=params, timeout=timeout_sec)
        response.raise_for_status()
        text = response.text
        if "Status=READY" in text:
            if "ThereAreHits=yes" in text or "ThereAreHits=no" in text:
                return
            return
        if "Status=FAILED" in text:
            raise RuntimeError("BLAST search failed.")
        if "Status=UNKNOWN" in text:
            raise RuntimeError("BLAST search expired or unknown RID.")
        time.sleep(poll_interval_sec)
        elapsed += poll_interval_sec
    raise TimeoutError("BLAST search did not complete in time.")


def _write_output(content: bytes, output: Optional[Path], format_type: str) -> Path:
    is_zip = content[:2] == b"PK"
    if output is None:
        if is_zip:
            suffix = ".zip"
        else:
            suffix = ".json" if format_type.upper().startswith("JSON") else ".txt"
        output = Path(f"blast_result{suffix}")
    output.write_bytes(content)
    return output


@app.command()
def run(
    program: str = typer.Option(
        "blastn",
        "--program",
        help="BLAST program (blastn, blastp, blastx, tblastn, tblastx).",
    ),
    database: str = typer.Option(
        "core_nt",
        "--database",
        help="BLAST database (e.g., core_nt, swissprot).",
    ),
    fasta: Optional[Path] = typer.Option(
        None,
        "--fasta",
        help="Path to FASTA file (single sequence only).",
    ),
    sequence: Optional[str] = typer.Option(
        None,
        "--sequence",
        help="Raw query sequence string.",
    ),
    format_type: str = typer.Option(
        "JSON2",
        "--format-type",
        help="Result format (JSON2, XML2, Text, HTML).",
    ),
    megablast: bool = typer.Option(
        False,
        "--megablast",
        help="Enable megablast (blastn only).",
    ),
    expect: Optional[float] = typer.Option(
        None,
        "--expect",
        help="E-value threshold (e.g., 10, 1e-3).",
    ),
    hitlist_size: Optional[int] = typer.Option(
        None,
        "--hitlist-size",
        help="Max number of hits to keep.",
    ),
    poll_interval_sec: int = typer.Option(
        10,
        "--poll-interval-sec",
        help="Polling interval in seconds (min 3).",
    ),
    max_wait_sec: int = typer.Option(
        600,
        "--max-wait-sec",
        help="Maximum wait time for completion.",
    ),
    timeout_sec: int = typer.Option(
        30,
        "--timeout-sec",
        help="HTTP request timeout in seconds.",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Output path for result file.",
    ),
) -> None:
    """
    Run BLAST via NCBI Common URL API.
    """
    if (fasta is None and sequence is None) or (fasta is not None and sequence is not None):
        raise typer.BadParameter("Provide exactly one of --fasta or --sequence.")
    if poll_interval_sec < MIN_POLL_INTERVAL_SEC:
        raise typer.BadParameter("poll-interval-sec must be >= 3.")

    if fasta is not None:
        query = _read_single_fasta_sequence(fasta)
    else:
        query = (sequence or "").strip()
        if not query:
            raise typer.BadParameter("sequence must be non-empty.")

    rid = _submit_query(
        query=query,
        program=program,
        database=database,
        megablast=megablast,
        expect=expect,
        hitlist_size=hitlist_size,
        timeout_sec=timeout_sec,
    )
    _poll_until_ready(
        rid=rid,
        poll_interval_sec=poll_interval_sec,
        max_wait_sec=max_wait_sec,
        timeout_sec=timeout_sec,
    )

    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": format_type,
    }
    response = requests.get(BLAST_URL, params=params, timeout=timeout_sec)
    response.raise_for_status()
    content = response.content
    output_path = _write_output(content, output, format_type)
    typer.echo(str(output_path))


if __name__ == "__main__":
    app()
