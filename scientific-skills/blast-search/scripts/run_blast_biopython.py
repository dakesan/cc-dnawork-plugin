#!/usr/bin/env python3
"""
Run NCBI BLAST using BioPython (NCBIWWW.qblast + NCBIXML parser).
"""
import json
from pathlib import Path
from typing import Optional

import typer
from Bio.Blast import NCBIWWW, NCBIXML

app = typer.Typer(add_completion=False)


def _read_single_fasta_sequence(path: Path) -> str:
    """Read a single FASTA sequence from file."""
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


def _parse_blast_results(blast_record) -> dict:
    """Convert BLAST record to JSON-serializable dict."""
    results = {
        "query": blast_record.query,
        "query_length": blast_record.query_length,
        "database": blast_record.database,
        "num_hits": len(blast_record.alignments),
        "hits": []
    }

    for i, alignment in enumerate(blast_record.alignments, 1):
        hsp = alignment.hsps[0]  # Best HSP for this alignment
        percent_identity = (hsp.identities / hsp.align_length) * 100

        hit_data = {
            "rank": i,
            "accession": alignment.accession,
            "title": alignment.title,
            "length": alignment.length,
            "e_value": float(hsp.expect),
            "bit_score": float(hsp.bits),
            "score": float(hsp.score),
            "identities": hsp.identities,
            "align_length": hsp.align_length,
            "percent_identity": round(percent_identity, 2),
            "gaps": hsp.gaps,
            "query_start": hsp.query_start,
            "query_end": hsp.query_end,
            "subject_start": hsp.sbjct_start,
            "subject_end": hsp.sbjct_end,
        }
        results["hits"].append(hit_data)

    return results


@app.command()
def run(
    program: str = typer.Option(
        "blastn",
        "--program",
        help="BLAST program (blastn, blastp, blastx, tblastn, tblastx).",
    ),
    database: str = typer.Option(
        "nt",
        "--database",
        help="BLAST database (nt, nr, refseq_rna, swissprot, etc.).",
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
    organism: Optional[str] = typer.Option(
        None,
        "--organism",
        help="Restrict search to organism (e.g., 'Homo sapiens', 'Escherichia coli').",
    ),
    expect: Optional[float] = typer.Option(
        0.001,
        "--expect",
        help="E-value threshold (default: 0.001).",
    ),
    hitlist_size: Optional[int] = typer.Option(
        10,
        "--hitlist-size",
        help="Maximum number of hits to return (default: 10).",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Output path for JSON results.",
    ),
) -> None:
    """
    Run BLAST via NCBI Web API using BioPython.

    Examples:
        # Search with FASTA file
        python run_blast_biopython.py --fasta query.fasta

        # Search with sequence string
        python run_blast_biopython.py --sequence ATGCGATCG...

        # Restrict to human sequences
        python run_blast_biopython.py --fasta query.fasta --organism "Homo sapiens"

        # Protein BLAST
        python run_blast_biopython.py --program blastp --database swissprot --sequence MTEYK...
    """
    # Validate input
    if (fasta is None and sequence is None) or (fasta is not None and sequence is not None):
        raise typer.BadParameter("Provide exactly one of --fasta or --sequence.")

    # Read sequence
    if fasta is not None:
        query_sequence = _read_single_fasta_sequence(fasta)
    else:
        query_sequence = (sequence or "").strip()
        if not query_sequence:
            raise typer.BadParameter("sequence must be non-empty.")

    # Prepare entrez_query for organism filtering
    entrez_query = None
    if organism:
        entrez_query = f"{organism}[Organism]"

    # Submit BLAST query
    typer.echo(f"Submitting BLAST search to NCBI...", err=True)
    typer.echo(f"Program: {program}, Database: {database}", err=True)
    if entrez_query:
        typer.echo(f"Organism filter: {organism}", err=True)

    result_handle = NCBIWWW.qblast(
        program=program,
        database=database,
        sequence=query_sequence,
        expect=expect,
        hitlist_size=hitlist_size,
        entrez_query=entrez_query,
    )

    # Parse results
    typer.echo("Parsing results...", err=True)
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    # Convert to JSON
    results = _parse_blast_results(blast_record)
    json_text = json.dumps(results, indent=2, ensure_ascii=False)

    # Output
    if output is not None:
        output.write_text(json_text, encoding="utf-8")
        typer.echo(f"Results saved to: {output}", err=True)

    typer.echo(json_text)


if __name__ == "__main__":
    app()
