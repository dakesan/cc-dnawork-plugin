#!/usr/bin/env python3
"""Extract reads from BAM file for specific genomic region.

This script extracts reads from a BAM/SAM/CRAM file for a specified genomic
region and outputs them as either BAM or JSON format.
"""

import json
import sys
from pathlib import Path
from typing import Optional

import pysam
import typer

app = typer.Typer()


def parse_region(region: str) -> tuple[str, int, int]:
    """Parse genomic region string into chromosome, start, and end.

    Args:
        region: Genomic region string (e.g., "chr1:1000-2000")

    Returns:
        Tuple of (chromosome, start, end)

    Raises:
        ValueError: If region format is invalid
    """
    try:
        if ":" not in region:
            raise ValueError("Region must contain ':'")

        chrom, positions = region.split(":")
        start, end = positions.split("-")

        return chrom, int(start), int(end)
    except (ValueError, AttributeError) as e:
        raise ValueError(
            f"Invalid region format: {region}. "
            "Expected format: chr1:1000-2000"
        ) from e


def read_to_dict(read: pysam.AlignedSegment) -> dict:
    """Convert pysam AlignedSegment to dictionary.

    Args:
        read: pysam AlignedSegment object

    Returns:
        Dictionary containing read information
    """
    return {
        "query_name": read.query_name,
        "reference_start": read.reference_start,
        "reference_end": read.reference_end,
        "sequence": read.query_sequence,
        "quality": "".join(chr(q + 33) for q in read.query_qualities) if read.query_qualities else None,
        "mapping_quality": read.mapping_quality,
        "is_reverse": read.is_reverse,
        "cigar": read.cigarstring,
        "is_proper_pair": read.is_proper_pair,
        "is_duplicate": read.is_duplicate,
        "is_paired": read.is_paired,
        "mate_is_reverse": read.mate_is_reverse if read.is_paired else None,
        "mate_reference_name": read.next_reference_name if read.is_paired else None,
        "mate_reference_start": read.next_reference_start if read.is_paired else None,
        "template_length": read.template_length if read.is_paired else None,
    }


@app.command()
def main(
    bam: Path = typer.Option(
        ...,
        "--bam",
        help="Input BAM/SAM/CRAM file path",
        exists=True,
        readable=True,
    ),
    region: str = typer.Option(
        ...,
        "--region",
        help="Genomic region (e.g., chr1:1000-2000)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Output file path (BAM or JSON based on extension)",
    ),
    format: str = typer.Option(
        "bam",
        "--format",
        help="Output format: 'bam' or 'json'",
    ),
    min_mapq: int = typer.Option(
        0,
        "--min-mapq",
        help="Minimum mapping quality",
    ),
    proper_pairs: bool = typer.Option(
        False,
        "--proper-pairs",
        help="Only properly paired reads",
    ),
    no_duplicates: bool = typer.Option(
        False,
        "--no-duplicates",
        help="Exclude duplicate reads",
    ),
) -> None:
    """Extract reads from BAM file for specific genomic region.

    Examples:
        # Extract reads as BAM
        python extract_reads.py --bam input.bam --region chr1:1000-2000 --output reads.bam

        # Extract reads as JSON
        python extract_reads.py --bam input.bam --region chr1:1000-2000 --format json --output reads.json

        # Extract high-quality properly-paired reads
        python extract_reads.py --bam input.bam --region chr1:1000-2000 --min-mapq 30 --proper-pairs --no-duplicates --output filtered.bam
    """
    # Validate format
    if format not in ["bam", "json"]:
        typer.echo(f"Error: Invalid format '{format}'. Must be 'bam' or 'json'.", err=True)
        raise typer.Exit(1)

    # Parse region
    try:
        chrom, start, end = parse_region(region)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    # Open BAM file
    try:
        bam_file = pysam.AlignmentFile(str(bam), "rb")
    except Exception as e:
        typer.echo(f"Error opening BAM file: {e}", err=True)
        raise typer.Exit(1)

    # Collect reads
    reads = []
    total_reads = 0

    try:
        for read in bam_file.fetch(chrom, start, end):
            total_reads += 1

            # Apply filters
            if read.mapping_quality < min_mapq:
                continue
            if proper_pairs and not read.is_proper_pair:
                continue
            if no_duplicates and read.is_duplicate:
                continue

            reads.append(read)
    except Exception as e:
        typer.echo(f"Error fetching reads: {e}", err=True)
        typer.echo(f"Make sure the region '{region}' exists and BAM file is indexed.", err=True)
        raise typer.Exit(1)

    # Output
    if format == "json":
        # JSON output
        output_data = {
            "region": region,
            "total_reads": len(reads),
            "filtered_from": total_reads,
            "filters": {
                "min_mapq": min_mapq,
                "proper_pairs_only": proper_pairs,
                "no_duplicates": no_duplicates,
            },
            "reads": [read_to_dict(read) for read in reads],
        }

        json_str = json.dumps(output_data, indent=2)

        if output:
            output.write_text(json_str)
            typer.echo(f"Extracted {len(reads)} reads (from {total_reads} total) to {output}")
        else:
            typer.echo(json_str)

    else:  # BAM output
        if not output:
            typer.echo("Error: --output is required for BAM format", err=True)
            raise typer.Exit(1)

        # Create output BAM file with same header
        try:
            out_bam = pysam.AlignmentFile(str(output), "wb", template=bam_file)

            for read in reads:
                out_bam.write(read)

            out_bam.close()
            typer.echo(f"Extracted {len(reads)} reads (from {total_reads} total) to {output}")
        except Exception as e:
            typer.echo(f"Error writing BAM file: {e}", err=True)
            raise typer.Exit(1)

    bam_file.close()


if __name__ == "__main__":
    app()
