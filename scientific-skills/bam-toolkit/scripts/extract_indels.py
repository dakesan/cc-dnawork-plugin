#!/usr/bin/env python3
"""Extract insertion and deletion sequences from BAM file.

This script extracts insertions and deletions from reads in a specified
genomic region by parsing CIGAR strings.
"""

import json
import sys
from collections import defaultdict
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


def extract_indels_from_read(
    read: pysam.AlignedSegment,
    min_indel_size: int = 1,
) -> tuple[list[dict], list[dict]]:
    """Extract insertions and deletions from a single read.

    Args:
        read: pysam AlignedSegment object
        min_indel_size: Minimum indel size in bp

    Returns:
        Tuple of (insertions list, deletions list)
    """
    insertions = []
    deletions = []

    if not read.cigartuples:
        return insertions, deletions

    # Track position in read and reference
    read_pos = 0
    ref_pos = read.reference_start

    for operation, length in read.cigartuples:
        # M (0): alignment match
        # I (1): insertion to reference
        # D (2): deletion from reference
        # N (3): skipped region from reference
        # S (4): soft clipping
        # H (5): hard clipping
        # P (6): padding
        # = (7): sequence match
        # X (8): sequence mismatch

        if operation == 0:  # M - match/mismatch
            read_pos += length
            ref_pos += length

        elif operation == 1:  # I - insertion
            if length >= min_indel_size:
                # Extract inserted sequence
                inserted_seq = read.query_sequence[read_pos:read_pos + length]
                insertions.append({
                    "position": ref_pos,
                    "size": length,
                    "sequence": inserted_seq,
                    "read_name": read.query_name,
                    "mapping_quality": read.mapping_quality,
                })
            read_pos += length

        elif operation == 2:  # D - deletion
            if length >= min_indel_size:
                deletions.append({
                    "position": ref_pos,
                    "size": length,
                    "read_name": read.query_name,
                    "mapping_quality": read.mapping_quality,
                })
            ref_pos += length

        elif operation == 4:  # S - soft clip
            read_pos += length

        elif operation == 3:  # N - skipped region
            ref_pos += length

        # H (5), P (6), = (7), X (8) don't consume query or reference

    return insertions, deletions


def group_indels(indels: list[dict], indel_type: str) -> list[dict]:
    """Group identical indels and count occurrences.

    Args:
        indels: List of indel dictionaries
        indel_type: "insertion" or "deletion"

    Returns:
        List of grouped indels with counts
    """
    # Group by position and sequence/size
    groups = defaultdict(list)

    for indel in indels:
        if indel_type == "insertion":
            key = (indel["position"], indel["sequence"])
        else:  # deletion
            key = (indel["position"], indel["size"])

        groups[key].append(indel)

    # Create grouped output
    grouped = []
    for key, group in groups.items():
        # Use first occurrence as representative
        representative = group[0].copy()
        representative["count"] = len(group)
        representative["supporting_reads"] = [r["read_name"] for r in group]
        grouped.append(representative)

    # Sort by position
    grouped.sort(key=lambda x: x["position"])

    return grouped


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
        help="JSON output file path (default: stdout)",
    ),
    min_mapq: int = typer.Option(
        20,
        "--min-mapq",
        help="Minimum mapping quality",
    ),
    min_indel_size: int = typer.Option(
        1,
        "--min-indel-size",
        help="Minimum indel size in bp",
    ),
) -> None:
    """Extract insertions and deletions from BAM file.

    Examples:
        # Extract all indels
        python extract_indels.py --bam input.bam --region chr1:1000-2000 --output indels.json

        # Extract large indels only (>= 5bp)
        python extract_indels.py --bam input.bam --region chr1:1000-2000 --min-indel-size 5 --output large_indels.json
    """
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

    # Collect indels
    all_insertions = []
    all_deletions = []
    total_reads = 0
    processed_reads = 0

    try:
        for read in bam_file.fetch(chrom, start, end):
            total_reads += 1

            # Filter by mapping quality
            if read.mapping_quality < min_mapq:
                continue

            processed_reads += 1

            # Extract indels from this read
            insertions, deletions = extract_indels_from_read(read, min_indel_size)
            all_insertions.extend(insertions)
            all_deletions.extend(deletions)

    except Exception as e:
        typer.echo(f"Error fetching reads: {e}", err=True)
        typer.echo(f"Make sure the region '{region}' exists and BAM file is indexed.", err=True)
        raise typer.Exit(1)

    bam_file.close()

    # Group identical indels
    grouped_insertions = group_indels(all_insertions, "insertion")
    grouped_deletions = group_indels(all_deletions, "deletion")

    # Create output
    output_data = {
        "region": region,
        "filters": {
            "min_mapq": min_mapq,
            "min_indel_size": min_indel_size,
        },
        "reads": {
            "total": total_reads,
            "processed": processed_reads,
        },
        "summary": {
            "total_insertions": len(all_insertions),
            "total_deletions": len(all_deletions),
            "unique_insertions": len(grouped_insertions),
            "unique_deletions": len(grouped_deletions),
        },
        "insertions": grouped_insertions,
        "deletions": grouped_deletions,
    }

    json_str = json.dumps(output_data, indent=2)

    if output:
        output.write_text(json_str)
        typer.echo(
            f"Extracted {len(grouped_insertions)} unique insertions "
            f"and {len(grouped_deletions)} unique deletions to {output}"
        )
    else:
        typer.echo(json_str)


if __name__ == "__main__":
    app()
