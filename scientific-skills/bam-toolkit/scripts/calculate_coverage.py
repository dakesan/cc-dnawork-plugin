#!/usr/bin/env python3
"""Calculate coverage statistics for genomic region.

This script calculates coverage statistics for a specified genomic region
using pileup operations.
"""

import json
import sys
from pathlib import Path
from statistics import mean, median
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
        0,
        "--min-mapq",
        help="Minimum mapping quality",
    ),
    min_baseq: int = typer.Option(
        0,
        "--min-baseq",
        help="Minimum base quality",
    ),
) -> None:
    """Calculate coverage statistics for genomic region.

    Examples:
        # Calculate coverage statistics
        python calculate_coverage.py --bam input.bam --region chr1:1000-2000 --output coverage.json

        # High-quality bases only
        python calculate_coverage.py --bam input.bam --region chr1:1000-2000 --min-mapq 30 --min-baseq 20 --output hq_coverage.json
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

    # Calculate coverage using pileup
    coverage_array = []

    try:
        # Initialize coverage array with zeros
        coverage_array = [0] * (end - start)

        # Pileup over the region
        for pileupcolumn in bam_file.pileup(
            chrom,
            start,
            end,
            truncate=True,
            min_mapping_quality=min_mapq,
            min_base_quality=min_baseq,
        ):
            # Get position relative to start
            pos = pileupcolumn.pos
            if start <= pos < end:
                coverage = pileupcolumn.n
                coverage_array[pos - start] = coverage

    except Exception as e:
        typer.echo(f"Error calculating coverage: {e}", err=True)
        typer.echo(f"Make sure the region '{region}' exists and BAM file is indexed.", err=True)
        raise typer.Exit(1)

    bam_file.close()

    # Calculate statistics
    total_bases = len(coverage_array)
    bases_with_coverage = sum(1 for c in coverage_array if c > 0)
    bases_without_coverage = total_bases - bases_with_coverage

    # Handle case where all bases have zero coverage
    if bases_with_coverage == 0:
        mean_cov = 0
        median_cov = 0
    else:
        mean_cov = mean(coverage_array)
        median_cov = median(coverage_array)

    min_cov = min(coverage_array) if coverage_array else 0
    max_cov = max(coverage_array) if coverage_array else 0
    percent_covered = (bases_with_coverage / total_bases * 100) if total_bases > 0 else 0

    # Create output
    output_data = {
        "region": region,
        "filters": {
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
        },
        "statistics": {
            "total_bases": total_bases,
            "mean_coverage": round(mean_cov, 2),
            "median_coverage": median_cov,
            "min_coverage": min_cov,
            "max_coverage": max_cov,
            "bases_with_coverage": bases_with_coverage,
            "bases_without_coverage": bases_without_coverage,
            "percent_covered": round(percent_covered, 2),
        },
    }

    json_str = json.dumps(output_data, indent=2)

    if output:
        output.write_text(json_str)
        typer.echo(
            f"Coverage statistics for {region}: "
            f"mean={mean_cov:.2f}, median={median_cov}, "
            f"coverage={percent_covered:.1f}%"
        )
    else:
        typer.echo(json_str)


if __name__ == "__main__":
    app()
