#!/usr/bin/env python3
"""
VCF Statistics Calculator

Calculate comprehensive statistics from VCF files and output as JSON.
"""

import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pysam
import typer

app = typer.Typer()


def calculate_stats(vcf_path: Path, chrom: Optional[str] = None, region: Optional[str] = None) -> dict:
    """Calculate comprehensive statistics from VCF file."""
    vcf = pysam.VariantFile(str(vcf_path))

    # Set up iterator based on region/chrom
    if region:
        iterator = vcf.fetch(region=region)
    elif chrom:
        iterator = vcf.fetch(chrom)
    else:
        iterator = vcf.fetch()

    # Initialize statistics collectors
    total_variants = 0
    filter_counts = defaultdict(int)
    quality_scores = []
    depths = []
    allele_frequencies = []
    variant_types = defaultdict(int)
    chrom_counts = defaultdict(int)

    # Collect statistics
    for variant in iterator:
        total_variants += 1

        # Filter counts
        for f in variant.filter:
            filter_counts[f] += 1

        # Quality score
        if variant.qual is not None:
            quality_scores.append(variant.qual)

        # Depth (INFO/DP)
        if "DP" in variant.info:
            depths.append(variant.info["DP"])

        # Allele frequency (INFO/AF)
        if "AF" in variant.info:
            af_values = variant.info["AF"]
            if isinstance(af_values, (list, tuple)):
                allele_frequencies.extend(af_values)
            else:
                allele_frequencies.append(af_values)

        # Variant type (SNP, insertion, deletion)
        ref_len = len(variant.ref)
        for alt in variant.alts:
            alt_len = len(alt)
            if ref_len == alt_len == 1:
                variant_types["SNP"] += 1
            elif alt_len > ref_len:
                variant_types["insertion"] += 1
            else:
                variant_types["deletion"] += 1

        # Chromosome counts
        chrom_counts[variant.chrom] += 1

    vcf.close()

    # Build result dictionary
    result = {
        "total_variants": total_variants,
        "filter_counts": dict(filter_counts),
        "variant_types": dict(variant_types),
        "chrom_counts": dict(chrom_counts),
    }

    # Quality score statistics
    if quality_scores:
        result["quality_stats"] = {
            "min": round(min(quality_scores), 2),
            "max": round(max(quality_scores), 2),
            "mean": round(statistics.mean(quality_scores), 2),
            "median": round(statistics.median(quality_scores), 2),
        }
    else:
        result["quality_stats"] = None

    # Depth statistics
    if depths:
        result["depth_stats"] = {
            "min": min(depths),
            "max": max(depths),
            "mean": round(statistics.mean(depths), 2),
            "median": statistics.median(depths),
        }
    else:
        result["depth_stats"] = None

    # Allele frequency statistics
    if allele_frequencies:
        result["allele_frequency_stats"] = {
            "min": round(min(allele_frequencies), 4),
            "max": round(max(allele_frequencies), 4),
            "mean": round(statistics.mean(allele_frequencies), 4),
            "median": round(statistics.median(allele_frequencies), 4),
        }
    else:
        result["allele_frequency_stats"] = None

    return result


@app.command()
def main(
    vcf: Path = typer.Option(..., "--vcf", help="Input VCF file path."),
    chrom: Optional[str] = typer.Option(None, "--chrom", help="Chromosome to analyze"),
    region: Optional[str] = typer.Option(
        None, "--region", help="Genomic region (e.g., chr1:1000-2000)"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", help="Output JSON file path (default: stdout)"
    ),
) -> None:
    """Calculate VCF statistics and output as JSON.

    Examples:
        # Statistics for all chromosomes
        python vcf_stats.py --vcf variants.vcf.gz --output stats.json

        # Statistics for chr1 only
        python vcf_stats.py --vcf variants.vcf.gz --chrom chr1

        # Statistics for specific region
        python vcf_stats.py --vcf variants.vcf.gz --region chr1:10000-20000
    """
    # Validate arguments
    if chrom and region:
        typer.echo("Error: Cannot specify both --chrom and --region.", err=True)
        raise typer.Exit(1)

    # Calculate statistics
    try:
        stats = calculate_stats(vcf, chrom, region)
    except Exception as e:
        typer.echo(f"Error calculating statistics: {e}", err=True)
        raise typer.Exit(1)

    # Output JSON
    json_output = json.dumps(stats, indent=2)

    if output:
        output.write_text(json_output)
        typer.echo(f"✓ Statistics written to {output}")
    else:
        typer.echo(json_output)


if __name__ == "__main__":
    app()
