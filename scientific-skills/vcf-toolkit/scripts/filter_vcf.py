#!/usr/bin/env python3
"""
VCF Filtering Script

Filter VCF files by quality, depth, allele frequency, and other criteria.
Output filtered variants as a new VCF file.
"""

from pathlib import Path
from typing import Optional

import pysam
import typer

app = typer.Typer()


def apply_filters(
    variant,
    min_qual: Optional[float] = None,
    min_dp: Optional[int] = None,
    min_af: Optional[float] = None,
    max_af: Optional[float] = None,
    pass_only: bool = False,
) -> bool:
    """Check if variant passes all filter criteria."""
    # Quality filter
    if min_qual is not None and (variant.qual is None or variant.qual < min_qual):
        return False

    # Depth filter
    if min_dp is not None:
        if "DP" not in variant.info or variant.info["DP"] < min_dp:
            return False

    # Allele frequency filters
    if min_af is not None or max_af is not None:
        if "AF" not in variant.info:
            return False
        af_values = variant.info["AF"]
        if isinstance(af_values, (list, tuple)):
            af = max(af_values)  # Use max AF for multi-allelic variants
        else:
            af = af_values
        if min_af is not None and af < min_af:
            return False
        if max_af is not None and af > max_af:
            return False

    # PASS filter
    if pass_only and "PASS" not in variant.filter:
        return False

    return True


@app.command()
def main(
    vcf: Path = typer.Option(..., "--vcf", help="Input VCF file path."),
    output: Path = typer.Option(..., "--output", help="Output VCF file path."),
    chrom: Optional[str] = typer.Option(None, "--chrom", help="Chromosome to filter"),
    region: Optional[str] = typer.Option(
        None, "--region", help="Genomic region (e.g., chr1:1000-2000)"
    ),
    min_qual: Optional[float] = typer.Option(None, "--min-qual", help="Minimum quality score"),
    min_dp: Optional[int] = typer.Option(None, "--min-dp", help="Minimum depth (INFO/DP)"),
    min_af: Optional[float] = typer.Option(None, "--min-af", help="Minimum allele frequency"),
    max_af: Optional[float] = typer.Option(None, "--max-af", help="Maximum allele frequency"),
    pass_only: bool = typer.Option(False, "--pass-only", help="Only include PASS variants"),
) -> None:
    """Filter VCF file and output filtered variants as a new VCF.

    Examples:
        # Filter by quality and depth
        python filter_vcf.py --vcf input.vcf.gz --output filtered.vcf --min-qual 30 --min-dp 10

        # Filter chr1 PASS variants only
        python filter_vcf.py --vcf input.vcf.gz --output chr1_pass.vcf --chrom chr1 --pass-only

        # Filter by allele frequency range
        python filter_vcf.py --vcf input.vcf.gz --output common.vcf --min-af 0.05 --max-af 0.95
    """
    # Validate arguments
    if chrom and region:
        typer.echo("Error: Cannot specify both --chrom and --region.", err=True)
        raise typer.Exit(1)

    # Open input VCF
    try:
        vcf_in = pysam.VariantFile(str(vcf))
    except Exception as e:
        typer.echo(f"Error opening input VCF: {e}", err=True)
        raise typer.Exit(1)

    # Create output VCF with same header
    try:
        vcf_out = pysam.VariantFile(str(output), "w", header=vcf_in.header)
    except Exception as e:
        vcf_in.close()
        typer.echo(f"Error creating output VCF: {e}", err=True)
        raise typer.Exit(1)

    # Set up iterator
    try:
        if region:
            iterator = vcf_in.fetch(region=region)
        elif chrom:
            iterator = vcf_in.fetch(chrom)
        else:
            iterator = vcf_in.fetch()
    except Exception as e:
        vcf_in.close()
        vcf_out.close()
        typer.echo(f"Error setting up iterator: {e}", err=True)
        raise typer.Exit(1)

    # Filter and write variants
    variant_count = 0
    passed_count = 0

    try:
        for variant in iterator:
            variant_count += 1
            if apply_filters(variant, min_qual, min_dp, min_af, max_af, pass_only):
                vcf_out.write(variant)
                passed_count += 1
    except Exception as e:
        typer.echo(f"Error processing variants: {e}", err=True)
        raise typer.Exit(1)
    finally:
        vcf_in.close()
        vcf_out.close()

    typer.echo(f"✓ Filtered {passed_count}/{variant_count} variants")
    typer.echo(f"✓ Output written to {output}")


if __name__ == "__main__":
    app()
