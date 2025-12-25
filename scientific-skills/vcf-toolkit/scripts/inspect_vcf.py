#!/usr/bin/env python3
"""
Inspect and filter VCF files, exporting results as JSON.
"""
import json
import sys
from pathlib import Path
from typing import Optional

import pysam
import typer

app = typer.Typer(add_completion=False)


def _parse_region(region: str) -> tuple[str, int, int]:
    """Parse region string like 'chr1:1000-2000' into components."""
    if ":" not in region:
        raise ValueError(f"Invalid region format: {region}. Expected 'chr:start-end'")

    chrom, pos_range = region.split(":", 1)
    if "-" not in pos_range:
        raise ValueError(f"Invalid region format: {region}. Expected 'chr:start-end'")

    start_str, end_str = pos_range.split("-", 1)
    start = int(start_str)
    end = int(end_str)

    return chrom, start, end


def _apply_filters(variant, min_qual, min_dp, min_af, max_af, pass_only) -> bool:
    """Apply filter conditions to a variant. Returns True if variant passes."""
    # Quality filter
    if min_qual is not None and variant.qual is not None:
        if variant.qual < min_qual:
            return False

    # PASS only filter
    if pass_only:
        if variant.filter.keys() and "PASS" not in variant.filter.keys():
            return False

    # INFO/DP filter
    if min_dp is not None:
        dp = variant.info.get("DP")
        if dp is None or dp < min_dp:
            return False

    # INFO/AF filter
    if min_af is not None or max_af is not None:
        af = variant.info.get("AF")
        if af is None:
            return False
        # AF can be a list for multi-allelic sites
        af_values = af if isinstance(af, (list, tuple)) else [af]
        for af_val in af_values:
            if min_af is not None and af_val < min_af:
                return False
            if max_af is not None and af_val > max_af:
                return False

    return True


def _variant_to_dict(variant) -> dict:
    """Convert pysam VariantRecord to dictionary."""
    # Convert INFO fields
    info_dict = {}
    for key in variant.info:
        value = variant.info[key]
        # Convert tuples to lists for JSON serialization
        if isinstance(value, tuple):
            value = list(value)
        info_dict[key] = value

    # Convert filter
    filter_list = list(variant.filter.keys()) if variant.filter.keys() else []

    # Convert sample genotypes
    samples_dict = {}
    for sample_name in variant.samples:
        sample = variant.samples[sample_name]
        sample_dict = {}
        for key in sample.keys():
            value = sample[key]
            # Convert tuples to lists
            if isinstance(value, tuple):
                value = list(value)
            sample_dict[key] = value
        samples_dict[sample_name] = sample_dict

    return {
        "chrom": variant.chrom,
        "pos": variant.pos,
        "id": variant.id if variant.id else None,
        "ref": variant.ref,
        "alts": list(variant.alts) if variant.alts else [],
        "qual": float(variant.qual) if variant.qual is not None else None,
        "filter": filter_list,
        "info": info_dict,
        "samples": samples_dict
    }


@app.command()
def main(
    vcf: Path = typer.Option(
        ...,
        "--vcf",
        help="Input VCF file path.",
    ),
    chrom: Optional[str] = typer.Option(
        None,
        "--chrom",
        help="Chromosome to analyze (e.g., chr1). Required unless --region is specified.",
    ),
    region: Optional[str] = typer.Option(
        None,
        "--region",
        help="Genomic region (e.g., chr1:1000-2000). Alternative to --chrom.",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Output JSON file path (default: stdout).",
    ),
    min_qual: Optional[float] = typer.Option(
        None,
        "--min-qual",
        help="Minimum quality score (QUAL >= X).",
    ),
    min_dp: Optional[int] = typer.Option(
        None,
        "--min-dp",
        help="Minimum depth (INFO/DP >= X).",
    ),
    min_af: Optional[float] = typer.Option(
        None,
        "--min-af",
        help="Minimum allele frequency (INFO/AF >= X).",
    ),
    max_af: Optional[float] = typer.Option(
        None,
        "--max-af",
        help="Maximum allele frequency (INFO/AF <= X).",
    ),
    pass_only: bool = typer.Option(
        True,
        "--pass-only/--all-filters",
        help="Only include FILTER=PASS variants (default: True). Use --all-filters to include all.",
    ),
    max_variants: int = typer.Option(
        100,
        "--max-variants",
        help="Maximum number of variants allowed (default: 100).",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Override variant limit (may produce very large JSON).",
    ),
) -> None:
    """
    Inspect and filter VCF files, exporting results as JSON.

    Examples:
        # Export all PASS variants from chr1 (up to 100)
        python inspect_vcf.py --vcf input.vcf --chrom chr1 --output chr1.json

        # Filter by quality and depth
        python inspect_vcf.py --vcf input.vcf --chrom chr1 --min-qual 30 --min-dp 10 --output filtered.json

        # Specify region
        python inspect_vcf.py --vcf input.vcf --region chr1:1000000-2000000 --output region.json

        # Include all filters (not just PASS)
        python inspect_vcf.py --vcf input.vcf --chrom chr1 --all-filters --output all.json

        # Override limit for large files
        python inspect_vcf.py --vcf large.vcf --chrom chr1 --force --output large.json
    """
    # Validate chrom or region
    if not chrom and not region:
        typer.echo("Error: Either --chrom or --region must be specified.", err=True)
        raise typer.Exit(1)

    if chrom and region:
        typer.echo("Error: Cannot specify both --chrom and --region. Use one or the other.", err=True)
        raise typer.Exit(1)

    # Validate input file
    if not vcf.exists():
        typer.echo(f"Error: VCF file not found: {vcf}", err=True)
        raise typer.Exit(1)

    # Open VCF file
    try:
        vcf_file = pysam.VariantFile(str(vcf))
    except Exception as e:
        typer.echo(f"Error: Failed to open VCF file: {e}", err=True)
        raise typer.Exit(1)

    # Get samples
    samples = list(vcf_file.header.samples)

    # Collect variants
    variants = []
    variant_count = 0

    # Determine iterator (with chrom or region)
    if region:
        try:
            chrom_name, start, end = _parse_region(region)
            iterator = vcf_file.fetch(chrom_name, start, end)
        except Exception as e:
            typer.echo(f"Error: Invalid region or fetch failed: {e}", err=True)
            raise typer.Exit(1)
    elif chrom:
        try:
            iterator = vcf_file.fetch(chrom)
        except Exception as e:
            typer.echo(f"Error: Failed to fetch chromosome {chrom}: {e}", err=True)
            raise typer.Exit(1)
    else:
        # This should not happen due to validation above
        iterator = vcf_file

    # Iterate and filter
    for variant in iterator:
        # Apply filters
        if not _apply_filters(variant, min_qual, min_dp, min_af, max_af, pass_only):
            continue

        variant_count += 1

        # Check limit
        if not force and variant_count > max_variants:
            vcf_file.close()

            # Build filter summary
            filter_summary = []
            if chrom:
                filter_summary.append(f"--chrom {chrom}")
            if region:
                filter_summary.append(f"--region {region}")
            if min_qual:
                filter_summary.append(f"--min-qual {min_qual}")
            if min_dp:
                filter_summary.append(f"--min-dp {min_dp}")
            if min_af:
                filter_summary.append(f"--min-af {min_af}")
            if max_af:
                filter_summary.append(f"--max-af {max_af}")
            if pass_only:
                filter_summary.append("--pass-only")
            else:
                filter_summary.append("--all-filters")

            filter_str = " ".join(filter_summary)

            typer.echo(f"\nError: VCF contains {variant_count}+ variants after filtering (limit: {max_variants}).\n", err=True)
            typer.echo("Suggestions:", err=True)
            typer.echo("  - Apply more restrictive filters: --min-qual, --min-dp, --pass-only", err=True)
            typer.echo("  - Specify a genomic region: --region chr1:1000-2000", err=True)
            typer.echo("  - Override limit with --force (warning: may produce very large JSON)", err=True)
            typer.echo("  - Use bcftools directly for large-scale processing\n", err=True)
            typer.echo(f"Current filter conditions:\n  {filter_str}\n", err=True)
            raise typer.Exit(1)

        # Convert to dict
        variants.append(_variant_to_dict(variant))

    vcf_file.close()

    # Build result
    result = {
        "num_variants": variant_count,
        "samples": samples,
        "variants": variants
    }

    # Convert to JSON
    json_str = json.dumps(result, indent=2, ensure_ascii=False)

    # Output
    if output:
        output.write_text(json_str, encoding="utf-8")

        # Calculate file size
        file_size_mb = output.stat().st_size / (1024 * 1024)

        if force and variant_count > max_variants:
            typer.echo(f"⚠ Warning: Exported {variant_count} variants (limit bypassed with --force)", err=True)
            typer.echo(f"✓ Successfully exported {variant_count} variants to {output} ({file_size_mb:.1f} MB)", err=True)
        else:
            typer.echo(f"✓ Successfully exported {variant_count} variants to {output}", err=True)
    else:
        # Output to stdout
        typer.echo(json_str)


if __name__ == "__main__":
    app()
