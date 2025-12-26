#!/usr/bin/env python3
"""Generate IGV snapshots for genomic regions with multiple BAM files.

This script generates an IGV batch script, executes it, and produces
PNG screenshots for specified genomic regions with multiple BAM tracks.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

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


def read_bed_file(bed_file: Path) -> List[tuple[str, int, int, Optional[str]]]:
    """Read regions from BED file.

    Args:
        bed_file: Path to BED file

    Returns:
        List of tuples (chromosome, start, end, name)

    Raises:
        ValueError: If BED file format is invalid
    """
    regions = []
    with open(bed_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"Invalid BED format at line {line_num}: {line}\n"
                    "Expected at least 3 columns: chrom, start, end"
                )

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else None

            regions.append((chrom, start, end, name))

    return regions


def generate_batch_script(
    genome: str,
    bam_files: List[Path],
    regions: List[tuple[str, int, int, Optional[str]]],
    output_dir: Path,
    max_panel_height: int = 500,
) -> str:
    """Generate IGV batch script content.

    Args:
        genome: Genome assembly ID (e.g., "hg38")
        bam_files: List of BAM file paths
        regions: List of (chrom, start, end, name) tuples
        output_dir: Directory for snapshot output
        max_panel_height: Maximum panel height in pixels

    Returns:
        IGV batch script content as string
    """
    lines = []

    # Initialize IGV session
    lines.append("new")
    lines.append(f"genome {genome}")

    # Load BAM files
    for bam_file in bam_files:
        lines.append(f"load {bam_file.absolute()}")

    # Set snapshot directory and panel height
    lines.append(f"snapshotDirectory {output_dir.absolute()}")
    lines.append(f"maxPanelHeight {max_panel_height}")

    # Add regions and snapshots
    for chrom, start, end, name in regions:
        region_str = f"{chrom}:{start}-{end}"
        lines.append(f"goto {region_str}")

        # Generate snapshot filename
        if name:
            snapshot_name = f"{name}.png"
        else:
            snapshot_name = f"{chrom}_{start}-{end}.png"

        lines.append(f"snapshot {snapshot_name}")

    # Exit IGV
    lines.append("exit")

    return "\n".join(lines)


def find_igv_command(igv_path: Optional[str] = None) -> List[str]:
    """Find IGV executable command.

    Args:
        igv_path: Optional path to IGV executable

    Returns:
        List of command components to run IGV

    Raises:
        FileNotFoundError: If IGV is not found
    """
    if igv_path:
        # User-specified path
        if not Path(igv_path).exists():
            raise FileNotFoundError(f"IGV not found at: {igv_path}")
        return [igv_path]

    # Try common IGV locations
    common_paths = [
        "igv.sh",
        "igv",
        "/Applications/IGV.app/Contents/MacOS/igv",
        str(Path.home() / "IGV" / "igv.sh"),
    ]

    for path in common_paths:
        try:
            result = subprocess.run(
                ["which", path] if "/" not in path else ["test", "-f", path],
                capture_output=True,
                check=False,
            )
            if result.returncode == 0:
                return [path]
        except Exception:
            continue

    raise FileNotFoundError(
        "IGV not found. Please install IGV or specify path with --igv-path.\n"
        "Install: brew install --cask igv (macOS) or download from "
        "https://software.broadinstitute.org/software/igv/download"
    )


@app.command()
def main(
    genome: str = typer.Option(
        ...,
        "--genome",
        help="Genome assembly ID (e.g., hg38, hg19, mm39)",
    ),
    bam: List[Path] = typer.Option(
        ...,
        "--bam",
        help="BAM file path(s). Can specify multiple files.",
        exists=True,
        readable=True,
    ),
    region: Optional[str] = typer.Option(
        None,
        "--region",
        help="Genomic region (e.g., chr1:1000-2000)",
    ),
    bed: Optional[Path] = typer.Option(
        None,
        "--bed",
        help="BED file with regions",
        exists=True,
        readable=True,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        help="Output directory for snapshots",
    ),
    igv_path: Optional[str] = typer.Option(
        None,
        "--igv-path",
        help="Path to IGV executable (default: auto-detect)",
    ),
    max_panel_height: int = typer.Option(
        500,
        "--max-panel-height",
        help="Maximum panel height in pixels",
    ),
    java_heap: str = typer.Option(
        "2g",
        "--java-heap",
        help="Java heap size (e.g., 2g, 4g)",
    ),
    save_batch_script: bool = typer.Option(
        False,
        "--save-batch-script",
        help="Save IGV batch script for debugging",
    ),
) -> None:
    """Generate IGV snapshots for genomic regions with multiple BAM files.

    Examples:
        # Single region, multiple BAMs
        python generate_igv_snapshots.py --genome hg38 --bam sample1.bam sample2.bam --region chr1:1000-2000 --output-dir ./snapshots

        # Multiple regions from BED file
        python generate_igv_snapshots.py --genome hg38 --bam sample1.bam sample2.bam --bed regions.bed --output-dir ./snapshots

        # Save batch script for debugging
        python generate_igv_snapshots.py --genome hg38 --bam sample.bam --region chr1:1000-2000 --output-dir ./snapshots --save-batch-script
    """
    # Validate input
    if not region and not bed:
        typer.echo("Error: Either --region or --bed must be specified", err=True)
        raise typer.Exit(1)

    if region and bed:
        typer.echo("Error: Cannot specify both --region and --bed", err=True)
        raise typer.Exit(1)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse regions
    regions = []
    if region:
        try:
            chrom, start, end = parse_region(region)
            regions.append((chrom, start, end, None))
        except ValueError as e:
            typer.echo(f"Error: {e}", err=True)
            raise typer.Exit(1)
    else:
        try:
            regions = read_bed_file(bed)
        except Exception as e:
            typer.echo(f"Error reading BED file: {e}", err=True)
            raise typer.Exit(1)

    # Generate batch script
    batch_content = generate_batch_script(
        genome=genome,
        bam_files=bam,
        regions=regions,
        output_dir=output_dir,
        max_panel_height=max_panel_height,
    )

    # Write batch script to file
    if save_batch_script:
        batch_file = output_dir / "igv_batch_script.txt"
        batch_file.write_text(batch_content)
        typer.echo(f"Batch script saved to: {batch_file}")
    else:
        # Use temporary file
        batch_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        )
        batch_file.write(batch_content)
        batch_file.close()
        batch_file = Path(batch_file.name)

    # Find IGV command
    try:
        igv_cmd = find_igv_command(igv_path)
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        if not save_batch_script:
            batch_file.unlink()
        raise typer.Exit(1)

    # Execute IGV
    typer.echo(f"Running IGV with {len(bam)} BAM file(s) and {len(regions)} region(s)...")

    try:
        # Run IGV in batch mode
        cmd = igv_cmd + ["-b", str(batch_file)]

        # Add Java heap size if using java command
        if "java" in igv_cmd[0]:
            cmd = ["java", f"-Xmx{java_heap}", "-jar"] + cmd

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            typer.echo(f"Error running IGV:\n{result.stderr}", err=True)
            if not save_batch_script:
                batch_file.unlink()
            raise typer.Exit(1)

    except Exception as e:
        typer.echo(f"Error executing IGV: {e}", err=True)
        if not save_batch_script:
            batch_file.unlink()
        raise typer.Exit(1)

    # Clean up temporary batch script
    if not save_batch_script:
        batch_file.unlink()

    # Report results
    typer.echo(f"\nSnapshots generated successfully in: {output_dir}")
    typer.echo(f"Total regions processed: {len(regions)}")

    # List generated snapshots
    snapshots = sorted(output_dir.glob("*.png"))
    if snapshots:
        typer.echo(f"\nGenerated {len(snapshots)} snapshot(s):")
        for snapshot in snapshots:
            typer.echo(f"  - {snapshot.name}")
    else:
        typer.echo("\nWarning: No snapshots found. Check IGV output for errors.")


if __name__ == "__main__":
    app()
