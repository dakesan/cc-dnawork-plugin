#!/usr/bin/env python3
"""Query COSMIC Cancer Gene Census for gene information.

This script queries the COSMIC Cancer Gene Census to retrieve information
about cancer-related genes.
"""

import json
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd
import typer

app = typer.Typer()


def get_default_gene_census_path() -> Path:
    """Get the default path to Cancer Gene Census file.

    Returns:
        Path to data/cancer_gene_census.csv relative to script location
    """
    script_dir = Path(__file__).parent
    skill_dir = script_dir.parent
    return skill_dir / "data" / "cancer_gene_census.csv"


def load_gene_census(file_path: Path) -> pd.DataFrame:
    """Load Cancer Gene Census CSV file.

    Args:
        file_path: Path to cancer_gene_census.csv

    Returns:
        DataFrame containing Cancer Gene Census data

    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If file format is invalid
    """
    if not file_path.exists():
        raise FileNotFoundError(
            f"Cancer Gene Census file not found at: {file_path}\n\n"
            "To use this tool, please download COSMIC data:\n\n"
            "1. Register for free academic access:\n"
            "   https://cancer.sanger.ac.uk/cosmic/register\n\n"
            "2. Download Cancer Gene Census:\n"
            "   https://cancer.sanger.ac.uk/cosmic/download\n"
            "   File: cancer_gene_census.csv (GRCh38)\n\n"
            "3. Place the file at:\n"
            f"   {file_path}\n\n"
            "For more information, see: cosmic-toolkit/data/README.md"
        )

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise ValueError(f"Error reading Cancer Gene Census file: {e}")

    # Verify required column exists
    if "Gene Symbol" not in df.columns:
        raise ValueError(
            "Invalid Cancer Gene Census format: 'Gene Symbol' column not found"
        )

    return df


def query_genes(
    gene_census: pd.DataFrame,
    genes: List[str]
) -> dict:
    """Query genes in Cancer Gene Census.

    Args:
        gene_census: DataFrame containing Cancer Gene Census data
        genes: List of gene symbols to query

    Returns:
        Dictionary with gene information
    """
    results = {}

    for gene in genes:
        # Search for gene (case-insensitive)
        gene_upper = gene.upper()
        matches = gene_census[
            gene_census["Gene Symbol"].str.upper() == gene_upper
        ]

        if len(matches) == 0:
            results[gene] = {"found": False}
        else:
            # Convert first match to dict (keep all columns)
            gene_info = matches.iloc[0].to_dict()

            # Convert NaN to None for JSON serialization
            gene_info = {
                k: (None if pd.isna(v) else v)
                for k, v in gene_info.items()
            }

            results[gene] = {
                "found": True,
                **gene_info
            }

    return results


@app.command()
def main(
    gene: Optional[str] = typer.Option(
        None,
        "--gene",
        help="Single gene symbol to query",
    ),
    genes: Optional[List[str]] = typer.Option(
        None,
        "--genes",
        help="Multiple gene symbols to query (space-separated)",
    ),
    gene_list: Optional[Path] = typer.Option(
        None,
        "--gene-list",
        help="File containing gene symbols (one per line)",
        exists=True,
        readable=True,
    ),
    gene_census: Optional[Path] = typer.Option(
        None,
        "--gene-census",
        help="Path to cancer_gene_census.csv (default: data/cancer_gene_census.csv)",
        exists=True,
        readable=True,
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Output JSON file path (default: stdout)",
    ),
) -> None:
    """Query COSMIC Cancer Gene Census for gene information.

    Examples:
        # Query single gene
        python query_cosmic_genes.py --gene TP53

        # Query multiple genes
        python query_cosmic_genes.py --genes TP53 BRCA1 EGFR

        # Query from file
        python query_cosmic_genes.py --gene-list genes.txt

        # Save to file
        python query_cosmic_genes.py --gene TP53 --output result.json

        # Use custom Cancer Gene Census file
        python query_cosmic_genes.py --gene TP53 --gene-census /path/to/cancer_gene_census.csv
    """
    # Validate input
    if not gene and not genes and not gene_list:
        typer.echo("Error: Must specify --gene, --genes, or --gene-list", err=True)
        raise typer.Exit(1)

    # Collect genes to query
    query_genes_list = []

    if gene:
        query_genes_list.append(gene)

    if genes:
        query_genes_list.extend(genes)

    if gene_list:
        try:
            with open(gene_list) as f:
                file_genes = [line.strip() for line in f if line.strip()]
                query_genes_list.extend(file_genes)
        except Exception as e:
            typer.echo(f"Error reading gene list file: {e}", err=True)
            raise typer.Exit(1)

    # Remove duplicates while preserving order
    seen = set()
    query_genes_list = [
        g for g in query_genes_list
        if not (g in seen or seen.add(g))
    ]

    # Determine gene census file path
    if gene_census is None:
        gene_census = get_default_gene_census_path()

    # Load Cancer Gene Census
    try:
        gene_census_df = load_gene_census(gene_census)
    except (FileNotFoundError, ValueError) as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(1)

    # Query genes
    results = query_genes(gene_census_df, query_genes_list)

    # Count results
    found_count = sum(1 for r in results.values() if r.get("found", False))
    not_found_count = len(results) - found_count

    # Add summary
    output_data = {
        "summary": {
            "total_genes": len(results),
            "found_in_cancer_census": found_count,
            "not_found": not_found_count,
        },
        "genes": results,
    }

    # Output
    json_str = json.dumps(output_data, indent=2)

    if output:
        output.write_text(json_str)
        typer.echo(
            f"Results saved to: {output}\n"
            f"Total genes: {len(results)}, "
            f"Found: {found_count}, "
            f"Not found: {not_found_count}"
        )
    else:
        typer.echo(json_str)


if __name__ == "__main__":
    app()
