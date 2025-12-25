# BLAT Integration Workflow Examples

Practical examples for WGS analysis and bioinformatics workflows using BLAT.

## Example 1: Simple Sequence Search

Search a single sequence against human genome (hg38).

```python
from blat_integration import BlatRunner

runner = BlatRunner(mode="api", verbose=True)

# Search sequence
sequence = "ATGCGTACGATCGATCGATCGATCGATCG"
results = runner.search(
    sequence=sequence,
    genome="hg38",
    output_format="dataframe"
)

# Display results
print(f"Found {len(results)} hits")
print(results[["q_name", "t_name", "t_start", "t_end", "identity"]])

# Save to CSV
results.to_csv("blat_results.csv", index=False)
```

## Example 2: Batch FASTA Search (WGS Reads)

Search multiple reads from a FASTA file.

```python
from blat_integration import BlatRunner
import pandas as pd

runner = BlatRunner(mode="auto")

# Search batch of reads
results_df = runner.search(
    fasta_file="wgs_reads.fasta",
    genome="hg38",
    output_format="dataframe",
    min_identity=0.95,  # Filter for >95% identity
    sort_by="score"
)

# Analyze results
print(f"Total reads searched: {len(results_df)}")
print(f"Mean identity: {results_df['identity'].mean():.2f}%")
print(f"Unique chromosomes: {results_df['t_name'].nunique()}")

# Save filtered results
results_df.to_excel("blat_filtered.xlsx", index=False)
```

## Example 3: Find Insertion Sites

Search for insertion sequences to find their genomic location.

```python
from blat_integration import BlatRunner
import pandas as pd

runner = BlatRunner(mode="api")

# Search insertion sequences
insertion_sequences_file = "insertions.fasta"
results = runner.search(
    fasta_file=insertion_sequences_file,
    genome="hg38",
    min_identity=0.98,  # High stringency
    output_format="dataframe",
    max_results=1  # Top hit only per insertion
)

# Extract insertion locations
insertions = []
for _, row in results.iterrows():
    insertions.append({
        "insertion_id": row["q_name"],
        "chromosome": row["t_name"],
        "position": (row["t_start"] + row["t_end"]) // 2,
        "strand": row["strand"],
        "identity": row["identity"]
    })

insertions_df = pd.DataFrame(insertions)
insertions_df.to_csv("insertion_locations.csv", index=False)
print(insertions_df)
```

## Example 4: Compare Multiple Genomes

Search same sequences against multiple genome assemblies.

```python
from blat_integration import BlatRunner
import pandas as pd

runner = BlatRunner(mode="api")

sequence = "ATGCGTACGATCGATCGATCGATCGATCG"
genomes = ["hg38", "mm39", "rn7"]  # Human, Mouse, Rat

results_by_genome = {}
for genome in genomes:
    print(f"Searching {genome}...")
    results = runner.search(
        sequence=sequence,
        genome=genome,
        output_format="dataframe"
    )
    results_by_genome[genome] = results
    print(f"  Found {len(results)} hits")

# Compare orthologous regions
comparison = []
for genome, df in results_by_genome.items():
    if len(df) > 0:
        top_hit = df.iloc[0]
        comparison.append({
            "genome": genome,
            "chromosome": top_hit["t_name"],
            "position": top_hit["t_start"],
            "identity": top_hit["identity"]
        })

comparison_df = pd.DataFrame(comparison)
print("\nCross-species comparison:")
print(comparison_df)
comparison_df.to_csv("cross_species_hits.csv", index=False)
```

## Example 5: WGS Analysis Pipeline

Complete pipeline for WGS data analysis with BLAT.

```python
from blat_integration import BlatRunner
import pandas as pd
from Bio import SeqIO

class WGSBlatAnalysis:
    def __init__(self, fasta_file, genome="hg38"):
        self.fasta_file = fasta_file
        self.genome = genome
        self.runner = BlatRunner(mode="auto", verbose=True)

    def run_analysis(self):
        # Step 1: Search all reads
        print("Step 1: Searching reads...")
        results = self.runner.search(
            fasta_file=self.fasta_file,
            genome=self.genome,
            output_format="dataframe",
            min_identity=0.90
        )

        # Step 2: Filter and categorize
        print("Step 2: Filtering results...")
        high_quality = results[results["identity"] >= 0.98].copy()
        medium_quality = results[(results["identity"] >= 0.90) &
                                (results["identity"] < 0.98)].copy()

        # Step 3: Analyze coverage
        print("Step 3: Analyzing coverage...")
        coverage_stats = self._compute_coverage(high_quality)

        # Step 4: Detect insertions
        print("Step 4: Detecting insertions...")
        insertions = self._find_insertions(results)

        # Step 5: Generate report
        print("Step 5: Generating report...")
        report = self._generate_report(
            results, high_quality, medium_quality, coverage_stats, insertions
        )

        return {
            "all_results": results,
            "high_quality": high_quality,
            "medium_quality": medium_quality,
            "coverage": coverage_stats,
            "insertions": insertions,
            "report": report
        }

    def _compute_coverage(self, results):
        """Compute coverage statistics per chromosome."""
        coverage = {}
        for chrom in results["t_name"].unique():
            chrom_results = results[results["t_name"] == chrom]
            coverage[chrom] = {
                "hit_count": len(chrom_results),
                "mean_identity": chrom_results["identity"].mean(),
                "regions": chrom_results[["t_start", "t_end"]].values.tolist()
            }
        return coverage

    def _find_insertions(self, results):
        """Identify potential insertions from gaps in alignment."""
        insertions = []
        for _, row in results.iterrows():
            # Check for gaps between query start and target start
            if row["q_start"] > 100 or row["t_base_insert"] > 50:
                insertions.append({
                    "query": row["q_name"],
                    "target": row["t_name"],
                    "position": row["t_start"],
                    "gap_size": row["t_base_insert"],
                    "identity": row["identity"]
                })
        return insertions

    def _generate_report(self, all_results, high, medium, coverage, insertions):
        """Generate analysis report."""
        report = f"""
# WGS BLAT Analysis Report

## Summary Statistics
- Total sequences searched: {len(all_results)}
- High quality hits (≥98% identity): {len(high)}
- Medium quality hits (90-98% identity): {len(medium)}
- Potential insertions found: {len(insertions)}

## Coverage by Chromosome
"""
        for chrom, stats in coverage.items():
            report += f"\n- {chrom}: {stats['hit_count']} hits, "
            report += f"{stats['mean_identity']:.2f}% avg identity"

        if insertions:
            report += "\n\n## Detected Insertions\n"
            for ins in insertions[:10]:  # Show first 10
                report += f"- {ins['query']}: {ins['target']} at {ins['position']} "
                report += f"(gap={ins['gap_size']}bp)\n"

        return report


# Usage
analysis = WGSBlatAnalysis("wgs_reads.fasta", genome="hg38")
results = analysis.run_analysis()

# Save results
results["all_results"].to_excel("wgs_blat_all.xlsx", index=False)
results["high_quality"].to_excel("wgs_blat_high_quality.xlsx", index=False)
results["insertions_df"] = pd.DataFrame(results["insertions"])
results["insertions_df"].to_csv("wgs_insertions.csv", index=False)

# Save report
with open("wgs_blat_report.txt", "w") as f:
    f.write(results["report"])

print("Analysis complete!")
print(results["report"])
```

## Example 6: Integration with Other Tools

Use BLAT results with downstream analysis.

```python
from blat_integration import BlatRunner
from Bio import SeqIO, Seq
import pandas as pd

# Run BLAT
runner = BlatRunner(mode="auto")
results = runner.search(
    fasta_file="sequences.fasta",
    genome="hg38",
    output_format="dataframe"
)

# Step 1: Filter for best hits
best_hits = results.groupby("q_name").apply(
    lambda x: x.loc[x["score"].idxmax()]
).reset_index(drop=True)

# Step 2: Create BED file for IGV visualization
bed_data = []
for _, row in best_hits.iterrows():
    bed_data.append({
        "chrom": row["t_name"],
        "chromStart": row["t_start"],
        "chromEnd": row["t_end"],
        "name": row["q_name"],
        "score": int(row["score"]),
        "strand": row["strand"]
    })

bed_df = pd.DataFrame(bed_data)
bed_df.to_csv("blat_results.bed", sep="\t", header=False, index=False)

# Step 3: Create GenBank-annotated sequences with BLAT locations
output_file = "annotated_sequences.gb"
with open(output_file, "w") as out:
    for record in SeqIO.parse("sequences.fasta", "fasta"):
        # Find matching BLAT hits
        hits = best_hits[best_hits["q_name"] == record.id]

        # Add BLAT annotations
        for _, hit in hits.iterrows():
            feature = {
                "type": "BLAT_hit",
                "location": f"{hit['t_start']}..{hit['t_end']}",
                "qualifiers": {
                    "target_genome": hit["t_name"],
                    "identity": f"{hit['identity']:.1f}%",
                    "strand": hit["strand"]
                }
            }
            record.annotations["features"].append(feature)

        SeqIO.write(record, out, "genbank")

print(f"Results saved to:")
print(f"  - BED file: blat_results.bed")
print(f"  - GenBank: annotated_sequences.gb")
```

## Troubleshooting

### Issue: API Rate Limit

```python
runner = BlatRunner(
    mode="api",
    retry_on_limit=True,
    retry_delay=15
)

# Or switch to local mode
runner = BlatRunner(mode="local")
```

### Issue: BLAT Not Found (Local Mode)

```bash
# macOS
brew install blat

# Linux
sudo apt-get install blat-suite

# Download from UCSC
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod +x blat
```

### Issue: Large Batch Searches

```python
# Split into chunks for API mode
from itertools import islice

runner = BlatRunner(mode="api")
fasta_records = SeqIO.parse("large_file.fasta", "fasta")

all_results = []
for chunk in islice(zip(*[fasta_records]*100), None):
    # Process 100 sequences at a time
    ...
```

## References

- [UCSC BLAT Home](https://genome.ucsc.edu/cgi-bin/hgBlat)
- [BLAT FAQ](https://genome.ucsc.edu/FAQ/FAQblat.html)
- [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html)
