---
name: blat-integration
description: "BLAT (BLAST-Like Alignment Tool) integration for sequence searching via UCSC REST API and local command-line tool. Support for both small-scale (online) and large-scale (local) genomic searches with automatic mode selection."
---

# BLAT Integration: Sequence Alignment Tool

## Overview

BLAT (BLAST-Like Alignment Tool) is a specialized alignment tool optimized for finding the longest perfect and near-perfect matches between DNA sequences. This skill provides both **UCSC REST API access** and **local command-line execution**, automatically selecting the optimal mode based on your search scale.

BLAT is particularly useful for:
- Rapid genomic sequence searches
- Finding similar genomic regions
- Validating sequence assembly
- Large-scale whole genome searches

## When to Use This Skill

Use this skill when:

- Searching DNA sequences against genomic databases
- Validating WGS (Whole Genome Sequencing) results
- Finding similar sequences in a reference genome
- Performing BLAT searches smaller than 5,000 hits/day (use API)
- Performing large-scale BLAT searches (use local tool)
- Comparing sequences across multiple genomic assemblies
- Identifying insertion sites, binding regions, or homologous sequences
- Automating genome-wide sequence mapping workflows

## Core Capabilities

### 1. **UCSC REST API Mode** (Small-scale searches)
- Direct REST API access to UCSC Genome Browser BLAT
- No installation required
- Built-in rate limiting (1 hit/15 sec, 5,000 hits/day)
- Suitable for validation and small batch queries
- Automatic retry with backoff on rate limit

### 2. **Local BLAT Mode** (Large-scale searches)
- Execute local `blat` command-line tool
- Unlimited searches (no rate limits)
- High-speed performance for large datasets
- Support for custom 2bit genome databases
- PSL output parsing and result conversion

### 3. **Automatic Mode Selection**
- Smart detection: API for <1000 sequences, local for larger
- Seamless fallback between modes
- Configurable threshold for mode selection

### 4. **Result Parsing & Format Conversion**
- PSL format parsing (BLAT standard output)
- Conversion to DataFrame (pandas), TSV, JSON, BED
- Statistical analysis (hit counts, identity %, E-values)
- Integration with BioPython for sequence analysis

### 5. **Quality Metrics & Filtering**
- Hit quality assessment (coverage, identity)
- Configurable filtering (min identity, max E-value)
- Sorting by score, coverage, or identity
- Result annotation with genomic coordinates

## Installation and Setup

### Prerequisites

- Python 3.9+
- pandas
- requests
- biopython

### Install Python Dependencies

```bash
uv pip install pandas requests biopython
```

### Optional: Install Local BLAT Tool

For large-scale searches, install the local BLAT tool:

**macOS (via Homebrew):**
```bash
brew install blat
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install blat-suite
```

**Manual Installation:**
Download from UCSC: https://genome.ucsc.edu/FAQ/FAQblat.html#blat

```bash
# Download BLAT binary
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat

# Make executable
chmod +x blat

# Move to PATH (e.g., /usr/local/bin)
sudo mv blat /usr/local/bin/
```

### Download Reference Genomes (2bit format)

For local BLAT searches, download 2bit genome files:

```bash
# Download hg38 (human)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

# Or other assemblies
# http://hgdownload.soe.ucsc.edu/goldenPath/[assembly]/bigZips/[assembly].2bit
```

## Usage Examples

### Example 1: Small-scale Search (API Mode)

```python
from blat_integration import BlatRunner

runner = BlatRunner()

# Single sequence search
sequence = "ATGCGTACGATCGATCG"
results = runner.search(sequence, genome="hg38", mode="api")

print(f"Found {len(results)} hits")
for hit in results:
    print(f"  {hit['chromosome']}:{hit['start']}-{hit['end']} "
          f"(identity: {hit['identity']:.1f}%)")

# Save results to CSV
results_df = runner.to_dataframe(results)
results_df.to_csv("blat_results.csv", index=False)
```

### Example 2: Large-scale Search (Local Mode)

```python
from blat_integration import BlatRunner

runner = BlatRunner()

# Large batch search using local BLAT
results = runner.search_batch(
    fasta_file="query_sequences.fasta",
    database="hg38.2bit",
    mode="local",
    output_format="psl"
)

# Parse and filter results
filtered = runner.filter_results(
    results,
    min_identity=0.95,
    min_coverage=0.90
)

print(f"Filtered results: {len(filtered)} hits")
```

### Example 3: Automatic Mode Selection

```python
from blat_integration import BlatRunner

runner = BlatRunner()

# Automatically selects mode based on input size
results = runner.search(
    fasta_file="sequences.fasta",
    genome="hg38",
    mode="auto"  # Selects API or local automatically
)
```

### Example 4: WGS Analysis Pipeline

```python
from blat_integration import BlatRunner, PSLParser
import pandas as pd

# Initialize BLAT runner
blat = BlatRunner(mode="auto", max_api_daily=5000)

# Step 1: Search reads against reference
results = blat.search(
    fasta_file="wgs_reads.fasta",
    database="hg38.2bit",
    output_format="psl"
)

# Step 2: Parse PSL results
psl_data = PSLParser.parse("output.psl")

# Step 3: Analyze coverage
coverage = blat.analyze_coverage(psl_data)
print(f"Coverage: {coverage['mean_identity']:.2f}% identity")

# Step 4: Find insertions
insertions = blat.find_insertions(psl_data, min_gap=50)

# Step 5: Export results
variants_df = pd.DataFrame(insertions)
variants_df.to_excel("insertion_variants.xlsx", index=False)

# Step 6: Generate report
report = blat.generate_report(
    results,
    output_file="blat_analysis_report.txt"
)
```

### Example 5: Comparing Multiple Genomes

```python
from blat_integration import BlatRunner
import pandas as pd

runner = BlatRunner()

sequence = "ATGCGTACGATCGATCG"
genomes = ["hg38", "mm39", "rn7"]

results = {}
for genome in genomes:
    results[genome] = runner.search(sequence, genome=genome, mode="api")

# Compare results across genomes
comparison = runner.compare_genomes(results)
comparison.to_csv("cross_genome_comparison.csv", index=False)
```

## Key Parameters

### Search Parameters

```python
runner.search(
    sequence=None,          # Single sequence string
    fasta_file=None,        # Path to FASTA file
    genome="hg38",          # Genome assembly (hg38, mm39, rn7, etc.)
    database=None,          # Path to local 2bit database
    mode="auto",            # "api", "local", or "auto"

    # Search parameters
    min_identity=0.0,       # Minimum % identity (0-100)
    max_gap_size=None,      # Maximum gap size
    min_tile_size=11,       # Minimum perfect match size

    # Output options
    output_format="psl",    # "psl", "json", "bed", "dataframe"
    output_file=None,       # Save to file

    # Filtering
    max_results=None,       # Limit number of results
    sort_by="score",        # "score", "identity", "coverage"
)
```

### API Rate Limiting

```python
runner = BlatRunner(
    mode="api",
    max_api_daily=5000,     # Daily hit limit (UCSC limit)
    retry_on_limit=True,    # Auto-retry on rate limit
    retry_delay=15          # Seconds between retries
)
```

## Installation and Troubleshooting

### Check BLAT Installation

```python
from blat_integration import BlatRunner

runner = BlatRunner()
runner.check_installation()
# Output: BLAT version, installation path, available databases
```

### Check API Connectivity

```python
from blat_integration import BlatRunner

runner = BlatRunner(mode="api")
runner.test_api_connection()
# Output: API status, rate limit status
```

### Common Issues

**Issue 1: "blat command not found" (local mode)**
- Solution: Install BLAT (see Installation section)
- Verify: `which blat`

**Issue 2: "API rate limit exceeded"**
- Solution: Wait 15 seconds, or use local mode for batch queries
- Alternative: Switch to `mode="local"`

**Issue 3: "Database file not found"**
- Solution: Download 2bit file for your genome
- Location: http://hgdownload.soe.ucsc.edu/goldenPath/

## Performance Characteristics

### API Mode
- **Speed**: ~2-3 seconds per sequence (with rate limit)
- **Throughput**: 5,000 hits/day (UCSC limit)
- **Suitable for**: Validation, small batches, exploration
- **Cost**: Free, within rate limits

### Local Mode
- **Speed**: <1 second per sequence (after database load)
- **Throughput**: Unlimited
- **Suitable for**: Large batches, WGS analysis, production pipelines
- **Cost**: Disk space for 2bit database (~1-3 GB per genome)

## References

For detailed information, see:
- `references/blat_installation.md` - Installation guide for all platforms
- `references/psl_format.md` - PSL output format specification
- `references/workflow_examples.md` - Additional workflow examples
- `references/ucsc_api.md` - UCSC REST API documentation

## Scripts

Pre-built utility scripts in `scripts/`:
- `blat_runner.py` - Main BlatRunner class
- `psl_parser.py` - PSL format parser
- `install_blat.sh` - Automated installation script
- `compare_blast_blat.py` - Comparison with BLAST results

## Related Skills

- **biopython** - Sequence manipulation and parsing
- **pysam** - VCF/SAM format handling for genomic coordinates
- **pubmed-database** - Find literature about detected sequences
- **ensembl-database** - Annotate BLAT hits with gene information
- **igv-integration** - Visualize BLAT hits in IGV

## Citation

If you use BLAT in your research, please cite:

```bibtex
@article{Kent2002BLAT,
  author = {Kent, W. James},
  title = {BLAT—the BLAST-like alignment tool},
  journal = {Genome Research},
  year = {2002},
  volume = {12},
  pages = {656--664}
}
```

Also cite UCSC Genome Browser:

```bibtex
@article{UCSC2023,
  author = {Haeussler, Max and others},
  title = {The UCSC Genome Browser database: 2023 update},
  journal = {Nucleic Acids Research},
  year = {2023},
  volume = {51},
  number = {D1}
}
```

---

For more information, visit:
- [UCSC BLAT Home](https://genome.ucsc.edu/cgi-bin/hgBlat)
- [UCSC REST API Documentation](https://genome.ucsc.edu/goldenPath/help/api.html)
- [BLAT FAQ](https://genome.ucsc.edu/FAQ/FAQblat.html)
