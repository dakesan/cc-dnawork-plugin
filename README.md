# cc-dnawork-plugin

A curated collection of **7 production-ready skills** for whole genome sequencing (WGS) analysis and bioinformatics workflows using Claude Code.

## Overview

This plugin provides essential tools for WGS/WES analysis pipelines:

- 🧬 DNA sequence I/O and manipulation
- 🔍 Sequence similarity search (BLAST, BLAT)
- 📊 Alignment file processing (BAM/SAM/CRAM)
- 🧪 Variant calling and annotation (VCF/BCF)
- 📸 Genomic visualization (IGV integration)
- 🎯 Cancer gene annotation (COSMIC database)

## Production Skills

### 1. sequence-io
Read, write, and manipulate biological sequences in multiple formats.

**Formats**: FASTA, GenBank, FASTQ
**Features**:
- Parse and write sequence files
- Format conversion
- Quality score handling (FASTQ)
- Sequence statistics and manipulation

**Use cases**: Sequence file conversion, quality filtering, format standardization

### 2. blast-search
NCBI BLAST sequence similarity search.

**Features**:
- BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX
- E-value filtering and alignment scoring
- Multiple database support
- Result parsing and annotation

**Use cases**: Homology search, sequence annotation, gene identification

### 3. blat-api-searching
Fast genome mapping with BLAT (BLAST-Like Alignment Tool).

**Features**:
- Rapid sequence-to-genome alignment
- UCSC Genome Browser integration
- High-speed mapping for large datasets
- Multiple genome builds support

**Use cases**: Quick genome mapping, PCR primer design, sequence localization

### 4. bam-toolkit
Comprehensive BAM/SAM/CRAM alignment file operations.

**Features**:
- Read/write BAM, SAM, CRAM files
- Read filtering and extraction
- Coverage calculation
- Alignment statistics
- Mate pair information analysis

**Use cases**: Post-alignment QC, read extraction, coverage analysis

### 5. vcf-toolkit
VCF/BCF variant file processing and analysis.

**Features**:
- Parse and write VCF/BCF files
- Variant filtering and annotation
- Multi-sample VCF operations
- Format conversion and validation

**Use cases**: Variant filtering, annotation, format conversion, quality control

### 6. igv-integration
Automated IGV (Integrative Genomics Viewer) snapshot generation.

**Features**:
- Batch IGV script generation
- Multiple BAM file visualization
- Single region or BED file input
- PNG screenshot output
- Automated IGV execution

**Use cases**: Visual variant validation, alignment inspection, publication figures

**Example**:
```bash
# Generate snapshots for multiple BAM files at specific region
python scripts/generate_igv_snapshots.py \
  --genome hg38 \
  --bam sample1.bam sample2.bam sample3.bam \
  --region chr17:7577001-7578000 \
  --output-dir ./snapshots
```

### 7. cosmic-toolkit
COSMIC Cancer Gene Census database annotation.

**Features**:
- Gene lookup in Cancer Gene Census
- Dynamic TSV-to-JSON conversion
- Batch gene queries
- All COSMIC metadata preserved

**Use cases**: Cancer gene validation, somatic mutation annotation, gene prioritization

**Example**:
```bash
# Query cancer genes
python scripts/query_cosmic_genes.py TP53 KRAS EGFR
```

**Data setup**: Download `cancer_gene_census.csv` from [COSMIC](https://cancer.sanger.ac.uk/cosmic) and place in `data/` directory.

## Installation

### Prerequisites
- **Claude Code**: Latest version
- **Python**: 3.9+ (3.12+ recommended)
- **uv**: Python package manager

### Install Claude Code

**macOS:**
```bash
curl -fsSL https://claude.ai/install.sh | bash
```

**Windows (PowerShell):**
```powershell
irm https://claude.ai/install.ps1 | iex
```

### Install uv

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### Install Plugin

1. Add marketplace:
   ```bash
   /plugin marketplace add dakesan/cc-dnawork-plugin
   ```

2. Install plugin:
   - Open Claude Code
   - Run `/plugin list`
   - Select **dnawork-skills**
   - Click **Install**

## Quick Start

### Example 1: WGS Variant Analysis Pipeline

```
Analyze this WGS dataset:
1. Check FASTQ quality and statistics (sequence-io)
2. Verify BAM alignment quality and coverage (bam-toolkit)
3. Parse and filter VCF variants (vcf-toolkit)
4. Annotate cancer-related genes (cosmic-toolkit)
5. Generate IGV snapshots for top variants (igv-integration)
6. Search homologous sequences (blast-search)
```

### Example 2: Cancer Somatic Variant Annotation

```
Process this tumor-normal VCF:
1. Extract somatic variants (vcf-toolkit)
2. Check genes in COSMIC Cancer Gene Census (cosmic-toolkit)
3. Generate IGV screenshots for driver mutations (igv-integration)
4. Create summary report with annotations
```

### Example 3: Sequence Validation

```
Validate this amplicon sequence:
1. Convert GenBank to FASTA (sequence-io)
2. Map to reference genome (blat-api-searching)
3. Verify alignment coverage (bam-toolkit)
4. Search for homologs in NCBI (blast-search)
```

## Architecture

```
cc-dnawork-plugin/
├── .claude-plugin/
│   └── marketplace.json          # Plugin metadata
├── scientific-skills/
│   ├── sequence-io/              # FASTA/GenBank/FASTQ I/O
│   ├── blast-search/             # NCBI BLAST
│   ├── blat-api-searching/       # BLAT genome mapping
│   ├── bam-toolkit/              # BAM/SAM/CRAM operations
│   ├── vcf-toolkit/              # VCF/BCF operations
│   ├── igv-integration/          # IGV automation
│   └── cosmic-toolkit/           # COSMIC annotation
├── STRUCTURE.md                  # Detailed structure
├── README.md                     # This file
└── LICENSE                       # MIT License
```

Each skill contains:
- `SKILL.md` - Complete documentation
- `references/` - Additional documentation
- `scripts/` - Python scripts (where applicable)
- `data/` - User-provided data directory (where applicable)

## Requirements

### System Requirements
- **OS**: macOS, Linux, or Windows with WSL2
- **Python**: 3.9+ (3.12+ recommended for best performance)
- **Memory**: 8GB+ RAM (16GB+ recommended for large datasets)
- **Disk**: 1GB+ for plugin + data storage

### Python Packages

Automatically installed when using skills:
- `biopython` - Sequence I/O and manipulation
- `pysam` - BAM/SAM/CRAM/VCF operations
- `pandas` - Data processing (COSMIC toolkit)
- `requests` - API access (BLAST, BLAT)

### External Tools

**Optional** (for full functionality):
- **IGV**: For igv-integration skill
  - Download from [IGV website](https://igv.org/doc/desktop/)
  - Ensure `igv.sh` (Linux/macOS) or `igv.bat` (Windows) is in PATH

### Database Access

**COSMIC toolkit** requires user-downloaded data:
1. Register at [COSMIC](https://cancer.sanger.ac.uk/cosmic)
2. Download Cancer Gene Census CSV
3. Place in `scientific-skills/cosmic-toolkit/data/cancer_gene_census.csv`

## Typical WGS Analysis Workflow

```
Raw Reads (FASTQ)
    ↓ [sequence-io]
Quality Check & Statistics
    ↓
Alignment (External: BWA/Bowtie2)
    ↓ [bam-toolkit]
BAM Quality Control & Coverage Analysis
    ↓
Variant Calling (External: GATK/FreeBayes)
    ↓ [vcf-toolkit]
VCF Filtering & Annotation
    ↓ [cosmic-toolkit]
Cancer Gene Annotation
    ↓ [igv-integration]
Visual Validation (IGV Screenshots)
    ↓ [blast-search]
Homology Search & Validation
    ↓
Final Report
```

## Troubleshooting

### Skills Not Loading
```bash
# Check installation
/plugin list

# Reinstall
/plugin marketplace remove dakesan/cc-dnawork-plugin
/plugin marketplace add dakesan/cc-dnawork-plugin
```

### Missing Python Dependencies
```bash
# Install specific package
uv pip install package-name

# List installed packages
uv pip list
```

### IGV Not Found
- Ensure IGV is installed and `igv.sh`/`igv.bat` is in PATH
- Test: `igv.sh --version` (Linux/macOS) or `igv.bat` (Windows)
- Alternative: Specify IGV path with `--igv-path` flag

### COSMIC Data Missing
- Download Cancer Gene Census from COSMIC website
- Place CSV file in `scientific-skills/cosmic-toolkit/data/`
- See `cosmic-toolkit/data/README.md` for detailed instructions

### Large File Processing
- Use BAM index files (`.bai`) for faster access
- Enable streaming mode for large VCF files
- Consider splitting large datasets into chunks

## Citation

If you use this plugin in your research:

```bibtex
@software{cc_dnawork_plugin_2025,
  author = {Odake, Hiroyuki},
  title = {cc-dnawork-plugin: WGS Analysis Skills for Claude Code},
  year = {2025},
  url = {https://github.com/dakesan/cc-dnawork-plugin}
}
```

## Contributing

Contributions are welcome! Please:
- Report bugs via GitHub Issues
- Submit pull requests for improvements
- Follow existing code style and documentation format

## Related Resources

- [Claude Code Documentation](https://docs.claude.com/)
- [IGV Documentation](https://igv.org/doc/desktop/)
- [COSMIC Database](https://cancer.sanger.ac.uk/cosmic)
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/)
- [UCSC BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat)

## License

MIT License - See LICENSE file for details.

---

**Designed for WGS/WES analysis workflows** 🧬
