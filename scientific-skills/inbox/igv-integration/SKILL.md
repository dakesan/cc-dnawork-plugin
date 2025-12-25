---
name: igv-integration
description: "IGV (Integrative Genomics Viewer) command-line integration with automatic batch script generation, headless screenshot capture, and Markdown report generation for WGS visualization."
---

# IGV Integration: Automated Genome Visualization & Reporting

## Overview

IGV Integration enables:
- **Batch script generation**: Automate IGV session creation
- **Headless screenshot capture**: Extract images without GUI
- **Region-based visualization**: Automated viewing of genomic regions
- **Report integration**: Generate Markdown reports with embedded visualizations
- **Multi-region processing**: Handle large genomic regions efficiently

## When to Use This Skill

Use this skill when:

- Visualizing WGS analysis results (BLAT hits, variants, insertions)
- Creating publication-ready genome browser screenshots
- Automating visualization of multiple genomic regions
- Generating comprehensive analysis reports with embedded images
- Validating sequence alignment in genomic context
- Creating presentations with IGV snapshots
- Batch processing many genomic regions

## Core Capabilities

### 1. **IGV Batch Script Generation**
- Create .bat/.sh scripts for IGV automation
- Specify genomic regions (chrom:start-end)
- Load BAM, VCF, BED, FASTA files
- Set zoom levels and color schemes

### 2. **Screenshot Capture**
- Headless IGV operation (no display needed)
- Automatic image export (PNG/SVG)
- Configurable resolution and size
- Automatic naming and organization

### 3. **Region Management**
- BED file support for region lists
- Coordinate validation
- Automatic coordinate conversion (0-based ↔ 1-based)
- Region annotation with metadata

### 4. **Markdown Integration**
- Generate Markdown reports with embedded images
- Table generation for region metadata
- Automatic image linking
- Multi-image page layouts

### 5. **Session Management**
- Save and load IGV sessions
- Reproducible visualizations
- Session templates for common analyses
- Project-based organization

## Installation and Setup

### Prerequisites

- Java Runtime Environment (JRE) 8+
- IGV application

### Install IGV

**macOS (Homebrew):**
```bash
brew install --cask igv
```

**Linux:**
```bash
# Download from Broad Institute
cd ~/tools
wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.0_withJava.zip
unzip IGV_Linux_2.16.0_withJava.zip
```

**Windows / Manual:**
- Download from: https://software.broadinstitute.org/software/igv/download
- Install Java if not present

### Verify Installation

```bash
igv.sh --version
# or
igv --version
```

## Usage Examples

### Example 1: Single Region Visualization

```python
from igv_integration import IGVBatchGenerator, ScreenshotProcessor

# Generate batch script
gen = IGVBatchGenerator()
gen.add_region("chr1", 1000000, 1100000)  # 100kb region
gen.set_zoom_level(5)

# Generate and run
gen.generate_batch_script("igv_session.txt")
IGVBatchGenerator.run_batch("igv_session.txt")

# Capture screenshot
processor = ScreenshotProcessor()
screenshots = processor.capture("igv_session.txt")

print(f"Captured {len(screenshots)} images")
for image in screenshots:
    print(f"  - {image}")
```

### Example 2: WGS Analysis Report

```python
from igv_integration import IGVBatchGenerator, ScreenshotProcessor, ReportGenerator
import pandas as pd

# Load BLAT results
blat_results = pd.read_csv("blat_results.csv")

# Create IGV visualizations for top hits
gen = IGVBatchGenerator()
gen.load_files({
    "bam": "alignment.bam",
    "vcf": "variants.vcf",
    "bed": "regions.bed"
})

# Add BLAT hit regions
for idx, row in blat_results.head(5).iterrows():
    gen.add_region(
        row["t_name"],
        int(row["t_start"]) - 5000,  # 5kb flanking
        int(row["t_end"]) + 5000,
        label=row["q_name"]
    )

# Generate batch script and capture
batch_file = "igv_blat_hits.txt"
gen.generate_batch_script(batch_file)
IGVBatchGenerator.run_batch(batch_file)

# Capture screenshots
processor = ScreenshotProcessor(resolution=(1920, 1080))
screenshots = processor.capture(batch_file)

# Generate report
report_gen = ReportGenerator()
report = report_gen.create_report(
    title="WGS BLAT Analysis - IGV Visualization",
    sections=[
        {
            "title": "Top BLAT Hits",
            "content": blat_results[["q_name", "t_name", "t_start", "t_end"]].to_html(),
            "images": screenshots
        }
    ]
)

report_gen.save(report, "wgs_igv_report.md")
```

### Example 3: Batch Region Processing

```python
from igv_integration import IGVBatchGenerator, ScreenshotProcessor
from Bio import SeqIO
import pandas as pd

# Process multiple insertion sequences
insertions_file = "insertions.fasta"
gen = IGVBatchGenerator()
gen.load_files({"bam": "wgs.bam"})

# Add insertion regions from BLAT results
blat_hits = pd.read_csv("blat_hits.csv")

for _, row in blat_hits.iterrows():
    gen.add_region(
        row["chromosome"],
        row["start"],
        row["end"],
        label=row["insertion_id"]
    )

# Generate batch script
gen.generate_batch_script("igv_insertions.txt")

# Run IGV and capture screenshots
IGVBatchGenerator.run_batch("igv_insertions.txt")
processor = ScreenshotProcessor()
screenshots = processor.capture("igv_insertions.txt")

# Organize output
for i, screenshot in enumerate(screenshots):
    insertion_id = blat_hits.iloc[i]["insertion_id"]
    print(f"Saved: {insertion_id} → {screenshot}")
```

### Example 4: Custom Styling & Colors

```python
from igv_integration import IGVBatchGenerator

gen = IGVBatchGenerator()

# Configure display settings
gen.set_zoom_level(3)
gen.set_genome_assembly("hg38")
gen.set_color_scheme({
    "insertion": "#FF0000",      # Red
    "deletion": "#0000FF",       # Blue
    "snp": "#FFA500"             # Orange
})

# Add regions with annotations
gen.add_region(
    "chr1",
    1000000,
    1100000,
    label="Insertion Site 1",
    color="#FF0000",
    feature_type="insertion"
)

gen.add_region(
    "chr2",
    2000000,
    2100000,
    label="Deletion Site 2",
    color="#0000FF",
    feature_type="deletion"
)

# Generate with custom viewport
gen.set_viewport_size(1920, 1080)
gen.generate_batch_script("custom_igv.txt")
```

## Key Parameters

```python
IGVBatchGenerator(
    genome_assembly="hg38",    # hg38, hg19, mm39, rn7, etc.
    java_heap="2g",            # Java heap size
    headless=True,             # Run without GUI
)

# Add regions
gen.add_region(
    chromosome,                # e.g., "chr1"
    start_position,           # 0-based or 1-based (auto-detect)
    end_position,
    label=None,               # Region label
    color=None,               # RGB or hex color
    feature_type=None         # Annotation type
)

# Screenshot parameters
processor = ScreenshotProcessor(
    resolution=(1920, 1080),  # Width x Height
    format="png",             # png or svg
    quality=90                # JPEG quality
)
```

## IGV Batch Script Format

```
# IGV Batch Script Format
genome hg38                          # Load genome assembly
load alignments.bam                  # Load BAM file
load variants.vcf                    # Load VCF file
goto chr1:1000000-1100000            # Navigate to region
snapshot chr1_1M.png                 # Take screenshot
goto chr2:2000000-2100000
snapshot chr2_2M.png
exit                                 # Exit IGV
```

## Performance Tips

1. **Large Alignments**: Pre-filter BAM files for regions of interest
2. **Memory**: Adjust Java heap size for large genomes
3. **Batch Processing**: Group regions geographically to minimize seeks
4. **Resolution**: Use lower resolution for speed, higher for publication

## Troubleshooting

### Problem: IGV not found
```bash
# Ensure IGV is installed and in PATH
which igv
# or
brew list igv
```

### Problem: Java error
```bash
# Update Java
java -version
# Should be 8 or newer

# Or specify Java path
export JAVA_HOME=/usr/libexec/java_home -v 11
```

### Problem: Memory errors
```python
gen = IGVBatchGenerator(java_heap="4g")  # Increase heap to 4GB
```

## Related Skills

- **blat-integration** - Find regions to visualize
- **msa-advanced** - Align sequences before viewing
- **clinical-reports** - Integrate IGV images into reports

## References

- [IGV User Guide](http://software.broadinstitute.org/software/igv/UserGuide)
- [IGV Batch Commands](http://software.broadinstitute.org/software/igv/batch)
- [IGV Download](https://software.broadinstitute.org/software/igv/download)

See `references/` for detailed documentation.
