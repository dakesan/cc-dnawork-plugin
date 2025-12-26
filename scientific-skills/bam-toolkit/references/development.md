# BAM Toolkit Skill Plan

## Purpose

BAM/SAM/CRAM alignment file operations for WGS/WES analysis. Extract reads from specific regions, identify insertions/deletions, and calculate coverage statistics.

## Scope (Responsibility)

### Included ✅

- **BAM/SAM/CRAM file operations** - Read alignment file processing
- **Read extraction** - Extract reads from specific genomic regions
- **Indel detection** - Identify and extract insertion/deletion sequences
- **Coverage calculation** - Calculate coverage statistics for specified regions

### Excluded ❌

- **VCF/BCF operations** - Handled by vcf-toolkit
- **FASTA/FASTQ operations** - Handled by sequence-io
- **Variant calling** - Use external tools (GATK, bcftools, etc.)
- **Alignment (mapping)** - Use external tools (bwa, bowtie2, etc.)

## Design Principles

1. **pysam foundation** - Use pysam for all BAM/SAM/CRAM operations
2. **Region-based processing** - All operations require region specification
3. **JSON output** - Structured output for downstream analysis
4. **Typer CLI** - Consistent CLI pattern with vcf-toolkit

## Scripts Specification

### 1. extract_reads.py - Read Extraction

**Purpose**: Extract reads from BAM file for specific genomic region and output as BAM or JSON.

**Arguments**:
```
# Required
--bam PATH              Input BAM file path
--region TEXT           Genomic region (e.g., chr1:1000-2000)

# Output
--output PATH           Output file path (BAM or JSON based on extension)
--format TEXT           Output format: "bam" or "json" (default: bam)

# Filters (optional)
--min-mapq INT          Minimum mapping quality (default: 0)
--proper-pairs          Only properly paired reads
--no-duplicates         Exclude duplicate reads
```

**Output (BAM)**:
- Filtered BAM file with reads in specified region

**Output (JSON)**:
```json
{
  "region": "chr1:1000-2000",
  "total_reads": 150,
  "reads": [
    {
      "query_name": "read1",
      "reference_start": 1050,
      "reference_end": 1150,
      "sequence": "ATCG...",
      "quality": "IIII...",
      "mapping_quality": 60,
      "is_reverse": false,
      "cigar": "100M"
    }
  ]
}
```

**Use cases**:
```bash
# Extract reads as BAM
python extract_reads.py --bam input.bam --region chr1:1000-2000 --output reads.bam

# Extract reads as JSON
python extract_reads.py --bam input.bam --region chr1:1000-2000 --output reads.json --format json

# Extract high-quality properly-paired reads
python extract_reads.py \
  --bam input.bam \
  --region chr1:1000-2000 \
  --min-mapq 30 \
  --proper-pairs \
  --no-duplicates \
  --output filtered.bam
```

### 2. extract_indels.py - Insertion/Deletion Extraction

**Purpose**: Extract insertion and deletion sequences from reads in specified region.

**Arguments**:
```
# Required
--bam PATH              Input BAM file path
--region TEXT           Genomic region (e.g., chr1:1000-2000)

# Output
--output PATH           JSON output file path (default: stdout)

# Filters (optional)
--min-mapq INT          Minimum mapping quality (default: 20)
--min-indel-size INT    Minimum indel size in bp (default: 1)
```

**Output (JSON)**:
```json
{
  "region": "chr1:1000-2000",
  "summary": {
    "total_insertions": 5,
    "total_deletions": 3,
    "unique_insertions": 2,
    "unique_deletions": 2
  },
  "insertions": [
    {
      "position": 1050,
      "size": 3,
      "sequence": "ATG",
      "read_name": "read1",
      "mapping_quality": 60,
      "count": 3
    }
  ],
  "deletions": [
    {
      "position": 1100,
      "size": 2,
      "deleted_bases": "CG",
      "read_name": "read2",
      "mapping_quality": 55,
      "count": 2
    }
  ]
}
```

**Implementation notes**:
- Parse CIGAR string to identify I (insertion) and D (deletion) operations
- Extract sequences for insertions
- Track deleted bases from reference
- Group identical indels and count occurrences
- Filter by mapping quality and indel size

**Use cases**:
```bash
# Extract all indels
python extract_indels.py --bam input.bam --region chr1:1000-2000 --output indels.json

# Extract large indels only (>= 5bp)
python extract_indels.py \
  --bam input.bam \
  --region chr1:1000-2000 \
  --min-indel-size 5 \
  --output large_indels.json
```

### 3. calculate_coverage.py - Coverage Calculation

**Purpose**: Calculate coverage statistics for specified genomic region.

**Arguments**:
```
# Required
--bam PATH              Input BAM file path
--region TEXT           Genomic region (e.g., chr1:1000-2000)

# Output
--output PATH           JSON output file path (default: stdout)

# Options
--per-base              Include per-base coverage (default: False)
--min-mapq INT          Minimum mapping quality (default: 0)
--min-baseq INT         Minimum base quality (default: 0)
```

**Output (JSON)**:
```json
{
  "region": "chr1:1000-2000",
  "statistics": {
    "total_bases": 1000,
    "mean_coverage": 45.3,
    "median_coverage": 48,
    "min_coverage": 0,
    "max_coverage": 120,
    "bases_with_coverage": 995,
    "bases_without_coverage": 5,
    "percent_covered": 99.5
  },
  "per_base_coverage": [
    {"position": 1000, "coverage": 45},
    {"position": 1001, "coverage": 47},
    ...
  ]
}
```

**Implementation notes**:
- Use pysam pileup for per-base coverage
- Calculate summary statistics (mean, median, min, max)
- Optional per-base coverage output (can be large for long regions)
- Filter by mapping quality and base quality

**Use cases**:
```bash
# Calculate summary statistics only
python calculate_coverage.py --bam input.bam --region chr1:1000-2000 --output coverage.json

# Include per-base coverage
python calculate_coverage.py \
  --bam input.bam \
  --region chr1:1000-2000 \
  --per-base \
  --output detailed_coverage.json

# High-quality bases only
python calculate_coverage.py \
  --bam input.bam \
  --region chr1:1000-2000 \
  --min-mapq 30 \
  --min-baseq 20 \
  --output hq_coverage.json
```

## Implementation Plan

1. ✅ Create directory structure
2. ✅ Write plan.md (this file)
3. ✅ Implement scripts
   - [x] extract_reads.py (with paired-end mate information)
   - [x] extract_indels.py
   - [x] calculate_coverage.py (statistics only, no per-base)
4. ✅ Create SKILL.md
5. ✅ Test with sample BAM files
6. ⏳ Integrate with project (move to production)

## Testing Strategy

### Test Data Requirements
- Small test BAM file (< 1MB)
- BAI index file
- Known region with reads
- Reference genome (optional, for validation)

### Test Cases

**extract_reads.py:**
- Test 1: Extract reads as BAM
- Test 2: Extract reads as JSON
- Test 3: Apply quality filters
- Test 4: Proper pairs only

**extract_indels.py:**
- Test 1: Extract all indels
- Test 2: Filter by size
- Test 3: Verify insertion sequences
- Test 4: Verify deletion positions

**calculate_coverage.py:**
- Test 1: Summary statistics only
- Test 2: Quality filters (--min-mapq, --min-baseq)
- Test 3: Zero coverage regions
- Test 4: Large regions (long-read safe)

## Integration with Existing Skills

| Skill | Responsibility |
|-------|---------------|
| **bam-toolkit** (new) | BAM/SAM/CRAM read operations |
| **vcf-toolkit** | VCF/BCF variant operations |
| **sequence-io** | FASTA/FASTQ sequence operations |
| **blast-search** | BLAST sequence search |
| **blat-api-searching** | BLAT genome mapping |

## Questions / Decisions

1. ⏳ Should extract_reads.py support multiple output formats (BAM, SAM, CRAM, JSON)?
   - **Decision**: BAM and JSON only for simplicity

2. ✅ Should we include paired-end mate information?
   - **Decision**: Yes. All mate information included (mate_reference_name, mate_reference_start, mate_is_reverse, template_length)

3. ✅ Per-base coverage support?
   - **Decision**: Removed. Statistics only (mean, median, etc.) for long-read safety. Per-base output can be gigantic for long reads.

4. ⏳ Should we support multiple regions in a single call?
   - **Decision**: Single region per call, use shell loops for multiple regions

---

**Note**: This plan will be updated as implementation progresses.
