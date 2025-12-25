# PSL Format Specification

PSL (PSLX/PSL eXtended) is the output format for BLAT. This document describes the PSL format fields.

## Overview

PSL format is tab-separated with 21+ fields per line. The first 5 lines are headers.

## Header Lines

The first 5 lines contain:
1. Browser position header (compatibility)
2. Browser position details
3. Blank line
4. Column header with field names
5. Blank line

## PSL Fields (21 columns)

| Field | Type | Description |
|-------|------|-------------|
| **match** | int | Number of bases that match |
| **mismatch** | int | Number of bases that don't match |
| **rep_match** | int | Number of bases that match but are in repeats |
| **n_count** | int | Number of 'N' bases |
| **q_num_insert** | int | Number of inserts in query |
| **q_base_insert** | int | Number of bases inserted in query |
| **t_num_insert** | int | Number of inserts in target |
| **t_base_insert** | int | Number of bases inserted in target |
| **strand** | char | '+' or '-' for strand |
| **q_name** | string | Query sequence name |
| **q_size** | int | Query sequence size |
| **q_start** | int | Alignment start position in query |
| **q_end** | int | Alignment end position in query |
| **t_name** | string | Target sequence name |
| **t_size** | int | Target sequence size |
| **t_start** | int | Alignment start position in target |
| **t_end** | int | Alignment end position in target |
| **block_count** | int | Number of blocks in the alignment |
| **blockSizes** | string | Comma-separated list of block sizes |
| **qStarts** | string | Comma-separated start positions in query |
| **tStarts** | string | Comma-separated start positions in target |

## Quality Metrics Derived from PSL

### Identity Percentage
```
identity (%) = 100 * match / (match + mismatch)
```

### Coverage (Query)
```
coverage (%) = 100 * (q_end - q_start) / q_size
```

### Coverage (Target)
```
coverage (%) = 100 * (t_end - t_start) / t_size
```

### UCSC Score
```
score = match - mismatch - q_base_insert - t_base_insert
```

### Bit Score (approximate)
```
bit_score = (match - mismatch) * log2(E-value) - t_base_insert - q_base_insert
```

## Example PSL Line

```
50  0  0  0  0  0  0  0  +  seq1  100  0  50  chr1  248956422  1000  1050  1  50,  0,  1000,
```

Interpretation:
- 50 matching bases
- 0 mismatches
- No repeats, no N's, no insertions
- Plus strand
- Query: seq1 (100bp total, aligned 0-50)
- Target: chr1 (aligned 1000-1050)

## PSL Variants

### Standard PSL
- 21 fields
- No sequence information

### PSLX (Extended)
- 23+ fields
- Includes actual sequences
- Fields 22-23: query sequence, target sequence

## Tips for PSL Analysis

1. **Filtering by quality:**
   - High identity: identity > 95%
   - High coverage: coverage > 80%
   - High score: score > 50

2. **Handling multiple hits:**
   - Sort by score (descending)
   - Select top hit per query
   - Filter by E-value equivalent

3. **Chain detection:**
   - Use multiple blocks for spliced alignment
   - Check blockCount > 1
   - Verify block continuity

## References

- [BLAT FAQ - Output Format](https://genome.ucsc.edu/FAQ/FAQblat.html#output)
- [UCSC Wiki - PSL Format](http://genome.ucsc.edu/FAQ/FAQformat.html)
