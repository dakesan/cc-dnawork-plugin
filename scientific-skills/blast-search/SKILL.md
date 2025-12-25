---
name: blast-search
description: "NCBI BLAST sequence similarity search. Web API (qblast) for standard searches, gget for quick queries. Find homologs, identify unknown sequences, validate primers. For genome mapping, use blat-search instead."
---

# BLAST Search

NCBI BLAST (Basic Local Alignment Search Tool) による配列類似性検索。

## When to Use This Skill

This skill should be used when:

- Finding similar sequences in databases (homolog search)
- Identifying unknown sequences
- Validating primer specificity
- Finding orthologs across species
- Checking sequence conservation

## When NOT to Use This Skill

- **Genome mapping / position identification** → Use `blat-search`
- **Reading/writing sequence files** → Use `sequence-io`
- **Multiple sequence alignment** → Use `msa-advanced`

## Tool Selection Guide

| Task | Tool | Speed |
|------|------|-------|
| Quick search (single sequence) | `gget.blast` | Fast |
| Standard search with options | `Bio.Blast.NCBIWWW` | Medium |
| Batch search / custom parameters | `Bio.Blast.NCBIWWW` | Medium |
| Large-scale local search | Local BLAST+ | TBD |

## Quick Start

### Installation

```bash
uv pip install biopython gget
```

### Quick BLAST with gget

```python
import gget

# Simple nucleotide BLAST
results = gget.blast("ATGCGATCGATCGATCG")

# Protein BLAST
results = gget.blast("MRHILKQWERTY", program="blastp")
```

### Standard BLAST with Biopython

```python
from Bio.Blast import NCBIWWW, NCBIXML

# Run BLAST
result_handle = NCBIWWW.qblast("blastn", "nt", "ATGCGATCGATCG")
blast_record = NCBIXML.read(result_handle)

# Display top hits
for alignment in blast_record.alignments[:5]:
    hsp = alignment.hsps[0]
    print(f"{alignment.title[:60]}...")
    print(f"  E-value: {hsp.expect}, Identity: {hsp.identities}/{hsp.align_length}")
```

## BLAST Programs

| Program | Query | Database | Use Case |
|---------|-------|----------|----------|
| `blastn` | DNA | DNA | Nucleotide similarity |
| `blastp` | Protein | Protein | Protein similarity |
| `blastx` | DNA (translated) | Protein | Find protein homologs |
| `tblastn` | Protein | DNA (translated) | Search genomic DNA |
| `tblastx` | DNA (translated) | DNA (translated) | Compare coding regions |

## Common Databases

### Nucleotide

| Database | Description |
|----------|-------------|
| `nt` | All GenBank + EMBL + DDBJ + PDB |
| `refseq_rna` | RefSeq RNA sequences |

### Protein

| Database | Description |
|----------|-------------|
| `nr` | Non-redundant protein sequences |
| `refseq_protein` | RefSeq protein sequences |
| `swissprot` | Curated UniProtKB/Swiss-Prot |
| `pdb` | Protein Data Bank sequences |

## Reference Documentation

### `references/ncbi_web.md`

- `NCBIWWW.qblast()` parameters
- Rate limiting and best practices
- Organism-specific searches (`entrez_query`)
- Saving and caching results

### `references/result_parsing.md`

- XML parsing with `NCBIXML`
- Extracting alignments and HSPs
- Calculating percent identity
- Filtering by E-value
- Tabular output format (outfmt 6)

## Common Patterns

### Filter by E-value

```python
E_VALUE_THRESH = 1e-10

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print(f"{alignment.accession}: E={hsp.expect}")
```

### Restrict to Organism

```python
result_handle = NCBIWWW.qblast(
    "blastn", "nt", sequence,
    entrez_query="Homo sapiens[Organism]"
)
```

### Get Percent Identity

```python
def percent_identity(hsp):
    return (hsp.identities / hsp.align_length) * 100
```

## Local BLAST (TBD)

Large-scale local BLAST support is planned. For now, use Web API.

## Best Practices

1. **Save results** - Don't re-run searches unnecessarily
2. **Set E-value threshold** - Default 10 is too permissive; use 0.001-0.01
3. **Use gget for quick searches** - Simpler API for single sequences
4. **Cache parsed data** - Avoid re-parsing large XML files
5. **Handle rate limits** - NCBI limits request frequency

## BLAST vs BLAT

| Aspect | BLAST | BLAT |
|--------|-------|------|
| Purpose | Similarity search | Genome mapping |
| Sensitivity | High | Medium |
| Speed | Medium | Very fast |
| Best for | Homolog search | Position finding |
