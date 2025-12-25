# NCBI Web BLAST (Bio.Blast.NCBIWWW)

## Overview

`Bio.Blast.NCBIWWW.qblast()` submits sequences to NCBI's online BLAST service. No installation required, but subject to rate limits.

## Basic Usage

```python
from Bio.Blast import NCBIWWW

# Run BLAST search
result_handle = NCBIWWW.qblast(
    program="blastn",           # BLAST program
    database="nt",              # Database to search
    sequence="ATGCGATCGATCG"    # Query sequence
)

# Save results to file
with open("blast_results.xml", "w") as f:
    f.write(result_handle.read())
result_handle.close()
```

## qblast() Parameters

### Required Parameters

```python
NCBIWWW.qblast(
    program,    # "blastn", "blastp", "blastx", "tblastn", "tblastx"
    database,   # "nt", "nr", "refseq_rna", "swissprot", etc.
    sequence    # Query sequence string, FASTA, or accession
)
```

### Common Optional Parameters

```python
NCBIWWW.qblast(
    program="blastn",
    database="nt",
    sequence=query,

    # Filtering
    expect=0.001,              # E-value threshold (default: 10)
    hitlist_size=50,           # Max hits to return (default: 50)

    # Alignment options
    word_size=11,              # Word size for initial match
    gapcosts="5 2",            # Gap costs (open, extend)

    # Output
    format_type="XML",         # Output format (default: XML)
    alignments=50,             # Number of alignments to show
    descriptions=50,           # Number of descriptions to show
)
```

### Organism-Specific Search

```python
# Restrict to specific organism
result_handle = NCBIWWW.qblast(
    "blastn", "nt", sequence,
    entrez_query="Homo sapiens[Organism]"
)

# Exclude organism
result_handle = NCBIWWW.qblast(
    "blastn", "nt", sequence,
    entrez_query="NOT Homo sapiens[Organism]"
)

# Multiple organisms
result_handle = NCBIWWW.qblast(
    "blastn", "nt", sequence,
    entrez_query="(Homo sapiens[Organism] OR Mus musculus[Organism])"
)
```

## Input Formats

### Sequence String

```python
result_handle = NCBIWWW.qblast("blastn", "nt", "ATGCGATCGATCG")
```

### FASTA String

```python
fasta_string = ">query\nATGCGATCGATCG"
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)
```

### From SeqRecord

```python
from Bio import SeqIO

record = SeqIO.read("sequence.fasta", "fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", str(record.seq))
```

### Accession Number

```python
# GenBank accession
result_handle = NCBIWWW.qblast("blastn", "nt", "EU490707")
```

## gget Alternative

For quick searches, `gget.blast` provides a simpler interface:

```python
import gget

# Nucleotide BLAST
results = gget.blast("ATGCGATCGATCG")

# Protein BLAST
results = gget.blast("MRHILKQWERTY", program="blastp")
```

Returns pandas DataFrame with: accession, description, e_value, percent_identity, etc.

## Rate Limiting

NCBI limits request frequency.

```python
import time

sequences = ["ATGC...", "GCTA...", "TACG..."]

for i, seq in enumerate(sequences):
    result = NCBIWWW.qblast("blastn", "nt", seq)
    with open(f"result_{i}.xml", "w") as f:
        f.write(result.read())
    result.close()

    if i < len(sequences) - 1:
        time.sleep(3)  # Wait between requests
```

## Error Handling

```python
from Bio.Blast import NCBIWWW
from urllib.error import HTTPError

try:
    result = NCBIWWW.qblast("blastn", "nt", sequence)
except HTTPError as e:
    if e.code == 429:
        print("Rate limited. Wait and retry.")
    else:
        raise
```

## Caching Results

```python
import os

cache_file = "blast_cache.xml"

if os.path.exists(cache_file):
    # Use cached result
    with open(cache_file) as f:
        from Bio.Blast import NCBIXML
        blast_record = NCBIXML.read(f)
else:
    # Run BLAST and cache
    result = NCBIWWW.qblast("blastn", "nt", sequence)
    with open(cache_file, "w") as f:
        f.write(result.read())
```

## Database Reference

### Nucleotide

| Database | Description |
|----------|-------------|
| `nt` | All GenBank + EMBL + DDBJ + PDB |
| `refseq_rna` | RefSeq RNA sequences |

### Protein

| Database | Description |
|----------|-------------|
| `nr` | Non-redundant protein |
| `refseq_protein` | RefSeq protein |
| `swissprot` | Curated Swiss-Prot |
| `pdb` | Protein Data Bank |
