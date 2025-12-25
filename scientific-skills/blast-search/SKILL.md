---
name: blast-search
description: "NCBI BLAST sequence similarity search via the Common URL API. Use when a user wants to run BLAST programmatically with blastn/blastp and retrieve results in JSON/XML/Text."
---

# BLAST Search

NCBI BLAST (Basic Local Alignment Search Tool) を Common URL API で実行するスキルです。

## Quick Start

### Install

```bash
uv pip install requests typer
```

### Run with FASTA

```bash
python scripts/run_blast_url_api.py run --program blastn --database core_nt --fasta path/to/query.fasta
```

### Run with raw sequence

```bash
python scripts/run_blast_url_api.py run --program blastp --database swissprot --sequence MTEYKLVVVG...
```

### Save output

```bash
python scripts/run_blast_url_api.py run --program blastn --database core_nt --sequence ACTG... --format-type JSON2 --output blast.json
```
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
