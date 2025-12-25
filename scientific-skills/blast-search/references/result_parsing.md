# BLAST Result Parsing

## XML Parsing with NCBIXML

### Single Result

```python
from Bio.Blast import NCBIXML

with open("blast_results.xml") as f:
    blast_record = NCBIXML.read(f)
```

### Multiple Results

```python
from Bio.Blast import NCBIXML

with open("batch_results.xml") as f:
    for blast_record in NCBIXML.parse(f):
        print(f"Query: {blast_record.query}")
        print(f"Hits: {len(blast_record.alignments)}")
```

## Blast Record Structure

```python
# Query information
blast_record.query           # Query sequence ID
blast_record.query_length    # Query length
blast_record.database        # Database searched

# Alignments (hits)
for alignment in blast_record.alignments:
    alignment.title          # Full hit description
    alignment.accession      # Accession number
    alignment.length         # Subject sequence length

    # HSPs (High-Scoring Pairs)
    for hsp in alignment.hsps:
        hsp.expect           # E-value
        hsp.score            # Raw score
        hsp.bits             # Bit score
        hsp.identities       # Number of identical matches
        hsp.align_length     # Alignment length
        hsp.gaps             # Number of gaps
        hsp.query            # Aligned query sequence
        hsp.match            # Match line (| for match)
        hsp.sbjct            # Aligned subject sequence
        hsp.query_start      # Query start position
        hsp.query_end        # Query end position
        hsp.sbjct_start      # Subject start position
        hsp.sbjct_end        # Subject end position
```

## Common Patterns

### Display Top Hits

```python
for alignment in blast_record.alignments[:5]:
    hsp = alignment.hsps[0]
    print(f"{alignment.accession}: E={hsp.expect:.2e}")
```

### Filter by E-value

```python
E_VALUE_THRESH = 1e-10

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print(f"{alignment.accession}: E={hsp.expect}")
```

### Calculate Percent Identity

```python
def percent_identity(hsp):
    return (hsp.identities / hsp.align_length) * 100

for alignment in blast_record.alignments[:5]:
    hsp = alignment.hsps[0]
    pct = percent_identity(hsp)
    print(f"{alignment.accession}: {pct:.1f}% identity")
```

### Extract Best Hits as List

```python
def get_best_hits(blast_record, max_hits=10, e_thresh=0.001):
    """Extract best hits as list of dicts."""
    hits = []
    for alignment in blast_record.alignments[:max_hits]:
        hsp = alignment.hsps[0]
        if hsp.expect < e_thresh:
            hits.append({
                'accession': alignment.accession,
                'title': alignment.title,
                'e_value': hsp.expect,
                'bit_score': hsp.bits,
                'identity': hsp.identities / hsp.align_length * 100,
                'align_length': hsp.align_length,
                'query_start': hsp.query_start,
                'query_end': hsp.query_end,
            })
    return hits
```

### Convert to DataFrame

```python
import pandas as pd

def blast_to_dataframe(blast_record, e_thresh=0.001):
    """Convert BLAST results to pandas DataFrame."""
    rows = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_thresh:
                rows.append({
                    'accession': alignment.accession,
                    'description': alignment.title,
                    'e_value': hsp.expect,
                    'bit_score': hsp.bits,
                    'identities': hsp.identities,
                    'align_length': hsp.align_length,
                    'percent_identity': hsp.identities / hsp.align_length * 100,
                    'gaps': hsp.gaps,
                })
    return pd.DataFrame(rows)
```

## Tabular Output (outfmt 6)

For local BLAST with tabular output:

```python
# Default columns for outfmt 6:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

with open("results.txt") as f:
    for line in f:
        fields = line.strip().split('\t')
        query_id = fields[0]
        subject_id = fields[1]
        percent_identity = float(fields[2])
        align_length = int(fields[3])
        e_value = float(fields[10])
        bit_score = float(fields[11])
```

### Parse to DataFrame

```python
import pandas as pd

columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
           'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

df = pd.read_csv("results.txt", sep='\t', names=columns)
```

## Handling No Hits

```python
if not blast_record.alignments:
    print("No hits found")
else:
    print(f"Found {len(blast_record.alignments)} hits")
```

## Fetching Hit Sequences

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

def fetch_hit_sequence(accession):
    """Fetch sequence for a BLAST hit."""
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    )
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record

# Fetch top 5 hits
for alignment in blast_record.alignments[:5]:
    seq = fetch_hit_sequence(alignment.accession)
    print(f">{seq.id}\n{seq.seq[:50]}...")
```
