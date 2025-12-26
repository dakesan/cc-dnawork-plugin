# COSMIC Data Directory

This directory contains COSMIC database files required for the cosmic-toolkit skill.

## Required Files

### cancer_gene_census.csv

**Description**: Expert-curated list of ~700+ cancer genes with substantial evidence of cancer involvement.

**Download Instructions**:

1. **Register for COSMIC account** (free for academic use):
   - Visit: https://cancer.sanger.ac.uk/cosmic/register
   - Fill in your details with institutional email
   - Verify your email address

2. **Login to COSMIC**:
   - Visit: https://cancer.sanger.ac.uk/cosmic
   - Login with your credentials

3. **Download Cancer Gene Census**:
   - Navigate to: https://cancer.sanger.ac.uk/cosmic/download
   - Or direct link: https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/latest/cancer_gene_census.csv
   - Select "Cancer Gene Census" (GRCh38, latest version)
   - Download the CSV file

4. **Place the file in this directory**:
   ```bash
   # Move downloaded file to this directory
   mv ~/Downloads/cancer_gene_census.csv /path/to/cosmic-toolkit/data/
   ```

### File Structure

After setup, this directory should contain:

```
data/
├── .gitkeep                    # Git tracking (do not delete)
├── README.md                   # This file
└── cancer_gene_census.csv      # COSMIC Cancer Gene Census (user downloads)
```

## Usage

Once `cancer_gene_census.csv` is in place, you can use the cosmic-toolkit scripts:

```bash
# Query a single gene
python scripts/query_cosmic_genes.py --gene TP53

# Query multiple genes
python scripts/query_cosmic_genes.py --genes TP53 BRCA1 EGFR

# Query from file
python scripts/query_cosmic_genes.py --gene-list my_genes.txt
```

## File Format

The `cancer_gene_census.csv` file contains columns such as:
- Gene Symbol
- Name
- Entrez GeneId
- Genome Location
- Tier (1 or 2)
- Hallmark (Yes/No)
- Role in Cancer (TSG, oncogene, fusion)
- Tumour Types (Somatic)
- Tumour Types (Germline)
- Cancer Syndrome
- And more...

**Note**: The script automatically reads all columns and outputs them in JSON format, so it adapts to COSMIC format updates.

## Data Updates

COSMIC is updated quarterly (approximately every 3 months). To get the latest data:

1. Download the newest version of `cancer_gene_census.csv`
2. Replace the file in this directory
3. Check the COSMIC release notes: https://cancer.sanger.ac.uk/cosmic/release_notes

## License

**Academic Use**: Free with registration
**Commercial Use**: License required through QIAGEN

Contact: cosmic-translation@sanger.ac.uk

## Citation

When using COSMIC data, please cite:

Tate JG, Bamford S, Jubb HC, et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Research. 2019;47(D1):D941-D947.

## Troubleshooting

### File Not Found Error

If you see:
```
Cancer Gene Census file not found at: data/cancer_gene_census.csv
```

Follow the download instructions above to obtain the file.

### Authentication Error

COSMIC requires authentication. Make sure you:
- Have registered an account
- Use institutional email for academic access
- Contact QIAGEN for commercial licensing if needed

### File Format Error

If the CSV format has changed:
- The script should still work (reads all columns dynamically)
- If issues persist, check COSMIC documentation for format changes
- Report issues to the skill maintainer
