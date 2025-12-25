# cc-dnawork-plugin

A curated collection of **68 scientific skills** for DNA work, molecular biology, genomics, and bioinformatics research using Claude Code.

## Overview

This plugin is a specialized subset of [K-Dense-AI/claude-scientific-skills](https://github.com/K-Dense-AI/claude-scientific-skills), specifically curated for DNA-related research and bioinformatics workflows. It includes all essential tools for:

- 🧬 DNA sequence analysis and manipulation
- 🔬 Genomic data processing and analysis
- 🧪 Molecular biology and chemistry
- 📊 Data visualization and analysis
- 💾 Access to major scientific databases
- 🤖 Advanced machine learning and AI tools
- 📚 Scientific communication and literature review
- 🧪 Laboratory automation workflows

## Skills Included

### DNA Sequence Analysis (6 skills)
- **BioPython** - Comprehensive Python toolkit for molecular biology
- **pysam** - SAM/BAM/VCF file processing
- **scikit-bio** - Advanced biological sequence operations
- **bioservices** - Access to biological web services
- **gget** - Fast genome information retrieval
- **gtars** - Genomic tools and resources

### Single-Cell & RNA-seq Analysis (4 skills)
- **Scanpy** - Single-cell RNA-seq analysis
- **CellXGene Census** - Large-scale single-cell data integration
- **PyDESeq2** - Differential expression analysis
- **Arboreto** - Gene regulatory network inference

### Genomic Tools (4 skills)
- **ETE Toolkit** - Phylogenetic analysis and tree manipulation
- **DeepTools** - Genomic signal processing
- **Geniml** - Machine learning for genomics
- **ESM** - Protein language models

### Chemistry & Molecular Design (6 skills)
- **RDKit** - Cheminformatics and molecular manipulation
- **Datamol** - Molecular data processing
- **DeepChem** - Deep learning for chemistry
- **DiffDock** - Molecular docking
- **MedChem** - Drug-likeness assessment
- **Molfeat** - Molecular feature computation

### Genomic Databases (14 skills)
- AlphaFold DB, Ensembl, NCBI Gene, UniProt, PDB
- PubMed, ClinVar, COSMIC, ChEMBL
- PubChem, ZINC, DrugBank
- KEGG, Reactome, STRING

### Visualization & Analysis (4 skills)
- **Matplotlib** - Publication-quality figures
- **Seaborn** - Statistical data visualization
- **Plotly** - Interactive visualizations
- **NetworkX** - Network analysis and visualization

### Scientific Communication (10 skills)
- Literature review, Scientific writing, Citation management
- Research lookup, Hypothesis generation
- Scientific visualization, Clinical reports
- And more...

### Laboratory Integration (7 skills)
- Benchling, DNAnexus, LatchBio, OMERO
- Opentrons, Protocols.io, LabArchives

### Plus 13 additional supporting skills

## Getting Started

### Installation

1. **Install Claude Code** (if not already installed)
   ```bash
   # macOS
   curl -fsSL https://claude.ai/install.sh | bash
   ```

2. **Add the Marketplace**
   ```bash
   /plugin marketplace add dakesan/cc-dnawork-plugin
   ```

3. **Install the Plugin**
   - Open Claude Code
   - Run `/plugin list`
   - Select **dnawork-skills** from the available plugins
   - Click **Install**

### Quick Start

Once installed, you can use Claude to execute complex bioinformatics workflows:

```
Analyze this DNA sequence for CpG islands, predict regulatory elements,
and search PubMed for related studies. Use NCBI Gene, Ensembl, and
UniProt databases to find homologous genes and their functions.
```

```
Process this RNA-seq dataset with Scanpy, identify differentially
expressed genes with PyDESeq2, infer gene regulatory networks with
Arboreto, and map results to KEGG pathways.
```

```
Design CRISPR sgRNAs for these target genes. Analyze off-target
binding with RDKit, check ClinVar for disease associations,
and find relevant clinical trials.
```

## Skill Categories

### By Functionality

**Data Processing**
- BioPython, pysam, scikit-bio, Scanpy, PyDESeq2

**Database Access**
- 14 scientific databases (Ensembl, PubMed, UniProt, etc.)

**Analysis & Visualization**
- Matplotlib, Seaborn, Plotly, NetworkX, ETE Toolkit

**Chemistry & Molecular Design**
- RDKit, DeepChem, DiffDock, Datamol, MedChem, Molfeat

**Machine Learning**
- Deep learning frameworks, protein language models (ESM)

**Scientific Communication**
- Literature review, scientific writing, visualization

**Laboratory Automation**
- Integration with Benchling, Opentrons, LatchBio, DNAnexus

## Architecture

```
cc-dnawork-plugin/
├── .claude-plugin/
│   └── marketplace.json       # Plugin configuration
├── scientific-skills/         # 68 individual skills
│   ├── biopython/
│   ├── pysam/
│   ├── scanpy/
│   ├── ... (65 more skills)
│   └── document-skills/
├── README.md                  # This file
└── LICENSE                    # MIT License
```

Each skill contains:
- `SKILL.md` - Comprehensive documentation
- `references/` - Additional documentation and examples
- `scripts/` - Executable code (if applicable)
- `assets/` - Templates and resources (if applicable)

## Usage Examples

### Example 1: Sequence Analysis Pipeline
```
Use BioPython to parse this FASTA file, calculate sequence statistics,
identify ORFs, and perform sequence alignment against UniProt database.
```

### Example 2: RNA-seq Data Integration
```
Load RNA-seq data with Scanpy, identify cell types, analyze differential
expression with PyDESeq2, visualize with Seaborn, and find pathways
in KEGG database.
```

### Example 3: Variant Annotation
```
Parse this VCF file with pysam, annotate variants with Ensembl,
check pathogenicity in ClinVar, find cancer mutations in COSMIC,
and search PubMed for related studies.
```

### Example 4: Drug Discovery
```
Query ChEMBL for kinase inhibitors, analyze SAR with RDKit,
perform virtual docking with DiffDock, and check drug-likeness
with MedChem.
```

## Requirements

- **Python**: 3.9+ (3.12+ recommended)
- **uv**: Python package manager
  ```bash
  # Install uv
  curl -LsSf https://astral.sh/uv/install.sh | sh
  ```
- **Claude Code**: Latest version
- **System**: macOS, Linux, or Windows with WSL2

## Citation

If you use these skills in your research, please cite:

```bibtex
@software{cc_dnawork_plugin,
  author = {Odake, Hiroyuki},
  title = {cc-dnawork-plugin: Curated DNA Work Skills for Claude Code},
  year = {2025},
  url = {https://github.com/dakesan/cc-dnawork-plugin},
  note = {Based on K-Dense Scientific Skills}
}
```

Also cite the original K-Dense Scientific Skills:
```bibtex
@software{claude_scientific_skills_2025,
  author = {{K-Dense Inc.}},
  title = {Claude Scientific Skills},
  year = {2025},
  url = {https://github.com/K-Dense-AI/claude-scientific-skills}
}
```

## Troubleshooting

### Skills Not Loading
- Ensure Claude Code is up to date
- Try: `/plugin marketplace remove` and add again
- Check: `ls ~/.claude/plugins/` for installation

### Missing Dependencies
- Check individual `SKILL.md` files for requirements
- Install: `uv pip install package-name`
- Verify: `uv pip list`

### API Rate Limits
- Most databases have rate limits
- Review individual skill documentation
- Implement caching or batch requests

### Authentication Issues
- Some databases require API keys
- Check `SKILL.md` for authentication setup
- Verify credentials and permissions

## Contributing

Contributions are welcome! If you:
- Find bugs or have suggestions
- Want to add new DNA-related skills
- Have improvements to existing skills

Please open an issue or submit a pull request.

## Related Resources

- [K-Dense Scientific Skills](https://github.com/K-Dense-AI/claude-scientific-skills) - Parent project with 125+ skills
- [Claude Code Documentation](https://docs.claude.com/)
- [Claude Agent SDK](https://github.com/anthropics/anthropic-sdk-python)

## License

MIT License - See LICENSE file for details

---

**Created for DNA work and bioinformatics research with Claude Code** 🧬
