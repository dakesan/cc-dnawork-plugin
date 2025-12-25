---
name: msa-advanced
description: "Advanced Multiple Sequence Alignment (MSA) with support for ClustalW, Muscle, and MAFFT. Includes quality metrics, consensus extraction, phylogenetic analysis, and publication-quality visualization."
---

# MSA Advanced: Multiple Sequence Alignment & Analysis

## Overview

MSA Advanced provides comprehensive Multiple Sequence Alignment capabilities with:
- **Alignment execution**: ClustalW, Muscle, MAFFT support
- **Quality analysis**: Identity, gap statistics, consensus scoring
- **Phylogenetic trees**: Automatic tree building and visualization
- **Visualization**: Heatmaps, sequence logos, alignment plots
- **Result parsing**: Parse and analyze alignment results

## When to Use This Skill

Use this skill when:

- Aligning DNA, RNA, or protein sequences
- Comparing multiple insertion sequences (Insertions from WGS)
- Analyzing BIND domain sequences
- Building phylogenetic trees
- Computing consensus sequences
- Creating publication-quality alignment visualizations
- Extracting sequence motifs from alignments
- Computing sequence identity matrices

## Core Capabilities

### 1. **MSA Execution**
- ClustalW (traditional, robust)
- Muscle (fast, accurate)
- MAFFT (very fast, good accuracy)
- Automatic tool selection based on alignment size

### 2. **Quality Metrics**
- Per-sequence identity (%)
- Overall alignment identity (%)
- Gap percentage analysis
- Conservation scoring
- Consensus quality metrics

### 3. **Consensus Analysis**
- Consensus sequence extraction
- Weighted consensus (with scoring)
- Ambiguity code generation
- High-confidence region identification

### 4. **Phylogenetic Analysis**
- Phylogenetic tree construction
- Distance matrix generation
- Tree visualization (Newick format)
- Clade identification

### 5. **Visualization**
- Sequence identity heatmaps
- Multiple alignment plots
- Sequence logos
- Conservation histograms
- Tree diagrams

### 6. **Result Management**
- Multi-format output (Fasta, Clustal, Phylip, Stockholm)
- Result statistics and summaries
- Block analysis for conserved regions

## Installation and Setup

### Prerequisites

```bash
uv pip install biopython pandas matplotlib seaborn scipy
```

### Install Alignment Tools

**macOS (Homebrew):**
```bash
brew install clustal-w muscle mafft
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install clustalw muscle mafft
```

**Manual Installation:**
- [ClustalW](http://www.clustal.org/)
- [Muscle](https://drive5.com/muscle/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)

## Usage Examples

### Example 1: Basic MSA

```python
from msa_advanced import MSARunner

# Initialize runner
msa = MSARunner(method="muscle", verbose=True)

# Run alignment
alignment = msa.align("sequences.fasta")

# Analyze
stats = msa.get_statistics(alignment)
print(f"Alignment length: {alignment.get_alignment_length()}")
print(f"Identity: {stats['identity']:.2f}%")

# Save
alignment.save("aligned.fasta", "fasta")
```

### Example 2: WGS Insertion Analysis

```python
from msa_advanced import MSARunner, MSAVisualizer
import pandas as pd

# Align multiple insertion sequences
msa = MSARunner(method="muscle")
insertions = ["insertion_1.fasta", "insertion_2.fasta", "insertion_3.fasta"]

all_results = []
for insertion_file in insertions:
    alignment = msa.align(insertion_file)
    consensus = msa.get_consensus(alignment)
    stats = msa.get_statistics(alignment)

    all_results.append({
        "insertion": insertion_file,
        "alignment_length": alignment.get_alignment_length(),
        "identity": stats["identity"],
        "consensus_seq": consensus
    })

# Save results
df = pd.DataFrame(all_results)
df.to_csv("insertion_analysis.csv", index=False)

# Visualize
visualizer = MSAVisualizer()
visualizer.plot_identity_heatmap(alignment, output="heatmap.png")
```

### Example 3: Phylogenetic Analysis

```python
from msa_advanced import MSARunner, PhyloBuilder

# Create alignment
msa = MSARunner(method="mafft")
alignment = msa.align("sequences.fasta")

# Build tree
phylo = PhyloBuilder()
tree = phylo.build_from_alignment(alignment)

# Visualize
phylo.visualize_tree(tree, output="phylo_tree.png")

# Save tree
with open("sequences.nwk", "w") as f:
    f.write(tree.format("newick"))
```

## Key Parameters

```python
MSARunner(
    method="muscle",           # muscle, clustalw, or mafft
    num_threads=4,             # Number of CPUs to use
    max_iterations=100,        # Max alignment iterations
    gap_open=10.0,            # Gap opening penalty
    gap_extend=0.1,           # Gap extension penalty
    verbose=False
)

# Alignment
alignment = msa.align(
    fasta_file,
    output_file=None,
    format="fasta",
    auto_select_method=True    # Auto-select tool based on size
)
```

## Performance

| Method | Speed | Accuracy | Memory | Best For |
|--------|-------|----------|--------|----------|
| MAFFT | Very Fast | Excellent | Low | Large datasets |
| Muscle | Fast | Good | Medium | Balanced |
| ClustalW | Slow | Good | High | Small-medium |

## Related Skills

- **biopython** - Sequence manipulation
- **etetoolkit** - Advanced phylogenetics
- **matplotlib/seaborn** - Advanced visualization
- **blat-integration** - Find insertion sequences first

## Citation

```bibtex
@article{Edgar2004Muscle,
  author = {Edgar, Robert C.},
  title = {MUSCLE: multiple sequence alignment with high accuracy and high throughput},
  journal = {Nucleic Acids Research},
  year = {2004},
  volume = {32},
  number = {5},
  pages = {1792--1797}
}

@article{Katoh2013MAFFT,
  author = {Katoh, Kazutaka and Standley, Daron M.},
  title = {MAFFT multiple sequence alignment software version 7: improvements in performance and usability},
  journal = {Molecular Biology and Evolution},
  year = {2013},
  volume = {30},
  pages = {772--780}
}
```

See `references/` for detailed documentation.
