# BLAT Installation Guide

Installation guide for BLAT tool on different platforms.

## Overview

BLAT (BLAST-Like Alignment Tool) has two usage modes:
- **UCSC REST API**: No installation needed (rate-limited)
- **Local command-line**: Requires BLAT binary installation

This guide covers local installation.

## macOS Installation

### Using Homebrew (Recommended)

```bash
# Install BLAT
brew install blat

# Verify installation
blat
```

### Manual Installation

```bash
# Download BLAT for macOS
cd ~/Downloads
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat

# Make executable
chmod +x blat

# Move to PATH
sudo mv blat /usr/local/bin/

# Verify
blat
```

## Linux Installation

### Ubuntu/Debian

```bash
# Using package manager
sudo apt-get update
sudo apt-get install blat-suite

# Verify
blat
```

### RedHat/CentOS/Fedora

```bash
# Install from source or download binary
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat

# Make executable and install
chmod +x blat
sudo mv blat /usr/local/bin/

# Verify
blat
```

### Manual Installation (All Linux)

```bash
# Create installation directory
mkdir -p ~/tools/blat
cd ~/tools/blat

# Download BLAT binary for your architecture
# For 64-bit Linux x86:
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat

# For other architectures, see:
# http://hgdownload.soe.ucsc.edu/admin/exe/

# Make executable
chmod +x blat

# Add to PATH
echo 'export PATH=$PATH:~/tools/blat' >> ~/.bashrc
source ~/.bashrc

# Verify
blat
```

## Windows Installation

Windows requires Windows Subsystem for Linux (WSL2) or manual setup.

### Using Windows Subsystem for Linux (WSL2) - Recommended

```bash
# In WSL2 terminal
sudo apt-get update
sudo apt-get install blat-suite

# Verify
blat
```

### Manual Installation in WSL2

```bash
# Download Linux version
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat

# Make executable
chmod +x blat

# Install
sudo mv blat /usr/local/bin/

# Verify
blat
```

## Download Reference Genome Files (2bit format)

For local BLAT searches, you need reference genome files in 2bit format.

### List of Available Genomes

```bash
# Browse available genomes:
# http://hgdownload.soe.ucsc.edu/goldenPath/
```

### Download Common Genomes

```bash
# Create directory for genomes
mkdir -p ~/genomes
cd ~/genomes

# Human (hg38)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

# Human (hg19 - older assembly)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

# Mouse (mm39)
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit

# Mouse (mm10 - older assembly)
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit

# Rat (rn7)
wget http://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.2bit

# Yeast (sacCer3)
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit

# Zebrafish (danRer11)
wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.2bit

# Arabidopsis (araPort11)
wget http://hgdownload.soe.ucsc.edu/goldenPath/araPort11/bigZips/araPort11.2bit
```

### File Size Reference

```
hg38.2bit     ~720 MB
hg19.2bit     ~720 MB
mm39.2bit     ~700 MB
mm10.2bit     ~700 MB
rn7.2bit      ~650 MB
sacCer3.2bit  ~12 MB
danRer11.2bit ~370 MB
```

## Verify Installation

### Check BLAT Version

```bash
blat
# Output:
# BLAT - Blast-like alignment tool by Jim Kent
# ...
# blat - BLAT v. 36x2
```

### Test BLAT Search

```bash
# Create test sequence file
cat > test.fa << 'EOF'
>test_seq
ATGCGTACGATCGATCGATCGATCGATCG
EOF

# Test with human genome (must have hg38.2bit first)
blat ~/genomes/hg38.2bit test.fa output.psl

# Check output
cat output.psl | head -10
```

### Test Python Integration

```python
from blat_integration import BlatRunner

runner = BlatRunner()

# Check installation
print(runner.check_installation())
# Output should show: blat_installed: True

# Test API (no local installation required)
print(runner.test_api_connection())
# Output should show: api_available: True
```

## Troubleshooting

### Problem: "blat: command not found"

**Solution 1**: Install BLAT
```bash
# macOS
brew install blat

# Linux
sudo apt-get install blat-suite

# Or download manually (see above)
```

**Solution 2**: Add to PATH
```bash
# If BLAT is installed but not in PATH
export PATH=$PATH:~/tools/blat
echo 'export PATH=$PATH:~/tools/blat' >> ~/.bashrc
```

### Problem: "cannot find database file"

**Solution**: Download 2bit file
```bash
mkdir -p ~/genomes
cd ~/genomes
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
```

### Problem: "Segmentation fault"

**Solution**: Check file paths and permissions
```bash
# Verify 2bit file exists
ls -lh ~/genomes/hg38.2bit

# Check BLAT binary
which blat
ls -l /usr/local/bin/blat

# Try small test
echo ">test\nATGC" | blat ~/genomes/hg38.2bit /dev/stdin output.psl
```

### Problem: Slow BLAT searches

**Solution 1**: Check hardware
- BLAT loads entire genome into memory
- Requires ~1.5-3 GB RAM for large genomes
- Searches should be fast once loaded

**Solution 2**: Check file system
- Use local SSD for best performance
- Network drives are much slower
- Copy 2bit file to local disk if needed

### Problem: API rate limit exceeded

**Solution**: Use local BLAT
```python
runner = BlatRunner(mode="local")
# Requires local BLAT installation (see above)
```

## Platform-Specific Tips

### macOS with Apple Silicon (M1/M2)

BLAT may not have native Apple Silicon binaries. Use Homebrew:
```bash
brew install blat
# Homebrew handles architecture automatically
```

### Docker Installation

If you prefer containerization:
```dockerfile
FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y blat-suite python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Download genomes
RUN mkdir -p /data/genomes && \
    cd /data/genomes && \
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

WORKDIR /work
```

Run with:
```bash
docker build -t blat-integration .
docker run -v $(pwd):/work blat-integration python3 your_script.py
```

## Environment Configuration

### Set Default Genome Location

```python
# In your script
import os

BLAT_DB_PATH = os.environ.get('BLAT_DB', os.path.expanduser('~/genomes'))

runner = BlatRunner()
results = runner.search(
    sequence="ATGC...",
    database=f"{BLAT_DB_PATH}/hg38.2bit",
    mode="local"
)
```

### Or set environment variable

```bash
export BLAT_DB=$HOME/genomes
# Now in Python:
# database = f"{os.environ['BLAT_DB']}/hg38.2bit"
```

## Performance Optimization

### Batch Searches

```bash
# For large FASTA files, BLAT is optimized
# Process all sequences in one run:
blat hg38.2bit large_file.fasta output.psl

# Not in a loop:
for seq in *.fa; do
    blat hg38.2bit $seq ${seq}.psl  # Slow!
done
```

### Memory Optimization

BLAT loads genome into RAM. Monitor usage:
```bash
# macOS
top -l 1 | grep blat

# Linux
ps aux | grep blat
top -p $(pgrep blat)
```

### Parallel Searches

```python
from multiprocessing import Pool
import subprocess

def blat_search(args):
    fasta, db, output = args
    subprocess.run(['blat', db, fasta, output])

# Search multiple files in parallel
queries = [('query1.fa', 'hg38.2bit', 'out1.psl'),
           ('query2.fa', 'hg38.2bit', 'out2.psl')]

with Pool(4) as p:
    p.map(blat_search, queries)
```

## See Also

- [UCSC BLAT Home](https://genome.ucsc.edu/cgi-bin/hgBlat)
- [BLAT FAQ](https://genome.ucsc.edu/FAQ/FAQblat.html)
- [Download BLAT Binaries](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [Available Genomes](http://hgdownload.soe.ucsc.edu/goldenPath/)
