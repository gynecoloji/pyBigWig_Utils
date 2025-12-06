# BigWig_Utils

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive Python toolkit for BigWig file analysis in computational biology and epigenomics research.

## 🔬 Overview

pyBigWig_Utils provides a robust, production-ready suite of tools for analyzing ChIP-seq, ATAC-seq, Cut&Run-seq and other genomic signal data stored in BigWig format. Designed with computational biologists in mind, it emphasizes clean code, comprehensive documentation, and biologically relevant analyses.

### Key Features

- **Signal Extraction**: Extract and quantify signals across genomic regions with flexible binning
- **Multi-Sample Analysis**: Correlation analysis across multiple histone marks or experimental conditions
- **Differential Analysis**: Identify regions with significant signal changes between conditions
- **Visualization**: Production-ready plots including heatmaps, metaprofiles, volcano plots, and genomic tracks

## 📦 Installation

### From Source

```bash
git clone https://github.com/gynecoloji/pyBigWig_Utils.git
cd pyBigWig_Utils
```

### Dependencies

```
pyBigWig>=0.3.18
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
```

## 🚀 Quick Start

### Basic Signal Extraction

```python
# Define genomic regions
regions = [
    ("chr1", 10000, 20000),
    ("chr1", 50000, 60000),
    ("chr2", 30000, 40000)
]

# Extract signal from a BigWig file
signal_matrix = extract_signal_regions(
    regions=regions,
    bw_file="path/to/H3K4me3.bigWig",
    bin_size=100
)

print(f"Signal matrix shape: {signal_matrix.shape}")
# Output: (3, 100) - 3 regions × 100 bins
```

### Multi-Mark Analysis

```python
# Analyze multiple histone marks simultaneously
bw_files = [
    "H3K4me1.bigWig",
    "H3K4me3.bigWig",
    "H3K27ac.bigWig",
    "H3K27me3.bigWig"
]

multi_signal = extract_signals_multi_marks(
    regions=regions,
    bw_files=bw_files,
    bin_size=100
)

print(f"Multi-mark shape: {multi_signal.shape}")
# Output: (3, 4, 100) - 3 regions × 4 marks × 100 bins
```

### Correlation Analysis

```python
# Calculate correlation between samples
corr_matrix = calculate_correlation_matrix(
    bw_files=["sample1.bigWig", "sample2.bigWig", "sample3.bigWig"],
    regions=regions,
    sample_names=["Control", "Treatment_4h", "Treatment_24h"],
    method='pearson'
)

# Visualize
fig = plot_correlation_heatmap(corr_matrix, title="Sample Correlation")
fig.savefig("correlation_heatmap.png", dpi=300)
```

### Differential Analysis

```python
# Compare treatment vs control
diff_results = find_differential_regions(
    treat_file="treatment.bigWig",
    ctrl_file="control.bigWig",
    regions=regions,
    log2fc_threshold=1.0,
    min_signal=2.0
)

# Save results
save_results_to_csv(diff_results, "differential_regions.csv")

# Visualize
fig = plot_volcano(diff_results, title="Treatment vs Control")
fig.savefig("volcano_plot.png", dpi=300)
```

## 📖 Documentation

### Core Functions

- `extract_signal_regions()`: Extract binned signals from genomic regions
- `extract_signals_multi_marks()`: Extract signals from multiple BigWig files
- `compare_treatment_control()`: Compare two conditions across regions
- `batch_extract_to_dataframe()`: Create tidy dataframes from multiple files

### Analysis Functions

- `calculate_correlation_matrix()`: Compute sample-sample correlations
- `pairwise_mark_correlation()`: Correlate different histone marks
- `calculate_enrichment_score()`: Compute enrichment over background
- `find_differential_regions()`: Identify differentially enriched regions
- `compute_frip_score()`: Calculate Fraction of Reads in Peaks

### Visualization Functions

- `plot_signal_heatmap()`: Heatmap of signals across regions
- `plot_correlation_heatmap()`: Correlation matrix visualization
- `plot_metaprofile()`: Average signal profiles
- `plot_volcano()`: Volcano plot for differential analysis
- `plot_genomic_track()`: Browser-style genomic tracks


## 👤 Author

**gynecoloji**
- GitHub: [@gynecoloji](https://github.com/gynecoloji)
- Email: your.email@example.com

## 🙏 Acknowledgments

- Built with [pyBigWig](https://github.com/deeptools/pyBigWig)
- Inspired by common workflows in epigenomics research

## 📚 Citation

If you use BigWig_Utils in your research, please cite:

```
[Your Name]. (2024). BigWig_Utils: A comprehensive toolkit for BigWig file analysis.
GitHub: https://github.com/gynecoloji/BigWig_Utils
```

## 🔗 Related Projects

- [deepTools](https://deeptools.readthedocs.io/) - Tools for exploring deep sequencing data
- [MACS](https://github.com/macs3-project/MACS) - Peak calling for ChIP-seq data
- [pyBigWig](https://github.com/deeptools/pyBigWig) - Python interface to BigWig files
