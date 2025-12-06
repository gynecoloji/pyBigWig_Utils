"""
Visualization functions for BigWig data analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple, Optional, Dict
import pyBigWig


def plot_signal_heatmap(
    signal_matrix: np.ndarray,
    region_labels: Optional[List[str]] = None,
    mark_labels: Optional[List[str]] = None,
    title: str = "Signal Heatmap",
    cmap: str = "YlOrRd",
    figsize: Tuple[int, int] = (10, 8)
):
    """
    Plot heatmap of signal across regions and marks.
    
    Parameters:
    -----------
    signal_matrix : np.ndarray
        2D array of signals (regions x marks or regions x bins)
    region_labels : list of str, optional
        Labels for regions
    mark_labels : list of str, optional
        Labels for marks/bins
    title : str
        Plot title
    cmap : str
        Colormap name
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(
        signal_matrix,
        cmap=cmap,
        xticklabels=mark_labels,
        yticklabels=region_labels,
        cbar_kws={'label': 'Signal Intensity'},
        ax=ax
    )
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('Position/Mark', fontsize=12)
    ax.set_ylabel('Region', fontsize=12)
    
    plt.tight_layout()
    return fig


def plot_correlation_heatmap(
    corr_matrix: pd.DataFrame,
    title: str = "Sample Correlation",
    cmap: str = "coolwarm",
    figsize: Tuple[int, int] = (8, 7)
):
    """
    Plot correlation matrix heatmap.
    
    Parameters:
    -----------
    corr_matrix : pd.DataFrame
        Correlation matrix from calculate_correlation_matrix
    title : str
        Plot title
    cmap : str
        Colormap
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(
        corr_matrix,
        annot=True,
        fmt='.3f',
        cmap=cmap,
        vmin=-1,
        vmax=1,
        center=0,
        square=True,
        linewidths=0.5,
        cbar_kws={'label': 'Correlation'},
        ax=ax
    )
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    return fig


def plot_metaprofile(
    signal_matrix: np.ndarray,
    labels: Optional[List[str]] = None,
    title: str = "Meta-Profile",
    xlabel: str = "Position (bins)",
    ylabel: str = "Average Signal",
    figsize: Tuple[int, int] = (10, 6),
    colors: Optional[List[str]] = None
):
    """
    Plot average signal profile across regions.
    
    Parameters:
    -----------
    signal_matrix : np.ndarray
        Matrix of signals (regions x bins) or (marks x bins)
    labels : list of str, optional
        Labels for each profile
    title : str
        Plot title
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    figsize : tuple
        Figure size
    colors : list of str, optional
        Colors for each profile
    
    Returns:
    --------
    matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if signal_matrix.ndim == 1:
        signal_matrix = signal_matrix.reshape(1, -1)
    
    n_profiles = signal_matrix.shape[0]
    
    if colors is None:
        colors = plt.cm.Set2(np.linspace(0, 1, n_profiles))
    
    for i in range(n_profiles):
        label = labels[i] if labels else f"Profile {i+1}"
        mean_signal = np.mean(signal_matrix[i, :])
        ax.plot(signal_matrix[i, :], label=label, color=colors[i], linewidth=2)
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def plot_volcano(
    diff_df: pd.DataFrame,
    log2fc_col: str = 'log2_fold_change',
    pval_col: Optional[str] = None,
    log2fc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    title: str = "Volcano Plot",
    figsize: Tuple[int, int] = (10, 8)
):
    """
    Create volcano plot for differential analysis.
    
    Parameters:
    -----------
    diff_df : pd.DataFrame
        Differential analysis results
    log2fc_col : str
        Column name for log2 fold-change
    pval_col : str, optional
        Column name for p-values (if None, uses signal intensity)
    log2fc_threshold : float
        Threshold for significance
    pval_threshold : float
        P-value threshold
    title : str
        Plot title
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # If no p-values, use negative treatment signal as proxy
    if pval_col is None:
        y_vals = -np.log10(1 / (diff_df['treatment'] + 1))
        ylabel = '-log10(1/signal)'
    else:
        y_vals = -np.log10(diff_df[pval_col] + 1e-300)
        ylabel = '-log10(p-value)'
    
    x_vals = diff_df[log2fc_col]
    
    # Color by significance
    colors = []
    for i, (x, y) in enumerate(zip(x_vals, y_vals)):
        if abs(x) >= log2fc_threshold:
            if x > 0:
                colors.append('red')
            else:
                colors.append('blue')
        else:
            colors.append('gray')
    
    ax.scatter(x_vals, y_vals, c=colors, alpha=0.6, s=30)
    
    # Add threshold lines
    ax.axvline(-log2fc_threshold, color='black', linestyle='--', alpha=0.5)
    ax.axvline(log2fc_threshold, color='black', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('log2(Fold Change)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def plot_genomic_track(
    bw_files: List[str],
    chrom: str,
    start: int,
    end: int,
    labels: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (14, 8)
):
    """
    Plot genomic tracks for multiple BigWig files.
    
    Parameters:
    -----------
    bw_files : list of str
        Paths to BigWig files
    chrom : str
        Chromosome
    start : int
        Start position
    end : int
        End position
    labels : list of str, optional
        Track labels
    colors : list of str, optional
        Track colors
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib figure
    """
    n_tracks = len(bw_files)
    fig, axes = plt.subplots(n_tracks, 1, figsize=figsize, sharex=True)
    
    if n_tracks == 1:
        axes = [axes]
    
    if labels is None:
        labels = [f"Track {i+1}" for i in range(n_tracks)]
    
    if colors is None:
        colors = plt.cm.Set2(np.linspace(0, 1, n_tracks))
    
    positions = np.linspace(start, end, end - start)
    
    for i, (bw_file, label, color) in enumerate(zip(bw_files, labels, colors)):
        with pyBigWig.open(bw_file) as bw:
            # Get values at single-base resolution
            try:
                values = bw.values(chrom, start, end)
                values = [v if v is not None else 0 for v in values]
            except:
                values = [0] * (end - start)
            
            axes[i].fill_between(range(len(values)), values, color=color, alpha=0.7)
            axes[i].set_ylabel(label, fontsize=10)
            axes[i].grid(True, alpha=0.3)
            axes[i].set_xlim(0, len(values))
    
    axes[-1].set_xlabel(f'{chrom}:{start:,}-{end:,}', fontsize=12)
    axes[0].set_title('Genomic Track View', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    return fig


def plot_enrichment_barplot(
    enrichment_df: pd.DataFrame,
    top_n: int = 20,
    score_col: str = 'enrichment_score',
    figsize: Tuple[int, int] = (10, 8)
):
    """
    Plot top enriched regions as bar plot.
    
    Parameters:
    -----------
    enrichment_df : pd.DataFrame
        Enrichment results
    top_n : int
        Number of top regions to show
    score_col : str
        Column with enrichment scores
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib figure
    """
    df_sorted = enrichment_df.nlargest(top_n, score_col).copy()
    df_sorted['region_label'] = df_sorted.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
    )
    
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = plt.cm.RdYlGn(np.linspace(0.3, 0.9, top_n))
    
    ax.barh(range(top_n), df_sorted[score_col].values, color=colors)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels(df_sorted['region_label'].values, fontsize=8)
    ax.set_xlabel('Enrichment Score', fontsize=12)
    ax.set_title(f'Top {top_n} Enriched Regions', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    return fig
