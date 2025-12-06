"""
Advanced analysis functions for BigWig data including correlation,
statistical tests, and quality metrics.
"""

import numpy as np
import pandas as pd
import pyBigWig
from typing import List, Tuple, Dict, Optional
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform


def calculate_correlation_matrix(
    bw_files: List[str],
    regions: List[Tuple[str, int, int]],
    sample_names: Optional[List[str]] = None,
    method: str = 'pearson',
    bin_size: int = 100
) -> pd.DataFrame:
    """
    Calculate correlation matrix between multiple samples.
    
    Parameters:
    -----------
    bw_files : list of str
        Paths to BigWig files
    regions : list of tuples
        Genomic regions to use for correlation
    sample_names : list of str, optional
        Sample names
    method : str
        'pearson' or 'spearman'
    bin_size : int
        Bin size for signal extraction
    
    Returns:
    --------
    pd.DataFrame
        Correlation matrix
    """
    from .core import extract_signal_regions
    
    if sample_names is None:
        sample_names = [f"sample_{i}" for i in range(len(bw_files))]
    
    # Extract signals for all samples
    all_signals = []
    for bw_file in bw_files:
        matrix = extract_signal_regions(regions, bw_file, bin_size=bin_size)
        # Flatten to 1D array
        all_signals.append(matrix.flatten())
    
    # Calculate correlation
    n_samples = len(all_signals)
    corr_matrix = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(n_samples):
            if method == 'pearson':
                corr, _ = pearsonr(all_signals[i], all_signals[j])
            elif method == 'spearman':
                corr, _ = spearmanr(all_signals[i], all_signals[j])
            else:
                raise ValueError("Method must be 'pearson' or 'spearman'")
            corr_matrix[i, j] = corr
    
    return pd.DataFrame(corr_matrix, index=sample_names, columns=sample_names)


def pairwise_mark_correlation(
    bw_files: List[str],
    regions: List[Tuple[str, int, int]],
    mark_names: Optional[List[str]] = None,
    bin_size: int = 100
) -> Dict[str, float]:
    """
    Calculate pairwise correlations between histone marks.
    
    Parameters:
    -----------
    bw_files : list of str
        Paths to BigWig files for different marks
    regions : list of tuples
        Regions to analyze
    mark_names : list of str, optional
        Names of histone marks
    bin_size : int
        Bin size
    
    Returns:
    --------
    dict
        Pairwise correlation coefficients
    """
    from .core import extract_signal_regions
    
    if mark_names is None:
        mark_names = [f"mark_{i}" for i in range(len(bw_files))]
    
    # Extract signals
    signals = {}
    for bw_file, mark in zip(bw_files, mark_names):
        matrix = extract_signal_regions(regions, bw_file, bin_size=bin_size)
        signals[mark] = matrix.flatten()
    
    # Calculate all pairwise correlations
    correlations = {}
    marks = list(signals.keys())
    for i in range(len(marks)):
        for j in range(i + 1, len(marks)):
            mark1, mark2 = marks[i], marks[j]
            corr, pval = pearsonr(signals[mark1], signals[mark2])
            correlations[f"{mark1}_vs_{mark2}"] = {
                'correlation': corr,
                'p_value': pval
            }
    
    return correlations


def calculate_enrichment_score(
    treatment_file: str,
    control_file: str,
    regions: List[Tuple[str, int, int]],
    background_regions: Optional[List[Tuple[str, int, int]]] = None
) -> pd.DataFrame:
    """
    Calculate enrichment scores comparing signal in regions vs background.
    
    Parameters:
    -----------
    treatment_file : str
        Treatment BigWig file
    control_file : str
        Control BigWig file
    regions : list of tuples
        Target regions (e.g., peaks)
    background_regions : list of tuples, optional
        Background regions (if None, uses whole genome stats)
    
    Returns:
    --------
    pd.DataFrame
        Enrichment scores per region
    """
    results = []
    
    with pyBigWig.open(treatment_file) as treat, \
         pyBigWig.open(control_file) as ctrl:
        
        # Calculate background if provided
        if background_regions:
            bg_treat = np.mean([
                treat.stats(c, s, e, type="mean")[0] or 0
                for c, s, e in background_regions
            ])
            bg_ctrl = np.mean([
                ctrl.stats(c, s, e, type="mean")[0] or 0
                for c, s, e in background_regions
            ])
        else:
            bg_treat = 1.0
            bg_ctrl = 1.0
        
        for chrom, start, end in regions:
            treat_val = treat.stats(chrom, start, end, type="mean")[0] or 0
            ctrl_val = ctrl.stats(chrom, start, end, type="mean")[0] or 0
            
            # Enrichment = (signal / background)
            enrichment = (treat_val / bg_treat) / max(ctrl_val / bg_ctrl, 0.01)
            
            results.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'treatment_signal': treat_val,
                'control_signal': ctrl_val,
                'enrichment_score': enrichment,
                'log2_enrichment': np.log2(enrichment + 0.01)
            })
    
    return pd.DataFrame(results)


def calculate_signal_to_noise(
    bw_file: str,
    peak_regions: List[Tuple[str, int, int]],
    background_regions: List[Tuple[str, int, int]]
) -> Dict[str, float]:
    """
    Calculate signal-to-noise ratio (SNR) for quality assessment.
    
    Parameters:
    -----------
    bw_file : str
        BigWig file to analyze
    peak_regions : list of tuples
        Regions with expected signal (peaks)
    background_regions : list of tuples
        Background regions
    
    Returns:
    --------
    dict
        SNR metrics
    """
    with pyBigWig.open(bw_file) as bw:
        # Calculate signal in peaks
        peak_signals = [
            bw.stats(c, s, e, type="mean")[0] or 0
            for c, s, e in peak_regions
        ]
        
        # Calculate background
        bg_signals = [
            bw.stats(c, s, e, type="mean")[0] or 0
            for c, s, e in background_regions
        ]
        
        peak_mean = np.mean(peak_signals)
        bg_mean = np.mean(bg_signals)
        bg_std = np.std(bg_signals)
        
        snr = (peak_mean - bg_mean) / (bg_std + 1e-10)
        
        return {
            'signal_mean': peak_mean,
            'background_mean': bg_mean,
            'background_std': bg_std,
            'SNR': snr,
            'fold_enrichment': peak_mean / (bg_mean + 1e-10)
        }


def find_differential_regions(
    treat_file: str,
    ctrl_file: str,
    regions: List[Tuple[str, int, int]],
    log2fc_threshold: float = 1.0,
    min_signal: float = 1.0
) -> pd.DataFrame:
    """
    Identify regions with significant differential signal.
    
    Parameters:
    -----------
    treat_file : str
        Treatment BigWig
    ctrl_file : str
        Control BigWig
    regions : list of tuples
        Candidate regions
    log2fc_threshold : float
        Minimum log2 fold-change
    min_signal : float
        Minimum signal in treatment
    
    Returns:
    --------
    pd.DataFrame
        Differential regions with statistics
    """
    from .core import normalize_to_control
    
    normalized = normalize_to_control(treat_file, ctrl_file, regions)
    df = pd.DataFrame(normalized)
    
    # Filter based on thresholds
    df['is_differential'] = (
        (np.abs(df['log2_fold_change']) >= log2fc_threshold) &
        (df['treatment'] >= min_signal)
    )
    
    df['direction'] = df['log2_fold_change'].apply(
        lambda x: 'up' if x > 0 else 'down'
    )
    
    return df.sort_values('log2_fold_change', ascending=False)
