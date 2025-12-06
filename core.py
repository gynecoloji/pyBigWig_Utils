"""
Core functions for BigWig file manipulation and signal extraction.
"""

import pyBigWig
import numpy as np
import pandas as pd
import math
from typing import List, Tuple, Optional, Dict


def extract_signal_regions(
    regions: List[Tuple[str, int, int]], 
    bw_file: str, 
    bin_size: int = 25
) -> np.ndarray:
    """
    Extract binned signal from specified genomic regions.
    
    Parameters:
    -----------
    regions : list of tuples
        Each tuple is (chromosome, start, end)
    bw_file : str
        Path to BigWig file
    bin_size : int
        Size of bins in base pairs (default: 25)
    
    Returns:
    --------
    np.ndarray
        Matrix of shape (n_regions, n_bins) containing signal values
    """
    with pyBigWig.open(bw_file) as bw:
        matrix = []
        for chrom, start, end in regions:
            n_bins = (end - start) // bin_size
            binned = bw.stats(chrom, start, end, nBins=n_bins)
            binned = [v if v is not None else 0 for v in binned]
            matrix.append(binned)
    
    return np.array(matrix)


def extract_signals_multi_marks(
    regions: List[Tuple[str, int, int]], 
    bw_files: List[str], 
    bin_size: int = 25,
    verbose: bool = True
) -> np.ndarray:
    """
    Extract signals from multiple histone marks/samples.
    
    Parameters:
    -----------
    regions : list of tuples
        Genomic regions to analyze
    bw_files : list of str
        Paths to multiple BigWig files
    bin_size : int
        Bin size for signal extraction
    verbose : bool
        Print progress messages
    
    Returns:
    --------
    np.ndarray
        3D array of shape (n_regions, n_marks, n_bins)
    """
    all_marks = []
    for bw_path in bw_files:
        if verbose:
            print(f"Processing {bw_path}...")
        matrix = extract_signal_regions(regions, bw_path, bin_size=bin_size)
        all_marks.append(matrix)
    
    stacked = np.stack(all_marks, axis=0)  # (n_marks, n_regions, n_bins)
    return np.transpose(stacked, (1, 0, 2))  # (n_regions, n_marks, n_bins)


def compare_treatment_control(
    treat_file: str,
    ctrl_file: str,
    regions: List[Tuple[str, int, int]],
    bin_size: int = 25
) -> List[Dict]:
    """
    Compare treatment vs control signals across regions.
    
    Parameters:
    -----------
    treat_file : str
        Path to treatment BigWig file
    ctrl_file : str
        Path to control BigWig file
    regions : list of tuples
        Genomic regions to compare
    bin_size : int
        Bin size for analysis
    
    Returns:
    --------
    list of dict
        Each dict contains region info, mean values, and fold-change
    """
    with pyBigWig.open(treat_file) as bw_treat, \
         pyBigWig.open(ctrl_file) as bw_ctrl:
        
        results = []
        for chrom, start, end in regions:
            n_bins = (end - start) // bin_size
            treat_mean = bw_treat.stats(chrom, start, end, type="mean", nBins=n_bins)[0]
            ctrl_mean = bw_ctrl.stats(chrom, start, end, type="mean", nBins=n_bins)[0]
            
            treat_mean = treat_mean if treat_mean else 0
            ctrl_mean = ctrl_mean if ctrl_mean else 0
            
            diff = treat_mean - ctrl_mean
            fold_change = math.log2((treat_mean + 1) / (ctrl_mean + 1))
            
            results.append({
                'region': f"{chrom}:{start}-{end}",
                'chrom': chrom,
                'start': start,
                'end': end,
                'treatment': treat_mean,
                'control': ctrl_mean,
                'difference': diff,
                'log2_FC': fold_change
            })
    
    return results


def batch_extract_to_dataframe(
    bw_files: List[str],
    regions: List[Tuple[str, int, int]],
    sample_names: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Extract data from multiple files into a tidy DataFrame.
    
    Parameters:
    -----------
    bw_files : list of str
        Paths to BigWig files
    regions : list of tuples
        Genomic regions
    sample_names : list of str, optional
        Names for samples (auto-generated if None)
    
    Returns:
    --------
    pd.DataFrame
        Tidy dataframe with columns: sample, chrom, start, end, mean, max
    """
    if sample_names is None:
        sample_names = [f"sample_{i}" for i in range(len(bw_files))]
    
    all_data = []
    
    for bw_file, sample_name in zip(bw_files, sample_names):
        with pyBigWig.open(bw_file) as bw:
            for chrom, start, end in regions:
                mean_val = bw.stats(chrom, start, end, type="mean")[0]
                max_val = bw.stats(chrom, start, end, type="max")[0]
                
                all_data.append({
                    'sample': sample_name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'mean': mean_val if mean_val else 0,
                    'max': max_val if max_val else 0
                })
    
    return pd.DataFrame(all_data)


def normalize_to_control(
    treat_file: str,
    control_file: str,
    regions: List[Tuple[str, int, int]],
    pseudocount: float = 0.1
) -> List[Dict]:
    """
    Normalize treatment signal to control (log2 fold-change).
    
    Parameters:
    -----------
    treat_file : str
        Treatment BigWig file
    control_file : str
        Control BigWig file
    regions : list of tuples
        Genomic regions
    pseudocount : float
        Value to add to avoid log(0)
    
    Returns:
    --------
    list of dict
        Normalized values for each region
    """
    with pyBigWig.open(treat_file) as treat, \
         pyBigWig.open(control_file) as control:
        
        normalized = []
        
        for chrom, start, end in regions:
            treat_val = treat.stats(chrom, start, end, type="mean")[0]
            ctrl_val = control.stats(chrom, start, end, type="mean")[0]
            
            treat_val = treat_val if treat_val else pseudocount
            ctrl_val = ctrl_val if ctrl_val else pseudocount
            
            log2fc = np.log2(treat_val / ctrl_val)
            
            normalized.append({
                'region': f"{chrom}:{start}-{end}",
                'chrom': chrom,
                'start': start,
                'end': end,
                'log2_fold_change': log2fc,
                'treatment': treat_val,
                'control': ctrl_val
            })
        
        return normalized
