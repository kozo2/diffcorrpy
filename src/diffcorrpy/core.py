# -*- coding: utf-8 -*-

"""Core functions for Differential Correlation Analysis."""

import numpy as np
from scipy import stats

def cor2_test(n, r, method="pearson"):
    """
    Calculates the p-value for a correlation coefficient.

    Python translation of the R function `cor2.test`.

    Args:
        n (int): Number of samples.
        r (float or np.ndarray): Correlation coefficient(s).
        method (str): Method used ("pearson", "spearman", or "kendall").
                      Currently, only "pearson" and "spearman" are implemented
                      similarly to the R version's primary logic.

    Returns:
        float or np.ndarray: The calculated p-value(s).

    References:
        http://aoki2.si.gunma-u.ac.jp/R/cor2.html
    """
    r = np.asarray(r)
    if method in ["pearson", "spearman"]:
        # Ensure r is within valid range to avoid NaNs in sqrt
        r_clipped = np.clip(r, -0.999999, 0.999999)
        t_stat = np.abs(r_clipped) * np.sqrt((n - 2) / (1 - r_clipped**2))
        df = n - 2
        # Use sf (survival function, 1 - cdf) for potentially better precision in the upper tail
        p_values = stats.t.sf(t_stat, df) * 2
        # Handle potential NaN results if input r was exactly 1 or -1 before clipping
        p_values[np.abs(r) >= 1.0] = 0.0
        return p_values
    elif method == "kendall":
        # The R code uses a normal approximation for Kendall's tau p-value.
        # Reference: http://aoki2.si.gunma-u.ac.jp/R/cor2.html
        # Variance of Kendall's tau estimate: (4n + 10) / (9n(n-1))
        # z = tau / sqrt(variance)
        variance = (4 * n + 10) / (9 * n * (n - 1))
        z_stat = np.abs(r) / np.sqrt(variance)
        p_values = stats.norm.sf(z_stat) * 2
        return p_values
    else:
        raise ValueError("Method must be 'pearson', 'spearman', or 'kendall'")


def compcorr(n1, r1, n2, r2):
    """
    Compares two correlation coefficients using Fisher's Z-transformation.

    Python translation of the R function `compcorr`.

    Args:
        n1 (int): Sample size under condition 1.
        r1 (float or np.ndarray): Correlation coefficient(s) under condition 1.
        n2 (int): Sample size under condition 2.
        r2 (float or np.ndarray): Correlation coefficient(s) under condition 2.

    Returns:
        tuple: A tuple containing:
            - diff (float or np.ndarray): Z-score difference.
            - pval (float or np.ndarray): p-value for the difference.

    References:
        http://www.fon.hum.uva.nl/Service/Statistics/Two_Correlations.html
    """
    r1 = np.asarray(r1)
    r2 = np.asarray(r2)

    # Fisher's Z-transformation
    # Clip values close to +/- 1 to avoid infinity in atanh
    r1_clipped = np.clip(r1, -0.999999, 0.999999)
    r2_clipped = np.clip(r2, -0.999999, 0.999999)

    z1 = np.arctanh(r1_clipped)
    z2 = np.arctanh(r2_clipped)

    # Difference of Z-scores, standardized
    diff_z = (z1 - z2) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))

    # p-value from the standard normal distribution
    pval = 2 * stats.norm.sf(np.abs(diff_z)) # sf is 1 - cdf

    return diff_z, pval

import pandas as pd
from statsmodels.stats.multitest import multipletests
import warnings

# Note on get_lfdr:
# The R version uses the 'fdrtool' package for the 'local' p.adjust.method.
# This implements a specific empirical Bayes approach (local FDR).
# A direct equivalent isn't readily available in standard Python stats libraries.
# We will use methods from statsmodels.stats.multitest for other adjustments
# and raise an error or warning for 'local'.

def comp_2_cc_fdr(data1, data2, method="pearson",
                  p_adjust_method="fdr_bh", # Changed default from 'local'
                  threshold=0.05,
                  output_file=None):
    """
    Calculates and exports differential correlations between two conditions.

    Python translation of the R function `comp.2.cc.fdr`.

    Args:
        data1 (pd.DataFrame): Data matrix under condition 1 (rows=molecules, cols=samples).
        data2 (pd.DataFrame): Data matrix under condition 2 (rows=molecules, cols=samples).
        method (str): Correlation method ("pearson", "spearman", "kendall").
        p_adjust_method (str): Method for p-value adjustment. Common options include
                               "bonferroni", "holm", "hochberg", "hommel", "fdr_bh" (Benjamini/Hochberg),
                               "fdr_by" (Benjamini/Yekutieli). 'local' from the R version
                               (using fdrtool) is not directly supported.
        threshold (float): Significance threshold for the adjusted p-value of the difference.
        output_file (str, optional): Path to save the results as a TSV file. If None,
                                     results are not saved to disk.

    Returns:
        pd.DataFrame: DataFrame containing significant differential correlations.
                      Columns include molecule pairs, correlations (r1, r2),
                      p-values (p1, p2, p_diff), correlation difference (r1-r2),
                      and adjusted p-values (p1_adj, p2_adj, p_diff_adj).

    Raises:
        ValueError: If an unsupported p_adjust_method is provided.
        NotImplementedError: If p_adjust_method='local' is requested.
    """
    if not isinstance(data1, pd.DataFrame) or not isinstance(data2, pd.DataFrame):
        raise TypeError("data1 and data2 must be pandas DataFrames.")
    if not data1.index.equals(data2.index):
        raise ValueError("Row indices (molecule names) must be identical for data1 and data2.")

    n1 = data1.shape[1]
    n2 = data2.shape[1]
    n_features = data1.shape[0]
    mol_names = data1.index

    if n_features < 2:
        warnings.warn("Need at least 2 features (rows) to calculate correlations.")
        return pd.DataFrame() # Return empty DataFrame

    # Calculate correlation matrices (pandas calculates column-wise, so transpose first)
    cc1 = data1.transpose().corr(method=method)
    cc2 = data2.transpose().corr(method=method)

    # Get lower triangle indices (excluding diagonal)
    lower_tri_indices = np.tril_indices(n_features, k=-1)
    ccc1 = cc1.values[lower_tri_indices]
    ccc2 = cc2.values[lower_tri_indices]

    # Calculate p-values for individual correlations
    p1 = cor2_test(n1, ccc1, method=method)
    p2 = cor2_test(n2, ccc2, method=method)

    # Calculate difference in correlations and its p-value
    _, pdiff = compcorr(n1, ccc1, n2, ccc2)
    diff = ccc1 - ccc2

    # Handle potential NaNs in pdiff (e.g., if n1 or n2 <= 3)
    pdiff = np.nan_to_num(pdiff, nan=1.0)

    # P-value adjustment
    valid_adjust_methods = ["bonferroni", "holm", "hochberg", "hommel", "fdr_bh", "fdr_by", "none"]
    if p_adjust_method == "local":
         # raise NotImplementedError("The 'local' FDR method (fdrtool) is not directly implemented in Python. "
         #                           "Consider using 'fdr_bh' or other methods from statsmodels.")
         warnings.warn("The 'local' FDR method is not directly implemented. Falling back to 'fdr_bh'.")
         p_adjust_method = 'fdr_bh' # Fallback to BH
    elif p_adjust_method not in valid_adjust_methods:
        raise ValueError(f"p_adjust_method '{p_adjust_method}' not recognized. "
                         f"Choose from: {valid_adjust_methods} or 'local' (uses fdr_bh).")

    if p_adjust_method != "none":
        _, p1_adj, _, _ = multipletests(p1, method=p_adjust_method)
        _, p2_adj, _, _ = multipletests(p2, method=p_adjust_method)
        _, pdiff_adj, _, _ = multipletests(pdiff, method=p_adjust_method)
    else:
        p1_adj = p1 # Or np.full(p1.shape, np.nan) or "not adjusted" string? R uses string. Let's use NaN.
        p2_adj = p2
        pdiff_adj = pdiff
        # R version uses string "not adjusted", let's stick to numeric for easier filtering
        # p1_adj = np.full(p1.shape, "not adjusted")
        # p2_adj = np.full(p2.shape, "not adjusted")
        # pdiff_adj = np.full(pdiff.shape, "not adjusted")


    # Generate molecule name pairs for the lower triangle
    mol_names1 = mol_names[lower_tri_indices[0]]
    mol_names2 = mol_names[lower_tri_indices[1]]

    # Filter results based on the adjusted difference p-value threshold
    significant_indices = np.where(pdiff_adj < threshold)[0]

    if significant_indices.size == 0:
        warnings.warn(f"No significant differential correlations found at threshold {threshold} "
                      f"with adjustment method '{p_adjust_method}'.")
        return pd.DataFrame() # Return empty DataFrame

    # Construct the results DataFrame
    results_data = {
        "molecule_X": mol_names1[significant_indices],
        "molecule_Y": mol_names2[significant_indices],
        "r1": ccc1[significant_indices],
        "p1": p1[significant_indices],
        "r2": ccc2[significant_indices],
        "p2": p2[significant_indices],
        "p_diff": pdiff[significant_indices],
        "r_diff": diff[significant_indices], # R version calls this (r1-r2)
        "p1_adj": p1_adj[significant_indices], # R version calls this lfdr (in cond. 1)
        "p2_adj": p2_adj[significant_indices], # R version calls this lfdr (in cond. 2)
        "p_diff_adj": pdiff_adj[significant_indices] # R version calls this lfdr (difference)
    }
    results_df = pd.DataFrame(results_data)

    # Reorder columns to match R output more closely if desired
    column_order = [
        "molecule_X", "molecule_Y", "r1", "p1", "r2", "p2",
        "p_diff", "r_diff", "p1_adj", "p2_adj", "p_diff_adj"
    ]
    results_df = results_df[column_order]

    # Save to file if requested
    if output_file:
        try:
            # R version saves without header and index, using tab separation
            results_df.to_csv(output_file, sep="\t", index=False, header=False, mode='w')
            print(f"Results saved to {output_file}")
        except Exception as e:
            warnings.warn(f"Could not save results to {output_file}: {e}")

    return results_df


def scaling_methods(data, method="auto"):
    """
    Applies various scaling methods to the data.

    Python translation of the R function `scalingMethods`. Assumes input data
    has features (molecules) as rows and samples as columns.

    Args:
        data (pd.DataFrame): Input data matrix (rows=features, cols=samples).
        method (str): Scaling method to use. Options:
                      "auto": Autoscaling (mean-center and divide by std dev).
                      "range": Range scaling (mean-center and divide by range).
                      "pareto": Pareto scaling (mean-center and divide by sqrt of std dev).
                      "vast": VAST scaling (Variance Stabilization Transformation).
                      "level": Level scaling (mean-center and divide by mean).
                      "power": Power scaling (sqrt transform and mean-center).

    Returns:
        pd.DataFrame: Scaled data matrix.

    Raises:
        ValueError: If an invalid method is specified or data format is incorrect.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if data.shape[1] <= 1:
         warnings.warn("Scaling requires more than one sample (column). Returning original data.")
         return data.copy() # Return a copy to avoid modifying original

    # Ensure numeric data, handling potential non-numeric entries gracefully
    # data_numeric = data.apply(pd.to_numeric, errors='coerce') # This might be too slow for large data
    # A faster approach might be needed if performance is critical
    data_numeric = data.copy() # Assume numeric for now, add checks if needed

    # Calculate stats row-wise (axis=1)
    mean_row = data_numeric.mean(axis=1, skipna=True)
    std_row = data_numeric.std(axis=1, skipna=True, ddof=1) # ddof=1 for sample std dev like R's sd()
    range_row = data_numeric.max(axis=1, skipna=True) - data_numeric.min(axis=1, skipna=True)

    # Avoid division by zero or near-zero
    std_row_safe = std_row.replace(0, np.nan) # Replace 0 std dev with NaN to avoid division error
    range_row_safe = range_row.replace(0, np.nan)
    mean_row_safe = mean_row.replace(0, np.nan)

    scaled_data = data_numeric.copy() # Start with a copy

    if method == "auto":
        # (x - mean) / std
        scaled_data = data_numeric.subtract(mean_row, axis=0).divide(std_row_safe, axis=0)
    elif method == "range":
        # (x - mean) / (max - min)
        scaled_data = data_numeric.subtract(mean_row, axis=0).divide(range_row_safe, axis=0)
    elif method == "pareto":
        # (x - mean) / sqrt(std)
        scaled_data = data_numeric.subtract(mean_row, axis=0).divide(np.sqrt(std_row_safe), axis=0)
    elif method == "vast":
        # mean * (x - mean) / (std^2)
        # Need to handle potential division by zero in std^2
        std_sq_safe = std_row_safe**2
        scaled_data = data_numeric.subtract(mean_row, axis=0).multiply(mean_row, axis=0).divide(std_sq_safe, axis=0)
    elif method == "level":
        # (x - mean) / mean
        scaled_data = data_numeric.subtract(mean_row, axis=0).divide(mean_row_safe, axis=0)
    elif method == "power":
        # sqrt(x) - mean(sqrt(x))
        # Apply sqrt element-wise, then calculate row mean of sqrt(x)
        sqrt_data = np.sqrt(data_numeric.clip(lower=0)) # Ensure non-negative before sqrt
        mean_sqrt_row = sqrt_data.mean(axis=1, skipna=True)
        scaled_data = sqrt_data.subtract(mean_sqrt_row, axis=0)
    else:
        raise ValueError(f"Unknown scaling method: {method}. Choose from 'auto', 'range', 'pareto', 'vast', 'level', 'power'.")

    # Handle rows where scaling resulted in NaNs (e.g., due to zero std dev/range/mean)
    # Option 1: Fill NaNs with 0 (or another value)
    # scaled_data = scaled_data.fillna(0)
    # Option 2: Keep NaNs (caller needs to handle them) - Let's keep NaNs for now

    return scaled_data


def cor_dist(data, method="pearson", absolute=False):
    """
    Calculates correlation distance (1-r or 1-|r|).

    Assumes input data has features (molecules) as rows and samples as columns.
    Calculates correlation between rows (features).

    Python translation of the R function `cor.dist`.

    Args:
        data (pd.DataFrame): Input data matrix (rows=features, cols=samples).
        method (str): Correlation method ("pearson", "spearman", "kendall").
        absolute (bool): If True, uses absolute correlation (1 - |r|).
                         If False (default), uses 1 - r.

    Returns:
        pd.DataFrame: A square DataFrame representing the correlation distance matrix
                      between features.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")

    # Calculate correlation between rows (features)
    # pandas .corr() calculates column-wise, so transpose first
    corr_matrix = data.transpose().corr(method=method)

    if absolute:
        dist_matrix = 1 - np.abs(corr_matrix)
    else:
        dist_matrix = 1 - corr_matrix

    # Ensure the diagonal is zero
    np.fill_diagonal(dist_matrix.values, 0)

    return dist_matrix


# Add other translated functions here...
# Note: Functions like uncent.cor2dist, uncent.cordist, cluster.molecule,
# get.eigen.molecule, generate_g, plotClusterMolecules, etc., are not
# translated yet. They might require additional dependencies (e.g., scikit-learn,
# networkx, matplotlib) and more complex implementation.
