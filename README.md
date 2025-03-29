# diffcorrpy

[![PyPI version](https://badge.fury.io/py/diffcorrpy.svg)](https://badge.fury.io/py/diffcorrpy) <!-- Placeholder badge -->

A Python package for Differential Correlation Analysis, translated from the R package [DiffCorr](https://github.com/afukushima/DiffCorr) by Atsushi Fukushima.

## Overview

`diffcorrpy` provides tools to analyze and visualize differential correlations in biological networks, particularly for large-scale "omics" data like transcriptomics and metabolomics. It helps identify changes in correlation patterns between two experimental conditions.

This package ports the core functionalities of the original R `DiffCorr` package to Python, focusing on the calculation of differential correlations using methods like Fisher's Z-test for comparing correlation coefficients.

**Note:** This is an initial translation focusing on the core statistical functions (`comp_2_cc_fdr`, `scaling_methods`, `cor2_test`, `compcorr`, `cor_dist`). Functions related to clustering, PCA, network generation, and plotting from the original R package are not included in this version, as similar functionalities are available in standard Python libraries (`scipy`, `scikit-learn`, `networkx`, `matplotlib`, `seaborn`).

## Installation

You can install `diffcorrpy` using pip (once published):

```bash
pip install diffcorrpy
```

Alternatively, install directly from the source directory after cloning:

```bash
pip install .
```

**Dependencies:**

*   Python >= 3.8
*   numpy
*   scipy
*   pandas
*   statsmodels

## Basic Usage

```python
import pandas as pd
import numpy as np
import diffcorrpy as dcp

# --- 1. Load or create your data ---
# Example: Create dummy data for two conditions
# Rows = features (e.g., genes, metabolites), Columns = samples
n_features = 50
n_samples1 = 15
n_samples2 = 20

# Condition 1 data
data1 = pd.DataFrame(np.random.rand(n_features, n_samples1),
                     index=[f'feature_{i+1}' for i in range(n_features)],
                     columns=[f'cond1_sample_{j+1}' for j in range(n_samples1)])

# Condition 2 data (introduce some correlation changes)
data2_base = np.random.rand(n_features, n_samples2)
# Make feature_1 and feature_2 correlated in cond2
data2_base[1, :] = data2_base[0, :] + np.random.normal(0, 0.1, n_samples2)
data2 = pd.DataFrame(data2_base,
                     index=[f'feature_{i+1}' for i in range(n_features)],
                     columns=[f'cond2_sample_{j+1}' for j in range(n_samples2)])

# --- 2. (Optional) Scale the data ---
# Example using Pareto scaling
scaled_data1 = dcp.scaling_methods(data1, method='pareto')
scaled_data2 = dcp.scaling_methods(data2, method='pareto')

# --- 3. Calculate differential correlations ---
# Use scaled or original data
results_df = dcp.comp_2_cc_fdr(
    data1=scaled_data1,
    data2=scaled_data2,
    method="pearson",          # Correlation method ('pearson', 'spearman', 'kendall')
    p_adjust_method="fdr_bh",  # P-value adjustment method
    threshold=0.05,            # Significance threshold for adjusted p-value of difference
    output_file="diffcorr_results.tsv" # Optional: Save results to a file
)

# --- 4. Inspect results ---
print(f"Found {len(results_df)} significant differential correlations.")
if not results_df.empty:
    print(results_df.head())

# --- Example: Calculate correlation distance ---
# Using original data1
dist_matrix = dcp.cor_dist(data1, method="pearson", absolute=False)
print("\nCorrelation distance matrix (first 5x5):")
print(dist_matrix.iloc[:5, :5])

```

## License

This package is licensed under the GNU General Public License v3 or later (GPLv3+), consistent with the original R package. See the [LICENSE](LICENSE) file for details. (Note: LICENSE file needs to be added).

## Acknowledgements

This package is a Python translation of the original R package [DiffCorr](https://github.com/afukushima/DiffCorr) developed by Atsushi Fukushima. Please cite the original work if you use this package in your research:

*   Fukushima, A. (2013). DiffCorr: An R package to analyze and visualize differential correlations in biological networks. *Gene*, *518*(1), 209-214. [https://doi.org/10.1016/j.gene.2012.12.048](https://doi.org/10.1016/j.gene.2012.12.048)
