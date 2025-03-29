# -*- coding: utf-8 -*-

"""
diffcorrpy: A Python package for Differential Correlation Analysis.

This package is a Python translation of the R package DiffCorr.
"""

__author__ = """Your Name"""  # Placeholder
__email__ = 'your.email@example.com'  # Placeholder
__version__ = '0.1.0'

# Import functions/classes here to make them available at the package level
from .core import comp_2_cc_fdr, scaling_methods, cor2_test, compcorr, cor_dist

__all__ = [
    'comp_2_cc_fdr',
    'scaling_methods',
    'cor2_test',
    'compcorr',
    'cor_dist',
]
