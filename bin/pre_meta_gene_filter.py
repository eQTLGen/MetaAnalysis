"""
Created:      03/05/2024
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2024 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import os
import sys
import argparse

import pandas as pd
import numpy as np

# Metadata
__program__ = "CNV-caller"
__author__ = "C.A. (Robert) Warmerdam"
__email__ = "c.a.warmerdam@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


# Constants

# Classes

# Functions

# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # Process input
    parser = argparse.ArgumentParser()
    parser.add_argument('--expression', nargs='+')
    parser.add_argument('--gene-inclusion', nargs='+')
    parser.add_argument('--output')
    args = parser.parse_args(argv)
    # Perform method

    gene_inclusion_list = list()
    total_sample_size = 0
    n_cohorts = len(args.expression)

    for i, expression_folder in enumerate(args.expression):
        expression_path = os.path.join(expression_folder, os.listdir(expression_folder)[0])
        gene_inclusion = pd.read_csv(args.gene_inclusion[i], sep="\t", header=0, names=['ID'])

        sample_size = len(pd.read_csv(expression_path, sep="\t", header=0).index)
        total_sample_size += sample_size

        gene_inclusion['sample_size'] = sample_size
        gene_inclusion['cohort_number'] = 1
        print(gene_inclusion.head())

        gene_inclusion_list.append(gene_inclusion)

    gene_inclusion_total = pd.concat(gene_inclusion_list, axis=0).groupby('ID').agg({'sample_size': 'sum', 'cohort_number': 'sum'}).reset_index()

    genes_to_analyse = (gene_inclusion_total.loc[
        np.logical_and(gene_inclusion_total['sample_size'] >= total_sample_size * 0.5,
        gene_inclusion_total['cohort_number'] >= n_cohorts * 0.5), 'ID'])

    genes_to_analyse.to_csv(args.output, sep="\t", index=False, header=True)
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
