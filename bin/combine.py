"""
Created:      29/11/2022
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2022 C.A. Warmerdam

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
import glob

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


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
    parser.add_argument('--path')
    parser.add_argument('--phenotypes')

    args = parser.parse_args(argv)
    # Perform method

    phenotypes = list()
    per_cohort = False

    with open(args.phenotypes) as opened:
        for line in opened:
            phenotypes.append(line.strip())

    for phenotype in phenotypes:

        print("Reading phenotype:", phenotype)
        parquet_dataset = pq.ParquetDataset(args.path, filters=[("phenotype", "=", phenotype)]))

        print("Writing dataset")
        pq.write_to_dataset(
            table=parquet_dataset,
            root_path=args.out,
            partition_cols=["phenotype", "cohort"] if per_cohort else ["phenotype"],
            max_rows_per_group=524288)

    return 0


if __name__ == "__main__":
    sys.exit(main())
