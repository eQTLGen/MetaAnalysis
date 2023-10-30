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
from __future__ import print_function

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
    parser.add_argument('--path', nargs='+')
    parser.add_argument('--phenotypes')
    parser.add_argument('--out')
    parser.add_argument('--ref')

    args = parser.parse_args(argv)
    # Perform method

    variant_reference = (
        pd.read_csv(args.ref, compression = 'gzip', sep = ' ')
        .drop(["allele1", "allele2", "str_allele1", "str_allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp", "CHR": "chromosome"}, axis=1)
        .set_index("variant"))

    variant_table = pa.Table.from_pandas(variant_reference)

    print(variant_reference.head())
    print(variant_reference.dtypes)

    phenotypes = list()

    with open(args.phenotypes) as opened:
        for line in opened:
            phenotypes.append(line.strip())

    for phenotype in phenotypes:
        print(phenotype)

        results_list_phen = list()
        for chunk_index, chunk_path in enumerate(args.path):
            phenotype_partition_glob = os.path.join(chunk_path, "phenotype_{}".format(phenotype), "*.parquet")
            print("(chunk {}/{})".format(chunk_index, len(args.path)))

            for i, (file_name) in enumerate(glob.glob(phenotype_partition_glob)):
                print("Reading file " + file_name)
                results_list_phen.append(pq.ParquetDataset(parquet_dataset, filters=[("phenotype".decode("utf8"), "=", )]).read())

        print("Concatenating datasets")
        results_concatenated = pa.concat_tables(results_list_phen)

        print("Is data per cohort?")
        per_cohort = "cohort" in results_concatenated.columns

        print("Writing dataset")
        pq.write_to_dataset(
            table=results_with_chromosome,
            root_path=args.out,
            partition_cols=["phenotype", "cohort"] if per_cohort else ["phenotype"],
            max_rows_per_group=524288)

    return 0


if __name__ == "__main__":
    sys.exit(main())
