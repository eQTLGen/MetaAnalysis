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
import pyarrow.feather as ft


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
PYARROW_SCHEMA_META = pa.schema([
    ("phenotype", pa.string()),
    ("chromosome", pa.int8()),
    ("variant", pa.string()),
    ("beta", pa.float64()),
    ("standard_error", pa.float64()),
    ("i_squared", pa.float64()),
    ("sample_size", pa.float64())])


PYARROW_SCHEMA_COHORT = pa.schema([
    ("phenotype", pa.string()),
    ("chromosome", pa.int8()),
    ("variant", pa.string()),
    ("beta", pa.float64()),
    ("standard_error", pa.float64()),
    ("sample_size", pa.float64()),
    ("cohort", pa.string())])


# Classes
class FeatherDataset:
    def __init__(self, path, filters):
        self.filters = filters
        self.path = path

        partitioning = list()

        for col_name, operator, value in filters:
            partitioning.append("{}={}".format(col_name, value))

        requested_files = os.path.join(self.path, os.path.join(*partitioning), "*.feather")
        self.matched_files = glob.glob(requested_files)

    def read(self):
        return pd.concat([ft.read_feather(f) for f in self.matched_files])


# Functions
def _combine(filters, partition_cols, path, out, schema):
    parquet_dataset = pq.ParquetDataset(path, schema=schema, validate_schema=False, filters=filters)
    if len(parquet_dataset.pieces) == 0:
        print("length of pieces equal to 0!")
    else:
        results_dataset = parquet_dataset.read()
        print("Writing dataset")
        pq.write_to_dataset(
            table=results_dataset,
            root_path=out,
            partition_cols=partition_cols,
            row_group_size=524288)


def _combine_feather(filters, partition_cols, path, out, schema):
    feather_dataset = FeatherDataset(path, filters=filters)
    if len(feather_dataset.matched_files) == 0:
        print("length of pieces equal to 0!")
    else:
        results_dataset = pa.Table.from_pandas(feather_dataset.read(), schema)
        print("Writing dataset")
        pq.write_to_dataset(
            table=results_dataset,
            root_path=out,
            partition_cols=partition_cols,
            row_group_size=524288)


def combine_per_cohort(path, out, phenotypes, cohorts):
    partition_cols = ["phenotype".decode("utf8"), "cohort".decode("utf8")]

    schema = PYARROW_SCHEMA_COHORT

    for cohort in cohorts:
        print(cohort)
        for phenotype in phenotypes:
            print(phenotype)
            filters = [("cohort".decode("utf8"), "=", cohort.decode("utf8")),
                       ("phenotype".decode("utf8"), "=", phenotype.decode("utf8"))]
            _combine_feather(filters, partition_cols, path, out, schema)


def combine_meta(path, out, phenotypes):
    partition_cols = ["phenotype".decode("utf8")]

    schema = PYARROW_SCHEMA_META

    for phenotype in phenotypes:
        print(phenotype)
        filters = [("phenotype".decode("utf8"), "=", phenotype.decode("utf8"))]
        _combine_feather(filters, partition_cols, path, out, schema)


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # Process input
    parser = argparse.ArgumentParser()
    parser.add_argument('--path')
    parser.add_argument('--out')
    parser.add_argument('--phenotypes')
    parser.add_argument('--cohorts', nargs="+", required=False, default=None)
    args = parser.parse_args(argv)
    # Perform method
    phenotypes_file = pd.read_csv(args.phenotypes)
    if args.cohorts is not None:
        combine_per_cohort(args.path, args.out, phenotypes_file.ID, args.cohorts)
    else:
        combine_meta(args.path, args.out, phenotypes_file.ID)
    return 0


if __name__ == "__main__":
    sys.exit(main())
