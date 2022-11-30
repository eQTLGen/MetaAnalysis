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
    parser.add_argument('--path')
    parser.add_argument('--out')

    args = parser.parse_args(argv)
    # Perform method

    file_names_parquet = glob.glob(os.path.join(args.path, "*.parquet"))
    for file_name in file_names_parquet:
        results = pq.ParquetFile(file_name)

        pyarrow_schema = pa.schema(
            [("variant", pa.string()),
             ("phenotype", pa.string()),
             ("beta", pa.float64()),
             ("standard_error", pa.float64()),
             ("i_squared", pa.float64()),
             ("sample_size", pa.float64())])

        # Output

        pq.write_to_dataset(
            table=pa.Table.from_pandas(results, pyarrow_schema),
            root_path=args.out,
            partition_cols=["phenotype"]
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())