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
import cPickle as pkl


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
MAX_SIZE = 4*10**8

# Classes

# Functions
def write_results(results_list, out):
    print("Writing results step")
    results = pd.concat(results_list)
    for index, (phenotype, phenotype_results) in enumerate(results.groupby(["phenotype"])):
        if index % 100 == 0:
            print(index, phenotype, end="\r")
        with open(os.path.join(out, 'phenotype_{}.pkl').format(phenotype), 'a') as f:
            pkl.dump(phenotype_results.drop('phenotype', inplace=False, axis=1), f)
        #(phenotype_results
        #    .drop('phenotype', inplace=False, axis=1)
        #    .to_csv(
        #        os.path.join(out, 'phenotype_{}.csv').format(phenotype),
        #        index=False, mode='a'))
    print("Finished writing step!")


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # Process input
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', nargs='+')
    parser.add_argument('--out')

    args = parser.parse_args(argv)
    # Perform method

    path_base_name = os.path.dirname(args.path[0])

    file_struct = dict()
    file_struct["file"] = list()
    file_struct["chunk"] = list()
    file_struct["iteration"] = list()

    for path in args.path:
        file_name = os.path.basename(path)
        splitted = file_name.split("_")
        file_struct['file'].append(path)
        file_struct['chunk'].append(int(splitted[2]))
        file_struct['iteration'].append(int(splitted[4]))

    file_df = pd.DataFrame(file_struct)
    file_df["rank"] = file_df.groupby("chunk")["iteration"].rank(method="dense")
    print(file_df)

    results_list = list()

#    for file_name in glob.glob(os.path.join(args.path, "*.parquet")):
    for phen_chunk in file_df.groupby(["rank"]):
        for i, (index, row) in enumerate(phen_chunk.iterrows()):
            file_name = row['file']
            print("Reading file " + file_name)
            print("(file {}/{})".format(i+1, len(phen_chunk)))
            results_list.append(pq.ParquetFile(file_name).read())

            # Output
            sum1 = sum([len(results) for results in results_list])
            print(sum1)
            if sum1 > MAX_SIZE:
                #write_results(results_list, args.out)
                pq.write_to_dataset(
                    table=pa.concat_tables(results_list),
                    root_path=args.out,
                    partition_cols=["phenotype"]
                )
                results_list = list()

        #write_results(results_list, args.out)

        pq.write_to_dataset(
            table=pa.concat_tables(results_list),
            root_path=args.out,
            partition_cols=["phenotype"]
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
