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
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.feather as feather


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
MAX_SIZE = 2.5*10**8


SCHEMA_DEFAULT = pa.schema([("variant", pa.string()), ("beta", pa.float64()),
                             ("standard_error", pa.float64()), ("i_squared", pa.float64()),
                             ("sample_size", pa.float64())])


SCHEMA_PER_COHORT = pa.schema(["cohort", pa.string()), ("variant", pa.string()), ("beta", pa.float64()),
                             ("standard_error", pa.float64()),
                             ("sample_size", pa.float64())])


def write_results(results_list, out, out_id, phenotype_list):
    print("Writing results step")

    results = pd.concat(results_list)
    results.columns = [col.decode('utf-8') for col in results.columns]

    results['variant'] = results['variant'].astype(pd.StringDtype())
    results['phenotype'] = results['phenotype'].astype(pd.StringDtype())

    schema = SCHEMA_DEFAULT

    per_cohort = "cohort" in results.columns
    if per_cohort:
        results['cohort'] = results['cohort'].astype(pd.StringDtype())
        schema = SCHEMA_PER_COHORT

    for index, (phenotype, phenotype_results) in enumerate(results.groupby(["phenotype"])):
        phenotype_list.append(phenotype)

        if index % 100 == 0:
            print(index, phenotype)

        partition_dir = os.path.join(out, 'phenotype_{}'.format(phenotype))
        if not os.path.exists(partition_dir):
            os.makedirs(partition_dir)

        pq.write_table(pa.Table.from_pandas(
            phenotype_results.drop("phenotype", inplace=False, axis=1), schema),
            os.path.join(partition_dir, 'results_{}{}.parquet').format(out_id, "_percohort" if per_cohort else ""))

    print("Finished writing step!")


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # Process input
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', nargs='+')
    parser.add_argument('--out')
    parser.add_argument('--out-list')

    args = parser.parse_args(argv)
    # Perform method

    path_base_name = os.path.dirname(args.path[0])

    file_struct = dict()
    file_struct["file"] = list()
    file_struct["chunk"] = list()
    file_struct["iteration"] = list()

    phenotype_list = list()

    for path in args.path:
        file_name = os.path.basename(path)
        splitted = file_name.split("_")
        file_struct['file'].append(path)
        file_struct['chunk'].append(int(splitted[2]))
        file_struct['iteration'].append(int(splitted[4]))

    file_df = pd.DataFrame(file_struct)
    file_df["rank"] = (
        file_df.sort_values(by=['iteration'])
        .groupby("chunk")["iteration"].rank(method="dense"))
    file_df = file_df.astype({'rank':'int'})
    print(file_df)

    out_index = 0

    #    for file_name in glob.glob(os.path.join(args.path, "*.parquet")):
    for group_id, phen_chunk in file_df.groupby(["rank"]):

        results_list = list()

        for i, (index, row) in enumerate(phen_chunk.iterrows()):
            file_name = row['file']
            print("Reading file " + file_name)
            print("(file {}/{})".format(i+1, len(phen_chunk)))
            results_list.append(pq.ParquetFile(file_name).read().to_pandas())

            # Output
            total_sum_lines = sum([len(results) for results in results_list])
            print(total_sum_lines)
            if total_sum_lines > MAX_SIZE:
                write_results(results_list, args.out, out_index, phenotype_list)
                out_index += 1
                results_list = list()

        # Output
        write_results(results_list, args.out, out_index, phenotype_list)
        out_index += 1

    with open(args.out_list, 'w') as opened_out_list:
        for line in phenotype_list:
            opened_out_list.write("{}\n".format(line))
    return 0


if __name__ == "__main__":
    sys.exit(main())
