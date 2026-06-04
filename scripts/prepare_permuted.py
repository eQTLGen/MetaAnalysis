"""
Created:      10/04/2024
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
def prepare_permuted(base_path, target_path):
    individuals = os.path.join(base_path, "permuted", "EncodedGenotypeData", "individuals")
    genotypes = os.path.join(base_path, "permuted", "EncodedGenotypeData", "genotype")
    probes = os.path.join(base_path, "empirical", "EncodedGenotypeData", "probes")
    if not os.path.exists(target_path):
        os.makedirs(target_path)
    os.symlink(individuals, os.path.join(target_path, "individuals"))
    os.symlink(genotypes, os.path.join(target_path, "genotype"))
    os.symlink(probes, os.path.join(target_path, "probes"))



# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    permuted_path = "/gpfs/space/GI/eQTLGen/PermutedDatasets"

    help_table = pd.read_csv(args.help_table, sep="\t")
    permuted_table = help_table.copy().rename({"dataset_path": "cohort_path"})

    inclusion_dir = "${inclusion_dir}"

    permuted_table["genotype"] = permuted_table.apply(
        lambda row: os.path.join(permuted_path, row["cohort_new_name"], "permuted", "EncodedGenotypeData"), axis=1)
    permuted_table["expression"] = permuted_table.apply(
        lambda row: os.path.join(row["cohort_path"], "permuted", "EncodedPhenotypeData"), axis=1)
    permuted_table["pd"] = permuted_table.apply(
        lambda row: os.path.join(row["cohort_path"], "permuted", "pd_shared"), axis=1)

    permuted_table["encoded"] = 1
    permuted_table["snp_inclusion"] = permuted_table.apply(
        lambda row: os.path.join(inclusion_dir, "{}_SnpsToInclude.txt".format(row["cohort_new_name"])), axis=1)
    permuted_table["gene_inclusion"] = permuted_table.apply(
        lambda row: os.path.join(inclusion_dir, "{}_GenesToInclude.txt".format(row["cohort_new_name"])), axis=1)

    for i, row in permuted_table.iterrows():
        prepare_permuted(row["cohort_path"], row["genotype"])






    return 0


if __name__ == "__main__":
    sys.exit(main())
