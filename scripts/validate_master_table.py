"""
Created:      30/01/2024
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
from __future__ import print_function
import os
import sys
import argparse
import shutil
import re

import numpy as np
import pandas as pd

from bs4 import BeautifulSoup


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

def replace_content_in_file(file_path, old_string, new_string, replace_func):
    replace_func(file_path, old_string, new_string)

def copy_and_rename_files_and_directories(src, dst, old_string, new_string, path_function_dict):
    for root, dirs, files in os.walk(src):
        # Calculate the corresponding destination path
        relative_path = os.path.relpath(root, src)
        dest_path = os.path.join(dst, relative_path)

        # Rename directories in the path
        new_dest_path = dest_path

        if not os.path.exists(new_dest_path):
            os.makedirs(new_dest_path)
            print('Created directory: {}'.format(new_dest_path))

        # Rename and copy files
        for filename in files:
            old_file_path = os.path.join(root, filename)
            new_filename = filename.replace(old_string, new_string)
            new_file_path = os.path.join(new_dest_path, new_filename)

            shutil.copy2(old_file_path, new_file_path)
            print('Copied file: {} -> {}'.format(old_file_path, new_file_path))

            # Check if file path matches any pattern in the dictionary
            for path_regex, func in path_function_dict.items():
                path_regex_comp = path_regex.replace('<new_string>', new_string)

                if re.search(path_regex_comp, new_file_path):
                    replace_content_in_file(new_file_path, old_string, new_string, func)
                    print('Processed file: {}'.format(new_file_path))

        # Rename and copy directories
        for dirname in dirs:
            new_dirname = dirname.replace(old_string, new_string)
            new_dir_path = os.path.join(new_dest_path, new_dirname)

            if not os.path.exists(new_dir_path):
                os.makedirs(new_dir_path)
                print('Created directory: {}'.format(new_dir_path))


# Example usage
path_function_dict = {
    r'EncodedGenotypeData/individuals/<new_string>\.h5$': replace_in_individuals_file,
    r'EncodedPhenotypeData/\d+_<new_string>\.csv$': replace_in_text_file,
    r'pd_shared/<new_string>_metadata\.npy$': replace_in_metadata_file
}

def replace_in_metadata_file(file_path, old_string, new_string):
    metadata = np.load(file_path).item()
    metadata["names"] = [name.replace(old_string, new_string) for name in metadata["names"]]
    metadata["id"] = np.char.replace(metadata["id"], old_string, new_string)
    np.save(file_path, metadata)

# Example replacement function for text files
def replace_in_text_file(file_path, old_string, new_string):
    temp_file_path = file_path + '.tmp'
    with open(file_path, 'r') as infile, open(temp_file_path, 'w') as outfile:
        for line in infile:
            new_line = line.replace(old_string, new_string)
            outfile.write(new_line)
    os.remove(file_path)
    os.rename(temp_file_path, file_path)

# Example replacement function for HDF5 files (assuming h5py is available)
def replace_in_individuals_file(file_path, old_string, new_string):
    import h5py
    with h5py.File(file_path, 'r+') as h5file:
        # This is a placeholder for the actual logic to modify the HDF5 file
        # Replace this with the actual code needed to update the HDF5 content
        dataset = h5file["individuals"]["table"]
        modified = dataset[...]
        for i, row in enumerate(modified[:]):
            modified[i][1] = np.char.replace(modified[i][1], old_string, new_string)
        dataset[...] = modified

def _prune_covariates(master_table, selected_covariates):
    pruned_covariates = dict()
    explained_variance = dict()
    for index, row in master_table.iterrows():
        cohort_new_name = row["cohort_new_name"]
        cohort = row["cohort"]
        print(cohort_new_name)
        metadata = np.load(os.path.join(row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"]))).item()
        a_cov = np.load(os.path.join(row["partial_derivatives"], "{}_a_cov.npy".format(row["cohort"])))
        encoded_covariates = np.array(metadata['names'])
        print(selected_covariates[cohort])
        covariate_indices = np.where(np.in1d(encoded_covariates, selected_covariates[cohort]))[0]
        assert len(selected_covariates[cohort]) == len(covariate_indices)
        a_selected = a_cov[np.ix_(covariate_indices, covariate_indices)]
        covariates_indices_mask = np.full_like(covariate_indices, True)
        print(np.linalg.det(a_selected))
        orig_eigenvalues, _ = np.linalg.eig(a_selected)
        print("condition number:", np.linalg.cond(a_selected))
        while np.linalg.det(a_selected) == 0:
            eigenvalues = np.full_like(covariate_indices, np.NAN, dtype=np.double)
            eigenvalues[covariates_indices_mask==1], _ = np.linalg.eig(a_selected)
            print(eigenvalues)
            if any(eigenvalues[covariates_indices_mask==1] == 0):
                covariates_indices_mask[np.where(eigenvalues[covariates_indices_mask==1] == 0)] = False
            else:
                covariates_indices_mask[np.nanargmin(np.abs(eigenvalues))] = False
            print(covariate_indices[covariates_indices_mask==0])
            print(encoded_covariates[covariate_indices[covariates_indices_mask==1]])
            a_selected = a_cov[np.ix_(
                covariate_indices[covariates_indices_mask==1],
                covariate_indices[covariates_indices_mask==1])]
            print("condition number:", np.linalg.cond(a_selected))
        pruned = list(encoded_covariates[covariate_indices[covariates_indices_mask==1]])
        print(pruned)
        pruned_covariates[cohort] = pruned
        explained_variance[cohort] = sum(orig_eigenvalues[covariates_indices_mask==0] / sum(orig_eigenvalues))
    return pruned_covariates

def prune_covariates(args):
    master_table = pd.read_csv(mastertable_path, sep="\t")
    covariates_table = pd.read_table("covariates_long_extended_2024-08-08.txt")
    analysis_columns = [col for col in covariates_table.columns if col.startswith('analysis')]
    for analysis in analysis_columns:
        selected_covariates = covariates_table[covariates_table[analysis] == True].groupby('cohort')['covariate'].apply(list).to_dict()
        pruned_covariates = _prune_covariates(master_table, selected_covariates)
        with open("covariates_curated_{}_2024-08-08.txt".format(analysis), 'w') as opened:
            for key, covars in pruned_covariates.items():
                print("{}\t{}".format(key, "\t".join(covars)), file=opened)




def align_mapper_paths(args):
    master_table = pd.read_csv(args.mastertable, sep="\t")
    mapper_path = ""
    # Perform method
    master_table_fixed = master_table.copy()

    mapper_file_dict = dict()
    mapper_files = os.listdir(mapper_path)
    for mapper_file in mapper_files:
        if mapper_file.startswith("values"):
            print(mapper_file)
            dataset = re.match("values_1000G-30x_ref_(.+)\.npy", mapper_file).group(1)
            print(dataset)
            mapper_file_path = os.path.dirname(os.path.realpath(os.path.join(mapper_path, mapper_file)))
            mapper_file_dict[dataset] = mapper_file_path

    for index, row in master_table.iterrows():
        cohort_new_name = row["cohort_new_name"]
        cohort = row["cohort"]
        path = row["base"]


    merged_table = pd.concat([master_table.set_index('cohort'), pd.Series(mapper_file_dict, name="mapper_path")], axis=1, join='inner')

    upload_path = "/gpfs/space/GI/eQTLGen/UploadedDatasets"
    fixed_path = "/gpfs/space/GI/eQTLGen/FixedDatasets"

    for index, row in master_table.iterrows():
        print(row["cohort_new_name"])



def fix_ensg_version_ids(args):
    master_table = pd.read_csv(args.mastertable, sep="\t")
    # Perform method
    master_table_fixed = master_table.copy()

    upload_path = "/gpfs/space/GI/eQTLGen/UploadedDatasets"
    fixed_path = "/gpfs/space/GI/eQTLGen/FixedDatasets"

    for index, row in master_table.iterrows():
        print(row["cohort_new_name"])
        fixed_row = master_table_fixed.iloc[index].copy()
        metadata = np.load(os.path.join(row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"]))).item()
        phenotypes = metadata['phenotype']
        expression_files = os.listdir(row["expression"])
        gene_inclusion_file = row['gene_inclusion']
        if not np.all(np.array(["." in phen for phen in phenotypes])):
            print("No dots in phenotype identifiers")
            invalid = [phen for phen in phenotypes if not re.match("^ENSG\d{11}$", phen)]
            print("Number of non-matching gene identifiers {}".format(len(invalid)))
            validate_uniqueness(phenotypes)
        else:
            print("Dots in phenotype identifiers, fixing files...")
            if not row["expression"].startswith(upload_path):
                print("WARNING: skipping. data from this cohort already is processed")
                continue
            phenotype_fixed_dict = {phen: phen.split('.')[0] for phen in phenotypes}
            new_metadata = dict(metadata)
            new_phenotypes = [phen.split('.')[0] for phen in phenotypes]
            new_metadata['phenotype'] = new_phenotypes
            print("Replacing gene ids as follows:")
            for key, value in list(phenotype_fixed_dict.items())[:10]:
                print("{:>12} -> {}".format(key, value))
            print("{:>12} -> {} (and {} more)".format("...", "...", len(phenotype_fixed_dict)-10))
            validate_uniqueness(new_phenotypes)
            fixed_row["expression"] = row["expression"].replace(upload_path, fixed_path)
            fixed_row["partial_derivatives"] = row["partial_derivatives"].replace(upload_path, fixed_path)
            fixed_row["gene_inclusion"] = gene_inclusion_file.replace("_GenesToInclude", "_GenesToInclude_Fixed")
            # Replace metadata
            if not os.path.exists(fixed_row["partial_derivatives"]):
                os.makedirs(fixed_row["partial_derivatives"])
            if not os.path.exists(fixed_row["expression"]):
                os.makedirs(fixed_row["expression"])
            for pd_file in os.listdir(row["partial_derivatives"]):
                if pd_file.endswith("metadata.npy") or os.path.exists(os.path.join(fixed_row["partial_derivatives"], pd_file)):
                    continue
                os.symlink(os.path.join(row["partial_derivatives"], pd_file),
                           os.path.join(fixed_row["partial_derivatives"], pd_file))
            print("Writing file {} to {} ...".format(row["partial_derivatives"], fixed_row["partial_derivatives"]))
            np.save(os.path.join(fixed_row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"])),
                    new_metadata)
            if os.path.exists(gene_inclusion_file):
                inclusion_list = pd.read_csv(gene_inclusion_file, sep="\t", header=0)
                inclusion_list.replace(to_replace=phenotype_fixed_dict, inplace=True)
                inclusion_list.to_csv(os.path.join(fixed_row["gene_inclusion"]), sep="\t", header=True, index=False)
            else:
                print("Gene inclusion file does not exist!", gene_inclusion_file)
            print("Path to fixed expression data: \n{}".format(fixed_row["expression"]))
            for expression_file in expression_files:
                expression = pd.read_csv(os.path.join(row["expression"], expression_file), sep="\t")
                expression.rename(columns=phenotype_fixed_dict, inplace=True)
                print("Writing file {} to {} ...".format(expression_file, fixed_row["expression"]))
                expression.to_csv(os.path.join(fixed_row["expression"], expression_file), sep="\t", index=False)
            master_table_fixed.iloc[index] = fixed_row

    master_table_fixed.to_csv("mastertable_empirical_2024-02-28_fixedSampleNames_fixedGeneIds.txt", sep="\t", index=False)


def validate_uniqueness(phenotypes):
    all_unique = np.all(np.unique(np.array(phenotypes), return_counts=True)[1] == 1)
    print("Number of keys {}, number of unique values {}".format(
        len(phenotypes), len(np.unique(np.array(phenotypes)))))
    if all_unique:
        print("All phenotypes identifiers are unique")
    else:
        raise ValueError("Not all phenotype identifiers are unique")


def extract_explained_variance(report_html_path):
    with open(report_html_path) as file:
        soup = BeautifulSoup(file, 'html.parser')
    tables = pd.read_html(str(soup))
    return tables[2]


def explained_variance_tables():
    master_table = pd.read_csv(args.mastertable, sep="\t")

    explained_variance_tables = dict()

    for index, row in master_table.iterrows():
        print(row["cohort"], row["cohort_new_name"])
        report_path = glob.glob(os.path.join(row['base'], "Report_DataQc_*.html"))[0]
        print(report_path)
        explained_variance_tables[row["cohort"]] = extract_explained_variance(report_path)

    combined_explained_variance_table = pd.concat(explained_variance_tables, names=["cohort"]).reset_index(level=0)
    combined_explained_variance_table['covariate'] = (
        combined_explained_variance_table.apply(lambda row: "_Exp".join([row.cohort, row.PC]), axis=1))
    combined_explained_variance_table = (
        combined_explained_variance_table.merge(master_table[["cohort", "cohort_new_name"]])
        .sort_values(by='explained_variance', ascending=False))
    combined_explained_variance_table["cumulative_explained_variance"] = (
        combined_explained_variance_table.groupby("cohort").cumsum())
    combined_covariate_table = covariates_selected.merge(combined_explained_variance_table, on=["cohort_new_name", "covariate"], how="outer")


def export_covariates_table():
    master_table = pd.read_csv(args.mastertable, sep="\t")

    covariates = dict()

    with open("covariates_export.txt", 'w') as opened:
        for index, row in master_table.iterrows():
            print(row["cohort"], row["cohort_new_name"])
            metadata = np.load(os.path.join(row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"]))).item()
            covariates[row["cohort"]] = metadata['names']
            print("{}\t{}".format(row["cohort"], "\t".join(metadata['names'])), file=opened)

    # After writing this file, use the following in R
    # dat <- fread("eQTLGen/freeze3/eqtl_mapping/input/covariates_export.txt", fill=T) %>%
    #   pivot_longer(cols=-V1) %>% select(c("cohort"="V1", "covariate"="value")) %>%
    #   mutate(analysis1_4GenPCNoExpPC = case_when(str_ends(covariate, "GenPC[1234]|intercept|GenSex") ~ T, TRUE ~ F)) %>%
    #   filter(covariate!="")
    # fwrite(dat, "eQTLGen/freeze3/eqtl_mapping/input/covariates_long.txt", col.names=T, row.names=F, quote=F, sep="\t")

    # After curation in Excel, execute the following in R
    # table <- read_excel("covariates_long_2024-04-05.xlsx")
    # fwrite(table %>% filter(analysis2_4GenPCNoExpPC) %>%
    # group_by(cohort) %>%
    # summarise(covariate_list=list(covariate)) %>%
    # unnest_wider(covariate_list),
    # "covariates_curated_2024-04-16.txt", sep="\t", quote=F, row.names=F, col.names=F)

class MasterTableGenerator:
    def __init__(self, base_table, inclusion_dir, permuted=None):
        self.base_table = base_table.copy()
        self.inclusion_dir = inclusion_dir
        self.type = "permuted" if permuted else "empirical"
        self.permuted = permuted
    def _generate(self, row):
        genotype_base_path = row.dataset_path if self.permuted is None else os.path.join(self.permuted, row.dataset_name_new)
        master_table_row = dict()
        master_table_row["cohort"] = row.dataset_name
        master_table_row["cohort_new_name"] = row.dataset_name_new
        master_table_row["encoded"] = 1
        master_table_row["genotype"] = os.path.join(genotype_base_path, self.type, "EncodedGenotypeData")
        master_table_row["expression"] = os.path.join(row.dataset_path, self.type, "EncodedPhenotypeData")
        master_table_row["partial_derivatives"] = os.path.join(row.dataset_path, self.type, "pd_shared")
        master_table_row["snp_inclusion"] = os.path.join(self.inclusion_dir, "".join([row.dataset_name_new, "_SnpsToInclude.txt"]))
        master_table_row["gene_inclusion"] = os.path.join(self.inclusion_dir, "".join([row.dataset_name_new, "_GenesToInclude.txt"]))
        master_table_row["base"] = row.dataset_path
        if self.type == "permuted":
            print(master_table_row["genotype"])
            self.prepare_permuted(master_table_row["base"], master_table_row["genotype"])
        return master_table_row
    def prepare_permuted(self, base_path, target_path):
        individuals = os.path.join(base_path, "permuted", "EncodedGenotypeData", "individuals")
        genotypes = os.path.join(base_path, "permuted", "EncodedGenotypeData", "genotype")
        probes = os.path.join(base_path, "empirical", "EncodedGenotypeData", "probes")
        if not os.path.exists(target_path):
            os.makedirs(target_path)
        if not os.path.exists(os.path.join(target_path, "individuals")):
            os.symlink(individuals, os.path.join(target_path, "individuals"))
        if not os.path.exists(os.path.join(target_path, "genotype")):
            os.symlink(genotypes, os.path.join(target_path, "genotype"))
        if not os.path.exists(os.path.join(target_path, "probes")):
            os.symlink(probes, os.path.join(target_path, "probes"))
    def run(self):
        return pd.DataFrame.from_records(self.base_table.apply(self._generate, axis=1))



class MasterTableValidator:
    def __init__(self, master_table, output_path):
        self.master_table = master_table.copy()
        self.output_path = output_path

    def validate(self):
        self.master_table.apply(self._validate, axis=1)

    def _validate(self, row):
        fixed_row = row.copy()
        metadata = np.load(os.path.join(row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"]))).item()
        self._validate_gene_ids(fixed_row, metadata)
        self._validate_sample_ids(fixed_row, metadata)

    def _validate_sample_ids(self, row, metadata):
        sample_ids = metadata['id']
        expression_files = os.listdir(row["expression"])
        expression = pd.read_csv(os.path.join(row["expression"], expression_files[0]), sep="\t")
        if np.all(expression.id == sample_ids):
            print("Sample IDs match, continuing")
        else:
            print("Mismatches found, fixing files...")
            fixed_row["expression"] = row["expression"].replace(upload_path, fixed_path)
            print("Path to fixed expression data: \n{}".format(fixed_row["expression"]))
            try:
                os.makedirs(fixed_row["expression"])
            except OSError as err:
                print(err.message)
            for file in expression_files:
                expression = pd.read_csv(os.path.join(row["expression"], file), sep="\t")
                expression.id = sample_ids
                print("Writing file {} to {} ...".format(file, fixed_row["expression"]))
                expression.to_csv(os.path.join(fixed_row["expression"], file), sep="\t", index=False)

    def _validate_gene_ids(self, row, metadata):
        print(row["cohort_new_name"])
        metadata = np.load(os.path.join(row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"]))).item()
        phenotypes = metadata['phenotype']
        if not np.all(np.array(["." in phen for phen in phenotypes])):
            print("No dots in phenotype identifiers")
            validate_uniqueness(phenotypes)
        else:
            print("Dots in phenotype identifiers, fixing files...")


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv
    # Process input
    parser = argparse.ArgumentParser()
    parser.add_argument('--mastertable')
    args = parser.parse_args(argv)

    upload_path = "/gpfs/space/GI/eQTLGen/UploadedDatasets"
    fixed_path = "/gpfs/space/GI/eQTLGen/PermutedDatasets"

    master_table = pd.read_csv(args.mastertable, sep="\t")
    # Perform method
    master_table_fixed = master_table.copy()

    for index, row in master_table_fixed.iterrows():
        print(row["cohort_new_name"])
        #fixed_row = master_table_fixed.iloc[index]
        sample_ids = np.load(os.path.join(row["partial_derivatives"], "{}_metadata.npy".format(row["cohort"])),
                             allow_pickle=True, encoding='latin1').item()['id']
        print(sample_ids[0], len(sample_ids), len(np.unique(sample_ids)))
        expression_files = os.listdir(row["expression"])
        expression = pd.read_csv(os.path.join(row["expression"], expression_files[0]), sep="\t")
        print(expression.id[0], len(expression.id))
        if np.all(expression.id == sample_ids):
            print("Sample IDs match, continuing")
        else:
            print("Mismatches found, fixing files...")
            if not row["expression"].startswith(upload_path):
                print("WARNING: skipping. data from this cohort already is processed")
                continue
            fixed_row["expression"] = row["expression"].replace(upload_path, fixed_path)
            print("Path to fixed expression data: \n{}".format(fixed_row["expression"]))
            try:
                os.makedirs(fixed_row["expression"])
            except OSError as err:
                print(err.message)
            for file in expression_files:
                expression = pd.read_csv(os.path.join(row["expression"], file), sep="\t")
                expression.id = sample_ids
                print("Writing file {} to {} ...".format(file, fixed_row["expression"]))
                expression.to_csv(os.path.join(row["expression"], file), sep="\t", index=False)
            master_table_fixed.iloc[index] = fixed_row





    master_table_fixed.to_csv("mastertable_empirical_2024-01-31.txt", sep="\t")

    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
