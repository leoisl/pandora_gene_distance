from pathlib import Path
import pandas as pd


# ======================================================
# Helper functions
# ======================================================
def update_to_absolute_path_core(path_series):
    return path_series.apply(lambda path: str(Path(path).absolute()))
def update_to_absolute_path(df, columns):
    for column in columns:
        df[column] = update_to_absolute_path_core(df[column])
    return df


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

# ======================================================
# Global variables
# ======================================================
output_folder = config['output_folder']
pandora_vcfs = pd.read_csv(config['pandora_vcfs'])
pandora_vcfs = update_to_absolute_path(pandora_vcfs, ["fasta", "matrix"])
samples = pd.read_csv(config["samples"])
samples = update_to_absolute_path(samples, ["fasta"])
references = pd.read_csv(config["references"])
references = update_to_absolute_path(references, ["fasta"])
samples_and_refs = pd.concat([samples, references], ignore_index=True)
samples_and_refs_and_pandora_vcfs = pd.concat([samples, references, pandora_vcfs[["id", "fasta"]]], ignore_index=True)
nb_of_samples = range(1, len(samples)+1)


# ======================================================
# Set pandas indexes
# ======================================================
samples = samples.set_index(["id"], drop=False)
references = references.set_index(["id"], drop=False)
pandora_vcfs = pandora_vcfs.set_index(["id"], drop=False)
samples_and_refs = samples_and_refs.set_index(["id"], drop=False)
samples_and_refs_and_pandora_vcfs = samples_and_refs_and_pandora_vcfs.set_index(["id"], drop=False)


# ======================================================
# Pipeline files
# ======================================================
files = []

# pandora vcf mappings
pandora_vcf_mapping_files = []
for pandora_id in pandora_vcfs.index:
    for index, row in samples_and_refs.iterrows():
        sample_ref_id = row["id"]
        pandora_vcf_mapping_files.append(f"{output_folder}/map_pandora_vcf_ref_to_truth_or_ref/{pandora_id}/{sample_ref_id}.bowtie.sam")
files.extend(pandora_vcf_mapping_files)

# sequences of genes from pandora vcf from the truths/ref
truth_or_ref_gene_sequences = []
for pandora_id in pandora_vcfs.index:
    for index, row in samples_and_refs.iterrows():
        sample_ref_id = row["id"]
        truth_or_ref_gene_sequences.append(f"{output_folder}/genes_from_truth_or_ref/{pandora_id}/{sample_ref_id}.csv")
    files.extend(truth_or_ref_gene_sequences)

# edit distance files
edit_distances_files = {}
for pandora_id in pandora_vcfs.index:
    edit_distances_files_for_this_pandora_id = []
    for truth_index, row in samples.iterrows():
        truth_id = row["id"]
        for ref_index, row in references.iterrows():
            ref_id = row["id"]
            edit_distances_files_for_this_pandora_id.append(f"{output_folder}/edit_distances/{pandora_id}/{truth_id}~~~{ref_id}.edit_distance.csv")
        edit_distances_files_for_this_pandora_id.append(f"{output_folder}/edit_distances/{pandora_id}/{truth_id}~~~{pandora_id}.edit_distance.csv")
    edit_distances_files[pandora_id] = edit_distances_files_for_this_pandora_id
    files.extend(edit_distances_files_for_this_pandora_id)
all_edit_distance_files_concatenated = expand(f"{output_folder}/edit_distances/{{pandora_id}}/all_edit_distances.csv", pandora_id = pandora_vcfs.index)
files.extend(all_edit_distance_files_concatenated)

files.extend(expand(f"{output_folder}/gene_presence_matrix/{{pandora_id}}/gene_presence_matrix_based_on_bowtie2", pandora_id = pandora_vcfs.index))
files.extend(expand(f"{output_folder}/gene_presence_matrix/{{pandora_id}}/gene_length_matrix", pandora_id = pandora_vcfs.index))
files = list(set(files))


# fp genes data and plots
for pandora_id in pandora_vcfs.index:
    files.extend([
          f"{output_folder}/FP_genes/{pandora_id}/gene_and_nb_of_FPs_counted.csv",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification.csv",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification.png",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification_by_sample.csv",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification_by_sample.png",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification_by_gene_length.csv",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification_by_gene_length.png",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification_by_gene_length_normalised.csv",
          f"{output_folder}/FP_genes/{pandora_id}/gene_classification_by_gene_length_normalised.png"])

for pandora_id in pandora_vcfs.index:
    files.extend([
            f"{output_folder}/gene_distance_plots/{pandora_id}/distribution_of_genes_per_ed.csv",
            f"{output_folder}/gene_distance_plots/{pandora_id}/distribution_of_genes_per_ed_counts.png",
            f"{output_folder}/gene_distance_plots/{pandora_id}/distribution_of_genes_per_ed_proportion.png",
            f"{output_folder}/gene_distance_plots/{pandora_id}/distribution_of_genes_per_nb_of_samples.csv",
            [f"{output_folder}/gene_distance_plots/{pandora_id}/distribution_of_genes_per_nb_of_samples_{nb_of_sample}.count.png" for nb_of_sample in nb_of_samples],
            [f"{output_folder}/gene_distance_plots/{pandora_id}/distribution_of_genes_per_nb_of_samples_{nb_of_sample}.proportion.png" for nb_of_sample in nb_of_samples],
            f"{output_folder}/gene_distance_plots/{pandora_id}/gene_sample_ref_ED_nbsamples_zam.csv",
    ])



# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("rules/")
include: str(rules_dir / "indexing_and_mapping.smk")
include: str(rules_dir / "finding_distance_between_loci_in_assemblies_and_refs.smk")
include: str(rules_dir / "FP_genes.smk")
