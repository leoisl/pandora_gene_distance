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
pandora_vcf_ref_to_find_genes = config['pandora_vcf_ref_to_find_genes']
pandora_vcf_ref_to_find_genes = str(Path(pandora_vcf_ref_to_find_genes).absolute())
samples = pd.read_csv(config["samples"])
samples = update_to_absolute_path(samples, ["fasta"])
references = pd.read_csv(config["references"])
references = update_to_absolute_path(references, ["fasta"])
samples_and_refs = pd.concat([samples, references], ignore_index=True)
pandora_multisample_matrices = pd.read_csv(config["pandora_multisample_matrices"])
pandora_multisample_matrices = update_to_absolute_path(pandora_multisample_matrices, ["matrix"])
nb_of_samples = len(samples)


# ======================================================
# Set pandas indexes
# ======================================================
samples = samples.set_index(["id"], drop=False)
references = references.set_index(["id"], drop=False)
samples_and_refs = samples_and_refs.set_index(["id"], drop=False)
pandora_multisample_matrices = pandora_multisample_matrices.set_index(["id"], drop=False)


# ======================================================
# Pipeline files
# ======================================================
files = []

# pandora vcf mappings
pandora_vcf_mapping_files = []
for index, row in samples_and_refs.iterrows():
    id = row["id"]
    pandora_vcf_mapping_files.append(f"{output_folder}/map_pandora_vcf_ref_to_truth_or_ref/{id}.bowtie.sam")
files.extend(pandora_vcf_mapping_files)


# sequences of genes from pandora vcf from the truths/ref
truth_or_ref_gene_sequences = []
for index, row in samples_and_refs.iterrows():
    id = row["id"]
    truth_or_ref_gene_sequences.append(f"{output_folder}/genes_from_truth_or_ref/{id}.csv")
files.extend(truth_or_ref_gene_sequences)


# edit distance files
edit_distances_files = []
for truth_index, row in samples.iterrows():
    truth_id = row["id"]
    for ref_index, row in references.iterrows():
        ref_id = row["id"]
        edit_distances_files.append(f"{output_folder}/edit_distances/{truth_id}~~~{ref_id}.edit_distance.csv")
files.extend(edit_distances_files)
all_edit_distance_files_concatenated = f"{output_folder}/edit_distances/all_edit_distances.csv"
files.append(all_edit_distance_files_concatenated)

files.append(f"{output_folder}/gene_presence_matrix/gene_presence_matrix_based_on_bowtie2")
files.append(f"{output_folder}/gene_presence_matrix/gene_length_matrix")
files = list(set(files))


# fp genes data and plots
for id in [reference_id for reference_id in references["id"] if reference_id.startswith("pandora")]:
    files.extend([
          f"{output_folder}/FP_genes/{id}/gene_and_nb_of_FPs_counted.csv",
          f"{output_folder}/FP_genes/{id}/gene_classification.csv",
          f"{output_folder}/FP_genes/{id}/gene_classification.png",
          f"{output_folder}/FP_genes/{id}/gene_classification_by_sample.csv",
          f"{output_folder}/FP_genes/{id}/gene_classification_by_sample.png",
          f"{output_folder}/FP_genes/{id}/gene_classification_by_gene_length.csv",
          f"{output_folder}/FP_genes/{id}/gene_classification_by_gene_length.png",
          f"{output_folder}/FP_genes/{id}/gene_classification_by_gene_length_normalised.csv",
          f"{output_folder}/FP_genes/{id}/gene_classification_by_gene_length_normalised.png"])

files.extend([
        f"{output_folder}/gene_distance_plots/distribution_of_genes_per_ed.csv",
        f"{output_folder}/gene_distance_plots/distribution_of_genes_per_ed_counts.png",
        f"{output_folder}/gene_distance_plots/distribution_of_genes_per_ed_proportion.png",
        f"{output_folder}/gene_distance_plots/distribution_of_genes_per_nb_of_samples.csv",
        expand(f"{output_folder}/gene_distance_plots/distribution_of_genes_per_nb_of_samples_{{nb_of_sample}}.count.png", nb_of_sample=nb_of_samples),
        expand(f"{output_folder}/gene_distance_plots/distribution_of_genes_per_nb_of_samples_{{nb_of_sample}}.proportion.png", nb_of_sample=nb_of_samples),
        f"{output_folder}/gene_distance_plots/gene_sample_ref_ED_nbsamples_zam.csv",
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
