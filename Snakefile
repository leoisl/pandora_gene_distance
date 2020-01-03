from pathlib import Path
from snakemake.utils import min_version
import pandas as pd
import pysam

min_version("5.1.0")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

# ======================================================
# Global variables
# ======================================================
output_folder = config['output_folder']
pandora_vcf_ref = config['pandora_vcf_ref']
truth_assemblies = pd.read_csv(config["truth_assemblies"])
references = pd.read_csv(config["references"])
assemblies_and_refs = pd.concat([truth_assemblies, references], ignore_index=True)
precision_reports = f"{config['pandora_eval_output_folder']}/precision/reports_from_probe_mappings"
recall_reports = f"{config['pandora_eval_output_folder']}/recall/reports_from_probe_mappings"



# ======================================================
# Set pandas indexes
# ======================================================
truth_assemblies = truth_assemblies.set_index(["id"], drop=False)
references = references.set_index(["id"], drop=False)
assemblies_and_refs = assemblies_and_refs.set_index(["id"], drop=False)


# ======================================================
# Pipeline files
# ======================================================
files = []

# add genes in vcf_ref
all_genes_in_vcf_ref = []
all_genes_filepaths_in_vcf_ref = []
with pysam.FastaFile(pandora_vcf_ref) as fasta_file:
    for gene in fasta_file.references:
        all_genes_in_vcf_ref.append(gene)
        all_genes_filepaths_in_vcf_ref.append(f"{output_folder}/genes_from_vcf_ref/{gene}.fa")
# files.extend(all_genes_filepaths_in_vcf_ref)

# gene mappings
gene_mapping_files = []
for index, row in assemblies_and_refs.iterrows():
    id = row["id"]
    for gene in all_genes_in_vcf_ref:
        gene_mapping_files.append(f"{output_folder}/map_gene_from_vcf_ref_to_truth_or_ref/{gene}~~~{id}.bowtie.sam")
files.extend(gene_mapping_files)


# genes in truth/ref
genes_in_truth_or_ref = []
for index, row in assemblies_and_refs.iterrows():
    id = row["id"]
    for gene in all_genes_in_vcf_ref:
        genes_in_truth_or_ref.append(f"{output_folder}/genes_from_truth_or_ref/{gene}~~~{id}.csv")
files.extend(genes_in_truth_or_ref)


# edit distance files
edit_distances_files = []
for gene in all_genes_in_vcf_ref:
    for truth_index, row in truth_assemblies.iterrows():
        truth_id = row["id"]
        for ref_index, row in references.iterrows():
            ref_id = row["id"]
            edit_distances_files.append(f"{output_folder}/edit_distances/{gene}~~~{truth_id}~~~{ref_id}.edit_distance.csv")
files.extend(edit_distances_files)
all_edit_distance_files_concatenated = f"{output_folder}/edit_distances/all_edit_distances.csv"
files.append(all_edit_distance_files_concatenated)

# variant distances files
variant_distance_files = []
for truth_index, row in truth_assemblies.iterrows():
    truth_id = row["id"]
    for ref_index, row in references.iterrows():
        ref_id = row["id"]
        variant_distance_files.append(f"{output_folder}/get_variant_precision_score_distance_csv/{truth_id}~~~{ref_id}.get_variant_precision_score_distance.csv")
        variant_distance_files.append(f"{output_folder}/get_variant_precision_score_distance_csv/{truth_id}~~~{ref_id}.get_variant_precision_score_distance.unmapped.csv")
files.extend(variant_distance_files)

gene_truth_ref_precision_proportion_distance_files = []
for truth_index, row in truth_assemblies.iterrows():
    truth_id = row["id"]
    for ref_index, row in references.iterrows():
        ref_id = row["id"]
        gene_truth_ref_precision_proportion_distance_files.append(f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/{truth_id}~~~{ref_id}.gene_truth_ref_precision_proportion_distance.csv")
files.extend(gene_truth_ref_precision_proportion_distance_files)
files.append(f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv")


files.append(f"{output_folder}/gene_distance.pdf")

files = list(set(files))

# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("rules/")
include: str(rules_dir / "indexing_and_mapping.smk")
include: str(rules_dir / "finding_distance_between_loci_in_assemblies_and_refs.smk")
include: str(rules_dir / "generating_csv_variant_and_gene_with_distances.smk")
include: str(rules_dir / "plot.smk")