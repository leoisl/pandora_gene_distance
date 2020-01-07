from pathlib import Path
from snakemake.utils import min_version
import pandas as pd
import pysam
import itertools

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
recall_reports = f"{config['pandora_eval_output_folder']}/recall/reports"
truth_pairs = [(truth_1, truth_2) for truth_1, truth_2 in itertools.combinations(sorted(truth_assemblies["id"]), r=2)]


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

# pandora vcf mappings
pandora_vcf_mapping_files = []
for index, row in assemblies_and_refs.iterrows():
    id = row["id"]
    pandora_vcf_mapping_files.append(f"{output_folder}/map_pandora_vcf_ref_to_truth_or_ref/{id}.bowtie.sam")
files.extend(pandora_vcf_mapping_files)


# sequences of genes from pandora vcf from the truths/ref
truth_or_ref_gene_sequences = []
for index, row in assemblies_and_refs.iterrows():
    id = row["id"]
    truth_or_ref_gene_sequences.append(f"{output_folder}/genes_from_truth_or_ref/{id}.csv")
files.extend(truth_or_ref_gene_sequences)


# edit distance files
edit_distances_files = []
for truth_index, row in truth_assemblies.iterrows():
    truth_id = row["id"]
    for ref_index, row in references.iterrows():
        ref_id = row["id"]
        edit_distances_files.append(f"{output_folder}/edit_distances/{truth_id}~~~{ref_id}.edit_distance.csv")
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

        for truth_1, truth_2 in [pair for pair in truth_pairs if truth_id in pair]:
            variant_distance_files.append(f"{output_folder}/get_variant_recall_score_distance_csv/{truth_id}~~~{ref_id}/{truth_1}_and_{truth_2}.get_variant_recall_score_distance.csv")
            variant_distance_files.append(f"{output_folder}/get_variant_recall_score_distance_csv/{truth_id}~~~{ref_id}/{truth_1}_and_{truth_2}.get_variant_recall_score_distance.unmapped.csv")
files.extend(variant_distance_files)

gene_truth_ref_precision_proportion_distance_files = []
gene_truth_ref_recall_proportion_distance_files = []
for truth_index, row in truth_assemblies.iterrows():
    truth_id = row["id"]
    for ref_index, row in references.iterrows():
        ref_id = row["id"]
        gene_truth_ref_precision_proportion_distance_files.append(f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/{truth_id}~~~{ref_id}.gene_truth_ref_precision_proportion_distance.csv")
        for truth_1, truth_2 in [pair for pair in truth_pairs if truth_id in pair]:
            gene_truth_ref_recall_proportion_distance_files.append(f"{output_folder}/get_gene_truth_ref_recall_proportion_distance/{truth_id}~~~{ref_id}/{truth_1}_and_{truth_2}.gene_truth_ref_recall_proportion_distance.csv")
files.extend(gene_truth_ref_precision_proportion_distance_files)
files.extend(gene_truth_ref_recall_proportion_distance_files)
files.append(f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv")
files.append(f"{output_folder}/get_gene_truth_ref_recall_proportion_distance/all_gene_truth_ref_recall_proportion_distance.csv")


files.append(f"{output_folder}/gene_distance_precision.png")
files.append(f"{output_folder}/gene_distance_precision.filtered.png")
files.append(f"{output_folder}/gene_distance_recall.png")
files.append(f"{output_folder}/gene_distance_recall.filtered.png")

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