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
all_assemblies = pd.concat([truth_assemblies, references], ignore_index=True)
# pandora_evaluation_folder = config["pandora_evaluation_folder"]




# ======================================================
# Set pandas indexes
# ======================================================
truth_assemblies = truth_assemblies.set_index(["id"], drop=False)
references = references.set_index(["id"], drop=False)
all_assemblies = all_assemblies.set_index(["id"], drop=False)


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
def get_gene_mapping_files(df, all_genes_in_vcf_ref):
    gene_mapping_files = []
    for index, row in df.iterrows():
        id = row["id"]
        for gene in all_genes_in_vcf_ref:
            gene_mapping_files.append(f"{output_folder}/map_gene_from_vcf_ref_to_truth_or_ref/{gene}~~~{id}.sam")
    return gene_mapping_files

files.extend(get_gene_mapping_files(truth_assemblies, all_genes_in_vcf_ref))
files.extend(get_gene_mapping_files(references, all_genes_in_vcf_ref))


print(f"Required files: {files}")




# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("rules/")
include: str(rules_dir / "finding_distance_between_loci_in_assemblies_and_refs.smk")