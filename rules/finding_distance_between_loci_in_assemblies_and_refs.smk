import editdistance
import pandas as pd
from Bio.Seq import Seq

rule get_truth_or_ref_gene_sequences:
    input:
        truth_or_ref = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}",
        truth_or_ref_index = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}.fai",
        pandora_vcf_ref_mapped_to_truth_or_ref_sam_file = rules.map_pandora_vcf_ref_to_truth_or_ref_using_bowtie.output.sam_file
    output:
         truth_or_ref_gene_sequences = f"{output_folder}/genes_from_truth_or_ref/{{id}}.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/get_truth_or_ref_gene_sequences/{id}.log"
    script:
        "../scripts/get_truth_or_ref_gene_sequences.py"


def reverse_complement(sequence):
    biopython_seq = Seq(sequence)
    return str(biopython_seq.reverse_complement())
def get_edit_distance(row):
    if row["status_truth_gene"] == "Mapped" and row["status_ref_gene"] == "Mapped":
            truth_gene_seq = row["sequence_truth_gene"]
            ref_gene_seq = row["sequence_ref_gene"]
            truth_gene_seq_rc = reverse_complement(truth_gene_seq)
            edit_distance_fw = editdistance.eval(truth_gene_seq, ref_gene_seq) / max(len(truth_gene_seq), len(ref_gene_seq))
            edit_distance_rc = editdistance.eval(truth_gene_seq_rc, ref_gene_seq) / max(len(truth_gene_seq_rc), len(ref_gene_seq))
            edit_distance = min(edit_distance_fw, edit_distance_rc)
            assert 0 <= edit_distance <= 1
    else:
        edit_distance = -1
    return edit_distance

rule get_edit_distance_between_genes_of_truth_assemblies_and_ref:
    input:
         all_mapped_truth_genes = lambda wildcards: f"{output_folder}/genes_from_truth_or_ref/{{truth_id}}.csv",
         all_mapped_ref_genes = lambda wildcards: f"{output_folder}/genes_from_truth_or_ref/{{ref_id}}.csv",
    output:
         edit_distance_between_genes_of_truth_assemblies_and_ref = f"{output_folder}/edit_distances/{{truth_id}}~~~{{ref_id}}.edit_distance.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/get_edit_distance_between_genes_of_truth_assemblies_and_ref/{{truth_id}}~~~{{ref_id}}.edit_distance.log"
    run:
        all_mapped_truth_genes_df = pd.read_csv(input.all_mapped_truth_genes, index_col="gene_name")
        all_mapped_ref_genes_df = pd.read_csv(input.all_mapped_ref_genes, index_col="gene_name")

        genes_mapping_to_both_truth_and_ref_df = all_mapped_truth_genes_df.join(all_mapped_ref_genes_df,
                                                how="inner", lsuffix="_truth_gene", rsuffix="_ref_gene")

        edit_distances = genes_mapping_to_both_truth_and_ref_df.apply(get_edit_distance, axis=1)
        genes_mapping_to_both_truth_and_ref_df["edit_distance"] = edit_distances
        genes_mapping_to_both_truth_and_ref_df.to_csv(output.edit_distance_between_genes_of_truth_assemblies_and_ref)

rule concatenate_edit_distance_files:
    input:
         edit_distances_files = edit_distances_files
    output:
         all_edit_distance_files_concatenated = f"{output_folder}/edit_distances/all_edit_distances.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/concatenate_edit_distance_files/all_edit_distances.log"
    run:
        dfs = [pd.read_csv(file) for file in input.edit_distances_files]
        concatenated_df = pd.concat(dfs, ignore_index=True)
        concatenated_df.to_csv(output.all_edit_distance_files_concatenated, index=False)
