import pysam
import editdistance
import pandas as pd


rule get_truth_or_ref_gene_sequences:
    input:
        truth_or_ref = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}",
        truth_or_ref_index = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}.fai",
        pandora_vcf_ref_mapped_to_truth_or_ref_sam_file = rules.map_pandora_vcf_ref_to_truth_or_ref_using_bowtie.output.sam_file
    output:
         truth_or_ref_gene_sequences = f"{output_folder}/genes_from_truth_or_ref/{{id}}.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/get_truth_or_ref_gene_sequences/{id}.log"
    run:
        gene_seq_infos = {
            "status": [],
            "gene_name": [],
            "ref_or_truth_id": [],
            "contig": [],
            "start": [],
            "stop": [],
            "sequence": []
        }

        mapped_genes = set()
        with pysam.AlignmentFile(input.pandora_vcf_ref_mapped_to_truth_or_ref_sam_file) as pandora_vcf_ref_mapped_to_truth_or_ref_sam_file, \
             pysam.FastaFile(filename = input.truth_or_ref, filepath_index=input.truth_or_ref_index) as truth_or_ref:
            ref_or_truth_id = wildcards.id
            for sam_record in pandora_vcf_ref_mapped_to_truth_or_ref_sam_file:
                if not sam_record.is_unmapped:
                    gene_seq_infos["status"].append("Mapped")

                    gene_name = sam_record.query_name
                    assert gene_name not in mapped_genes, f"Gene {gene_name} has 2+ mappings in {input.pandora_vcf_ref_mapped_to_truth_or_ref_sam_file}, this should not have happened"
                    mapped_genes.add(gene_name)
                    gene_seq_infos["gene_name"].append(gene_name)

                    gene_seq_infos["ref_or_truth_id"].append(ref_or_truth_id)

                    contig = sam_record.reference_name
                    gene_seq_infos["contig"].append(contig)

                    start = int(sam_record.reference_start)
                    gene_seq_infos["start"].append(start)

                    stop = int(sam_record.reference_end)
                    gene_seq_infos["stop"].append(stop)

                    truth_or_ref_contig_sequence = truth_or_ref.fetch(reference = contig)
                    sequence = truth_or_ref_contig_sequence[start:stop].upper()
                    gene_seq_infos["sequence"].append(sequence)
                else:
                    gene_seq_infos["status"].append("Unmapped")
                    gene_seq_infos["gene_name"].append(sam_record.query_name)
                    gene_seq_infos["ref_or_truth_id"].append(ref_or_truth_id)
                    gene_seq_infos["contig"].append("")
                    gene_seq_infos["start"].append(-1)
                    gene_seq_infos["stop"].append(-1)
                    gene_seq_infos["sequence"].append("")

        df = pd.DataFrame(gene_seq_infos)
        df.to_csv(output.truth_or_ref_gene_sequences, index=False)



def get_edit_distance(row):
    if row["status_truth_gene"] == "Mapped" and row["status_ref_gene"] == "Mapped":
            truth_gene_seq = row["sequence_truth_gene"]
            ref_gene_seq = row["sequence_ref_gene"]
            edit_distance = editdistance.eval(truth_gene_seq, ref_gene_seq) / max(len(truth_gene_seq), len(ref_gene_seq))
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
        mem_mb = lambda wildcards, attempt: 2000 * attempt
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
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/concatenate_edit_distance_files/all_edit_distances.log"
    run:
        dfs = [pd.read_csv(file) for file in input.edit_distances_files]
        concatenated_df = pd.concat(dfs, ignore_index=True)
        concatenated_df.to_csv(output.all_edit_distance_files_concatenated, index=False)
