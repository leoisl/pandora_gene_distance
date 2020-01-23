import pandas as pd
from scripts.dtypes import (
    variant_precision_score_distance_csv_dtype_dict,
    variant_recall_score_distance_csv_dtype_dict,
    gene_truth_ref_precision_dtype_dict,
    gene_truth_ref_recall_dtype_dict)

rule generate_variant_calls_for_precision_folder_if_it_does_not_exist:
    input:
        original_variant_calls = lambda wildcards: f"{precision_reports}/{wildcards.truth_id}/{coverage}/{tool}/coverage_filter_{coverage_filter}/strand_bias_filter_{strand_bias_filter}/gaps_filter_{gaps_filter}/variant_calls_probeset_report.tsv"
    output:
        variants_calls =  f"{precision_reports}/{{truth_id}}/{coverage}/{tool}_{{ref_id}}/coverage_filter_{coverage_filter}/strand_bias_filter_{strand_bias_filter}/gaps_filter_{gaps_filter}/variant_calls_probeset_report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * 2**(attempt-1)
    log:
        "logs/generate_variant_calls_for_precision_folder_if_it_does_not_exist/{truth_id}~~~{ref_id}.generate_variant_calls_for_precision_folder_if_it_does_not_exist.log"
    shell:
        "ln -s `realpath {input.original_variant_calls}` {output.variants_calls}"


rule get_variant_precision_score_distance_csv:
    input:
         edit_distance_csv = rules.concatenate_edit_distance_files.output.all_edit_distance_files_concatenated,
         variants_calls = rules.generate_variant_calls_for_precision_folder_if_it_does_not_exist.output.variants_calls,
    output:
         variant_precision_score_distance_file = f"{output_folder}/get_variant_precision_score_distance_csv/{{truth_id}}~~~{{ref_id}}.get_variant_precision_score_distance.csv",
         variant_precision_score_distance_unmapped_probes_file = f"{output_folder}/get_variant_precision_score_distance_csv/{{truth_id}}~~~{{ref_id}}.get_variant_precision_score_distance.unmapped.csv",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/get_variant_precision_score_distance_csv/{truth_id}~~~{ref_id}.get_variant_precision_score_distance.log"
    script:
        "../scripts/get_variant_precision_score_distance_csv.py"


rule generate_truth_calls_folder_if_it_does_not_exist:
    input:
        original_variant_calls = lambda wildcards: f"{recall_reports}/{wildcards.truth_id}/{coverage}/{tool}/coverage_filter_{coverage_filter}/strand_bias_filter_{strand_bias_filter}/gaps_filter_{gaps_filter}/{wildcards.sample_pair}.report.tsv",
    output:
        variants_calls =  f"{recall_reports}/{{truth_id}}/{coverage}/{tool}_{{ref_id}}/coverage_filter_{coverage_filter}/strand_bias_filter_{strand_bias_filter}/gaps_filter_{gaps_filter}/{{sample_pair}}.report.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * 2**(attempt-1)
    log:
        "logs/generate_truth_calls_folder_if_it_does_not_exist/{truth_id}~~~{ref_id}~~~{sample_pair}.generate_truth_calls_folder_if_it_does_not_exist.log"
    shell:
        "ln -s `realpath {input.original_variant_calls}` {output.variants_calls}"



rule get_variant_recall_score_distance_csv:
    input:
         edit_distance_csv = rules.concatenate_edit_distance_files.output.all_edit_distance_files_concatenated,
         variants_calls = rules.generate_truth_calls_folder_if_it_does_not_exist.output.variants_calls
    output:
         variant_recall_score_distance_file = f"{output_folder}/get_variant_recall_score_distance_csv/{{truth_id}}~~~{{ref_id}}/{{sample_pair}}.get_variant_recall_score_distance.csv",
         variant_recall_score_distance_unmapped_probes_file = f"{output_folder}/get_variant_recall_score_distance_csv/{{truth_id}}~~~{{ref_id}}/{{sample_pair}}.get_variant_recall_score_distance.unmapped.csv",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/get_variant_recall_score_distance_csv/{truth_id}~~~{ref_id}~~~{sample_pair}.get_variant_recall_score_distance.log"
    script:
        "../scripts/get_variant_recall_score_distance_csv.py"


rule get_gene_truth_ref_precision_proportion_distance:
    input:
         variant_precision_score_distance_csv = rules.get_variant_precision_score_distance_csv.output.variant_precision_score_distance_file
    output:
         gene_truth_ref_precision_proportion_distance_file = f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/{{truth_id}}~~~{{ref_id}}.gene_truth_ref_precision_proportion_distance.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/get_gene_truth_ref_precision_proportion_distance/{truth_id}~~~{ref_id}.get_gene_truth_ref_precision_proportion_distance.log"
    run:
        variant_precision_score_distance_csv = pd.read_csv(input.variant_precision_score_distance_csv,
                                                           dtype=variant_precision_score_distance_csv_dtype_dict)
        del variant_precision_score_distance_csv['variant']
        gene_truth_ref_precision_proportion_distance = variant_precision_score_distance_csv.groupby(["gene", "truth", "ref", "distance"])\
            .agg({
                "precision_score": [
                    ("max_precision", "count"),
                    ("observed_precision", "sum"),
                    ("precision_ratio", lambda values: values.sum() / values.count())]
        })
        gene_truth_ref_precision_proportion_distance.reset_index(inplace=True)
        gene_truth_ref_precision_proportion_distance.columns = ["gene", "truth", "ref", "distance", "max_precision", "observed_precision", "precision_ratio"]
        gene_truth_ref_precision_proportion_distance.to_csv(output.gene_truth_ref_precision_proportion_distance_file, index=False)



rule get_gene_truth_ref_recall_proportion_distance:
    input:
         variant_recall_score_distance_csv = rules.get_variant_recall_score_distance_csv.output.variant_recall_score_distance_file
    output:
         gene_truth_ref_recall_proportion_distance_file = f"{output_folder}/get_gene_truth_ref_recall_proportion_distance/{{truth_id}}~~~{{ref_id}}/{{sample_pair}}.gene_truth_ref_recall_proportion_distance.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/get_gene_truth_ref_recall_proportion_distance/{truth_id}~~~{ref_id}~~~{sample_pair}.get_gene_truth_ref_recall_proportion_distance.log"
    run:
        variant_recall_score_distance_csv = pd.read_csv(input.variant_recall_score_distance_csv,
                                                        dtype=variant_recall_score_distance_csv_dtype_dict)
        del variant_recall_score_distance_csv['variant']
        gene_truth_ref_recall_proportion_distance = variant_recall_score_distance_csv.groupby(["gene", "truth", "ref", "distance"])\
            .agg({
                "recall_score": [
                    ("max_recall", "count"),
                    ("observed_recall", "sum"),
                    ("recall_ratio", lambda values: values.sum() / values.count())]
        })
        gene_truth_ref_recall_proportion_distance.reset_index(inplace=True)
        gene_truth_ref_recall_proportion_distance.columns = ["gene", "truth", "ref", "distance", "max_recall", "observed_recall", "recall_ratio"]
        gene_truth_ref_recall_proportion_distance.to_csv(output.gene_truth_ref_recall_proportion_distance_file, index=False)


rule concat_gene_truth_ref_precision_proportion_distance_files:
    input:
         gene_truth_ref_precision_proportion_distance_files = gene_truth_ref_precision_proportion_distance_files
    output:
         gene_truth_ref_precision_proportion_distance_concatenated = f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/concat_gene_truth_ref_precision_proportion_distance_files/all_gene_truth_ref_precision_proportion_distance.log"
    run:
        dfs = [pd.read_csv(file, dtype=gene_truth_ref_precision_dtype_dict) for file in input.gene_truth_ref_precision_proportion_distance_files]
        concatenated_df = pd.concat(dfs)
        concatenated_df.set_index(["gene", "truth", "ref"], inplace=True)
        concatenated_df.to_csv(output.gene_truth_ref_precision_proportion_distance_concatenated)



rule concat_gene_truth_ref_recall_proportion_distance_files:
    input:
         gene_truth_ref_recall_proportion_distance_files = gene_truth_ref_recall_proportion_distance_files
    output:
         gene_truth_ref_recall_proportion_distance_concatenated = f"{output_folder}/get_gene_truth_ref_recall_proportion_distance/all_gene_truth_ref_recall_proportion_distance.csv",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/concat_gene_truth_ref_recall_proportion_distance_files/all_gene_truth_ref_recall_proportion_distance.log"
    run:
        dfs = [pd.read_csv(file, dtype=gene_truth_ref_recall_dtype_dict) for file in input.gene_truth_ref_recall_proportion_distance_files]
        concatenated_df = pd.concat(dfs)
        concatenated_df.set_index(["gene", "truth", "ref"], inplace=True)
        concatenated_df.to_csv(output.gene_truth_ref_recall_proportion_distance_concatenated)
