import pandas as pd
import numpy as np

rule get_variant_precision_score_distance_csv:
    input:
         edit_distance_csv = rules.concatenate_edit_distance_files.output.all_edit_distance_files_concatenated,
         variants_calls = lambda wildcards: f"{pandora_eval_precision_reports_from_probe_mappings_folder}/{wildcards.truth_id}/all/snippy_{wildcards.ref_id}/coverage_filter_Not_App/strand_bias_filter_Not_App/gaps_filter_Not_App/variant_calls_probeset_report.tsv"
    output:
         variant_precision_score_distance_file = f"{output_folder}/get_variant_precision_score_distance_csv/{{truth_id}}~~~{{ref_id}}.get_variant_precision_score_distance.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/get_variant_precision_score_distance_csv/{truth_id}~~~{ref_id}.get_variant_precision_score_distance.log"
    script:
        "../scripts/get_variant_precision_score_distance_csv.py"


rule get_gene_truth_ref_precision_proportion_distance:
    input:
         variant_precision_score_distance_csv = rules.get_variant_precision_score_distance_csv.output.variant_precision_score_distance_file
    output:
         gene_truth_ref_precision_proportion_distance_file = f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/{{truth_id}}~~~{{ref_id}}.gene_truth_ref_precision_proportion_distance.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/get_gene_truth_ref_precision_proportion_distance/{truth_id}~~~{ref_id}.get_gene_truth_ref_precision_proportion_distance.log"
    run:
        variant_precision_score_distance_csv = pd.read_csv(input.variant_precision_score_distance_csv)
        del variant_precision_score_distance_csv['variant']
        gene_truth_ref_precision_proportion_distance = variant_precision_score_distance_csv.groupby(["gene", "truth", "ref", "distance"])\
            .agg({
                "precision_score": [
                    ("max_precision", "count"),
                    ("observed_precision", "sum"),
                    ("precision_ratio", lambda values: values.sum() / values.count())]
        })
        gene_truth_ref_precision_proportion_distance.to_csv(output.gene_truth_ref_precision_proportion_distance_file)



rule concat_gene_truth_ref_precision_proportion_distance_files:
    input:
         gene_truth_ref_precision_proportion_distance_files = gene_truth_ref_precision_proportion_distance_files
    output:
         gene_truth_ref_precision_proportion_distance_concatenated = f"{output_folder}/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/concat_gene_truth_ref_precision_proportion_distance_files/all_gene_truth_ref_precision_proportion_distance.log"
    run:
        dfs = [pd.read_csv(file, header=[0,1,2]) for file in input.gene_truth_ref_precision_proportion_distance_files]
        concatenated_df = pd.concat(dfs)
        concatenated_df.columns = ["gene", "truth", "ref", "distance", "max_precision", "observed_precision", "precision_ratio"]
        concatenated_df.set_index(["gene", "truth", "ref"], inplace=True)
        concatenated_df.to_csv(output.gene_truth_ref_precision_proportion_distance_concatenated)