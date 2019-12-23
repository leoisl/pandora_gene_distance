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
