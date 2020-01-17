import pandas as pd
from scripts.plot import *


rule plot_gene_distance_plot_for_precision:
    input:
         gene_truth_ref_precision_proportion_distance_concatenated = rules.concat_gene_truth_ref_precision_proportion_distance_files.output.gene_truth_ref_precision_proportion_distance_concatenated
    output:
         gene_distance_lineplot_0001 = f"{output_folder}/gene_distance_precision.lineplot.0.001_bins.png",
         gene_distance_lineplot_001 = f"{output_folder}/gene_distance_precision.lineplot.0.01_bins.png",
         gene_distance_violinplot_001_area = f"{output_folder}/gene_distance_precision.violinplot.0.01_bins.area.png",
         gene_distance_violinplot_001_count = f"{output_folder}/gene_distance_precision.violinplot.0.01_bins.count.png",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/plot_gene_distance_plot_for_precision.log"
    run:
        df = pd.read_csv(input.gene_truth_ref_precision_proportion_distance_concatenated)
        df_with_step_0001 = get_df_with_edit_distance_labels_for_precision(df, step = 0.001, max_bin = 0.2)
        df_with_step_001 = get_df_with_edit_distance_labels_for_precision(df, step = 0.01, max_bin = 0.2)

        plot_lineplot_mean_precision_and_count_in_genes_in_several_bins(df_with_step_0001, step=0.001,
                                                                     output_filepath=output.gene_distance_lineplot_0001)
        plot_lineplot_mean_precision_and_count_in_genes_in_several_bins(df_with_step_001, step=0.01,
                                                                     output_filepath=output.gene_distance_lineplot_001)

        plot_violinplot_precision_in_genes_in_several_bins(df_with_step_001, step=0.01, scale="area",
                                                        output_filepath=output.gene_distance_violinplot_001_area)
        plot_violinplot_precision_in_genes_in_several_bins(df_with_step_001, step=0.01, scale="count",
                                                        output_filepath=output.gene_distance_violinplot_001_count)


rule plot_gene_distance_plot_for_recall:
    input:
         gene_truth_ref_recall_proportion_distance_concatenated = rules.concat_gene_truth_ref_recall_proportion_distance_files.output.gene_truth_ref_recall_proportion_distance_concatenated
    output:
         gene_distance_lineplot_0001 = f"{output_folder}/gene_distance_recall.lineplot.0.001_bins.png",
         gene_distance_lineplot_001 = f"{output_folder}/gene_distance_recall.lineplot.0.01_bins.png",
         gene_distance_violinplot_001_area = f"{output_folder}/gene_distance_recall.violinplot.0.01_bins.area.png",
         gene_distance_violinplot_001_count = f"{output_folder}/gene_distance_recall.violinplot.0.01_bins.count.png",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/plot_gene_distance_plot_for_recall.log"
    run:
        df = pd.read_csv(input.gene_truth_ref_recall_proportion_distance_concatenated)
        df_with_step_0001 = get_df_with_edit_distance_labels_for_recall(df, step = 0.001, max_bin = 0.2)
        df_with_step_001 = get_df_with_edit_distance_labels_for_recall(df, step = 0.01, max_bin = 0.2)

        plot_lineplot_mean_recall_and_count_in_genes_in_several_bins(df_with_step_0001, step=0.001,
                                                                     output_filepath=output.gene_distance_lineplot_0001)
        plot_lineplot_mean_recall_and_count_in_genes_in_several_bins(df_with_step_001, step=0.01,
                                                                     output_filepath=output.gene_distance_lineplot_001)

        plot_violinplot_recall_in_genes_in_several_bins(df_with_step_001, step=0.01, scale="area",
                                                        output_filepath=output.gene_distance_violinplot_001_area)
        plot_violinplot_recall_in_genes_in_several_bins(df_with_step_001, step=0.01, scale="count",
                                                        output_filepath=output.gene_distance_violinplot_001_count)