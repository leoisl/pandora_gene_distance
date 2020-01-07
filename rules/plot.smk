import seaborn as sns
import pandas as pd

rule plot_gene_distance_plot_for_precision:
    input:
         gene_truth_ref_precision_proportion_distance_concatenated = rules.concat_gene_truth_ref_precision_proportion_distance_files.output.gene_truth_ref_precision_proportion_distance_concatenated,
         gene_truth_ref_precision_proportion_distance_concatenated_filtered = rules.concat_gene_truth_ref_precision_proportion_distance_files.output.gene_truth_ref_precision_proportion_distance_concatenated_filtered
    output:
         gene_distance_plot = f"{output_folder}/gene_distance_precision.pdf",
         gene_distance_plot_filtered = f"{output_folder}/gene_distance_precision.filtered.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/plot_gene_distance_plot_for_precision.log"
    run:
        def plot(gene_truth_ref_precision_proportion_distance, output_filepath):
            df = pd.read_csv(gene_truth_ref_precision_proportion_distance)
            ax = sns.scatterplot(x=df["distance"], y=df["precision_ratio"], alpha=0.3, color="black", marker=".")
            figure = ax.get_figure()
            figure.savefig(output_filepath)
            figure.clf()

        plot(input.gene_truth_ref_precision_proportion_distance_concatenated, output.gene_distance_plot)
        plot(input.gene_truth_ref_precision_proportion_distance_concatenated_filtered, output.gene_distance_plot_filtered)



rule plot_gene_distance_plot_for_recall:
    input:
         gene_truth_ref_recall_proportion_distance_concatenated = rules.concat_gene_truth_ref_recall_proportion_distance_files.output.gene_truth_ref_recall_proportion_distance_concatenated,
         gene_truth_ref_recall_proportion_distance_concatenated_filtered = rules.concat_gene_truth_ref_recall_proportion_distance_files.output.gene_truth_ref_recall_proportion_distance_concatenated_filtered
    output:
         gene_distance_plot = f"{output_folder}/gene_distance_recall.pdf",
         gene_distance_plot_filtered = f"{output_folder}/gene_distance_recall.filtered.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/plot_gene_distance_plot_for_recall.log"
    run:
        def plot(gene_truth_ref_recall_proportion_distance, output_filepath):
            df = pd.read_csv(gene_truth_ref_recall_proportion_distance)
            ax = sns.scatterplot(x=df["distance"], y=df["recall_ratio"], alpha=0.3, color="black", marker=".")
            figure = ax.get_figure()
            figure.savefig(output_filepath)
            figure.clf()

        plot(input.gene_truth_ref_recall_proportion_distance_concatenated, output.gene_distance_plot)
        plot(input.gene_truth_ref_recall_proportion_distance_concatenated_filtered, output.gene_distance_plot_filtered)