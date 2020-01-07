import seaborn as sns
import pandas as pd

def plot(gene_truth_ref_measure_and_distance, precision_or_recall_ration, output_filepath):
    df = pd.read_csv(gene_truth_ref_measure_and_distance)
    ax = sns.scatterplot(x=df["distance"], y=df[precision_or_recall_ration], alpha=0.3, color="black", marker=".")
    figure = ax.get_figure()
    figure.savefig(output_filepath, format="png", dpi=1200)
    figure.clf()



rule plot_gene_distance_plot_for_precision:
    input:
         gene_truth_ref_precision_proportion_distance_concatenated = rules.concat_gene_truth_ref_precision_proportion_distance_files.output.gene_truth_ref_precision_proportion_distance_concatenated,
         gene_truth_ref_precision_proportion_distance_concatenated_filtered = rules.concat_gene_truth_ref_precision_proportion_distance_files.output.gene_truth_ref_precision_proportion_distance_concatenated_filtered
    output:
         gene_distance_plot = f"{output_folder}/gene_distance_precision.png",
         gene_distance_plot_filtered = f"{output_folder}/gene_distance_precision.filtered.png"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/plot_gene_distance_plot_for_precision.log"
    run:
        plot(input.gene_truth_ref_precision_proportion_distance_concatenated, "precision_ratio", output.gene_distance_plot)
        plot(input.gene_truth_ref_precision_proportion_distance_concatenated_filtered, "precision_ratio", output.gene_distance_plot_filtered)



rule plot_gene_distance_plot_for_recall:
    input:
         gene_truth_ref_recall_proportion_distance_concatenated = rules.concat_gene_truth_ref_recall_proportion_distance_files.output.gene_truth_ref_recall_proportion_distance_concatenated,
         gene_truth_ref_recall_proportion_distance_concatenated_filtered = rules.concat_gene_truth_ref_recall_proportion_distance_files.output.gene_truth_ref_recall_proportion_distance_concatenated_filtered
    output:
         gene_distance_plot = f"{output_folder}/gene_distance_recall.png",
         gene_distance_plot_filtered = f"{output_folder}/gene_distance_recall.filtered.png"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/plot_gene_distance_plot_for_recall.log"
    run:
        plot(input.gene_truth_ref_recall_proportion_distance_concatenated, "recall_ratio", output.gene_distance_plot)
        plot(input.gene_truth_ref_recall_proportion_distance_concatenated_filtered, "recall_ratio", output.gene_distance_plot_filtered)