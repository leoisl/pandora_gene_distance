import seaborn as sns

rule plot_gene_distance:
    input:
         gene_truth_ref_precision_proportion_distance_concatenated = rules.concat_gene_truth_ref_precision_proportion_distance_files.output.gene_truth_ref_precision_proportion_distance_concatenated
    output:
         gene_distance_plot = f"{output_folder}/gene_distance.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/plot_gene_distance/plot_gene_distance.log"
    run:
        df = pd.read_csv(input.gene_truth_ref_precision_proportion_distance_concatenated)
        ax = sns.scatterplot(x=df["distance"], y=df["precision_ratio"], alpha=0.3, color="black", marker=".")
        figure = ax.get_figure()
        figure.savefig(output.gene_distance_plot)