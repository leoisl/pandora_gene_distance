rule compute_gene_presence_matrix_based_on_bowtie2:
    input:
        truth_assemblies = expand(f"{output_folder}/genes_from_truth_or_ref/{{truth_assembly}}.csv", truth_assembly=truth_assemblies.index.tolist())
    output:
        gene_presence_matrix_based_on_bowtie2 = f"{output_folder}/gene_presence_matrix/gene_presence_matrix_based_on_bowtie2",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/compute_gene_presence_matrix_based_on_bowtie2.log"
    run:
        import pandas as pd
        from functools import reduce
        from pathlib import Path
        samples_names = [Path(truth_assembly).with_suffix("").name for truth_assembly in input.truth_assemblies]
        dfs = [pd.read_csv(truth_assembly) for truth_assembly in input.truth_assemblies]
        dfs_with_status_and_gene_name = [df[["status", "gene_name"]] for df in dfs]
        dfs_with_status_and_gene_name = [df.rename(columns={"status": sample_name}) for df, sample_name in zip(dfs_with_status_and_gene_name, samples_names)]
        merged_df = reduce(lambda left, right : left.merge(right, on="gene_name"), dfs_with_status_and_gene_name)
        binarized_merged_df = merged_df[samples_names].applymap(lambda value: 0 if value=="Unmapped" else 1)
        binarized_merged_df["gene_name"] = merged_df["gene_name"]
        binarized_merged_df.set_index(keys=["gene_name"], inplace=True)
        binarized_merged_df.to_csv(output.gene_presence_matrix_based_on_bowtie2, sep="\t")


rule get_gene_lengths:
    input:
         pandora_vcf_ref = pandora_vcf_ref
    output:
         gene_length_matrix = f"{output_folder}/gene_presence_matrix/gene_length_matrix"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * 2**(attempt-1)
    log:
        "logs/get_gene_lengths.log"
    run:
        import pysam
        import pandas as pd

        with pysam.FastxFile(input.pandora_vcf_ref) as pandora_vcf_ref_filehandler:
            gene_names = []
            gene_lengths = []
            for entry in pandora_vcf_ref_filehandler:
                gene_names.append(entry.name)
                gene_lengths.append(len(entry.sequence))

        df = pd.DataFrame(data = {"gene_name": gene_names, "gene_length": gene_lengths})
        df.to_csv(output.gene_length_matrix, sep="\t", index=False)


rule make_FP_genes_plot:
    input:
        pandora_multisample_matrix = pandora_multisample_matrix,
        gene_presence_matrix_based_on_bowtie2 = rules.compute_gene_presence_matrix_based_on_bowtie2.output.gene_presence_matrix_based_on_bowtie2,
        gene_length_matrix = rules.get_gene_lengths.output.gene_length_matrix,
    output:
        gene_and_nb_of_FPs_counted_data = f"{output_folder}/FP_genes/gene_and_nb_of_FPs_counted.csv",
        gene_classification_plot_data = f"{output_folder}/FP_genes/gene_classification.csv",
        gene_classification_plot = f"{output_folder}/FP_genes/gene_classification.png",
        gene_classification_by_sample_plot_data = f"{output_folder}/FP_genes/gene_classification_by_sample.csv",
        gene_classification_by_sample_plot = f"{output_folder}/FP_genes/gene_classification_by_sample.png",
        gene_classification_by_gene_length_plot_data = f"{output_folder}/FP_genes/gene_classification_by_gene_length.csv",
        gene_classification_by_gene_length_plot = f"{output_folder}/FP_genes/gene_classification_by_gene_length.png",
        gene_classification_by_gene_length_normalised_plot_data = f"{output_folder}/FP_genes/gene_classification_by_gene_length_normalised.csv",
        gene_classification_by_gene_length_normalised_plot = f"{output_folder}/FP_genes/gene_classification_by_gene_length_normalised.png",
    params:
        samples = truth_assemblies["id"].to_list(),
        title = pandora_run_that_was_done
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        notebook="logs/make_FP_genes_plot/make_FP_genes_plot.ipynb"
    notebook:
        "../notebooks/FP_genes/FP_genes.ipynb"


rule make_gene_distance_plots:
    input:
        all_edit_distances = rules.concatenate_edit_distance_files.output.all_edit_distance_files_concatenated,
    output:
        distribution_of_genes_per_ed_plot_data =                      f"{output_folder}/gene_distance_plots/distribution_of_genes_per_ed.csv",
        distribution_of_genes_per_ed_counts_plot =                    f"{output_folder}/gene_distance_plots/distribution_of_genes_per_ed_counts.png",
        distribution_of_genes_per_ed_proportion_plot =                f"{output_folder}/gene_distance_plots/distribution_of_genes_per_ed_proportion.png",
        distribution_of_genes_per_nb_of_samples_plot_data =           f"{output_folder}/gene_distance_plots/distribution_of_genes_per_nb_of_samples.csv",
        distribution_of_genes_nb_of_samples_count_plots =      expand(f"{output_folder}/gene_distance_plots/distribution_of_genes_per_nb_of_samples_{{nb_of_sample}}.count.png", nb_of_sample=nb_of_samples),
        distribution_of_genes_nb_of_samples_proportion_plots = expand(f"{output_folder}/gene_distance_plots/distribution_of_genes_per_nb_of_samples_{{nb_of_sample}}.proportion.png", nb_of_sample=nb_of_samples),
        gene_sample_ref_ED_nbsamples_zam =                            f"{output_folder}/gene_distance_plots/gene_sample_ref_ED_nbsamples_zam.csv",
    params:
        nb_of_samples = nb_of_samples,
        output_folder = output_folder
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        notebook="logs/make_gene_distance_plots/make_gene_distance_plots.ipynb"
    notebook:
        "../notebooks/gene_distance_plots/gene_distance_plots.ipynb"
