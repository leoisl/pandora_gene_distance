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
        samples_names = [truth_assembly[:-4] for truth_assembly in input.truth_assemblies]
        dfs = [pd.read_csv(truth_assembly) for truth_assembly in input.truth_assemblies]
        dfs_with_status_and_gene_name = [df[["status", "gene_name"]] for df in dfs]
        dfs_with_status_and_gene_name = [df.rename(columns={"status": sample_name}) for df, sample_name in zip(dfs_with_status_and_gene_name, samples_names)]
        merged_df = reduce(lambda left, right : left.merge(right, on="gene_name"), dfs_with_status_and_gene_name)
        binarized_merged_df = merged_df[samples_names].applymap(lambda value: 0 if value=="Unmapped" else 1)
        binarized_merged_df["gene_name"] = merged_df["gene_name"]
        binarized_merged_df.set_index(keys=["gene_name"], inplace=True)
        binarized_merged_df.to_csv(output.gene_presence_matrix_based_on_bowtie2, sep="\t")