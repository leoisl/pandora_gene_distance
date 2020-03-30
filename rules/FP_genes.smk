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
    shell:
        "touch {output.gene_presence_matrix_based_on_bowtie2}"