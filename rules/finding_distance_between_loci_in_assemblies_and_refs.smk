import pysam

rule split_genes_from_vcf_ref:
    input:
        pandora_vcf_ref = pandora_vcf_ref
    output:
        all_genes_filepaths_in_vcf_ref = all_genes_filepaths_in_vcf_ref
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/split_genes_from_vcf_ref.log"
    run:
        with pysam.FastxFile(input.pandora_vcf_ref) as fasta_file:
            for fasta_record in fasta_file:
                with open(f"{output_folder}/genes_from_vcf_ref/{fasta_record.name}.fa", "w") as fasta_record_file:
                    fasta_record_file.write(str(fasta_record))


rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.log"
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    shell:
        "bwa index {input.fasta} > {log} 2>&1"


rule map_gene_from_vcf_ref_to_truth_or_ref:
    input:
        gene = lambda wildcards: f"{output_folder}/genes_from_vcf_ref/{wildcards.gene}.fa",
        truth_or_ref = lambda wildcards: f"{all_assemblies.xs(wildcards.id)['fasta']}",
        truth_or_ref_index = lambda wildcards: f"{all_assemblies.xs(wildcards.id)['fasta']}.amb"
    output:
        sam_file = f"{output_folder}/map_gene_from_vcf_ref_to_truth_or_ref/{{gene}}~~~{{id}}.sam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/map_gene_from_vcf_ref_to_truth_or_ref/{{gene}}~~~{{id}}.log"
    shell:
        "bwa mem {input.truth_or_ref} {input.gene} -t {threads} -o {output.sam_file}"