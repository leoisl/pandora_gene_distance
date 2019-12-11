rule bwa_index:
    input:
        fasta = "{fasta}"
    output:
        indexed_fasta = "{fasta}.amb"
    threads: 1
    log: "{fasta}.bwa_index.log"
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    shell:
        "bwa index {input.fasta} > {log} 2>&1"


rule bowtie2_build:
    input:
        fasta = "{fasta}"
    output:
        "{fasta}.bowtie_index.1.bt2",
        "{fasta}.bowtie_index.2.bt2",
        "{fasta}.bowtie_index.3.bt2",
        "{fasta}.bowtie_index.4.bt2",
        "{fasta}.bowtie_index.rev.1.bt2",
        "{fasta}.bowtie_index.rev.2.bt2"
    threads: 1
    log: "{fasta}.bowtie2_build.log"
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    shell:
        "bowtie2-build {input.fasta} {input.fasta}.bowtie_index"


rule map_gene_from_vcf_ref_to_truth_or_ref_using_bowtie:
    input:
        gene = lambda wildcards: f"{output_folder}/genes_from_vcf_ref/{wildcards.gene}.fa",
        truth_or_ref = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}",
        truth_or_ref_index = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}.bowtie_index.1.bt2"
    output:
        sam_file = f"{output_folder}/map_gene_from_vcf_ref_to_truth_or_ref/{{gene}}~~~{{id}}.bowtie.sam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/map_gene_from_vcf_ref_to_truth_or_ref/{gene}~~~{id}.log"
    shell:
        "bowtie2 --very-sensitive --end-to-end -f -x {input.truth_or_ref}.bowtie_index -U {input.gene} -S {output.sam_file}"


rule map_gene_from_vcf_ref_to_truth_or_ref_using_bwamem:
    input:
        gene = lambda wildcards: f"{output_folder}/genes_from_vcf_ref/{wildcards.gene}.fa",
        truth_or_ref = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}",
        truth_or_ref_index = lambda wildcards: f"{assemblies_and_refs.xs(wildcards.id)['fasta']}.amb"
    output:
        sam_file = f"{output_folder}/map_gene_from_vcf_ref_to_truth_or_ref/{{gene}}~~~{{id}}.bwa.sam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/map_gene_from_vcf_ref_to_truth_or_ref/{gene}~~~{id}.log"
    shell:
        "bwa mem {input.truth_or_ref} {input.gene} -t {threads} -o {output.sam_file}"


rule index_fasta_file:
    input:
        fasta = "{fasta}"
    output:
        fasta_index = "{fasta}.fai"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    log: "{fasta}.index_fasta_file.log"
    shell:
        "samtools faidx {input.fasta}"
