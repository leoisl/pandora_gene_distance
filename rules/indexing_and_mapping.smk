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
        mem_mb = lambda wildcards, attempt: 4000 * 2**(attempt-1)
    singularity:
        "docker://leandroishilima/pandora_gene_distance_indexing_mapping:pandora_paper_tag1"
    shell:
        "bowtie2-build {input.fasta} {input.fasta}.bowtie_index > {log} 2>&1"


rule map_pandora_vcf_ref_to_truth_or_ref_using_bowtie:
    input:
        pandora_vcf_ref = lambda wildcards: f"{pandora_vcfs.xs(wildcards.pandora_id)['fasta']}",
        truth_or_ref = lambda wildcards: f"{samples_and_refs_and_pandora_vcfs.xs(wildcards.sample_ref_id)['fasta']}",
        truth_or_ref_index = lambda wildcards: f"{samples_and_refs_and_pandora_vcfs.xs(wildcards.sample_ref_id)['fasta']}.bowtie_index.1.bt2"
    output:
        sam_file = f"{output_folder}/map_pandora_vcf_ref_to_truth_or_ref/{{pandora_id}}/{{sample_ref_id}}.bowtie.sam"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * 2**(attempt-1)
    log:
        "logs/map_pandora_vcf_ref_to_truth_or_ref/{pandora_id}/{sample_ref_id}.bowtie.log"
    singularity:
        "docker://leandroishilima/pandora_gene_distance_indexing_mapping:pandora_paper_tag1"
    shell:
        "bowtie2 --very-sensitive --end-to-end -f -x {input.truth_or_ref}.bowtie_index -U {input.pandora_vcf_ref} -S {output.sam_file} > {log} 2>&1"


rule index_fasta_file:
    input:
        fasta = "{fasta}"
    output:
        fasta_index = "{fasta}.fai"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * 2**(attempt-1)
    log: "{fasta}.index_fasta_file.log"
    singularity:
        "docker://leandroishilima/pandora_gene_distance_indexing_mapping:pandora_paper_tag1"
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"
