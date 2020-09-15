# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora_gene_distance_indexing_mapping:pandora_paper_tag1

# This container has basic indexing and mapping softwares to be used by the pandora1_paper pipeline
# to build: sudo docker build . -t leandroishilima/pandora_gene_distance_indexing_mapping:pandora_paper_tag1
FROM continuumio/miniconda3:4.8.2
ENV PATH="/opt/conda/bin/:${PATH}"
RUN conda install -c conda-forge -c bioconda samtools==1.10 bwa==0.7.17 bowtie2==2.4.1
