import pysam
import pandas as pd

# snakemake vars dealing
input_pandora_vcf_ref_mapped_to_truth_or_ref_sam_file = snakemake.input.pandora_vcf_ref_mapped_to_truth_or_ref_sam_file
input_truth_or_ref = snakemake.input.truth_or_ref
input_truth_or_ref_index = snakemake.input.truth_or_ref_index
ref_or_truth_id = snakemake.wildcards.sample_ref_id
output_truth_or_ref_gene_sequences = snakemake.output.truth_or_ref_gene_sequences



def get_truth_or_ref_gene_sequences(pandora_vcf_ref_mapped_to_truth_or_ref_sam_file, truth_or_ref):
    # init
    gene_seq_infos = {
        "status": [],
        "gene_name": [],
        "ref_or_truth_id": [],
        "contig": [],
        "start": [],
        "stop": [],
        "sequence": [],
        "strand": []
    }

    # processing
    mapped_genes = set()
    for sam_record in pandora_vcf_ref_mapped_to_truth_or_ref_sam_file:
        # check for multimapping genes
        gene_name = sam_record.query_name
        assert gene_name not in mapped_genes, f"Gene {gene_name} has 2+ mappings in {input_pandora_vcf_ref_mapped_to_truth_or_ref_sam_file}, this should not have happened"
        mapped_genes.add(gene_name)

        # add common stuff
        gene_seq_infos["gene_name"].append(gene_name)
        gene_seq_infos["ref_or_truth_id"].append(ref_or_truth_id)

        if not sam_record.is_unmapped:
            gene_seq_infos["status"].append("Mapped")

            contig = sam_record.reference_name
            gene_seq_infos["contig"].append(contig)

            start = int(sam_record.reference_start)
            gene_seq_infos["start"].append(start)

            stop = int(sam_record.reference_end)
            gene_seq_infos["stop"].append(stop)

            is_reverse_complemented = sam_record.flag & 0x10
            if is_reverse_complemented:
                strand = "-"
            else:
                strand = "+"
            gene_seq_infos["strand"].append(strand)

            truth_or_ref_contig_sequence = truth_or_ref.fetch(reference=contig)
            sequence = truth_or_ref_contig_sequence[start:stop].upper()
            gene_seq_infos["sequence"].append(sequence)
        else:
            gene_seq_infos["status"].append("Unmapped")
            gene_seq_infos["contig"].append("")
            gene_seq_infos["start"].append(-1)
            gene_seq_infos["stop"].append(-1)
            gene_seq_infos["sequence"].append("")
            gene_seq_infos["strand"].append("")

    df = pd.DataFrame(gene_seq_infos)
    df.to_csv(output_truth_or_ref_gene_sequences, index=False)


def main():
    with pysam.AlignmentFile(input_pandora_vcf_ref_mapped_to_truth_or_ref_sam_file) as pandora_vcf_ref_mapped_to_truth_or_ref_sam_file, \
         pysam.FastaFile(filename=input_truth_or_ref, filepath_index=input_truth_or_ref_index) as truth_or_ref:
        get_truth_or_ref_gene_sequences(pandora_vcf_ref_mapped_to_truth_or_ref_sam_file, truth_or_ref)



main()
