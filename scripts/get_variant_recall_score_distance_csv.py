import pandas as pd
from GetVariantDistanceHelpers import GetVariantDistanceHelpersForRecall, ProbeDoesNotMapToAnyGene, ProbeMapsToSeveralGenes
import logging
log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)
import copy


def get_recall_score(classification):
    if classification == "primary_correct" or classification == "secondary_correct" or classification == "supplementary_correct":
        return 1.0
    elif classification == "unmapped" or classification == "partially_mapped" or classification == "primary_incorrect" \
            or classification == "secondary_incorrect" or classification == "supplementary_incorrect":
        return 0.0
    else:
        assert False, f"Fatal error in get_recall_score(): unknown classification: {classification}"

# functions
def get_variant_output_dict(edit_distance_df, variant_calls_df, truth_id, ref_id):
    output_dict = {
        "gene": [],
        "truth": [],
        "ref": [],
        "variant": [],
        "recall_score": [],
        "distance": []
    }

    output_dict_probes_do_no_map_to_any_gene = copy.deepcopy(output_dict)

    for index, row in variant_calls_df.iterrows():
        truth_probe_header = row["query_probe_header"]
        classification = row["classification"]
        recall_score = get_recall_score(classification)
        contig_truth_probe_originates_from = GetVariantDistanceHelpersForRecall.parse_field_from_header("CHROM",
                                                                                                         truth_probe_header)
        pos_truth_probe_originates_from = GetVariantDistanceHelpersForRecall.parse_field_from_header("POS",
                                                                                                      truth_probe_header,
                                                                                                      return_type=int)

        try:
            gene_name, edit_distance = GetVariantDistanceHelpersForRecall.get_gene_name_and_edit_distance_of_gene_this_truth_probe_maps_to(
                edit_distance_df, contig_truth_probe_originates_from, pos_truth_probe_originates_from)
            output_dict["gene"].append(gene_name)
            output_dict["truth"].append(truth_id)
            output_dict["ref"].append(ref_id)
            output_dict["variant"].append(truth_probe_header)
            output_dict["recall_score"].append(recall_score)
            output_dict["distance"].append(edit_distance)
        except ProbeDoesNotMapToAnyGene:
            output_dict_probes_do_no_map_to_any_gene["gene"].append(truth_probe_header)
            output_dict_probes_do_no_map_to_any_gene["truth"].append(truth_id)
            output_dict_probes_do_no_map_to_any_gene["ref"].append(ref_id)
            output_dict_probes_do_no_map_to_any_gene["variant"].append(truth_probe_header)
            output_dict_probes_do_no_map_to_any_gene["recall_score"].append(recall_score)
            output_dict_probes_do_no_map_to_any_gene["distance"].append(-1.0)
        except ProbeMapsToSeveralGenes:
            assert False, f"Probe {truth_probe_header} maps to several genes, which does not make sense (or it does: https://en.wikipedia.org/wiki/Overlapping_gene#Evolution)..."

    return output_dict, output_dict_probes_do_no_map_to_any_gene

# snakemake config
edit_distance_csv = snakemake.input.edit_distance_csv
truth_id = snakemake.wildcards.truth_id
ref_id = snakemake.wildcards.ref_id
variant_calls = snakemake.input.variants_calls
variant_recall_score_distance_file = snakemake.output.variant_recall_score_distance_file
variant_recall_score_distance_unmapped_probes_file = snakemake.output.variant_recall_score_distance_unmapped_probes_file

# API usage
edit_distance_df = pd.read_csv(edit_distance_csv)
edit_distance_df = GetVariantDistanceHelpersForRecall.get_edit_distance_df_given_truth_and_ref(edit_distance_df, truth_id, ref_id)
variant_calls_df = pd.read_csv(variant_calls, sep="\t")

output_dict, output_dict_probes_do_no_map_to_any_gene = get_variant_output_dict(edit_distance_df, variant_calls_df, truth_id, ref_id)

# output
output_df = pd.DataFrame(data=output_dict)
output_df.to_csv(variant_recall_score_distance_file)
output_df_probes_do_no_map_to_any_gene = pd.DataFrame(data=output_dict_probes_do_no_map_to_any_gene)
output_df_probes_do_no_map_to_any_gene.to_csv(variant_recall_score_distance_unmapped_probes_file)