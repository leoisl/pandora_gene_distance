from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))

import re

class ProbeDoesNotMapToAnyGene(Exception):
    pass

class GetVariantDistanceHelpers:
    @staticmethod
    def get_edit_distance_df_given_truth_and_ref(edit_distance_df, truth_id, ref_id):
        # restricts the edit_distance_df to this truth and ref
        edit_distance_df = edit_distance_df.query("ref_or_truth_id_truth_gene == @truth_id and "
                                                  "ref_or_truth_id_ref_gene == @ref_id")

        # restricts only to mapped pair of genes between truth and ref
        edit_distance_df = edit_distance_df.query("status_truth_gene == 'Mapped' and status_ref_gene == 'Mapped'")

        return edit_distance_df

    @staticmethod
    def get_gene_names_and_edit_distances_of_genes_this_probe_maps_to(edit_distance_df, contig_column, start_gene_column,
                                                                      stop_gene_column, contig_probe_originates_from,
                                                                      pos_probe_originates_from):
        # get all genes the vcf probe maps to
        df_where_probe_maps = edit_distance_df.query(f"{contig_column} == @contig_probe_originates_from and "
                                                     f"{start_gene_column} <= @pos_probe_originates_from and "
                                                     f"@pos_probe_originates_from < {stop_gene_column}")

        if len(df_where_probe_maps) == 0:
            raise ProbeDoesNotMapToAnyGene()
        else:
            gene_names = df_where_probe_maps["gene_name"].to_list()
            edit_distances = df_where_probe_maps["edit_distance"].to_list()
            return gene_names, edit_distances


    # NOTE: copied from https://github.com/iqbal-lab/pandora1_paper
    # NOTE: already tested, no need to unit test
    @staticmethod
    def parse_field_from_header(
        field: str, header: str, return_type = str, delim = ";"
    ):
        regex = re.compile(f"{field}=(.+?){delim}")
        match = regex.search(header)
        if match:
            return return_type(match.group(1))
        else:
            return return_type()


class GetVariantDistanceHelpersForPrecision(GetVariantDistanceHelpers):
    @staticmethod
    def get_gene_names_and_edit_distances_of_genes_this_vcf_probe_maps_to(edit_distance_df, contig_vcf_probe_originates_from,
                                                                          pos_vcf_probe_originates_from):
        return GetVariantDistanceHelpers.get_gene_names_and_edit_distances_of_genes_this_probe_maps_to\
            (edit_distance_df=edit_distance_df, contig_column="contig_ref_gene", start_gene_column="start_ref_gene",
             stop_gene_column="stop_ref_gene", contig_probe_originates_from=contig_vcf_probe_originates_from,
             pos_probe_originates_from=pos_vcf_probe_originates_from)


class GetVariantDistanceHelpersForRecall(GetVariantDistanceHelpers):
    @staticmethod
    def get_gene_names_and_edit_distances_of_genes_this_truth_probe_maps_to(edit_distance_df, contig_truth_probe_originates_from,
                                                                            pos_truth_probe_originates_from):
        return GetVariantDistanceHelpers.get_gene_names_and_edit_distances_of_genes_this_probe_maps_to\
            (edit_distance_df=edit_distance_df, contig_column="contig_truth_gene", start_gene_column="start_truth_gene",
             stop_gene_column="stop_truth_gene", contig_probe_originates_from=contig_truth_probe_originates_from,
             pos_probe_originates_from=pos_truth_probe_originates_from)