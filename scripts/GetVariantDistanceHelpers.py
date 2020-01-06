from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))

import re
from scripts.utils import get_first_value_of_series

class ProbeDoesNotMapToAnyGene(Exception):
    pass
class ProbeMapsToSeveralGenes(Exception):
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
    def get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(edit_distance_df, contig_column, start_gene_column,
                                                                       stop_gene_column, contig_probe_originates_from,
                                                                       pos_probe_originates_from):
        # get all genes the vcf probe maps to
        df_where_probe_maps = edit_distance_df.query(f"{contig_column} == @contig_probe_originates_from and "
                                                     f"{start_gene_column} <= @pos_probe_originates_from and "
                                                     f"@pos_probe_originates_from < {stop_gene_column}")

        if len(df_where_probe_maps) == 0:
            raise ProbeDoesNotMapToAnyGene()
        elif len(df_where_probe_maps) == 1:
            truth_gene_name = get_first_value_of_series(df_where_probe_maps["gene_name_truth_gene"])
            ref_gene_name = get_first_value_of_series(df_where_probe_maps["gene_name_ref_gene"])
            assert truth_gene_name == ref_gene_name
            edit_distance = get_first_value_of_series(df_where_probe_maps["edit_distance"])
            return truth_gene_name, edit_distance
        else:
            raise ProbeMapsToSeveralGenes()


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
    def get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(edit_distance_df, contig_vcf_probe_originates_from,
                                                                       pos_vcf_probe_originates_from):
        return GetVariantDistanceHelpers.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to\
            (edit_distance_df=edit_distance_df, contig_column="contig_ref_gene", start_gene_column="start_ref_gene",
             stop_gene_column="stop_ref_gene", contig_probe_originates_from=contig_vcf_probe_originates_from,
             pos_probe_originates_from=pos_vcf_probe_originates_from)