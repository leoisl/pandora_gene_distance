from unittest import TestCase
from scripts.GetVariantDistanceHelpers import GetVariantDistanceHelpers, GetVariantDistanceHelpersForPrecision, GetVariantDistanceHelpersForRecall, ProbeMapsToSeveralGenes, ProbeDoesNotMapToAnyGene
import pandas as pd


class TestGetVariantDistanceHelpersBaseClass(TestCase):
    def assert_df_are_equal(self, df1, df2):
        self.assertEqual(df1.to_dict("records"), df2.to_dict("records"))

    def setUp(self) -> None:
        self.empty_edit_distance_df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": [],
            "ref_or_truth_id_ref_gene": [],
            "status_truth_gene": [],
            "status_ref_gene": [],
            "contig_ref_gene": [],
            "start_ref_gene": [],
            "stop_ref_gene": [],
        })



class TestGetVariantDistanceHelpers(TestGetVariantDistanceHelpersBaseClass):
    def test___get_edit_distance_df_given_truth_and_ref___empty_df(self):
        df = self.empty_edit_distance_df

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = self.empty_edit_distance_df

        self.assert_df_are_equal(actual, expected)

    def test___get_edit_distance_df_given_truth_and_ref___one_record___all_fields_are_good(self):
        df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["truth"],
            "ref_or_truth_id_ref_gene": ["ref"],
            "status_truth_gene": ["Mapped"],
            "status_ref_gene": ["Mapped"],
        })

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = df

        self.assert_df_are_equal(actual, expected)

    def test___get_edit_distance_df_given_truth_and_ref___one_record___ref_or_truth_id_truth_gene_not_good(self):
        df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["bad_truth"],
            "ref_or_truth_id_ref_gene": ["ref"],
            "status_truth_gene": ["Mapped"],
            "status_ref_gene": ["Mapped"],
        })

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = self.empty_edit_distance_df

        self.assert_df_are_equal(actual, expected)

    def test___get_edit_distance_df_given_truth_and_ref___one_record___ref_or_truth_id_ref_gene_not_good(self):
        df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["truth"],
            "ref_or_truth_id_ref_gene": ["bad_ref"],
            "status_truth_gene": ["Mapped"],
            "status_ref_gene": ["Mapped"],
        })

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = self.empty_edit_distance_df

        self.assert_df_are_equal(actual, expected)

    def test___get_edit_distance_df_given_truth_and_ref___one_record___status_truth_gene_is_unmapped(self):
        df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["truth"],
            "ref_or_truth_id_ref_gene": ["ref"],
            "status_truth_gene": ["Unmapped"],
            "status_ref_gene": ["Mapped"],
        })

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = self.empty_edit_distance_df

        self.assert_df_are_equal(actual, expected)

    def test___get_edit_distance_df_given_truth_and_ref___one_record___status_ref_gene_is_unmapped(self):
        df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["truth"],
            "ref_or_truth_id_ref_gene": ["ref"],
            "status_truth_gene": ["Mapped"],
            "status_ref_gene": ["Unmapped"],
        })

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = self.empty_edit_distance_df

        self.assert_df_are_equal(actual, expected)


    def test___get_edit_distance_df_given_truth_and_ref___several_records___two_are_fine___four_are_not(self):
        df = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["truth",    "truth",   "truth",  "truth",    "bad_truth", "truth"],
            "ref_or_truth_id_ref_gene":   ["ref",      "bad_ref", "ref",    "ref",      "ref",       "ref"],
            "status_truth_gene":          ["Unmapped", "Mapped",  "Mapped", "Mapped",   "Mapped",    "Mapped"],
            "status_ref_gene":            ["Mapped",   "Mapped",  "Mapped", "Unmapped", "Mapped",    "Mapped"],
        })

        actual = GetVariantDistanceHelpers.get_edit_distance_df_given_truth_and_ref(df, "truth", "ref")
        expected = pd.DataFrame({
            "ref_or_truth_id_truth_gene": ["truth"]*2,
            "ref_or_truth_id_ref_gene":   ["ref"]*2,
            "status_truth_gene":          ["Mapped"]*2,
            "status_ref_gene":            ["Mapped"]*2,
        })

        self.assert_df_are_equal(actual, expected)



class TestGetVariantDistanceHelpersForPrecision(TestGetVariantDistanceHelpersBaseClass):
    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___empty_edit_distance_df(self):
        df = self.empty_edit_distance_df
        contig_vcf_probe_originates_from = "contig"
        pos_vcf_probe_originates_from = 10

        with self.assertRaises(ProbeDoesNotMapToAnyGene):
            GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(
                df, contig_vcf_probe_originates_from, pos_vcf_probe_originates_from)


    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___good_contig_and_pos_just_before_left_border(self):
        edit_distance = 200
        df = pd.DataFrame({
            "contig_ref_gene": ["contig"],
            "start_ref_gene": [10],
            "stop_ref_gene": [20],
            "edit_distance": [edit_distance],
        })
        contig_vcf_probe_originates_from = "contig"
        pos_vcf_probe_originates_from = 9

        with self.assertRaises(ProbeDoesNotMapToAnyGene):
            GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(df, contig_vcf_probe_originates_from, pos_vcf_probe_originates_from)


    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___good_contig_and_pos_just_on_left_border(self):
        df = pd.DataFrame({
            "contig_ref_gene": ["contig"],
            "start_ref_gene": [10],
            "stop_ref_gene": [20],
            "edit_distance": [200],
            "gene_name_truth_gene": ["gene"],
            "gene_name_ref_gene": ["gene"],
        })
        contig_vcf_probe_originates_from = "contig"
        pos_vcf_probe_originates_from = 10

        actual = GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(
                df, contig_vcf_probe_originates_from, pos_vcf_probe_originates_from)

        self.assertEqual(("gene", 200), actual)

    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___good_contig_and_pos_just_on_right_border(self):
        df = pd.DataFrame({
            "contig_ref_gene": ["contig"],
            "start_ref_gene": [10],
            "stop_ref_gene": [20],
            "edit_distance": [200],
            "gene_name_truth_gene": ["gene"],
            "gene_name_ref_gene": ["gene"],
        })
        contig_vcf_probe_originates_from = "contig"
        pos_vcf_probe_originates_from = 19

        actual = GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(
                df, contig_vcf_probe_originates_from, pos_vcf_probe_originates_from)

        self.assertEqual(("gene", 200), actual)

    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___good_contig_and_pos_just_after_right_border(self):
        edit_distance = 200
        df = pd.DataFrame({
            "contig_ref_gene": ["contig"],
            "start_ref_gene": [10],
            "stop_ref_gene": [20],
            "edit_distance": [edit_distance],
        })
        contig_vcf_probe_originates_from = "contig"
        pos_vcf_probe_originates_from = 20

        with self.assertRaises(ProbeDoesNotMapToAnyGene):
            GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(df,
                                                                                                                 contig_vcf_probe_originates_from,
                                                                                                                 pos_vcf_probe_originates_from)


    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___two_contigs_with_same_positions___only_maps_to_one(self):
        df = pd.DataFrame({
            "contig_ref_gene": ["contig_2", "contig_1"],
            "start_ref_gene": [10, 10],
            "stop_ref_gene": [20, 20],
            "edit_distance": [100, 185],
            "gene_name_truth_gene": ["gene_2", "gene_1"],
            "gene_name_ref_gene": ["gene_2", "gene_1"],
        })
        contig_vcf_probe_originates_from = "contig_1"
        pos_vcf_probe_originates_from = 15

        actual = GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(df,
                                                                                                                      contig_vcf_probe_originates_from,
                                                                                                                      pos_vcf_probe_originates_from)
        self.assertEqual(("gene_1", 185), actual)



    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___two_contigs_with_same_positions___maps_to_none(self):
        df = pd.DataFrame({
            "contig_ref_gene": ["contig_2", "contig"],
            "start_ref_gene": [10, 10],
            "stop_ref_gene": [20, 20],
            "edit_distance": [100, 185],
        })
        contig_vcf_probe_originates_from = "unmapped_contig"
        pos_vcf_probe_originates_from = 15

        with self.assertRaises(ProbeDoesNotMapToAnyGene):
            GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(df,
                                                                                                                 contig_vcf_probe_originates_from,
                                                                                                                 pos_vcf_probe_originates_from)


    def test___get_edit_distance_of_gene_this_vcf_probe_maps_to___position_mapping_to_two_genes(self):
        df = pd.DataFrame({
            "contig_ref_gene": ["contig", "contig"],
            "start_ref_gene": [10, 19],
            "stop_ref_gene": [20, 30],
            "edit_distance": [100, 185],
        })
        contig_vcf_probe_originates_from = "contig"
        pos_vcf_probe_originates_from = 19

        with self.assertRaises(ProbeMapsToSeveralGenes):
            GetVariantDistanceHelpersForPrecision.get_gene_name_and_edit_distance_of_gene_this_vcf_probe_maps_to(df,
                                                                                                                 contig_vcf_probe_originates_from,
                                                                                                                 pos_vcf_probe_originates_from)


class TestGetVariantDistanceHelpersForRecall(TestGetVariantDistanceHelpersBaseClass):
    def test___get_edit_distance_of_gene_this_truth_probe_maps_to___only_maps_to_one(self):
        df = pd.DataFrame({
            "contig_truth_gene": ["contig_1", "contig_1"],
            "gene_name_truth_gene": ["gene_1", "gene_2"],
            "gene_name_ref_gene": ["gene_1", "gene_2"],
            "edit_distance": [100, 185],
            "start_ref_gene": [10, 12],  # should not be used here
            "stop_ref_gene": [20, 30],  # should not be used here
            "start_truth_gene": [1, 11],  # test should check this is used
            "stop_truth_gene": [10, 100],  # test should check this is used
        })
        contig_truth_probe_originates_from = "contig_1"
        pos_truth_probe_originates_from = 15

        actual = GetVariantDistanceHelpersForRecall.get_gene_name_and_edit_distance_of_gene_this_truth_probe_maps_to\
            (df, contig_truth_probe_originates_from, pos_truth_probe_originates_from)
        self.assertEqual(("gene_2", 185), actual)
