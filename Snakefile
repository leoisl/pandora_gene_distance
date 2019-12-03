from pathlib import Path
import itertools
from snakemake.utils import min_version, validate
import pandas as pd
from collections import defaultdict

min_version("5.1.0")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

validate(config, "pipeline/schemas/config.schema.yaml")
samples = pd.read_csv(config["samples"])
validate(samples, "pipeline/schemas/samples.schema.yaml")
samples.rename(columns={"reference": "reference_assembly"}, inplace=True)

variant_calls = pd.read_csv(config["variant_calls"])
validate(variant_calls, "pipeline/schemas/variant_calls.schema.yaml")
variant_calls.rename(columns={"reference": "vcf_reference"}, inplace=True)


# ======================================================
# Global variables
# ======================================================
output_folder = config['output_folder']
data: pd.DataFrame = pd.merge(variant_calls, samples, on="sample_id")
data = data.set_index(["sample_id", "coverage", "tool"], drop=False)
samples = samples.set_index(["sample_id"], drop=False)
sample_pairs = [(sample1, sample2) for sample1, sample2 in itertools.combinations(sorted(samples["sample_id"]), r=2)]
number_of_points_in_ROC_curve = config['number_of_points_in_ROC_curve']



# ======================================================
# Helper functions
# ======================================================
def get_coverage_filters(tool):
    return [str(elem) for elem in config['coverage_filters']]

def get_strand_bias_filters(tool):
    if tool.startswith("pandora"):
        return [str(elem) for elem in config['strand_bias_filters']]
    else:
        return ["Not_App"]

def get_gaps_filters(tool):
    if tool.startswith("pandora"):
        return [str(elem) for elem in config['gaps_filters']]
    else:
        return ["Not_App"]


# ======================================================
# Pipeline files
# ======================================================
files = []

# Common files
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    files_with_filters = expand(f"{output_folder}/variant_calls_probesets/{sample_id}/{coverage}/{tool}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/variant_calls_probeset.fa", coverage_threshold = get_coverage_filters(tool), strand_bias_threshold = get_strand_bias_filters(tool), gaps_threshold = get_gaps_filters(tool))
    files.extend(files_with_filters)



# Precision files
all_precision_files=[]
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    files_with_filters = expand(f"{output_folder}/precision/variant_calls_probesets_mapped_to_refs/{sample_id}/{coverage}/{tool}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/variant_calls_probeset_mapped.sam", coverage_threshold = get_coverage_filters(tool), strand_bias_threshold = get_strand_bias_filters(tool), gaps_threshold = get_gaps_filters(tool))
    all_precision_files.extend(files_with_filters)

cov_tool_and_filters_to_precision_report_files = defaultdict(list)
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for coverage_threshold in get_coverage_filters(tool):
        for strand_bias_threshold in get_strand_bias_filters(tool):
            for gaps_threshold in get_gaps_filters(tool):
                report_file = f"{output_folder}/precision/reports_from_probe_mappings/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/variant_calls_probeset_report.tsv"
                all_precision_files.append(report_file)
                cov_tool_and_filters_to_precision_report_files[(coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold)].append(report_file)

for coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold in cov_tool_and_filters_to_precision_report_files:
    all_precision_files.append(f"{output_folder}/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv")

files.extend(all_precision_files)



# Recall files
all_recall_files=[]
for sample1, sample2 in sample_pairs:
    all_recall_files.extend(
        [
            f"{output_folder}/recall/truth_probesets/{sample1}/{sample1}_and_{sample2}.truth_probeset.fa",
            f"{output_folder}/recall/truth_probesets/{sample2}/{sample1}_and_{sample2}.truth_probeset.fa",
        ]
    )

for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for sample1, sample2 in [pair for pair in sample_pairs if sample_id in pair]:
        filename_prefix = f"{sample1}_and_{sample2}"
        files_with_filters = expand(f"{output_folder}/recall/map_probes/{sample_id}/{coverage}/{tool}/coverage_filter_{{coverage_threshold}}/strand_bias_filter_{{strand_bias_threshold}}/gaps_filter_{{gaps_threshold}}/{filename_prefix}.sam", coverage_threshold = get_coverage_filters(tool), strand_bias_threshold = get_strand_bias_filters(tool), gaps_threshold = get_gaps_filters(tool))
        all_recall_files.extend(files_with_filters)

cov_tool_and_filters_to_recall_report_files = defaultdict(list)
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    for sample1, sample2 in [pair for pair in sample_pairs if sample_id in pair]:
        filename_prefix = f"{sample1}_and_{sample2}"
        for coverage_threshold in get_coverage_filters(tool):
            for strand_bias_threshold in get_strand_bias_filters(tool):
                for gaps_threshold in get_gaps_filters(tool):
                    report_file = f"{output_folder}/recall/reports/{sample_id}/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/{filename_prefix}.report.tsv"
                    all_recall_files.append(report_file)
                    cov_tool_and_filters_to_recall_report_files[(coverage, tool, str(coverage_threshold), str(strand_bias_threshold), str(gaps_threshold))].append(report_file)

for coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold in cov_tool_and_filters_to_recall_report_files:
    all_recall_files.append(f"{output_folder}/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv")

files.extend(all_recall_files)



# Plot files
all_plot_data_intermediate_files = []
for coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold in cov_tool_and_filters_to_recall_report_files:
    all_plot_data_intermediate_files.append(f"{output_folder}/plot_data/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.tsv")
final_plot_data_file = f"{output_folder}/plot_data/ROC_data.tsv"

files.extend(all_plot_data_intermediate_files)
files.append(final_plot_data_file)



# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("pipeline/rules/")
include: str(rules_dir / "common.smk")
include: str(rules_dir / "recall.smk")
include: str(rules_dir / "precision.smk")
include: str(rules_dir / "plot.smk")
