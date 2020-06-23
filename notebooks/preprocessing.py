#!/usr/bin/env python
# coding: utf-8

# # Gene distance plots exploratory data analysis notebook
# 
# ## What are the gene distance plots?
# 
# An attempt to investigate if the precision and recall of variant callers is affected by the distance between samples and references. However here we focus on regions that are common between the samples and the refs, thus we do not address the issues of reference incompleteness (i.e. we don't analyse regions that exclusive to the sample or the ref).
# 
# ## How do we do this?
# 
# We first identify genes that are common between a pair `(sample, ref)`. For `snippy`, `ref` is each 226 references we use. For pandora, `ref` is simply the `pandora_multisample.vcf_ref.fa`. To identify these genes, we need a gene sequence that is as close as possible to the samples/refs. Thus we use the gene sequences in the `pandora_multisample.vcf_ref.fa`. We map the genes in `pandora_multisample.vcf_ref.fa` for each pair `(sample, ref)`, allowing us to know if the gene appears in the sample and its coordinates, in the ref and its coordinates, and the edit distance between them (in case they appear in both sample and ref).
# 
# Using the gene coordinates, we can infer how many variant calls lie on each identified gene, and if they are correct or incorrect (using the 4-way pipeline assessment), allowing us to know the precision ratio of the caller in that gene. The same for recall - we know how many truth variants lie on the gene, and how many the caller discovered. We plot these values by bins of genes split by its edit distance between the sample and the refs.
# 
# All of this is comprised in a snakemake pipeline, and this EDA notebook analyse the pipeline's output to check what is the best way to convey this data.

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from pathlib import Path
import plotly.express as px


# # Main configs

# In[2]:


testing = False
data_path = Path("/hps/nobackup/iqbal/leandro/snippy_calls_gene_distance")
tools = ["pandora", "snippy", "samtools"]
colors = ["blue", "red", "green"]

plt.rcParams['figure.dpi'] = 600
sns.set()


# # Data input
# 
# For each triple `(gene, truth, reference)`, we know the nb of truth variants in that gene (`max_recall`), the nb of truth variants the tool was able to call (`observed_recall`) and thus the recall performance (`recall_ratio = observed_recall / max_recall`). We also know the edit distance of the gene between the ref and the truth (`distance`).
# 
# We also have the same stuff for precision.

# In[3]:


print("Loading...")

df_pandora_precision = pd.read_csv(data_path / "gene_distance_plots_20_way_pandora_illumina_100x/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv")
df_pandora_precision["tool"] = "pandora"
df_snippy_precision = pd.read_csv(data_path / "gene_distance_plots_20_way_snippy_illumina_100x/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv")
df_snippy_precision["tool"] = "snippy"
df_samtools_precision = pd.read_csv(data_path / "gene_distance_plots_20_way_samtools_illumina_100x/get_gene_truth_ref_precision_proportion_distance/all_gene_truth_ref_precision_proportion_distance.csv")
df_samtools_precision["tool"] = "samtools"
df_precision = pd.concat([df_pandora_precision, df_snippy_precision, df_samtools_precision], ignore_index = True)
# display(df_precision)

df_pandora_recall = pd.read_csv(data_path / "gene_distance_plots_20_way_pandora_illumina_100x/get_gene_truth_ref_recall_proportion_distance/all_gene_truth_ref_recall_proportion_distance.csv")
df_pandora_recall["tool"] = "pandora"
df_snippy_recall = pd.read_csv(data_path / "gene_distance_plots_20_way_snippy_illumina_100x/get_gene_truth_ref_recall_proportion_distance/all_gene_truth_ref_recall_proportion_distance.csv")
df_snippy_recall["tool"] = "snippy"
df_samtools_recall = pd.read_csv(data_path / "gene_distance_plots_20_way_samtools_illumina_100x/get_gene_truth_ref_recall_proportion_distance/all_gene_truth_ref_recall_proportion_distance.csv")
df_samtools_recall["tool"] = "samtools"
df_recall = pd.concat([df_pandora_recall, df_snippy_recall, df_samtools_recall], ignore_index = True)
# display(df_recall)


# In[4]:


# test data
if testing:
    tools = ["snippy"]
    colors = ["red"]
    
#     df_pandora_precision = pd.read_csv("eda_data/all_gene_truth_ref_precision_proportion_distance.csv")
#     df_pandora_precision["tool"] = "pandora"
    df_snippy_precision = pd.read_csv("eda_data/all_gene_truth_ref_precision_proportion_distance.csv")
    df_snippy_precision["tool"] = "snippy"
#     df_samtools_precision = pd.read_csv("eda_data/all_gene_truth_ref_precision_proportion_distance.csv")
#     df_samtools_precision["tool"] = "samtools"
#     df_precision = pd.concat([df_pandora_precision, df_snippy_precision, df_samtools_precision], ignore_index = True)
    df_precision = df_snippy_precision
    display(df_precision)

#     df_pandora_recall = pd.read_csv("eda_data/all_gene_truth_ref_recall_proportion_distance.csv")
#     df_pandora_recall["tool"] = "pandora"
    df_snippy_recall = pd.read_csv("eda_data/all_gene_truth_ref_recall_proportion_distance.csv")
    df_snippy_recall["tool"] = "snippy"
#     df_samtools_recall = pd.read_csv("eda_data/all_gene_truth_ref_recall_proportion_distance.csv")
#     df_samtools_recall["tool"] = "samtools"
#     df_recall = pd.concat([df_pandora_recall, df_snippy_recall, df_samtools_recall], ignore_index = True)
    df_recall = df_snippy_recall
    display(df_recall)


# # Main helper functions (please skip, go direct to the Results)

# In[6]:


from matplotlib.lines import Line2D

def get_bounds(step, max_bin):
    bounds = [0.0]
    current_step = step
    while current_step <= max_bin:
        bounds.append(current_step)
        current_step += step
    assert len(bounds) > 0
    return bounds

def get_edit_distances(row, bounds):
    distance = row["distance"]
    for upper_bound in bounds[1:]:
        if distance < upper_bound:
            return upper_bound
    return bounds[-1]

def round_ed_label_and_return_str(row):
    return "%.2f"%round(row["edit_distance_labels"], 2)

def get_df_with_edit_distance_labels (df, max_precision_or_recall_column, having_at_least_n_variants, step, max_bin):
    df = df.query(f"{max_precision_or_recall_column} >= @having_at_least_n_variants")
    bounds = get_bounds(step, max_bin)
    edit_distance_labels = df.apply(get_edit_distances, axis=1, bounds=bounds)
    df["edit_distance_labels"] = edit_distance_labels
    df["edit_distance_labels_as_str"] = df.apply(round_ed_label_and_return_str, axis=1)
    return df

def get_df_with_edit_distance_labels_for_precision (df, step, max_bin):
    return get_df_with_edit_distance_labels(df, max_precision_or_recall_column="max_precision",
                                     having_at_least_n_variants=1, step=step, max_bin=max_bin)

def get_df_with_edit_distance_labels_for_recall (df, step, max_bin):
    return get_df_with_edit_distance_labels(df, max_precision_or_recall_column="max_recall",
                                     having_at_least_n_variants=1, step=step, max_bin=max_bin)

def init_plot(figsize=(20, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax

def save_figure(fig, output_filepath):
    fig.savefig(output_filepath, format="png", dpi=600)
    
def plot_lineplot_count_in_datapoints_in_several_bins(count_df, step, display_plot=False, display_dfs=False, output_filepath=None):
    fig, ax = init_plot()
    ax.set_xlabel(f"Edit distance (gene bins at each {step * 100}%)")
    ax.set_xlim([0.0, 0.2])
    ax.set_ylabel("Number of datapoints (gene, truth, ref) per bin")
    lineplot = sns.lineplot(x="edit_distance_labels", y="nb_of_datapoints", data=count_df, hue="tool",
                            palette=colors, alpha=0.3)
    lineplot.set_yscale("log")
    
    if output_filepath is not None:
        save_figure(fig, output_filepath)
        
    if display_plot:
        ax.plot()
    
    if display_dfs:
        display_full(count_df.query("tool == 'pandora'"))
        display_full(count_df.query("tool == 'snippy'"))
        display_full(count_df.query("tool == 'samtools'"))

def get_edit_distance_threshold_where_all_tools_have_at_least_the_min_nb_of_datapoints(count_df, min_nb_of_datapoints, nb_of_tools):
    edit_distance_and_nb_of_tools = count_df.query("nb_of_datapoints >= @min_nb_of_datapoints").groupby("edit_distance_labels").count()
    del edit_distance_and_nb_of_tools["nb_of_datapoints"]
    edit_distance_threshold_to_use = max(edit_distance_and_nb_of_tools.query("tool == @nb_of_tools").index)
    return edit_distance_threshold_to_use


def plot_violin_box_to_axes(df, tool, color, ax):
    df = df.query("tool == @tool")
    
    sns.violinplot(x="edit_distance_labels_as_str", y="recall_ratio", data=df,
                   scale="count", cut=0, color=color, width=0.9/len(tools), inner=None, linewidth=1.0,
                   ax=ax)
        
    sns.boxplot(x="edit_distance_labels_as_str", y="recall_ratio", data=df,
             color=color, width=0.02, linewidth=0.3, showfliers=False, ax=ax)
    
    plt.setp(ax.collections, alpha=.3)
    ax.set_ylim([-0.05, 1.05])
    
    
def plot_lines_to_axes(df, tool, color, ax):
    df = df.query("tool == @tool")
    
    sns.lineplot(x="edit_distance_labels_as_str", y="recall_ratio", data=df,
                 estimator=np.median, color=color, linewidth=1.0, ci=None, alpha=0.3, ax=ax)
    sns.lineplot(x="edit_distance_labels_as_str", y="recall_ratio", alpha=0.3, data=df,
                 estimator=np.mean, color=color, linewidth=1.0, ci=None, ax=ax)
    ax.lines[-1].set_linestyle(":")
    ax.set_ylim([-0.05, 1.05])


def plot_violin_for_recall_in_genes_in_several_bins(df, step, edit_distance_threshold, display_plot=False, output_filepath=None):
    fig, ax = init_plot((40, 3))
    df_with_edit_distance_threshold = df.query("edit_distance_labels <= @edit_distance_threshold")
    sorted_ed_labels = sorted(df_with_edit_distance_threshold["edit_distance_labels_as_str"].unique())
    
    
    plot_violin_box_to_axes(df_with_edit_distance_threshold, tools[0], colors[0], ax)
    tool_index=1
    for tool, color in zip(tools[1:], colors[1:]):
        ax_bbox = ax.get_position()
        ax_2_bbox = ax_bbox
        x_shift = (ax_bbox.x1-ax_bbox.x0)/len(sorted_ed_labels)/len(tools)
        ax_2_bbox.x0 += x_shift * tool_index
        ax_2_bbox.x1 += x_shift * tool_index
        ax_2 = fig.add_axes(ax_2_bbox, frameon=False)
        plot_violin_box_to_axes(df_with_edit_distance_threshold, tool, color, ax_2)
        ax_2.get_xaxis().set_visible(False)
        ax_2.get_yaxis().set_visible(False)
        tool_index += 1
    

    ax.set(xlabel=f"Edit distance (gene bins at each {step*100}%)", ylabel='Recall ratio per bin')    
    ax.xaxis.set(ticks=np.arange(1/len(sorted_ed_labels)/2+ax_bbox.x0, len(sorted_ed_labels)+ax_bbox.x0), ticklabels=sorted_ed_labels)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16) 
    ax.xaxis.label.set_fontsize(30) 
    ax.yaxis.label.set_fontsize(26) 
    
    legend_elements = [Line2D([0], [0], color=color, lw=4, label=tool)
                       for tool, color in zip(tools, colors)]
    ax.legend(handles=legend_elements, loc="lower left")
   
    if output_filepath is not None:
        save_figure(fig, output_filepath)
        
    if display_plot:
        plt.show()

        
def plot_violin_box_to_axes_for_a_tool(df, tool, ax):
    df = df.query("tool == @tool")
    
    sns.violinplot(x="edit_distance_labels_as_str", y="recall_ratio", data=df,
                   cut=0, inner=None, linewidth=1.0, ax=ax)
    ax.set_ylim([-0.05, 1.05])

def plot_violin_for_recall_in_genes_in_several_bins_for_a_tool(df, tool, step, edit_distance_threshold, display_plot=False, output_filepath=None):
    fig, ax = init_plot((20, 3))
    df_with_edit_distance_threshold = df.query("edit_distance_labels <= @edit_distance_threshold")
    sorted_ed_labels = sorted(df_with_edit_distance_threshold["edit_distance_labels_as_str"].unique())
    
    
    plot_violin_box_to_axes_for_a_tool(df_with_edit_distance_threshold, tool, ax)
    ax.set(xlabel=f"Edit distance (gene bins at each {step*100}%)", ylabel='Recall ratio per bin', title=tool)    
   
    if output_filepath is not None:
        save_figure(fig, output_filepath)
        
    if display_plot:
        plt.show()
        
        
def plot_line_in_genes_in_several_bins(df, measure, step, edit_distance_threshold, display_plot=False, output_filepath=None):
    fig, ax = init_plot((10, 3))
    df_with_edit_distance_threshold = df.query("edit_distance_labels <= @edit_distance_threshold")
    
    sns.lineplot(x="edit_distance_labels", y=measure, data=df_with_edit_distance_threshold, hue="tool",
                 estimator=np.mean, linewidth=1.0, ci=None, palette=colors, alpha=0.3, ax=ax)
#     ax_2 = fig.add_axes(ax.get_position(), frameon=False)
#     sns.lineplot(x="edit_distance_labels", y=measure, data=df_with_edit_distance_threshold, hue="tool",
#                  estimator=np.median, linewidth=1.0, ci=None, ax=ax_2, palette=colors, alpha=0.3, legend=False)
#     ax_2.lines[0].set_linestyle(":")
#     ax_2.lines[1].set_linestyle(":")
#     ax_2.get_xaxis().set_visible(False)
#     ax_2.get_yaxis().set_visible(False)

    ax.set(xlabel=f"Edit distance (gene bins at each {step*100}%)", ylabel=f'{measure} per bin')    
       
    if output_filepath is not None:
        save_figure(fig, output_filepath)
        
    if display_plot:
        plt.show()

        
def plot_line_for_precision_in_genes_in_several_bins(df, step, edit_distance_threshold, display_plot=False, output_filepath=None):
    plot_line_in_genes_in_several_bins(df, "precision_ratio", step, edit_distance_threshold, display_plot, output_filepath)

def plot_line_for_recall_in_genes_in_several_bins(df, step, edit_distance_threshold, display_plot=False, output_filepath=None):
    plot_line_in_genes_in_several_bins(df, "recall_ratio", step, edit_distance_threshold, display_plot, output_filepath)


# # Cached dfs for easier processing (please skip, go direct to the Results)

# In[7]:


def get_count_df(df):
    count_df = df.groupby(by=["edit_distance_labels", "tool"], as_index=False).count()
    count_df = count_df[["edit_distance_labels", "tool", "gene"]]
    count_df.columns = ["edit_distance_labels", "tool", "nb_of_datapoints"]
    return count_df


print("Caching...")
df_precision_with_step_001 = get_df_with_edit_distance_labels_for_precision(df_precision, step = 0.01, max_bin = 0.2)
count_df_precision_with_step_001 = get_count_df(df_precision_with_step_001)
df_recall_with_step_001 = get_df_with_edit_distance_labels_for_recall(df_recall, step = 0.01, max_bin = 0.2)
count_df_recall_with_step_001 = get_count_df(df_recall_with_step_001)


# # Serializing stuff

# In[8]:


df_precision.to_pickle("df_precision.pickle")
df_recall.to_pickle("df_recall.pickle")
df_precision_with_step_001.to_pickle("df_precision_with_step_001.pickle")
count_df_precision_with_step_001.to_pickle("count_df_precision_with_step_001.pickle")
df_recall_with_step_001.to_pickle("df_recall_with_step_001.pickle")
count_df_recall_with_step_001.to_pickle("count_df_recall_with_step_001.pickle")

