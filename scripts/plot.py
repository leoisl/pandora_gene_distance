import seaborn as sns
sns.set()
import pandas as pd
import matplotlib.pyplot as plt

def get_bounds(step, max_bin):
    bounds = [0.0]
    current_step = step
    while current_step <= max_bin:
        bounds.append(current_step)
        current_step += step
    assert len(bounds) > 0
    return bounds

def get_edit_distance_label(row, bounds):
    distance = row["distance"]
    for upper_bound in bounds[1:]:
        if distance < upper_bound:
            return upper_bound
    return bounds[-1]

def get_df_with_edit_distance_labels (df, having_at_least_n_variants, step, max_bin):
    df = df.query("max_recall >= @having_at_least_n_variants")
    bounds = get_bounds(step, max_bin)
    edit_distance_labels = df.apply(get_edit_distance_label, axis=1, bounds=bounds)
    df["edit_distance_labels"] = edit_distance_labels
    return df

def init_plot(figsize=(20, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax

def save_figure(fig, output_filepath):
    fig.savefig(output_filepath, format="png", dpi=600)

def plot_lineplot_mean_recall_and_count_in_genes_in_several_bins(df, step, output_filepath):
    fig, ax = init_plot()

    sns.lineplot(x="edit_distance_labels", y="recall_ratio", data=df, color='b')
    count_df = df.groupby(by="edit_distance_labels", as_index=False).count()
    ax2 = ax.twinx()
    sns.lineplot(x="edit_distance_labels", y="gene", data=count_df, ax=ax2, color='r')

    ax.set(xlabel=f"Edit distance (gene bins at each {step * 100}%)", ylabel='Mean recall score per bin')
    ax.set_ylim([0.0, 1.0])
    ax.set_xlim([0.0, 0.2])
    ax2.set(ylabel='Number of genes per bin')

    save_figure(fig, output_filepath)

def plot_violinplot_recall_in_genes_in_several_bins(df, step, scale, output_filepath):
    fig, ax = init_plot()
    sns.violinplot(x="edit_distance_labels", y="recall_ratio", data=df, scale=scale)
    ax.set(xlabel=f"Edit distance (gene bins at each {step*100}%)", ylabel='Recall ratio per bin')
    ax.set_ylim([-0.25, 1.25])
    ax.set_xlim([0.0, 0.2])
    save_figure(fig, output_filepath)

def plot_scatterplot_for_precision(gene_truth_ref_measure_and_distance, output_filepath):
    fig, ax = init_plot()
    df = pd.read_csv(gene_truth_ref_measure_and_distance)
    sns.scatterplot(x=df["distance"], y=df["precision_ratio"], alpha=0.3, color="black", marker=".")
    save_figure(fig, output_filepath)
