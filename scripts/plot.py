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

def get_edit_distances(row, bounds):
    distance = row["distance"]
    for upper_bound in bounds[1:]:
        if distance < upper_bound:
            return upper_bound
    return bounds[-1]

def get_df_with_edit_distance_labels (df, max_precision_or_recall_column, having_at_least_n_variants, step, max_bin):
    df = df.query(f"{max_precision_or_recall_column} >= @having_at_least_n_variants")
    bounds = get_bounds(step, max_bin)
    edit_distance_labels = df.apply(get_edit_distances, axis=1, bounds=bounds)
    df["edit_distance_labels"] = edit_distance_labels
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

def plot_lineplot_mean_precision_or_recall_and_count_in_genes_in_several_bins(df, recall_or_precision_ratio, step, output_filepath):
    fig, ax = init_plot()

    sns.lineplot(x="edit_distance_labels", y=recall_or_precision_ratio, data=df, color='b')
    count_df = df.groupby(by="edit_distance_labels", as_index=False).count()
    ax2 = ax.twinx()
    sns.lineplot(x="edit_distance_labels", y="gene", data=count_df, ax=ax2, color='r')

    ax.set(xlabel=f"Edit distance (gene bins at each {step * 100}%)", ylabel=f'Mean {recall_or_precision_ratio} per bin')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlim([0.0, 0.25])
    ax2.set(ylabel='Number of genes per bin')

    save_figure(fig, output_filepath)

def plot_lineplot_mean_precision_and_count_in_genes_in_several_bins(df, step, output_filepath):
    plot_lineplot_mean_precision_or_recall_and_count_in_genes_in_several_bins(df, "precision_ratio", step, output_filepath)

def plot_lineplot_mean_recall_and_count_in_genes_in_several_bins(df, step, output_filepath):
    plot_lineplot_mean_precision_or_recall_and_count_in_genes_in_several_bins(df, "recall_ratio", step, output_filepath)

def plot_violinplot_precision_or_recall_in_genes_in_several_bins(df, recall_or_precision_ratio, step, scale, output_filepath):
    fig, ax = init_plot()
    sns.violinplot(x="edit_distance_labels", y=recall_or_precision_ratio, data=df, scale=scale)
    ax.set(xlabel=f"Edit distance (gene bins at each {step*100}%)", ylabel=f'{recall_or_precision_ratio} per bin')
    ax.set_ylim([-0.25, 1.25])
    save_figure(fig, output_filepath)

def plot_violinplot_precision_in_genes_in_several_bins(df, step, scale, output_filepath):
    plot_violinplot_precision_or_recall_in_genes_in_several_bins(df, "precision_ratio", step, scale, output_filepath)

def plot_violinplot_recall_in_genes_in_several_bins(df, step, scale, output_filepath):
    plot_violinplot_precision_or_recall_in_genes_in_several_bins(df, "recall_ratio", step, scale, output_filepath)