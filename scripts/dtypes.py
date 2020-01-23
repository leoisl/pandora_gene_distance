edit_distance_dtype_dict = {
    "gene_name": str,
    "status_truth_gene": str,
    "ref_or_truth_id_truth_gene": str,
    "contig_truth_gene": str,
    "start_truth_gene": int,
    "stop_truth_gene": int,
    "sequence_truth_gene": str,
    "status_ref_gene": str,
    "ref_or_truth_id_ref_gene": str,
    "contig_ref_gene": str,
    "start_ref_gene": int,
    "stop_ref_gene": int,
    "sequence_ref_gene": str,
    "edit_distance": float
}

mapped_genes_dtype_dict = {
    "status": str,
    "gene_name": str,
    "ref_or_truth_id": str,
    "contig": str,
    "start": int,
    "stop": int,
    "sequence": str
}

variant_precision_score_distance_csv_dtype_dict = {
    "gene": str,
    "truth": str,
    "ref": str,
    "variant": str,
    "precision_score": float,
    "distance": float
}

variant_recall_score_distance_csv_dtype_dict = {
    "gene": str,
    "truth": str,
    "ref": str,
    "variant": str,
    "recall_score": float,
    "distance": float
}


gene_truth_ref_precision_dtype_dict = {
  "gene": str,
  "truth": str,
  "ref": str,
  "distance": float,
  "max_precision": float,
  "observed_precision": float,
  "precision_ratio": float,
}


gene_truth_ref_recall_dtype_dict = {
  "gene": str,
  "truth": str,
  "ref": str,
  "distance": float,
  "max_recall": float,
  "observed_recall": float,
  "recall_ratio": float,
}