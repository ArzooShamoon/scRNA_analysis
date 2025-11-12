# scRNA_analysis
The evaluate_cluster_resolutions() function provides a quantitative and visual framework to assess clustering quality across multiple resolution settings in Seurat single-cell RNA-seq analyses.

This function computes several clustering evaluation metrics — including Silhouette Score, Davies–Bouldin Index (DBI), and Root Mean Square Deviation (RMSD) — for all clustering resolutions stored in a Seurat object’s metadata (e.g., integrated_snn_res.* columns).

It helps identify the optimal clustering resolution that balances between compact (tight) and well-separated clusters.
Key Features

Automatically detects all resolution columns in integrated@meta.data

Computes:

Global Silhouette Score — higher is better

Mean Cluster Silhouette — per-cluster quality measure

RMSD (within clusters) — lower is better

RMSD (between clusters) — higher is better

Davies–Bouldin Index (DBI) — lower is better (cluster separation)

Generates two dual-axis ggplot visualizations:

Silhouette vs DBI (blue vs red)

RMSD Between vs Within (green vs orange)

Returns both the metric table and the ggplot objects for further customization.



**Usage Example**
out <- evaluate_cluster_resolutions(
  seurat_obj = integrated,
  reduction_use = "pca",
  dims_use = 1:15
)


# View results
View(out$results)

# Display plots
print(out$plot_sil_dbi)

print(out$plot_rmsd)
