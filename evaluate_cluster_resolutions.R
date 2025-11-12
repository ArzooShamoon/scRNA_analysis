library(dplyr)
library(purrr)
library(cluster)
library(ggplot2)
library(scales)

evaluate_cluster_resolutions <- function(
  seurat_obj,
  reduction_use = "pca",
  dims_use = 1:15,
  pattern = "^integrated_snn_res"
) {
  # Identify clustering resolution columns
  res_cols <- grep(pattern, colnames(seurat_obj@meta.data), value = TRUE)
  if (length(res_cols) == 0) {
    stop("No clustering resolution columns found matching pattern: ", pattern)
  }
  
  # Helper function: RMSD within a cluster
  calc_rmsd_within <- function(emb_matrix) {
    centroid <- colMeans(emb_matrix)
    sqrt(mean((emb_matrix - matrix(centroid,
                                   nrow = nrow(emb_matrix),
                                   ncol = ncol(emb_matrix),
                                   byrow = TRUE))^2))
  }
  
  # Helper function: Davies-Bouldin Index
  calc_dbi <- function(emb_matrix, clust_factor) {
    clusters <- levels(clust_factor)
    k <- length(clusters)
    
    centroids <- map_dfr(clusters, function(cl) {
      mat <- emb_matrix[clust_factor == cl, , drop = FALSE]
      tibble(cluster = cl, centroid = list(colMeans(mat)))
    })
    
    scatters <- map_dbl(clusters, function(cl) {
      mat <- emb_matrix[clust_factor == cl, , drop = FALSE]
      centroid <- colMeans(mat)
      sqrt(mean(rowSums((mat - matrix(centroid, nrow(mat), ncol(mat), byrow = TRUE))^2)))
    })
    
    db_vals <- c()
    for (i in seq_len(k)) {
      max_ratio <- 0
      for (j in seq_len(k)) {
        if (i != j) {
          dist_ij <- sqrt(sum((centroids$centroid[[i]] - centroids$centroid[[j]])^2))
          ratio <- (scatters[i] + scatters[j]) / dist_ij
          max_ratio <- max(max_ratio, ratio)
        }
      }
      db_vals[i] <- max_ratio
    }
    
    mean(db_vals)
  }
  
  # Precompute PCA (or chosen reduction) embedding
  emb_all <- Embeddings(seurat_obj, reduction = reduction_use)[, dims_use, drop = FALSE]
  
  # Compute metrics for each resolution
  results <- map_dfr(res_cols, function(res_col) {
    clust_factor <- factor(seurat_obj[[res_col, drop = TRUE]])
    
    # Silhouette analysis
    sil <- silhouette(as.integer(clust_factor), dist(emb_all))
    sil_df <- as.data.frame(unclass(sil))
    colnames(sil_df) <- c("cluster_num", "neighbor", "sil_width")
    sil_df$subcluster <- levels(clust_factor)[sil_df$cluster_num]
    
    sil_global <- mean(sil_df$sil_width)
    sil_mean_clusters <- sil_df %>%
      group_by(subcluster) %>%
      summarise(mean_sil = mean(sil_width), .groups = "drop") %>%
      summarise(mean_of_means = mean(mean_sil)) %>%
      pull(mean_of_means)
    
    # Between-cluster RMSD
    cluster_means <- aggregate(emb_all, by = list(cluster = clust_factor), mean)
    cluster_matrix <- as.matrix(cluster_means[, -1])
    rmsd_between <- sqrt(mean(dist(cluster_matrix)^2))
    
    # Within-cluster RMSD
    rmsd_within <- map_dbl(split(emb_all, clust_factor), function(mat) {
      mat <- matrix(as.numeric(mat), ncol = ncol(emb_all), byrow = FALSE)
      calc_rmsd_within(mat)
    })
    rmsd_within_avg <- mean(rmsd_within)
    
    # Davies-Bouldin Index
    dbi <- calc_dbi(emb_all, clust_factor)
    
    tibble(
      resolution = res_col,
      Silhouette_global = sil_global,
      Silhouette_meanCluster = sil_mean_clusters,
      RMSD_between = rmsd_between,
      RMSD_within = rmsd_within_avg,
      DBI = dbi
    )
  })
  
  # Sort results numerically by resolution
  results <- results %>%
    mutate(res_num = as.numeric(sub(".*res\\.", "", resolution))) %>%
    arrange(res_num)
  
  # --- Plot 1: Silhouette vs DBI ---
  max_sil <- max(results$Silhouette_global)
  max_dbi <- max(results$DBI)
  scale_factor <- max_sil / max_dbi
  
  p1 <- ggplot(results, aes(x = res_num)) +
    geom_line(aes(y = Silhouette_global), color = "steelblue", size = 1.2) +
    geom_point(aes(y = Silhouette_global), color = "steelblue", size = 2) +
    geom_line(aes(y = DBI * scale_factor), color = "firebrick", size = 1.2, linetype = "dashed") +
    geom_point(aes(y = DBI * scale_factor), color = "firebrick", size = 2) +
    scale_y_continuous(
      name = "Silhouette Score (Higher is Better)",
      sec.axis = sec_axis(~ . / scale_factor, name = "DBI (Lower is Better)")
    ) +
    labs(
      title = "Global Silhouette Score vs. DBI Across Resolutions",
      x = "Clustering Resolution"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(color = "steelblue"),
      axis.title.y.right = element_text(color = "firebrick")
    )
  
  # --- Plot 2: RMSD Between vs Within ---
  max_between <- max(results$RMSD_between)
  max_within <- max(results$RMSD_within)
  scale_factor_rmsd <- max_between / max_within
  
  p2 <- ggplot(results, aes(x = res_num)) +
    geom_line(aes(y = RMSD_between), color = "forestgreen", size = 1.2) +
    geom_point(aes(y = RMSD_between), color = "forestgreen", size = 2) +
    geom_line(aes(y = RMSD_within * scale_factor_rmsd),
              color = "orange", size = 1.2, linetype = "dashed") +
    geom_point(aes(y = RMSD_within * scale_factor_rmsd),
               color = "orange", size = 2) +
    scale_y_continuous(
      name = "RMSD Between Clusters (Higher is Better)",
      sec.axis = sec_axis(~ . / scale_factor_rmsd,
                          name = "RMSD Within Clusters (Lower is Better)")
    ) +
    labs(
      title = "RMSD Between vs Within Across Resolutions",
      x = "Clustering Resolution"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(color = "forestgreen"),
      axis.title.y.right = element_text(color = "orange")
    )
  
  # Return both results and plots
  return(list(
    results = results,
    plot_sil_dbi = p1,
    plot_rmsd = p2
  ))
}



#out <- evaluate_cluster_resolutions(integrated, reduction_use = "pca", dims_use = 1:15)

# View table
#View(out$results)

# Show plots
#print(out$plot_sil_dbi)
#print(out$plot_rmsd)
