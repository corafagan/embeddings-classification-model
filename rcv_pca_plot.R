# ===============================================================
# Script Name:     rcv_pca_plot.R
# Project:         Embedding Visualization (Control vs Treatment)
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This function projects high-dimensional document embeddings for
#   two groups (e.g., control vs. treatment) into 2D space and produces
#   a publication-ready scatterplot. It supports PCA (fast via irlba)
#   and UMAP projections, optional balancing/downsampling, centroid
#   visualization, and bootstrap confidence ellipses.
#
# Expected Inputs:
#   - emb_a:
#       • Numeric matrix of embeddings (rows = documents, cols = dims).
#       • Corresponds to group labels[1] (default = "control").
#   - emb_b:
#       • Numeric matrix of embeddings (rows = documents, cols = dims).
#       • Corresponds to group labels[2] (default = "treatment").
#
# Key Parameters:
#   - labels (default = c("control","treatment")):
#       • Group names (used in plotting legend and centroid labels).
#   - projection (default = "pca"):
#       • Either "pca" (irlba-based PCA) or "umap" (uwot package).
#   - balance_for_plot (default = TRUE):
#       • Subsample the larger group to match the smaller, preventing
#         visual dominance.
#   - sample_size (default = NULL):
#       • Optional cap on total points plotted; each group gets half.
#   - zoom_to_centroids (default = TRUE):
#       • If TRUE, axis limits zoom to include centroids with padding.
#   - ellipse_group (default = "none"):
#       • If "a" or "b", draws a bootstrap 95% CI ellipse around the
#         centroid of that group.
#   - B_boot (default = 500):
#       • Number of bootstrap resamples for the CI ellipse.
#   - stats (optional):
#       • List of summary statistics to annotate on the plot, e.g.:
#         list(diff_mean=..., p_diff=..., cent_cos=..., p_centroid=...).
#   - save_path (default = NULL):
#       • If non-NULL, saves plot to file (PNG, PDF, etc.).
#   - base_family (default = "Times New Roman"):
#       • Font family for plot text.
#   - base_size (default = 12):
#       • Base font size for plot.
#
# Outputs:
#   - Returns (invisibly) a ggplot2 object representing the scatterplot.
#   - Side effects:
#       • Optionally saves plot to disk if save_path is provided.
#       • Prints progress messages at each step.
#
# Notes:
#   - Groups are aligned to the minimum shared dimensionality.
#   - Centroid separation is visualized with a dashed line.
#   - Bootstrap CI ellipse uses multivariate normal approximation
#     (via eigen decomposition of bootstrap covariance).
#   - Default aesthetics:
#       • Group A (labels[1]) = circle marker
#       • Group B (labels[2]) = triangle marker
# ===============================================================

# load packages
library(dplyr)
library(ggplot2)
library(ggforce)
library(irlba)

plot_embeddings_2d <- function(
    emb_a, emb_b,
    labels = c("control","treatment"),
    projection = c("pca","umap"),
    balance_for_plot = TRUE,
    zoom_to_centroids = TRUE,
    ellipse_group = c("none","a","b"),   # bootstrap CI ellipse around centroid
    ellipse_level = 0.95,
    B_boot = 500,
    point_alpha = 0.1,
    base_family = "Times New Roman",
    base_size = 12,
    stats = NULL,                        # list(diff_mean=..., p_diff=..., cent_cos=..., p_centroid=...)
    stats_title = "Document Embeddings (2D)",
    stats_subtitle = "Shapes by group; dashed: centroid separation",
    save_path = NULL, width = 10, height = 8, dpi = 300,
    seed = 123, 
    sample_size = NULL
){
  
  # ---- Validate inputs & normalize args ---------------------------------------
  # Expect two numeric matrices of embeddings, one per group (A=labels[1], B=labels[2]).
  stopifnot(is.matrix(emb_a), is.matrix(emb_b))
  projection    <- match.arg(projection)     # "pca" or "umap"
  ellipse_group <- match.arg(ellipse_group)  # "none", "a", or "b"
  set.seed(seed)                             # reproducible sampling / projection
  
  # ---- 1) Align dimensions first ----------------------------------------------
  # Ensure both groups have the same number of columns (embedding dims).
  k <- min(ncol(emb_a), ncol(emb_b))
  emb_a <- emb_a[, seq_len(k), drop = FALSE]
  emb_b <- emb_b[, seq_len(k), drop = FALSE]
  
  print("step 1 complete")
  
  # ---- 2) Start from full sets for plotting -----------------------------------
  # Keep pristine copies (emb_a/emb_b) and create working copies (emb_a_plot/emb_b_plot)
  emb_a_plot <- emb_a
  emb_b_plot <- emb_b
  
  print("step 2 complete")
  
  # ---- 3) Optional: balance sizes (equalize row counts across groups) ---------
  # Prevent the larger group from visually dominating the scatter by randomly
  # subsampling the larger one down to the size of the smaller.
  if (isTRUE(balance_for_plot)) {
    n_bal <- min(nrow(emb_a_plot), nrow(emb_b_plot))
    emb_a_plot <- emb_a_plot[sample(seq_len(nrow(emb_a_plot)), n_bal), , drop = FALSE]
    emb_b_plot <- emb_b_plot[sample(seq_len(nrow(emb_b_plot)), n_bal), , drop = FALSE]
  }
  
  print("step 3 complete")
  
  # ---- 4) Optional: downsample to a total `sample_size` (split evenly) --------
  # If you just want a lighter plot, cap total points; each group gets half.
  if (!is.null(sample_size) && is.finite(sample_size) && sample_size > 0) {
    half <- max(1L, floor(sample_size / 2))
    emb_a_plot <- emb_a_plot[sample(seq_len(nrow(emb_a_plot)), min(half, nrow(emb_a_plot))), , drop = FALSE]
    emb_b_plot <- emb_b_plot[sample(seq_len(nrow(emb_b_plot)), min(half, nrow(emb_b_plot))), , drop = FALSE]
  }
  
  print("step 4 complete")
  
  # ---- 5) Combine for projection/plotting -------------------------------------
  X <- rbind(emb_a_plot, emb_b_plot)
  grp <- factor(c(rep(labels[1], nrow(emb_a_plot)),
                  rep(labels[2], nrow(emb_b_plot))),
                levels = labels)
  
  print("step 5 complete")
  
  # ---- 6) 2D projection --------------------------------------------------------
  # Choose fast PCA (irlba) or UMAP; both return two coordinates per row of X.
  if (projection == "pca") {
    # memory-efficient PCA for large matrices
    if (!requireNamespace("irlba", quietly = TRUE))
      stop("projection='pca' requires the 'irlba' package.")
    pc <- irlba::prcomp_irlba(as.matrix(X), n = 2, center = TRUE, scale. = FALSE)
    pc2 <- as.data.frame(pc$x)
    names(pc2) <- c("PC1", "PC2")
  } else {
    if (!requireNamespace("uwot", quietly = TRUE))
      stop("projection='umap' requires the 'uwot' package.")
    coords <- uwot::umap(X, n_neighbors = 15, min_dist = 0.1, n_components = 2, ret_model = FALSE)
    pc2 <- data.frame(PC1 = coords[,1], PC2 = coords[,2])
  }
  pc2$group <- grp
  
  # 3) 2D projection
  if (projection == "pca") {
    proj <- prcomp(X, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
  } else {
    if (!requireNamespace("uwot", quietly = TRUE))
      stop("UMAP selected but package 'uwot' not installed.")
    set.seed(seed)
    proj <- uwot::umap(X, n_neighbors = 15, min_dist = 0.1, n_components = 2)
  }
  pc <- prcomp_irlba(as.matrix(X), n = 2, center = TRUE, scale. = FALSE)
  pc2 <- as.data.frame(pc$x)
  names(pc2) <- c("PC1", "PC2")
  pc2$group <- grp
  
  print("step 6 complete")
  
  # Step 7
  cent <- pc2 |>
    group_by(group) |>
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  
  seg <- data.frame(
    x    = cent$PC1[cent$group == labels[1]],
    y    = cent$PC2[cent$group == labels[1]],
    xend = cent$PC1[cent$group == labels[2]],
    yend = cent$PC2[cent$group == labels[2]]
  )
  
  print("step 7 complete")
  
  # 8) Optional bootstrap CI ellipse around chosen centroid
  boot_df <- NULL
  if (ellipse_group != "none") {
    target_lab <- if (ellipse_group == "a") labels[1] else labels[2]
    dat_tgt <- dplyr::filter(pc2, group == target_lab)[, c("PC1","PC2")]
    set.seed(seed)
    boot_cent <- replicate(
      B_boot,
      colMeans(dat_tgt[sample(nrow(dat_tgt), replace = TRUE), , drop = FALSE]),
      simplify = TRUE
    )
    boot_df <- data.frame(PC1 = boot_cent[1, ], PC2 = boot_cent[2, ])
  }
  
  print("step 8 complete")
  
  # 9) Optional zoom around centroids
  xlim <- ylim <- NULL
  if (zoom_to_centroids) {
    xr <- range(cent$PC1); yr <- range(cent$PC2)
    pad_x <- diff(xr); if (pad_x == 0) pad_x <- 1
    pad_y <- diff(yr); if (pad_y == 0) pad_y <- 1
    xlim <- xr + c(-0.8, 0.8) * pad_x
    ylim <- yr + c(-0.8, 0.8) * pad_y
  }
  
  print("step 9 complete")
  
  # 7) Build figure
  g <- ggplot(pc2, aes(PC1, PC2)) +
    geom_point(aes(shape = group), alpha = point_alpha, size = 0.7) +
    geom_point(data = cent, aes(shape = group), size = 4, fill = "white", stroke = 1.2) +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, linetype = "dashed") +
    scale_shape_manual(values = c(21, 24), name = "Group") +
    labs(title = stats_title, subtitle = stats_subtitle, x = "PC1", y = "PC2") +
    theme_minimal(base_family = base_family, base_size = base_size) +
    theme(legend.position = "top")
  
  print("step 10 complete")
  
  # Bootstrap CI ellipse layer (if any)
  if (!is.null(boot_df)) {
    # draw a clean 95% ellipse around bootstrap centroids
    mu  <- colMeans(boot_df); S <- cov(boot_df)
    eig <- eigen(S); r95 <- sqrt(qchisq(ellipse_level, df = 2))
    a <- r95 * sqrt(eig$values[1]); b <- r95 * sqrt(eig$values[2])
    ang <- atan2(eig$vectors[2,1], eig$vectors[1,1])
    g <- g + ggforce::geom_ellipse(aes(x0 = mu[1], y0 = mu[2], a = a, b = b, angle = ang),
                                   colour = "black", fill = NA, linewidth = 1.2)
  }
  
  if (!is.null(xlim) && !is.null(ylim)) {
    g <- g + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  }
  
  print("step 11 complete")
  
  # Optional stats annotation (if provided)
  if (!is.null(stats)) {
    fmt <- function(x) if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)
    lab <- c(
      if (!is.null(stats$diff_mean))        paste0("Own–Other diff = ", fmt(stats$diff_mean)),
      if (!is.null(stats$p_diff))           paste0("Perm p (diff) = ", fmt(stats$p_diff)),
      if (!is.null(stats$cent_cos))         paste0("Centroid cosine = ", fmt(stats$cent_cos)),
      if (!is.null(stats$p_centroid))       paste0("Perm p (centroid) = ", fmt(stats$p_centroid))
    )
    g <- g + annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                      label = paste(lab, collapse = "\n"))
  }
  
  print("step 12 complete")
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = g, width = width, height = height, dpi = dpi)
  }
  
  print("step 13 complete")
  
  return(invisible(g))
}
