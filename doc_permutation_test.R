# ===============================================================
# Script Name:     doc_permutation_test.R
# Project:         Doc-Level Embedding Inference
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This script runs a permutation-based inference pipeline for
#   document-level embeddings comparing treatment vs. control groups.
#   It computes:
#     - Similarity of each document to its own vs. other group centroid
#     - Difference in within- vs. cross-group similarities
#     - Cosine similarity between treatment and control centroids
#   The analysis includes parametric OLS summaries, observed statistics,
#   permutation null distributions, and two-sided p-values.
#
# Expected Inputs:
#   - rcv_doc_embeddings:
#       • A data.frame or matrix of embeddings (rows = docs, cols = dims).
#       • Represents the treatment group (RCV).
#   - combined_control_embeddings:
#       • A data.frame or matrix of embeddings (rows = docs, cols = dims).
#       • Represents the control group.
#
# Parameters:
#   - n_perm (default 1000):
#       • Number of permutations for null distributions.
#   - seed (default 123):
#       • RNG seed for reproducibility.
#   - prog_every (default 500):
#       • Print progress message every N permutations.
#
# Outputs (returned as a list of class "embedding_perm_results"):
#   - n_perm: number of permutations used
#   - models:
#       • lm_diff_summary: OLS on sim_diff ~ group
#       • intra_model_summary: OLS on sim_intra ~ group
#   - observed:
#       • diff_means_simdiff: observed group difference in similarity difference
#       • centroid_cosine: observed cosine similarity between centroids
#   - permutation:
#       • diff_means_stats: vector of permuted diff-of-means stats
#       • centroid_cos_stats: vector of permuted centroid cosine stats
#   - p_values:
#       • diff_means_two_sided: two-sided permutation p-value for diff-of-means
#       • centroid_cos_two_sided: two-sided permutation p-value for centroid cosine
#   - data:
#       • combined_embeddings: tidy dataframe with doc_id, group, and similarity measures
#
# Notes:
#   - All embeddings are L2-normalized row-wise so that dot products
#     correspond to cosine similarities.
#   - The "sim_diff" metric = similarity to own centroid − similarity to other centroid.
#   - Permutation procedure resamples group labels while holding sample sizes fixed.
# ===============================================================

run_embedding_permutation_pipeline <- function(rcv_doc_embeddings,
                                               combined_control_embeddings,
                                               n_perm = 1000,
                                               seed   = 123,
                                               prog_every = 500) {
  stopifnot(nrow(rcv_doc_embeddings) > 0,
            nrow(combined_control_embeddings) > 0)
  
  # --- helpers ---
  cosine_sim <- function(a, b) {
    denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (denom == 0) return(NA_real_)
    sum(a * b) / denom
  }
  
  l2_normalize_rows <- function(M) {
    M <- as.matrix(M)
    nrms <- sqrt(rowSums(M^2))
    nz   <- nrms > 0
    M[nz, ] <- M[nz, , drop = FALSE] / nrms[nz]
    M
  }
  
  # --- 0) Coerce to matrices; add ids and group labels ---
  rcv_dense <- as.matrix(rcv_doc_embeddings)
  rownames(rcv_dense) <- as.character(seq_len(nrow(rcv_dense)))
  
  control_dense <- as.matrix(combined_control_embeddings)
  rownames(control_dense) <- as.character(seq_len(nrow(control_dense)))
  
  rcv_df <- as.data.frame(rcv_dense)
  rcv_df$doc_id <- rownames(rcv_dense)
  rcv_df$group  <- "treatment"
  
  control_df <- as.data.frame(control_dense)
  control_df$doc_id <- rownames(control_dense)
  control_df$group  <- "control"
  
  combined_embeddings <- rbind(rcv_df, control_df)
  combined_embeddings$group <- factor(combined_embeddings$group,
                                      levels = c("control", "treatment"))
  
  # Identify numeric embedding columns (robust to column names)
  embed_cols <- names(combined_embeddings)[vapply(combined_embeddings, is.numeric, TRUE)]
  embed_cols <- setdiff(embed_cols, c("doc_id"))  # safety
  
  # --- 1) Normalize (row-wise L2) so dot product = cosine ---
  embedding_matrix <- l2_normalize_rows(combined_embeddings[, embed_cols, drop = FALSE])
  
  # Split index
  is_treat <- combined_embeddings$group == "treatment"
  
  # --- 2) Group centroids (same space, normalized) ---
  treatment_mat <- embedding_matrix[is_treat, , drop = FALSE]
  control_mat   <- embedding_matrix[!is_treat, , drop = FALSE]
  
  centroid_treat <- colMeans(treatment_mat)
  centroid_ctrl  <- colMeans(control_mat)
  
  # Normalize centroids so dot is cosine similarity
  centroid_treat <- centroid_treat / sqrt(sum(centroid_treat^2))
  centroid_ctrl  <- centroid_ctrl  / sqrt(sum(centroid_ctrl^2))
  
  # --- 3) Similarity to each centroid for every doc ---
  sim_to_treat   <- as.vector(embedding_matrix %*% centroid_treat)
  sim_to_control <- as.vector(embedding_matrix %*% centroid_ctrl)
  
  # Intra vs inter; sim_diff is "own centroid minus other centroid"
  sim_intra <- ifelse(is_treat, sim_to_treat,   sim_to_control)
  sim_inter <- ifelse(is_treat, sim_to_control, sim_to_treat)
  sim_diff  <- sim_intra - sim_inter
  
  # Attach to data
  combined_embeddings$sim_intra <- sim_intra
  combined_embeddings$sim_inter <- sim_inter
  combined_embeddings$sim_diff  <- sim_diff
  
  # --- 4) Parametric OLS (optional but nice to report) ---
  lm_diff     <- lm(sim_diff  ~ group, data = combined_embeddings)
  intra_model <- lm(sim_intra ~ group, data = combined_embeddings)
  
  # Observed statistics
  obs_diff_means   <- mean(sim_diff[is_treat]) - mean(sim_diff[!is_treat])
  obs_centroid_cos <- cosine_sim(centroid_treat, centroid_ctrl)
  
  # --- 5) Permutation tests (both stats in one pass) ---
  run_joint_permutations <- function(n_perm,
                                     groups,
                                     embedding_matrix,
                                     sim_diff,
                                     prog_every = 500,
                                     seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    perm_stats_diff     <- numeric(n_perm)
    perm_stats_centroid <- numeric(n_perm)
    
    for (i in seq_len(n_perm)) {
      if (i %% prog_every == 0 || i == 1 || i == n_perm) {
        message("Permutation ", i, " of ", n_perm, " (diff means + centroid cosine)")
      }
      
      # Permute labels once, reuse for both stats
      perm_labels <- sample(groups)
      
      # ---- diff of means on sim_diff ----
      perm_stats_diff[i] <-
        mean(sim_diff[perm_labels == "treatment"]) -
        mean(sim_diff[perm_labels == "control"])
      
      # ---- centroid cosine ----
      t_mat <- embedding_matrix[perm_labels == "treatment", , drop = FALSE]
      c_mat <- embedding_matrix[perm_labels == "control",   , drop = FALSE]
      
      if (nrow(t_mat) == 0L || nrow(c_mat) == 0L) {
        perm_stats_centroid[i] <- NA_real_
      } else {
        ct <- colMeans(t_mat); cc <- colMeans(c_mat)
        ct <- ct / sqrt(sum(ct^2)); cc <- cc / sqrt(sum(cc^2))
        perm_stats_centroid[i] <- sum(ct * cc)  # cosine similarity
      }
    }
    
    list(diff_means_stats = perm_stats_diff,
         centroid_cos_stats = perm_stats_centroid)
  }
  
  perms <- run_joint_permutations(
    n_perm           = n_perm,
    groups           = combined_embeddings$group,
    embedding_matrix = embedding_matrix,
    sim_diff         = sim_diff,
    prog_every       = prog_every,
    seed             = seed
  )
  
  perm_stats_diff     <- perms$diff_means_stats
  perm_stats_centroid <- perms$centroid_cos_stats
  
  # Two-sided p for diff-of-means
  p_diff <- mean(abs(perm_stats_diff) >= abs(obs_diff_means))
  
  # Two-sided p for centroid cosine: center around permuted mean
  perm_mean_cos <- mean(perm_stats_centroid, na.rm = TRUE)
  obs_dev       <- abs(obs_centroid_cos - perm_mean_cos)
  p_centroid    <- mean(abs(perm_stats_centroid - perm_mean_cos) >= obs_dev, na.rm = TRUE)
  
  # --- 6) Package results ---
  out <- list(
    n_perm = n_perm,
    models = list(
      lm_diff_summary     = summary(lm_diff),
      intra_model_summary = summary(intra_model)
    ),
    observed = list(
      diff_means_simdiff = obs_diff_means,
      centroid_cosine    = obs_centroid_cos
    ),
    permutation = list(
      diff_means_stats   = perm_stats_diff,
      centroid_cos_stats = perm_stats_centroid
    ),
    p_values = list(
      diff_means_two_sided   = p_diff,
      centroid_cos_two_sided = p_centroid
    ),
    data = list(
      combined_embeddings = combined_embeddings  # handy for plots
    )
  )
  
  class(out) <- c("embedding_perm_results", class(out))
  return(out)
}

# -------- Example usage --------
res <- run_embedding_permutation_pipeline(
  rcv_doc_embeddings,
  combined_control_embeddings,
  n_perm = 1000,
  prog_every = 100  # print more frequently
)
 
# Quick peek:
res$p_values
res$observed
res$models$lm_diff_summary
