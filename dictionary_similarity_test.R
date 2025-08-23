# ===============================================================
# Script Name:     dictionary_similarity_test.R
# Project:         Micro-level Cohesion Analysis
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This script provides utility functions to measure *within-label
#   cohesion* of embeddings using cosine similarity. It computes:
#     - Mean within-cluster cosine similarity for each label
#     - Permutation-based null distributions (holding class sizes fixed)
#     - One-sided p-values testing whether observed cohesion exceeds null
#     - 95% null quantile intervals for context
#   The output can be rendered as a LaTeX table summarizing results
#   for pro- and anti-democracy dictionary terms.
#
# Expected Inputs:
#   - emb_mat:
#       • A numeric matrix of embeddings (rows = terms, cols = dimensions).
#       • Must be dense; no missing values allowed.
#   - labels:
#       • A factor or character vector of group labels.
#       • Length must equal nrow(emb_mat).
#       • Each entry corresponds to the row in emb_mat.
#       • Example: c("pro","pro","anti","anti",...)
#   - n_perm (default 1000):
#       • Integer, number of random permutations for null distributions.
#   - seed (default 123):
#       • Optional integer seed for reproducibility of permutations.
#
# Outputs (returned as a list by per_label_cohesion):
#   - counts: table of group sizes
#   - observed_within: observed mean within-label cosine similarities
#   - null_within: permutation null matrix (rows = permutations, cols = labels)
#   - p_greater: one-sided p-values (observed > null)
#   - null_ci95: 2.5% and 97.5% quantiles of null distribution per label
#
# Example Workflow:
#   1. Subset embeddings to matched pro/anti terms
#   2. Run per_label_cohesion(emb_mat, labels, n_perm = 1000)
#   3. Extract observed vs. null statistics
#   4. Build a publication-ready LaTeX table with kableExtra
#
# Notes:
#   - Uses coop::cosine for efficient similarity computation.
#   - P-values are one-sided: test whether cohesion is stronger
#     than expected under random labeling with fixed class sizes.
#   - The table builder adds stars for conventional significance
#     thresholds (*** < .001, ** < .01, * < .05, . < .1).
# ===============================================================

# ---- Utilities ----
# creates micro-separation table

# Mean within-cluster cosine for a given set of indices
.mean_within_idx <- function(S, idx) {
  if (length(idx) < 2) return(NA_real_)
  M <- S[idx, idx, drop = FALSE]
  mean(M[upper.tri(M)])
}

# Compute observed per-label within means + permutation nulls + p-values
per_label_cohesion <- function(emb_mat, labels, n_perm = 1000, seed = 123) {
  stopifnot(nrow(emb_mat) == length(labels))
  labs <- as.character(labels)
  classes <- unique(labs)
  
  # Precompute cosine similarity once
  S <- coop::cosine(t(emb_mat))
  
  # Observed within means per label
  obs <- setNames(
    vapply(classes, function(cl) .mean_within_idx(S, which(labs == cl)), numeric(1)),
    classes
  )
  
  # Permutation nulls per label (keep class sizes)
  if (!is.null(seed)) set.seed(seed)
  null_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(classes))
  colnames(null_mat) <- classes
  
  for (b in seq_len(n_perm)) {
    shuf <- sample(labs)
    for (j in seq_along(classes)) {
      cl <- classes[j]
      null_mat[b, j] <- .mean_within_idx(S, which(shuf == cl))
    }
  }
  
  # One-sided p (greater): are groups more cohesive than random?
  pvals <- setNames(
    vapply(seq_along(classes), function(j) {
      mean(null_mat[, j] >= obs[j], na.rm = TRUE)
    }, numeric(1)),
    classes
  )
  
  # Simple 95% null quantiles for context
  qlo <- apply(null_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  qhi <- apply(null_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  list(
    counts = table(labs),
    observed_within = obs,
    null_within = null_mat,            # rows = permutations, cols = labels
    p_greater = pvals,
    null_ci95 = rbind(`2.5%` = qlo, `97.5%` = qhi)
  )
}

# ---- Example usage (plug your objects in) ----
# Suppose you already have:
   terms_all <- c(matched_pro_terms, matched_anti_terms)
   emb_mat   <- rcv_embedding_matrix[terms_all, , drop = FALSE]
   labels    <- c(rep("pro", length(matched_pro_terms)),
                  rep("anti", length(matched_anti_terms)))

# Run:
 res_per <- per_label_cohesion(emb_mat, labels, n_perm = 1000, seed = 123)
 res_per$counts
 res_per$observed_within
 res_per$p_greater
 res_per$null_ci95
 
 # Keep a consistent label order
 lbls <- c("anti", "pro")
 lbls <- lbls[lbls %in% names(res_per$observed_within)]
 
 # Extract pieces in the same order
 N_vals       <- as.integer(res_per$counts[lbls])
 obs_within   <- as.numeric(res_per$observed_within[lbls])
 null_mean    <- colMeans(res_per$null_within[, lbls, drop = FALSE], na.rm = TRUE)
 ci_low       <- as.numeric(res_per$null_ci95["2.5%",  lbls, drop = TRUE])
 ci_high      <- as.numeric(res_per$null_ci95["97.5%", lbls, drop = TRUE])
 p_vals       <- as.numeric(res_per$p_greater[lbls])
 
 # Significance stars
 sig_stars <- cut(
   p_vals,
   breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
   labels = c("***", "**", "*", ".", ""),
   right = FALSE
 )
 
 # Build the data frame
 tbl_df <- data.frame(
   Label = lbls,
   N = N_vals,
   `Obs. within`       = round(obs_within, 3),
   `Null mean`         = round(null_mean, 3),
   `Null 95% CI low`   = round(ci_low, 3),
   `Null 95% CI high`  = round(ci_high, 3),
   `p`                 = ifelse(p_vals < 0.001, "< 0.001", sprintf("%.3f", p_vals)),
   `Sig.`              = as.character(sig_stars),
   check.names = FALSE
 )
 
 str(tbl_df)
 
 # Render LaTeX table
 kbl(tbl_df,
     format = "latex",
     booktabs = TRUE,
     linesep = "",
     caption = "Micro-level: Cohesion by label for pro- and anti-democracy dictionary terms. Pro-democracy terms exhibit significantly greater within-cluster similarity than expected under random labeling, whereas anti-democracy terms do not.",
     label = "tab:micro-separation"
 ) |>
   kable_styling(latex_options = c("hold_position", "scale_down")) |>
   add_header_above(c(" " = 2, "Observed" = 1, "Null (Permutation)" = 3, "Inference" = 2)) |>
   footnote(
     general = "One-sided p-values test whether within-label cohesion exceeds the permutation null; class sizes held fixed.",
     threeparttable = TRUE
   )
