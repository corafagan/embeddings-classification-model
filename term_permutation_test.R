# ===============================================================
# Script Name:     term_permutation_test.R
# Project:         Robust Within-Cluster Similarity (Micro-Cohesion)
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This script tests whether pro- and anti-democratic dictionary
#   terms form more internally cohesive clusters in embedding space
#   than expected under random label assignments. It computes:
#     - Mean within-cluster cosine similarity (weighted by pairs
#       or by clusters)
#     - Permutation-based null distribution and p-values
#     - Centroid cosine similarity between pro and anti groups
#   Outputs can be summarized as LaTeX tables or visualized with
#   null-density plots.
#
# Expected Inputs:
#   - rcv_embedding_matrix:
#       • Numeric matrix of embeddings (rows = terms, cols = dims).
#       • Must have rownames set to tokens/terms.
#   - rcv_sentiment_dict:
#       • List or quanteda::dictionary2 object containing dictionary
#         patterns (globs).
#       • By default assumes anti-democratic at index [[1]] and
#         pro-democratic at index [[2]].
#
# Key Parameters:
#   - pro_index / anti_index:
#       • Indices (or names) of dictionary entries corresponding to
#         pro- and anti-democratic patterns.
#   - match_mode:
#       • How glob patterns are applied:
#           "substring" (default), "anchored", or "word".
#   - ignore_case (default = TRUE):
#       • Case-insensitive matching of terms to dictionary patterns.
#   - n_perm (default = 1000):
#       • Number of permutations for null distribution.
#   - weight:
#       • "pairs" (each term pair weighted equally) or "equal"
#         (each cluster weighted equally).
#   - side:
#       • "greater" (test if observed > null mean),
#         "less", or "two.sided".
#   - seed (default = 123):
#       • RNG seed for reproducibility.
#
# Outputs (list of class "dict_cohesion_result"):
#   - counts: number of pro vs. anti terms matched
#   - terms: vector of matched terms
#   - observed_mean_within: observed mean within-cluster cosine similarity
#   - null_mean_within: vector of permuted mean-within values
#   - p_value: permutation p-value (according to `side`)
#   - centroid_cosine_pro_vs_anti: cosine similarity between pro vs anti centroids
#   - side / weight: test configuration used
#
# Visualization:
#   - plot_null_density(result): plots null distribution of within-cluster
#     similarity with observed value and p-value annotated.
#   - LaTeX table creation via kableExtra summarizing cohesion, separation,
#     and significance stars.
#
# Notes:
#   - If one group has fewer than 2 matched terms, within-group means
#     cannot be computed for that cluster.
#   - Centroid cosine provides a complementary measure of group separation.
#   - Significance codes: ** < .01, * < .05.
# ===============================================================

# load packages
library(kableExtra)
library(quanteda)   
library(coop)

# ---- Helpers ----

# Safe glob -> regex converter
glob_to_regex <- function(glob, mode = c("substring", "anchored", "word")) {
  mode <- match.arg(mode)
  # Escape regex metacharacters except * and ?
  esc <- gsub("([.\\^$+(){}|\\[\\]\\\\])", "\\\\\\1", glob, perl = TRUE)
  esc <- gsub("\\?", ".", esc, perl = TRUE)
  esc <- gsub("\\*", ".*", esc, perl = TRUE)
  if (mode == "anchored") {
    paste0("^", esc, "$")
  } else if (mode == "word") {
    # Word boundaries on both ends when appropriate (letters/digits/underscore)
    paste0("\\b", esc, "\\b")
  } else {
    # substring (anywhere)
    esc
  }
}

# Match terms to a set of glob patterns
match_terms_by_globs <- function(terms, patterns, ignore_case = TRUE,
                                 mode = c("substring", "anchored", "word")) {
  mode <- match.arg(mode)
  hits <- character(0)
  for (p in patterns) {
    rx <- glob_to_regex(p, mode = mode)
    m <- grep(rx, terms, value = TRUE, ignore.case = ignore_case, perl = TRUE)
    if (length(m)) hits <- union(hits, m)
  }
  hits
}

# Compute full cosine similarity between terms (rows) once
cosine_matrix_from_embeddings <- function(emb_row_mat) {
  # emb_row_mat: rows = terms, cols = embedding dims
  # coop::cosine expects columns as observations; so transpose.
  coop::cosine(t(emb_row_mat))
}

# Mean within-cluster similarity from a precomputed S
# weight = "pairs" (current behavior) or "equal" (each cluster equally weighted)
mean_within_from_S <- function(S, labels, weight = c("pairs", "equal"), verbose = FALSE) {
  weight <- match.arg(weight)
  labs <- unique(labels)
  cluster_means <- numeric(0)
  cluster_weights <- numeric(0)
  
  for (cl in labs) {
    idx <- which(labels == cl)
    if (length(idx) > 1) {
      m <- S[idx, idx, drop = FALSE]
      v <- m[upper.tri(m)]
      mu <- mean(v)
      if (verbose) {
        cat("Cluster:", cl, " n =", length(idx), " mean within =", sprintf("%.6f", mu), "\n")
      }
      cluster_means <- c(cluster_means, mu)
      if (weight == "pairs") {
        cluster_weights <- c(cluster_weights, length(idx) * (length(idx) - 1) / 2)
      } else {
        cluster_weights <- c(cluster_weights, 1)
      }
    } else if (verbose) {
      cat("Cluster:", cl, " has <2 terms; skipping.\n")
    }
  }
  
  if (!length(cluster_means)) return(NA_real_)
  weighted.mean(cluster_means, w = cluster_weights)
}

# Permutation test around the mean-within statistic
# side = "greater" (cohesion higher than null), "less", or "two.sided"
perm_test_within_similarity <- function(S, labels, n_perm = 1000, seed = NULL,
                                        weight = c("pairs", "equal"),
                                        side = c("greater", "less", "two.sided"),
                                        verbose = FALSE) {
  weight <- match.arg(weight)
  side <- match.arg(side)
  if (!is.null(seed)) set.seed(seed)
  
  observed <- mean_within_from_S(S, labels, weight = weight, verbose = verbose)
  
  null_vals <- replicate(n_perm, {
    mean_within_from_S(S, sample(labels), weight = weight, verbose = FALSE)
  })
  
  if (all(is.na(null_vals)) || is.na(observed)) {
    return(list(observed = observed, null = null_vals, p_value = NA_real_))
  }
  
  if (side == "greater") {
    p <- mean(null_vals >= observed, na.rm = TRUE)
  } else if (side == "less") {
    p <- mean(null_vals <= observed, na.rm = TRUE)
  } else {
    # two-sided: double the one-sided tail wrt null mean
    mu0 <- mean(null_vals, na.rm = TRUE)
    p <- mean(abs(null_vals - mu0) >= abs(observed - mu0), na.rm = TRUE)
  }
  
  list(observed = observed, null = null_vals, p_value = p)
}

# Centroid cosine between two label groups (optional diagnostic)
centroid_cosine <- function(emb_row_mat, labels, a, b) {
  A <- emb_row_mat[labels == a, , drop = FALSE]
  B <- emb_row_mat[labels == b, , drop = FALSE]
  if (nrow(A) < 1 || nrow(B) < 1) return(NA_real_)
  ca <- colMeans(A)
  cb <- colMeans(B)
  sum(ca * cb) / (sqrt(sum(ca^2)) * sqrt(sum(cb^2)))
}

# ---- Pipeline ----
# Inputs expected:
# - rcv_embedding_matrix: rows named by terms, columns are embedding dims
# - rcv_sentiment_dict: list with anti at [[1]] and pro at [[2]] (adjust if different)

run_dictionary_cohesion_test <- function(rcv_embedding_matrix,
                                         rcv_sentiment_dict,
                                         pro_index = 2,
                                         anti_index = 1,
                                         match_mode = c("substring", "anchored", "word"),
                                         ignore_case = TRUE,
                                         n_perm = 1000,
                                         weight = c("pairs", "equal"),
                                         side = c("greater", "less", "two.sided"),
                                         seed = 123,
                                         verbose = TRUE) {
  match_mode <- match.arg(match_mode)
  weight <- match.arg(weight)
  side <- match.arg(side)
  
  if (is.null(rownames(rcv_embedding_matrix))) {
    stop("Embedding matrix must have rownames (terms).")
  }
  terms <- rownames(rcv_embedding_matrix)
  
  pro_patterns  <- rcv_sentiment_dict[[pro_index]]
  anti_patterns <- rcv_sentiment_dict[[anti_index]]
  
  matched_pro  <- match_terms_by_globs(terms, pro_patterns,  ignore_case = ignore_case, mode = match_mode)
  matched_anti <- match_terms_by_globs(terms, anti_patterns, ignore_case = ignore_case, mode = match_mode)
  
  if (verbose) {
    cat("Matched pro terms:", length(matched_pro), "\n")
    cat("Matched anti terms:", length(matched_anti), "\n")
  }
  
  # Prepare matrices and labels
  terms_all <- c(matched_pro, matched_anti)
  if (length(terms_all) < 2) stop("Fewer than 2 matched terms overall; cannot compute similarities.")
  
  labels <- c(rep("pro", length(matched_pro)), rep("anti", length(matched_anti)))
  if (sum(labels == "pro") < 2 || sum(labels == "anti") < 2) {
    warning("One or both groups have < 2 terms. Within-group mean will be NA for that group; test may be uninformative.")
  }
  
  emb_mat <- rcv_embedding_matrix[terms_all, , drop = FALSE]
  
  # Precompute S once
  if (verbose) cat("Computing cosine similarity matrix once...\n")
  S <- cosine_matrix_from_embeddings(emb_mat)
  
  # Observed + permutation test
  if (verbose) cat("Running permutation test...\n")
  pt <- perm_test_within_similarity(S, labels,
                                    n_perm = n_perm, seed = seed,
                                    weight = weight, side = side,
                                    verbose = FALSE)
  
  # Optional diagnostics
  cos_centroids <- centroid_cosine(emb_mat, labels, "pro", "anti")
  
  # Tidy outputs
  out <- list(
    counts = table(labels),
    labels = labels,
    terms  = terms_all,
    observed_mean_within = pt$observed,
    null_mean_within = pt$null,
    p_value = pt$p_value,
    side = side,
    weight = weight,
    centroid_cosine_pro_vs_anti = cos_centroids
  )
  
  if (verbose) {
    cat("Observed mean within-cluster similarity:", sprintf("%.6f", pt$observed), "\n")
    cat("Centroid cosine (pro vs anti):", ifelse(is.na(cos_centroids), "NA", sprintf("%.6f", cos_centroids)), "\n")
    cat("Permutation p-value (", side, "): ", sprintf("%.4g", pt$p_value), "\n", sep = "")
  }
  
  class(out) <- c("dict_cohesion_result", class(out))
  out
}

# ---- Plotting ----

plot_null_density <- function(result) {
  if (is.null(result$null_mean_within) || all(is.na(result$null_mean_within))) {
    stop("Null distribution is empty or NA; cannot plot.")
  }
  df <- data.frame(similarity = result$null_mean_within)
  rng <- range(df$similarity, na.rm = TRUE)
  x_off <- diff(rng) * 0.02
  
  ggplot(df, aes(x = similarity)) +
    geom_density(size = 1.1) +
    geom_vline(xintercept = result$observed_mean_within, linetype = "dashed", size = 1) +
    annotate("text",
             x = result$observed_mean_within - x_off,
             y = max(stats::density(df$similarity, na.rm = TRUE)$y) * 0.9,
             label = paste0("Observed = ", round(result$observed_mean_within, 5),
                            "\nP-value = ", format(result$p_value, digits = 3)),
             hjust = 1, size = 4.5) +
    labs(x = "Mean within-cluster cosine similarity", y = "Density") +
    theme_minimal(base_size = 14)
}

# =========================
# Example usage (adjust indices if your dictionary order is different):
# =========================
 res <- run_dictionary_cohesion_test(
   rcv_embedding_matrix = rcv_embedding_matrix,
   rcv_sentiment_dict = rcv_sentiment_dict,
   pro_index = 2, anti_index = 1,          # matches your current structure
   match_mode = "substring",                # or "anchored" / "word"
   ignore_case = TRUE,
   n_perm = 1000,
   weight = "pairs",                        # or "equal"
   side = "greater",                        # or "less" / "two.sided"
   seed = 123,
   verbose = TRUE
 )
print(res$counts)

# simulated_cossim_plot <- plot_null_density(res)

# Extract values from res
n_pro   <- as.integer(res$counts["pro"])
n_anti  <- as.integer(res$counts["anti"])
null_mu <- mean(res$null_mean_within, na.rm = TRUE)

# If you want CI for the null distribution:
null_ci <- quantile(res$null_mean_within, c(.025, .975), na.rm = TRUE)

# Make table
df_overall <- tibble(
  `Group size`         = sprintf("anti = %d; pro = %d", n_anti, n_pro),
  `Obs. within`        = round(res$observed_mean_within, 5),
  `Null mean`          = round(null_mu, 5),
  `Null 95% CI low`    = round(null_ci[1], 5),
  `Null 95% CI high`   = round(null_ci[2], 5),
  `Centroid (pro–anti)`= round(res$centroid_cosine_pro_vs_anti, 3),
  `p`                  = ifelse(res$p_value < 0.001, "< 0.001",
                                sprintf("%.3f", res$p_value)),
  `Sig.`               = case_when(
    res$p_value < 0.01 ~ "**",
    res$p_value < 0.05 ~ "*",
    TRUE               ~ ""
  )
)

# Render as LaTeX table
kbl(
  df_overall,
  format = "latex",
  booktabs = TRUE,
  caption = "Micro-level: Overall cohesion and separation of pro- and anti-democracy dictionary terms in embedding space."
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  footnote(
    general = "Permutation test on mean within-cluster cosine similarity. Null based on label shuffles.",
    threeparttable = TRUE
  )
