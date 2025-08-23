
# ===============================================================
# Script Name:     rcv-vs-control-inference.R
# Project:         RCV vs Control Embedding Inference
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This script compares RCV (treatment) vs. control word embeddings
#   using sentiment dictionary centroids and tests whether the groups
#   differ on a pro-democratic vs. anti-democratic index. The analysis
#   includes cosine similarity computation, permutation testing, and
#   bootstrap confidence intervals, as well as effect size estimates.
#
# Expected Inputs:
#   - rcv_word_embeddings_df:
#       • Columns V1..Vk: numeric embedding dimensions (k ≥ 1)
#       • term: character token/word
#       • group: factor/character, must equal "treatment"
#   - control_word_embeddings_df:
#       • Columns V1..Km: numeric embedding dimensions (m ≥ 1)
#       • term: character token/word
#       • group: factor/character, must equal "control"
#   - rcv_sentiment_dict:
#       • quanteda::dictionary2 object with named keys:
#         - "anti_democratic": glob patterns for anti-democratic terms
#         - "pro_democratic":  glob patterns for pro-democratic terms
#
# Outputs:
#   - Observed difference (treatment - control) in democratic index
#   - Permutation p-value (two-sided)
#   - Bootstrap 95% CI for difference
#   - Effect size (Cohen’s d, Hedges’ g)
#   - Per-term scores (cosine similarities and index)
#
# Notes:
#   - Embedding dimensions between treatment and control are aligned
#     to the minimum shared k.
#   - Deduplication removes empty/duplicate terms.
#   - The “index” is defined as sim_pro - sim_anti; higher = more
#     pro-democratic orientation.
# ===============================================================
# ---------- 1) Align dimensions ----------
get_embed_cols <- function(df) grep("^V\\d+$", names(df), value = TRUE)
cols_rcv     <- get_embed_cols(rcv_word_embeddings_df)
cols_control <- get_embed_cols(control_word_embeddings_df)
k <- min(length(cols_rcv), length(cols_control))
stopifnot(k > 0)

keep_rcv     <- cols_rcv[seq_len(k)]
keep_control <- cols_control[seq_len(k)]

rcv <- rcv_word_embeddings_df[, c(keep_rcv, "term", "group"), drop = FALSE]
ctl <- control_word_embeddings_df[, c(keep_control, "term", "group"), drop = FALSE]

# ---------- 2) Basic cleaning ----------
dedup <- function(df) df[!is.na(df$term) & df$term != "" & !duplicated(df$term), , drop = FALSE]
rcv <- dedup(rcv)
ctl <- dedup(ctl)

# Combine
both <- rbind(rcv, ctl)
row.names(both) <- NULL
embed_cols <- get_embed_cols(both)
X <- as.matrix(both[, embed_cols, drop = FALSE])

# ---------- 3) Expand dictionary (glob → regex) and match ----------
dict_list <- quanteda::as.list(rcv_sentiment_dict)  # list with $anti_democratic, $pro_democratic; each is list(length=1) of character patterns

glob_or <- function(patterns) {
  # combine multiple globs into a single OR-regex
  rx <- vapply(patterns, utils::glob2rx, "", USE.NAMES = FALSE)
  paste0("(", paste(rx, collapse = ")|("), ")")
}
rx_anti <- glob_or(dict_list$anti_democratic[[1]])
rx_pro  <- glob_or(dict_list$pro_democratic[[1]])

terms_lower <- tolower(both$term)
m_anti <- grepl(rx_anti, terms_lower, ignore.case = TRUE)
m_pro  <- grepl(rx_pro,  terms_lower, ignore.case = TRUE)

if (!any(m_anti)) stop("No tokens matched anti_democratic patterns.")
if (!any(m_pro))  stop("No tokens matched pro_democratic patterns.")

# ---------- 4) Centroids ----------
centroid <- function(M) colMeans(M, na.rm = TRUE)
c_anti <- centroid(X[m_anti, , drop = FALSE])
c_pro  <- centroid(X[m_pro,  , drop = FALSE])

# ---------- 5) Cosine similarities & index ----------
row_norms <- sqrt(rowSums(X^2))
safe_div  <- function(num, den) { den[den == 0] <- NA_real_; num / den }

cosine_to <- function(c_vec) {
  c_norm <- sqrt(sum(c_vec^2))
  safe_div(as.numeric(X %*% c_vec), row_norms * c_norm)
}

both$sim_pro  <- cosine_to(c_pro)
both$sim_anti <- cosine_to(c_anti)
both$index    <- both$sim_pro - both$sim_anti  # higher => more "pro_democratic"

# ---------- 6) Inference: permutation test + bootstrap CI ----------
# Keep rows with finite index
ok <- is.finite(both$index)
both <- both[ok, , drop = FALSE]

mean_by_group <- tapply(both$index, both$group, mean, na.rm = TRUE)
diff_obs <- unname(mean_by_group["treatment"] - mean_by_group["control"])

set.seed(123)
B <- 10000L

# Permutation (two-sided)
n <- nrow(both)
n_treat <- sum(both$group == "treatment")
perm_null <- replicate(B, {
  idx <- sample.int(n)
  mt  <- mean(both$index[idx[1:n_treat]], na.rm = TRUE)
  mc  <- mean(both$index[idx[(n_treat + 1):n]], na.rm = TRUE)
  mt - mc
})
pval <- mean(abs(perm_null) >= abs(diff_obs))

# Bootstrap CI
boot_diff <- replicate(B, {
  mt <- mean(sample(both$index[both$group == "treatment"], replace = TRUE), na.rm = TRUE)
  mc <- mean(sample(both$index[both$group == "control"],    replace = TRUE), na.rm = TRUE)
  mt - mc
})
ci <- quantile(boot_diff, c(0.025, 0.975), na.rm = TRUE)

# ---------- 7) Outputs ----------
cat("Embedding dims used (k):", k, "\n")
cat(sprintf("Observed difference (treatment - control): %.4f\n", diff_obs))
cat(sprintf("Permutation p-value (two-sided): %.4g\n", pval))
cat(sprintf("Bootstrap 95%% CI: [%.4f, %.4f]\n", ci[1], ci[2]))

results <- list(
  k_used = k,
  group_means = mean_by_group,
  diff_obs = diff_obs,
  perm_p = pval,
  boot_ci = ci,
  scores = both[, c("term", "group", "sim_pro", "sim_anti", "index"), drop = FALSE]
)

# Example: inspect top tokens by index in each group (optional)
# head(results$scores[results$scores$group=="treatment"][order(-results$scores$index), ], 20)
# head(results$scores[results$scores$group=="control"][order(-results$scores$index), ], 20)

results

ix_t <- results$scores$index[results$scores$group == "treatment"]
ix_c <- results$scores$index[results$scores$group == "control"]

sd_pooled <- sqrt(((length(ix_t)-1)*var(ix_t) + (length(ix_c)-1)*var(ix_c)) /
                    (length(ix_t) + length(ix_c) - 2))
d <- (mean(ix_t) - mean(ix_c)) / sd_pooled

# Hedges' g small-sample correction:
J <- 1 - 3/(4*(length(ix_t)+length(ix_c))-9)
g <- J * d
c(d = d, g = g)

