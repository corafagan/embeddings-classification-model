# ===============================================================
# Script Name:     targeted-permutation.R
# Project:         Dictionary-Anchored Projection (RCV vs Control)
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This script compares treatment vs. control embeddings using a
#   dictionary-anchored projection axis. It computes per-term indices:
#       index = cos(term, pro-centroid) − cos(term, anti-centroid)
#   and evaluates whether the average index differs between groups.
#   The analysis includes permutation tests, bootstrap confidence
#   intervals, and effect size estimates, with results formatted as
#   a LaTeX table for reporting.
#
# Expected Inputs:
#   - combined_shared_embeddings:
#       • Data frame with columns:
#           term   : character token
#           V1..Vk : numeric embedding dimensions
#           group  : factor/character ("control" or "treatment")
#   - rcv_sentiment_dict:
#       • quanteda::dictionary2 object with entries:
#           $pro_democratic  : list of glob patterns for pro terms
#           $anti_democratic : list of glob patterns for anti terms
#
# Parameters:
#   - B (default = 10,000): number of permutations / bootstrap resamples
#   - Seed is set to 123 for reproducibility
#
# Outputs:
#   - Printed summary:
#       • Mean index difference (treatment − control)
#       • Permutation p-value (two-sided, Monte Carlo)
#       • Bootstrap 95% CI
#       • Effect sizes: Cohen’s d and Hedges’ g
#   - Data frame `df_results` with formatted statistics
#   - LaTeX table `latex_code` (via kableExtra) summarizing results
#
# Notes:
#   - Dictionary terms are pooled across groups to compute centroids.
#   - Positive index values indicate greater pro-democratic alignment.
#   - Permutation procedure holds group sizes fixed.
#   - If permutation p-value is 0, it is displayed as < 1/B.
# ===============================================================

# load packages
library(dplyr)
library(stringr)
library(kableExtra)

# --- assume you already created `combined_shared_embeddings` with term, V1:V10, group ---
X <- as.matrix(combined_shared_embeddings %>% dplyr::select(dplyr::starts_with("V")))

terms <- tolower(combined_shared_embeddings$term)
groups <- combined_shared_embeddings$group

# --- expand pro/anti patterns to regex ORs ---
dict_list <- quanteda::as.list(rcv_sentiment_dict)  # $anti_democratic, $pro_democratic
glob_or <- function(patterns) paste(vapply(patterns, utils::glob2rx, ""), collapse="|")
rx_pro  <- glob_or(dict_list$pro_democratic[[1]])
rx_anti <- glob_or(dict_list$anti_democratic[[1]])

m_pro  <- grepl(rx_pro,  terms, ignore.case = TRUE)
m_anti <- grepl(rx_anti, terms, ignore.case = TRUE)

stopifnot(any(m_pro), any(m_anti))

# --- (1) centroids for the axis, pooled across groups ---
c_pro  <- colMeans(X[m_pro, , drop = FALSE])
c_anti <- colMeans(X[m_anti, , drop = FALSE])

# --- helpers ---
row_norms <- sqrt(rowSums(X^2))
safe_cos <- function(mat, v){
  vn <- sqrt(sum(v^2))
  as.numeric(mat %*% v) / (row_norms * vn)
}

# --- (2) per-token index score ---
sim_pro  <- safe_cos(X, c_pro)
sim_anti <- safe_cos(X, c_anti)
index    <- sim_pro - sim_anti   # higher => more pro-democratic alignment

df_scores <- tibble(term = terms, group = groups,
                    sim_pro = sim_pro, sim_anti = sim_anti, index = index) %>%
  filter(is.finite(index))

# --- (3) group comparison: mean Index (RCV - control) ---
mean_by_group <- tapply(df_scores$index, df_scores$group, mean, na.rm = TRUE)
diff_obs <- unname(mean_by_group["treatment"] - mean_by_group["control"])

# permutation test (two-sided)
set.seed(123)
B <- 10000L
n  <- nrow(df_scores)
n_treat <- sum(df_scores$group == "treatment")
perm_null <- replicate(B, {
  idx <- sample.int(n)
  mt  <- mean(df_scores$index[idx[1:n_treat]], na.rm = TRUE)
  mc  <- mean(df_scores$index[idx[(n_treat + 1):n]], na.rm = TRUE)
  mt - mc
})
p_two_sided <- mean(abs(perm_null - mean(perm_null)) >= abs(diff_obs - mean(perm_null)))

# bootstrap CI for diff of means
boot_diff <- replicate(B, {
  mt <- mean(sample(df_scores$index[df_scores$group=="treatment"], replace=TRUE), na.rm=TRUE)
  mc <- mean(sample(df_scores$index[df_scores$group=="control"],    replace=TRUE), na.rm=TRUE)
  mt - mc
})
ci <- quantile(boot_diff, c(.025, .975), na.rm = TRUE)

# effect sizes
ix_t <- df_scores$index[df_scores$group=="treatment"]
ix_c <- df_scores$index[df_scores$group=="control"]
sd_pooled <- sqrt(((length(ix_t)-1)*var(ix_t) + (length(ix_c)-1)*var(ix_c)) /
                    (length(ix_t)+length(ix_c)-2))
d <- (mean(ix_t) - mean(ix_c)) / sd_pooled
J <- 1 - 3/(4*(length(ix_t)+length(ix_c)) - 9)
g <- J * d

cat(sprintf("Mean Index (treat - control): %.4f\n", diff_obs))
cat(sprintf("Permutation p (two-sided): %.4g\n", p_two_sided))
cat(sprintf("Bootstrap 95%% CI: [%.4f, %.4f]\n", ci[1], ci[2]))
cat(sprintf("Effect sizes: d = %.3f, g = %.3f\n", d, g))

perm_p_display <- if (p_two_sided == 0) {
  sprintf("< %.4f", 1 / B)
} else {
  sprintf("%.3f", p_two_sided)
}

df_results <- data.frame(
  check.names = FALSE,  # keep your header text as-is
  `Mean Index (Treat − Control)` = sprintf("%.4f", diff_obs),
  `95% CI Low`                = sprintf("%.4f", ci[1]),
  `95% CI High`               = sprintf("%.4f", ci[2]),
  `Permutation p`             = perm_p_display,
  `Cohen’s d`                 = sprintf("%.3f", d),
  `Hedges g`                 = sprintf("%.3f", g)
)

latex_code <- kbl(
  df_results,
  format  = "latex",
  booktabs = TRUE,
  caption = "Dictionary-anchored projection: difference in mean index (treatment minus control)",
  label   = "tab:dict_index",
  align   = "lrrrrr",
  linesep = ""
) |>
  kable_styling(latex_options = c("hold_position", "scale_down")) |>
  add_header_above(c(" " = 1, "Statistic" = 5)) |>
  footnote(
    general = "Index = cos(term, pro-centroid) − cos(term, anti-centroid). Positive values indicate greater similarity to the pro-democracy centroid. Permutation p reports two-sided Monte Carlo p-value with B = 10,000.",
    threeparttable = TRUE
  ) |>
  as.character()

cat(latex_code)
