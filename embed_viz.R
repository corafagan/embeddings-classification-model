# ===============================================================
# Script Name:     embed_viz.R
# Project:         Democratic Legitimacy Dictionary Projection
# Author:          Cora Fagan
# Date Created:    2025-08-23
#
# Description:
#   This script visualizes dictionary terms (pro- vs. anti-democratic)
#   in embedding space. It computes each term’s projection index:
#       index = cos(term, pro-centroid) − cos(term, anti-centroid)
#   and produces a lollipop plot of the top ±K terms ranked by |index|.
#   Positive values indicate stronger alignment with the pro-democratic
#   centroid; negative values indicate stronger alignment with the
#   anti-democratic centroid.
#
# Expected Inputs:
#   - control_word_embeddings_df:
#       • Data frame of embeddings (rows = terms, cols = V1..Vk).
#       • Must contain:
#           V1..Vk: numeric embedding dimensions
#           term:   character token
#           group:  factor/character (ignored here, but preserved)
#   - matched_pro_terms:
#       • Character vector of pro-democratic dictionary terms.
#   - matched_anti_terms:
#       • Character vector of anti-democratic dictionary terms.
#
# Parameters:
#   - vec_cols: embedding column names, e.g. V1..V10
#   - top_k (default = 10): number of top positive/negative terms to plot
#
# Outputs:
#   - Console message: counts of pro/anti terms found in embedding vocab
#   - Data frame of term-level indices (term, label, index)
#   - Lollipop plot (saved to disk as PNG):
#       "Top Terms by Democratic-Legitimacy Projection"
#
# Notes:
#   - Terms not present in the embedding vocabulary are ignored.
#   - Centroids are computed as the mean embedding of matched terms
#     within each dictionary (pro vs. anti).
#   - Plot aesthetics:
#       • Pro terms: solid line, circle marker
#       • Anti terms: dashed line, triangle marker
#       • Theme: minimal, Times New Roman font, legend at bottom
# ===============================================================

# load packages

library(dplyr)
library(ggplot2)
library(stringr)

# ---- inputs you provide ----
emb <- control_word_embeddings_df          # columns V1:V10, term, group

pro_terms  <- tolower(matched_pro_terms)
anti_terms <- tolower(matched_anti_terms)

# ensure consistent casing and drop duplicates
emb <- emb %>% mutate(term = tolower(term)) %>% distinct(term, .keep_all = TRUE)

# helper: cosine
cos_sim <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))

# pull embedding matrix
vec_cols <- paste0("V", 1:10)   # adjust if you have a different dimensionality
M <- as.matrix(emb[, vec_cols])

# look up embeddings for dictionary terms that exist in the space
pro_tbl  <- emb %>% filter(term %in% pro_terms)
anti_tbl <- emb %>% filter(term %in% anti_terms)

# guardrails
message("Found ", nrow(pro_tbl),  " pro terms, and ", nrow(anti_tbl), " anti terms in the embedding vocab.")

# centroids
c_pro  <- colMeans(pro_tbl[, vec_cols, drop = FALSE])
c_anti <- colMeans(anti_tbl[, vec_cols, drop = FALSE])

# compute index for DICTIONARY terms only (clear story for the figure)
dict_tbl <- bind_rows(
  pro_tbl  %>% mutate(label = "Pro"),
  anti_tbl %>% mutate(label = "Anti")
)

term_index <- dict_tbl %>%
  rowwise() %>%
  mutate(
    cos_pro  = cos_sim(c_across(all_of(vec_cols)), c_pro),
    cos_anti = cos_sim(c_across(all_of(vec_cols)), c_anti),
    index    = cos_pro - cos_anti
  ) %>%
  ungroup() %>%
  dplyr::select(term, label, index)

# === Lollipop of top ±K by |index| ===
top_k <- 10
plot_df <- term_index %>%
  arrange(desc(index)) %>% slice_head(n = top_k) %>%
  bind_rows(term_index %>% arrange(index) %>% slice_head(n = top_k)) %>%
  mutate(term = factor(term, levels = unique(term[order(index)])),
         side = if_else(index >= 0, "Pro", "Anti"))

ggplot(plot_df, aes(x = term, y = index)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_segment(aes(xend = term, y = 0, yend = index, linetype = side),
               linewidth = 0.6, color = "black") +
  geom_point(aes(shape = side, fill = side), size = 2.8, color = "black") +
  scale_linetype_manual(values = c("Pro" = "solid", "Anti" = "dashed")) +
  scale_shape_manual(values   = c("Pro" = 21,     "Anti" = 24)) +
  scale_fill_manual(values    = c("Pro" = "white","Anti" = "grey70")) +
  coord_flip() +
  labs(
    title = "Top Terms by Democratic-Legitimacy Projection",
    subtitle = "Index = cos(term, pro) − cos(term, anti) in the embedding space",
    x = NULL, y = "Projection index", linetype = "Alignment",
    shape = "Alignment", fill = "Alignment"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"))
