# embeddings-classification-model

# Media Coverage & Election Reform (NLP)

This repository contains a set of R scripts for analyzing **media coverage of election reform** using word/document embeddings and permutation-based inference. The workflow combines dictionary-based approaches, embedding visualizations, and predictive modeling.

---

## Repository Structure

### Core Analysis Scripts

- **`dictionary_similarity_test.R`**  
  Compares pro- vs. anti-democratic dictionary terms in embedding space.  
  - Computes centroids for dictionary-based categories.  
  - Tests whether terms cluster more tightly than expected under random labeling.  
  - Outputs: within-cluster similarity stats, permutation p-values, and optional LaTeX tables.

- **`word_embedding_dict_inf_test.R`**  
  Anchors embeddings to pro-/anti-democratic dictionary centroids and computes an **index** for each token.  
  - Evaluates group differences in average index (treatment vs. control).  
  - Outputs: observed difference, permutation p-values, bootstrap confidence intervals, and effect sizes (Cohen’s *d*, Hedges’ *g*).

- **`term_permutation_test.R`**  
  Tests whether specific subsets of terms (e.g., pro-democratic vs. anti-democratic) differ in embedding alignment.  
  - Uses permutation tests to assess robustness.  
  - Produces micro-level tables of cohesion and separation.

- **`doc_permutation_test.R`**  
  Runs permutation tests on **document-level embeddings**.  
  - Computes centroid similarity and difference in within-group similarity.  
  - Packages results with parametric OLS summaries, permutation distributions, and p-values.  
  - Useful for testing whether treatment vs. control document sets are separable in embedding space.

- **`targeted-permutation.R`**  
  Implements a more **fine-grained permutation test** for targeted subsets of tokens or documents.  
  - Allows hypothesis testing on pre-specified slices of the data.

### Visualization Scripts

- **`embed_viz.R`**  
  Provides visualization utilities for embeddings.  
  - Functions for plotting document embeddings in 2D space (PCA/UMAP).  
  - Options for balancing group sizes, drawing centroid separation, and bootstrap CI ellipses.

- **`rcv_pca_plot.R`**  
  Specialized PCA plotting script for comparing RCV (treatment) vs. control embeddings.  
  - Produces 2D projections with centroids and annotated statistics.

### Predictive Modeling

- **`pred_model.R`**  
  Fits predictive models (e.g., logistic regression with regularization) using embeddings as features.  
  - Allows classification of treatment vs. control groups.  
  - Reports accuracy, AUC, confusion matrices, and selected features.

### Supporting Files

- **`README.md`**  
  You are here. Describes repository structure and usage.  

---

## Expected Data Inputs

Most scripts expect data frames or matrices of embeddings with the following structure:

- **Embedding columns:** `V1 ... Vk` (numeric, k dimensions).  
- **term:** character token (for term-level analyses).  
- **doc_id:** unique document ID (for doc-level analyses).  
- **group:** factor/character with at least two levels: `"control"`, `"treatment"`.  

Additionally:

- **`rcv_sentiment_dict`:** A `quanteda::dictionary2` object with two entries:  
  - `"pro_democratic"` → glob patterns for pro-democratic terms  
  - `"anti_democratic"` → glob patterns for anti-democratic terms  

---

## Outputs

Across scripts, outputs include:

- **Permutation test results** (observed stats vs. null distributions)  
- **Bootstrap confidence intervals** for robustness  
- **Effect sizes** (Cohen’s *d*, Hedges’ *g*)  
- **Visualizations** (2D projections, lollipop plots, null distributions)  
- **Publication-ready LaTeX tables** (via `kableExtra`)  

---

## Getting Started

### Packages required

- quanteda  
- tidyverse  
- tidytext  
- topicmodels  
- stm  
- textdata  
- quanteda.textstats  
- quanteda.textplots  
- MatchIt  
- bit64  
- glmnet  
- umap  
- uwot  
- text2vec  
- Matrix  
- stargazer  
- cobalt  
- kableExtra  
- qs  
- Hotelling  
- reshape2  
- irr  
- ggpattern  
- Rtsne  
- lsa  
- ggrepel  
- plotly  
- coop  
- caret  
- DiagrammeR  
- DiagrammeRsvg  
- rsvg  
- irlba  
- ggforce  

### Prepare your embeddings

- **`combined_embeddings`** — for document-level analyses  
- **`rcv_word_embeddings_df` / `control_word_embeddings_df`** — for term-level analyses  
- **`rcv_sentiment_dict`** — for dictionary-based tests  

You can run the scripts independently, or source them into a single workflow.
