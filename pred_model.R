##### LASSO Regression for Classification Between Treatment and Control Groups #####
# Author: Cora Fagan
# Date: August 2025

# This code takes document-level embeddings (using GloVe) at 50 dimensions and classifies them using LASSO Regression 

# Load packages
library(dplyr)
library(caret)
library(glmnet)
library(pROC)

set.seed(123) # for reproducability

# Build full frame (50 dims + group)
full_df <- combined_embeddings %>%
  mutate(group = factor(group, levels = c("control","treatment")))

# Split once; everything downstream uses ONLY train to choose features
idx <- createDataPartition(full_df$group, p = 0.8, list = FALSE)
train_df <- full_df[idx, ]
test_df  <- full_df[-idx, ]

# ---- Feature selection on training data ----
X_train_full <- as.matrix(train_df %>% dplyr::select(starts_with("V")))
y_train      <- train_df$group

pvals <- sapply(1:ncol(X_train_full), function(i) {
  x <- X_train_full[, i]
  t.test(x[y_train=="treatment"], x[y_train=="control"])$p.value
})
pvals_adj <- p.adjust(pvals, method = "BH")

top_n   <- 10
top_ix  <- order(pvals_adj)[1:top_n]
top_cols <- colnames(X_train_full)[top_ix]

# Build train/test matrices with the SAME columns; standardize with train stats
X_train <- scale(X_train_full[, top_ix, drop=FALSE])
ctr    <- attr(X_train, "scaled:center"); scl <- attr(X_train, "scaled:scale")

X_test  <- as.matrix(test_df %>% dplyr::select(all_of(top_cols)))
X_test  <- sweep(sweep(X_test, 2, ctr, "-"), 2, scl, "/")
y_test  <- test_df$group

# ---- LASSO logistic with CV on training data  ----
cv_fit <- cv.glmnet(x = X_train, y = y_train, family = "binomial", alpha = 1)

# Predict with lambda.min and lambda.1se
for (lam in c("lambda.min","lambda.1se")) {
  prob_hat <- predict(cv_fit, newx = X_test, s = lam, type = "response")
  pred_cls <- factor(ifelse(prob_hat >= 0.5, "treatment", "control"),
                     levels = levels(y_test))
  
  cm  <- caret::confusionMatrix(pred_cls, y_test, positive = "treatment")
  roc <- pROC::roc(y_test, as.numeric(prob_hat), levels = c("control","treatment"))
  
  cat("\nLASSO (", lam, "), top", top_n, "dims\n",
      "Accuracy:", round(cm$overall["Accuracy"], 3), "\n",
      "AUC:", round(as.numeric(pROC::auc(roc)), 3), "\n",
      "Selected λ:", cv_fit[[lam]], "\n", sep = "")
  print(cm)
  
  # Non-zero coefficients
  nz <- coef(cv_fit, s = lam)
  cat("Non-zero coefficients:\n")
  print(nz[nz[,1]!=0, , drop=FALSE])
}
