# 05_metrics.R

## classification metrics from predicted P(class 1)
classification_metrics <- function(prob, y_true) {  # y_true numeric 0/1
  if (anyNA(prob) || !all(is.finite(prob)))
    return(list(acc = NA, auc = NA, f1 = NA))
  pred <- as.integer(prob >= 0.5)
  acc  <- mean(pred == y_true)
  auc  <- tryCatch(as.numeric(pROC::roc(y_true, prob, quiet = TRUE)$auc),
                   error = function(e) NA)
  tp <- sum(pred == 1 & y_true == 1); fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  prec <- if (tp + fp > 0) tp / (tp + fp) else 0
  rec  <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  list(acc = acc, auc = auc, f1 = f1)
}

## selection metric
# selected, truth: character vectors of feature names. Returns TPR/FPR/F1 plus size. `p_total` is needed for FPR (number of true noise features).
selection_metrics <- function(selected, signal_names, p_total) {
  s  <- length(signal_names); noise <- p_total - s
  tp <- length(intersect(selected, signal_names))
  fp <- length(setdiff(selected, signal_names))
  fn <- s - tp
  tpr <- if (s > 0) tp / s else NA                       # sensitivity
  fpr <- if (noise > 0) fp / noise else NA                # false-positive rate
  prec <- if (tp + fp > 0) tp / (tp + fp) else 0
  f1  <- if (prec + tpr > 0 && !is.na(tpr)) 2 * prec * tpr / (prec + tpr) else 0
  list(n_selected = length(selected), tpr = tpr, fpr = fpr, f1 = f1)
}

## timing helper
# Runs expr, returns its value plus elapsed seconds. Use one timer around
# fit+predict so the reported time is the end-to-end cost of a method.
timed <- function(expr) {
  t0  <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, seconds = proc.time()[["elapsed"]] - t0)
}
