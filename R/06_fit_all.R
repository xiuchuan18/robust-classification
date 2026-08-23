# 06_fit_all.R ---------------------------------------------------------------
# Fit every classifier on one train/test split, time each end-to-end (fit + predict, including each method's own tuning/CV), and return a tidy
# data.frame of AUC / accuracy / F1 / seconds and, where meaningful, selection accuracy. Shared by the simulation and real-data runners.

MODEL_ORDER <- c("RoLLR", "S-RoLLR", "MOKE", "GMM-NB", "L1-RL", "SVM", "RF", "LightGBM")

run_all_models <- function(X_train, y_train, X_test, y_test,
                           signal_names = NULL, tune = TRUE) {
  X_train <- as.matrix(X_train); X_test <- as.matrix(X_test)
  if (is.null(colnames(X_train))) {
    colnames(X_train) <- paste0("V", seq_len(ncol(X_train)))
    colnames(X_test)  <- colnames(X_train)
  }
  yq <- as.numeric(as.character(y_train))
  yt <- as.numeric(as.character(y_test))
  p  <- ncol(X_train)
  yf <- factor(as.character(y_train), levels = c("0", "1"))

  blank <- function(model) data.frame(
    model = model, auc = NA, acc = NA, f1 = NA, seconds = NA,
    n_selected = NA, sel_tpr = NA, sel_fpr = NA, sel_f1 = NA,
    stringsAsFactors = FALSE)

  row <- function(model, prob, seconds, selected = NULL) {
    m <- classification_metrics(prob, yt)
    r <- blank(model)
    r$auc <- m$auc; r$acc <- m$acc; r$f1 <- m$f1; r$seconds <- seconds
    if (!is.null(selected)) {
      r$n_selected <- length(selected)
      if (!is.null(signal_names)) {
        sm <- selection_metrics(selected, signal_names, p)
        r$sel_tpr <- sm$tpr; r$sel_fpr <- sm$fpr; r$sel_f1 <- sm$f1
      }
    }
    r
  }

  out <- list()

  ## RoLLR
  out$RoLLR <- tryCatch({
    tm <- timed({ f <- fit_rollr(X_train, yq); predict_rollr_prob(f, X_test) })
    row("RoLLR", tm$value, tm$seconds)
  }, error = function(e) blank("RoLLR"))

  ## S-RoLLR
  out$`S-RoLLR` <- tryCatch({
    f <- NULL                     
    tm <- timed({ f <- fit_srollr(X_train, yq, parallel = FALSE)
                  predict_srollr(f, X_test) })
    row("S-RoLLR", tm$value, tm$seconds, selected = srollr_selected(f))
  }, error = function(e) blank("S-RoLLR"))

  ## MOKE
  out$MOKE <- tryCatch({
    tm <- timed({ f <- fit_moke(X_train, yf); predict_moke_prob(f, X_test) })
    row("MOKE", tm$value, tm$seconds)
  }, error = function(e) blank("MOKE"))

  ## GMM-NB
  out$`GMM-NB` <- tryCatch({
    tm <- timed({ f <- fit_gmm_nb(X_train, yq); predict_gmm_nb_prob(f, X_test) })
    row("GMM-NB", tm$value, tm$seconds)
  }, error = function(e) blank("GMM-NB"))

  ## L1-RL
    out$`L1-RL` <- tryCatch({
    f <- NULL
    tm <- timed({ f <- cv.glmnet(X_train, yq, family = "binomial", alpha = 1, parallel = FALSE)
                  as.numeric(predict(f, newx = X_test, s = "lambda.min", type = "response")) })
    co  <- as.matrix(coef(f, s = "lambda.min"))
    sel <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
    row("L1-RL", tm$value, tm$seconds, selected = sel)
  }, error = function(e) blank("L1-RL"))

  ## SVM (RBF)
  out$SVM <- tryCatch({
    tm <- timed({
      if (tune) {
        tr <- e1071::tune(svm, train.x = X_train, train.y = yf, kernel = "radial",
                          ranges = list(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1)),
                          tunecontrol = tune.control(cross = 3, best.model = FALSE))
        cost <- tr$best.parameters$cost; gam <- tr$best.parameters$gamma
      } else { cost <- 1; gam <- 1 / p }
      f  <- svm(x = X_train, y = yf, kernel = "radial", cost = cost, gamma = gam,
                probability = TRUE)
      pm <- attr(predict(f, X_test, probability = TRUE), "probabilities")
      if ("1" %in% colnames(pm)) pm[, "1"] else pm[, 2]
    })
    row("SVM", tm$value, tm$seconds)
  }, error = function(e) blank("SVM"))

  ## Random Forest
  out$RF <- tryCatch({
    tm <- timed({
      mtry <- if (tune) {
        tf <- randomForest::tuneRF(X_train, yf, stepFactor = 1.5, improve = 0.01,
                                   ntreeTry = 200, plot = FALSE, trace = FALSE, doBest = FALSE)
        tf[which.min(tf[, "OOBError"]), "mtry"]
      } else floor(sqrt(p))
      f <- randomForest(x = X_train, y = yf, ntree = 200, mtry = mtry)
      pm <- predict(f, X_test, type = "prob")
      if ("1" %in% colnames(pm)) pm[, "1"] else pm[, 2]
    })
    row("RF", tm$value, tm$seconds)
  }, error = function(e) blank("RF"))

  ## LightGBM
  out$LightGBM <- tryCatch({
    tm <- timed({
      dtrain <- lgb.Dataset(data = X_train, label = yq)
      f <- lgb.train(list(objective = "binary", metric = "auc", verbose = -1,
                          learning_rate = 0.1, num_leaves = 31,
                          min_data_in_leaf = 20, feature_fraction = 0.8,
                          bagging_fraction = 0.8, bagging_freq = 5,
                          force_row_wise = TRUE),
                     dtrain, nrounds = 100, verbose = -1)
      predict(f, X_test)
    })
    row("LightGBM", tm$value, tm$seconds)
  }, error = function(e) blank("LightGBM"))

  res <- do.call(rbind, out[MODEL_ORDER])
  rownames(res) <- NULL
  res
}
