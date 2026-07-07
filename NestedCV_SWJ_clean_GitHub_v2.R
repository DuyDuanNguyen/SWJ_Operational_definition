# ============================================================
# Nested Cross-Validation for SWJ Diagnostic Performance
# Clean GitHub version
#
# Purpose
# -------
# This script evaluates the diagnostic performance of square-wave jerk (SWJ)
# rate for distinguishing PSP from comparator groups. It uses nested
# cross-validation:
#   1. Inner loop: select the SWJ criterion set with the highest validation AUC.
#   2. Outer loop: fit a logistic model on the outer-training set and predict
#      held-out outer-test participants.
#   3. Pooled OOF analysis: average repeated out-of-fold PSP probabilities at
#      participant level.
#   4. Bootstrap: estimate participant-level confidence intervals.
#
# Important assumption
# --------------------
# This version assumes that all required variables are complete in the final
# analysis dataset. Therefore, missing-value handling, NA fallbacks, and
# tryCatch-based silent failures have been removed for clarity.
#
# Required input columns
# ----------------------
#   Diagnosis
#   File                         participant ID
#   SWJ_Rate_per_min             continuous SWJ score
#   AmpMin, AmpMax, ISImin,
#   ISImax, DirTol_deg, SimThr   SWJ criterion-set identifiers
#
# Outputs
# -------
# Per comparison:
#   - fold-level outer-test performance
#   - inner-loop AUC results
#   - selected criteria per outer fold
#   - selected-criteria frequency
#   - row-level out-of-fold predictions
#   - participant-level pooled out-of-fold predictions
#   - bootstrap confidence intervals
#
# Combined outputs:
#   - pooled OOF diagnostic-performance summary for all comparisons
#   - outer-fold performance distribution for all comparisons
#   - selected-criteria frequency for all comparisons
# ============================================================


# =========================
# Packages
# =========================

library(readr)
library(dplyr)
library(purrr)
library(rsample)
library(pROC)


# =========================
# User settings
# =========================

project_dir <- "/home"

# Replace this with the exact CSV filename used for the nested-CV input.
input_file <- file.path(project_dir, "YOUR_SWJ_GRID_RESULTS.csv")

output_dir <- file.path(project_dir, "nested_cv_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

thr_cols <- c("AmpMin", "AmpMax", "ISImin", "ISImax", "DirTol_deg", "SimThr")

group_var <- "File"              # participant ID variable
pos_label <- "PSP"               # positive class
score_var <- "SWJ_Rate_per_min"  # continuous SWJ score
# or "SWJ_amp_per_min"

# Cross-validation settings
v_outer <- 5
repeats_outer <- 40

v_inner <- 5
repeats_inner <- 5

# Minimum numbers per class required inside each model/evaluation fold.
# These are checks for cross-validation design, not missing-value handling.
min_per_group_inner <- 4
min_per_group_outer <- 5

n_boot <- 2000
seed <- 123

prefix <- "NestedCV_SCORE_PooledOOF"


# =========================
# Load and prepare data
# =========================

dat0 <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(Diagnosis %in% c("PSP", "CONT", "PD", "AD")) %>%
  mutate(Diagnosis = as.character(Diagnosis))

# Keep baseline visit only when the file contains a Visit column.
if ("Visit" %in% names(dat0)) {
  dat0 <- dat0 %>% filter(Visit == 1)
}

required_cols <- c("Diagnosis", group_var, score_var, thr_cols)
stopifnot(all(required_cols %in% names(dat0)))


# =========================
# Basic helper functions
# =========================

confusion_metrics <- function(truth, pred, pos_label, neg_label) {
  truth <- factor(truth, levels = c(neg_label, pos_label))
  pred  <- factor(pred,  levels = c(neg_label, pos_label))

  tab <- table(truth, pred)

  TN <- as.integer(tab[neg_label, neg_label])
  FP <- as.integer(tab[neg_label, pos_label])
  FN <- as.integer(tab[pos_label, neg_label])
  TP <- as.integer(tab[pos_label, pos_label])

  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  ppv  <- TP / (TP + FP)
  npv  <- TN / (TN + FN)

  tibble(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Sens = sens,
    Spec = spec,
    PPV = ppv,
    NPV = npv,
    LR_pos = sens / (1 - spec),
    LR_neg = (1 - sens) / spec
  )
}




auc_from_score <- function(d, pos_label, neg_label, score_var, min_per_group = 4) {
  d <- d %>% filter(Diagnosis %in% c(pos_label, neg_label))

  n_pos <- sum(d$Diagnosis == pos_label)
  n_neg <- sum(d$Diagnosis == neg_label)

  if (n_pos < min_per_group || n_neg < min_per_group) {
    stop(
      "Insufficient class counts for ROC/AUC: ",
      pos_label, " = ", n_pos, ", ",
      neg_label, " = ", n_neg,
      ". Reduce the number of folds or the minimum class-size setting."
    )
  }

  roc_obj <- pROC::roc(
    response = factor(d$Diagnosis, levels = c(neg_label, pos_label)),
    predictor = d[[score_var]],
    quiet = TRUE,
    direction = "auto"
  )

  as.numeric(pROC::auc(roc_obj))
}


youden_cutpoint <- function(roc_obj) {
  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "youden"),
    transpose = FALSE
  ) %>%
    as.data.frame() %>%
    arrange(desc(youden), desc(sensitivity), desc(specificity), threshold) %>%
    slice(1)

  list(
    threshold = as.numeric(best$threshold)[1],
    sens = as.numeric(best$sensitivity)[1],
    spec = as.numeric(best$specificity)[1],
    direction = roc_obj$direction
  )
}


filter_by_criteria <- function(d, criteria_row, thr_cols) {
  criteria_key <- criteria_row %>% select(all_of(thr_cols))
  semi_join(d, criteria_key, by = thr_cols)
}


# =========================
# Inner-loop criterion selection
# =========================

evaluate_criteria_auc <- function(d,
                                  pos_label,
                                  neg_label,
                                  thr_cols,
                                  score_var,
                                  min_per_group = 4) {
  d %>%
    filter(Diagnosis %in% c(pos_label, neg_label)) %>%
    group_by(across(all_of(thr_cols))) %>%
    group_modify(~ {
      tibble(
        AUC = auc_from_score(
          d = .x,
          pos_label = pos_label,
          neg_label = neg_label,
          score_var = score_var,
          min_per_group = min_per_group
        )
      )
    }) %>%
    ungroup()
}


select_best_criteria_inner_cv <- function(train_outer,
                                          pos_label,
                                          neg_label,
                                          thr_cols,
                                          group_var,
                                          score_var,
                                          v_inner = 5,
                                          repeats_inner = 5,
                                          min_per_group = 4,
                                          seed = 123) {
  set.seed(seed)

  inner_folds <- rsample::group_vfold_cv(
    train_outer,
    group = !!sym(group_var),
    v = v_inner,
    repeats = repeats_inner,
    strata = Diagnosis
  )

  inner_eval <- map_dfr(seq_along(inner_folds$splits), function(j) {
    inner_valid <- rsample::assessment(inner_folds$splits[[j]])

    evaluate_criteria_auc(
      d = inner_valid,
      pos_label = pos_label,
      neg_label = neg_label,
      thr_cols = thr_cols,
      score_var = score_var,
      min_per_group = min_per_group
    ) %>%
      mutate(InnerFold = j)
  })

  inner_summary <- inner_eval %>%
    group_by(across(all_of(thr_cols))) %>%
    summarise(
      n_valid_inner = n(),
      mean_inner_AUC = mean(AUC),
      median_inner_AUC = median(AUC),
      sd_inner_AUC = sd(AUC),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_inner_AUC), desc(median_inner_AUC), desc(n_valid_inner))

  best_criteria <- inner_summary %>%
    slice(1) %>%
    select(all_of(thr_cols))

  list(
    inner_eval = inner_eval,
    inner_summary = inner_summary,
    best_criteria = best_criteria
  )
}


# =========================
# Outer-fold probability model
# =========================

fit_outer_probability_model <- function(train_df,
                                        test_df,
                                        pos_label,
                                        neg_label,
                                        score_var,
                                        min_per_group_inner = 4,
                                        min_per_group_outer = 5) {
  train_df <- train_df %>%
    filter(Diagnosis %in% c(pos_label, neg_label)) %>%
    mutate(
      y = ifelse(Diagnosis == pos_label, 1, 0),
      Score_model = .data[[score_var]]
    )

  test_df <- test_df %>%
    filter(Diagnosis %in% c(pos_label, neg_label)) %>%
    mutate(
      y = ifelse(Diagnosis == pos_label, 1, 0),
      Score_model = .data[[score_var]]
    )

  n_train_pos <- sum(train_df$Diagnosis == pos_label)
  n_train_neg <- sum(train_df$Diagnosis == neg_label)
  n_test_pos  <- sum(test_df$Diagnosis == pos_label)
  n_test_neg  <- sum(test_df$Diagnosis == neg_label)

  if (n_train_pos < min_per_group_inner || n_train_neg < min_per_group_inner ||
      n_test_pos < min_per_group_outer || n_test_neg < min_per_group_outer) {
    stop(
      "Insufficient class counts in outer fold: ",
      "training ", pos_label, " = ", n_train_pos, ", training ", neg_label, " = ", n_train_neg,
      "; test ", pos_label, " = ", n_test_pos, ", test ", neg_label, " = ", n_test_neg,
      ". Reduce the number of folds or the minimum class-size setting."
    )
  }

  fit <- glm(y ~ Score_model, data = train_df, family = binomial())

  train_df <- train_df %>%
    mutate(OOF_Prob_PSP = predict(fit, newdata = train_df, type = "response"))

  test_df <- test_df %>%
    mutate(OOF_Prob_PSP = predict(fit, newdata = test_df, type = "response"))

  roc_train <- pROC::roc(
    response = factor(train_df$Diagnosis, levels = c(neg_label, pos_label)),
    predictor = train_df$OOF_Prob_PSP,
    quiet = TRUE,
    direction = "<"
  )

  cutpoint <- youden_cutpoint(roc_train)$threshold

  test_df <- test_df %>%
    mutate(Pred_fold = ifelse(OOF_Prob_PSP >= cutpoint, pos_label, neg_label))

  roc_test <- pROC::roc(
    response = factor(test_df$Diagnosis, levels = c(neg_label, pos_label)),
    predictor = test_df$OOF_Prob_PSP,
    quiet = TRUE,
    direction = "<"
  )

  conf <- confusion_metrics(
    truth = test_df$Diagnosis,
    pred = test_df$Pred_fold,
    pos_label = pos_label,
    neg_label = neg_label
  )

  list(
    fit = fit,
    train_df = train_df,
    test_df = test_df,
    auc_test = as.numeric(pROC::auc(roc_test)),
    cutpoint = cutpoint,
    conf = conf
  )
}


# =========================
# Participant-level metrics and bootstrap
# =========================

participant_metrics <- function(df,
                                pos_label,
                                neg_label,
                                score_col = "Mean_OOF_Prob_PSP") {
  df <- df %>% filter(Diagnosis %in% c(pos_label, neg_label))

  roc_obj <- pROC::roc(
    response = factor(df$Diagnosis, levels = c(neg_label, pos_label)),
    predictor = df[[score_col]],
    quiet = TRUE,
    direction = "<"
  )

  cutpoint <- youden_cutpoint(roc_obj)$threshold
  pred <- ifelse(df[[score_col]] >= cutpoint, pos_label, neg_label)

  conf <- confusion_metrics(
    truth = df$Diagnosis,
    pred = pred,
    pos_label = pos_label,
    neg_label = neg_label
  )

  tibble(
    AUC = as.numeric(pROC::auc(roc_obj)),
    Cutpoint = cutpoint,
    TP = conf$TP,
    FP = conf$FP,
    TN = conf$TN,
    FN = conf$FN,
    Sens = conf$Sens,
    Spec = conf$Spec,
    PPV = conf$PPV,
    NPV = conf$NPV,
    LR_pos = conf$LR_pos,
    LR_neg = conf$LR_neg
  )
}


bootstrap_participant_metrics <- function(participant_oof,
                                          pos_label,
                                          neg_label,
                                          group_var,
                                          score_col = "Mean_OOF_Prob_PSP",
                                          n_boot = 2000,
                                          seed = 123) {
  set.seed(seed)

  ids <- participant_oof %>%
    distinct(.data[[group_var]], Diagnosis, Neg)

  boot_res <- map_dfr(seq_len(n_boot), function(b) {
    sampled_ids <- ids %>%
      group_by(Diagnosis) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup() %>%
      mutate(BootRow = row_number())

    boot_df <- sampled_ids %>%
      left_join(participant_oof, by = c(group_var, "Diagnosis", "Neg")) %>%
      mutate(Boot_ID = paste0(.data[[group_var]], "_boot", BootRow))

    participant_metrics(
      df = boot_df,
      pos_label = pos_label,
      neg_label = neg_label,
      score_col = score_col
    ) %>%
      mutate(Boot = b)
  })

  boot_ci <- boot_res %>%
    summarise(
      across(
        c(AUC, Sens, Spec, PPV, NPV, LR_pos, LR_neg),
        list(
          low = ~ quantile(.x, 0.025),
          mid = ~ median(.x),
          high = ~ quantile(.x, 0.975)
        ),
        .names = "{.col}_{.fn}"
      )
    )

  list(
    boot_res = boot_res,
    boot_ci = boot_ci
  )
}


summarise_outer_distribution <- function(fold_summary) {
  fold_summary %>%
    summarise(
      n_outer_iterations = n(),

      AUC_median = median(AUC_test),
      AUC_low = quantile(AUC_test, 0.025),
      AUC_high = quantile(AUC_test, 0.975),

      Sens_median = median(Sens),
      Sens_low = quantile(Sens, 0.025),
      Sens_high = quantile(Sens, 0.975),

      Spec_median = median(Spec),
      Spec_low = quantile(Spec, 0.025),
      Spec_high = quantile(Spec, 0.975),

      PPV_median = median(PPV),
      PPV_low = quantile(PPV, 0.025),
      PPV_high = quantile(PPV, 0.975),

      NPV_median = median(NPV),
      NPV_low = quantile(NPV, 0.025),
      NPV_high = quantile(NPV, 0.975),

      LR_pos_median = median(LR_pos),
      LR_pos_low = quantile(LR_pos, 0.025),
      LR_pos_high = quantile(LR_pos, 0.975),

      LR_neg_median = median(LR_neg),
      LR_neg_low = quantile(LR_neg, 0.025),
      LR_neg_high = quantile(LR_neg, 0.975)
    )
}


summarise_selected_criteria_frequency <- function(inner_selected, thr_cols) {
  inner_selected %>%
    group_by(across(all_of(thr_cols))) %>%
    summarise(
      n_selected = n(),
      mean_selected_inner_AUC = mean(mean_inner_AUC),
      median_selected_inner_AUC = median(median_inner_AUC),
      .groups = "drop"
    ) %>%
    mutate(
      denominator = nrow(inner_selected),
      selected_percent = 100 * n_selected / denominator
    ) %>%
    arrange(desc(n_selected), desc(mean_selected_inner_AUC))
}


# =========================
# Nested CV for one diagnostic comparison
# =========================

nested_cv_one_comparison <- function(dat,
                                     neg_labels,
                                     neg_name = NULL,
                                     pos_label = "PSP",
                                     thr_cols,
                                     group_var,
                                     score_var,
                                     v_outer = 5,
                                     repeats_outer = 40,
                                     v_inner = 5,
                                     repeats_inner = 5,
                                     min_per_group_inner = 4,
                                     min_per_group_outer = 5,
                                     n_boot = 2000,
                                     seed = 123) {
  set.seed(seed)

  if (is.null(neg_name)) {
    neg_name <- if (length(neg_labels) == 1) neg_labels else "OTHERS"
  }

  dat_pair <- dat %>%
    filter(Diagnosis %in% c(pos_label, neg_labels)) %>%
    mutate(Diagnosis = ifelse(Diagnosis == pos_label, pos_label, neg_name))

  outer_folds <- rsample::group_vfold_cv(
    dat_pair,
    group = !!sym(group_var),
    v = v_outer,
    repeats = repeats_outer,
    strata = Diagnosis
  )

  fold_summary_list <- vector("list", length(outer_folds$splits))
  inner_eval_list <- vector("list", length(outer_folds$splits))
  inner_summary_list <- vector("list", length(outer_folds$splits))
  inner_selected_list <- vector("list", length(outer_folds$splits))
  oof_list <- vector("list", length(outer_folds$splits))

  for (i in seq_along(outer_folds$splits)) {
    message("Running ", pos_label, " vs ", neg_name,
            " | outer fold ", i, " / ", length(outer_folds$splits))

    outer_train <- rsample::analysis(outer_folds$splits[[i]])
    outer_test  <- rsample::assessment(outer_folds$splits[[i]])

    inner_out <- select_best_criteria_inner_cv(
      train_outer = outer_train,
      pos_label = pos_label,
      neg_label = neg_name,
      thr_cols = thr_cols,
      group_var = group_var,
      score_var = score_var,
      v_inner = v_inner,
      repeats_inner = repeats_inner,
      min_per_group = min_per_group_inner,
      seed = seed + i
    )

    inner_eval <- inner_out$inner_eval %>%
      mutate(OuterFold = i, Neg = neg_name)

    inner_summary <- inner_out$inner_summary %>%
      mutate(OuterFold = i, Neg = neg_name)

    selected_row <- inner_summary %>%
      slice(1) %>%
      select(
        all_of(thr_cols),
        n_valid_inner,
        mean_inner_AUC,
        median_inner_AUC,
        sd_inner_AUC,
        OuterFold,
        Neg
      )

    best_criteria <- inner_out$best_criteria

    train_best <- filter_by_criteria(outer_train, best_criteria, thr_cols)
    test_best  <- filter_by_criteria(outer_test,  best_criteria, thr_cols)

    outer_out <- fit_outer_probability_model(
      train_df = train_best,
      test_df = test_best,
      pos_label = pos_label,
      neg_label = neg_name,
      score_var = score_var,
      min_per_group_inner = min_per_group_inner,
      min_per_group_outer = min_per_group_outer
    )

    met <- outer_out$conf

    fold_summary_list[[i]] <- tibble(
      Fold = i,
      Neg = neg_name,
      AUC_test = outer_out$auc_test,
      Cutpoint_train = outer_out$cutpoint,

      TP = met$TP,
      FP = met$FP,
      TN = met$TN,
      FN = met$FN,

      Sens = met$Sens,
      Spec = met$Spec,
      PPV = met$PPV,
      NPV = met$NPV,
      LR_pos = met$LR_pos,
      LR_neg = met$LR_neg,

      N_train_PSP = sum(train_best$Diagnosis == pos_label),
      N_train_NEG = sum(train_best$Diagnosis == neg_name),
      N_test_PSP = sum(test_best$Diagnosis == pos_label),
      N_test_NEG = sum(test_best$Diagnosis == neg_name)
    ) %>%
      bind_cols(best_criteria)

    oof_list[[i]] <- outer_out$test_df %>%
      transmute(
        Fold = i,
        Neg = neg_name,
        !!group_var := .data[[group_var]],
        Diagnosis = Diagnosis,
        Score_raw = .data[[score_var]],
        OOF_Prob_PSP = OOF_Prob_PSP,
        Cutpoint_train = outer_out$cutpoint,
        Pred_fold = Pred_fold
      )

    inner_eval_list[[i]] <- inner_eval
    inner_summary_list[[i]] <- inner_summary
    inner_selected_list[[i]] <- selected_row
  }

  fold_summary <- bind_rows(fold_summary_list)
  inner_eval <- bind_rows(inner_eval_list)
  inner_summary <- bind_rows(inner_summary_list)
  inner_selected <- bind_rows(inner_selected_list)
  oof <- bind_rows(oof_list)

  participant_oof <- oof %>%
    group_by(across(all_of(c(group_var, "Diagnosis", "Neg")))) %>%
    summarise(
      Mean_OOF_Prob_PSP = mean(OOF_Prob_PSP),
      Median_OOF_Prob_PSP = median(OOF_Prob_PSP),
      SD_OOF_Prob_PSP = sd(OOF_Prob_PSP),
      N_OOF_predictions = n(),
      Mean_raw_score = mean(Score_raw),
      .groups = "drop"
    )

  pooled_point <- participant_metrics(
    df = participant_oof,
    pos_label = pos_label,
    neg_label = neg_name,
    score_col = "Mean_OOF_Prob_PSP"
  )

  boot_out <- bootstrap_participant_metrics(
    participant_oof = participant_oof,
    pos_label = pos_label,
    neg_label = neg_name,
    group_var = group_var,
    score_col = "Mean_OOF_Prob_PSP",
    n_boot = n_boot,
    seed = seed
  )

  boot_ci <- boot_out$boot_ci

  pooled_metrics <- tibble(
    Comparison = paste0(pos_label, " vs ", neg_name),

    Pooled_OOF_AUC = pooled_point$AUC,
    AUC_boot_low = boot_ci$AUC_low,
    AUC_boot_mid = boot_ci$AUC_mid,
    AUC_boot_high = boot_ci$AUC_high,

    Final_Cutpoint_on_OOF_probability = pooled_point$Cutpoint,

    TP = pooled_point$TP,
    FP = pooled_point$FP,
    TN = pooled_point$TN,
    FN = pooled_point$FN,

    Sens = pooled_point$Sens,
    Sens_boot_low = boot_ci$Sens_low,
    Sens_boot_mid = boot_ci$Sens_mid,
    Sens_boot_high = boot_ci$Sens_high,

    Spec = pooled_point$Spec,
    Spec_boot_low = boot_ci$Spec_low,
    Spec_boot_mid = boot_ci$Spec_mid,
    Spec_boot_high = boot_ci$Spec_high,

    PPV = pooled_point$PPV,
    PPV_boot_low = boot_ci$PPV_low,
    PPV_boot_mid = boot_ci$PPV_mid,
    PPV_boot_high = boot_ci$PPV_high,

    NPV = pooled_point$NPV,
    NPV_boot_low = boot_ci$NPV_low,
    NPV_boot_mid = boot_ci$NPV_mid,
    NPV_boot_high = boot_ci$NPV_high,

    LR_pos = pooled_point$LR_pos,
    LR_pos_boot_low = boot_ci$LR_pos_low,
    LR_pos_boot_mid = boot_ci$LR_pos_mid,
    LR_pos_boot_high = boot_ci$LR_pos_high,

    LR_neg = pooled_point$LR_neg,
    LR_neg_boot_low = boot_ci$LR_neg_low,
    LR_neg_boot_mid = boot_ci$LR_neg_mid,
    LR_neg_boot_high = boot_ci$LR_neg_high,

    N_PSP = sum(participant_oof$Diagnosis == pos_label),
    N_NEG = sum(participant_oof$Diagnosis == neg_name),
    N_total = nrow(participant_oof)
  )

  outer_distribution <- summarise_outer_distribution(fold_summary) %>%
    mutate(Comparison = paste0(pos_label, " vs ", neg_name), .before = 1)

  selected_frequency <- summarise_selected_criteria_frequency(
    inner_selected = inner_selected,
    thr_cols = thr_cols
  ) %>%
    mutate(Comparison = paste0(pos_label, " vs ", neg_name), .before = 1)

  list(
    fold_summary = fold_summary,
    inner_eval = inner_eval,
    inner_summary = inner_summary,
    inner_selected = inner_selected,
    selected_frequency = selected_frequency,
    oof = oof,
    participant_oof = participant_oof,
    pooled_metrics = pooled_metrics,
    boot_res = boot_out$boot_res,
    boot_ci = boot_ci,
    outer_distribution = outer_distribution
  )
}


# =========================
# Run all diagnostic comparisons
# =========================

run_all_comparisons <- function() {
  list(
    CONT = nested_cv_one_comparison(
      dat = dat0,
      neg_labels = "CONT",
      neg_name = "CONT",
      pos_label = pos_label,
      thr_cols = thr_cols,
      group_var = group_var,
      score_var = score_var,
      v_outer = v_outer,
      repeats_outer = repeats_outer,
      v_inner = v_inner,
      repeats_inner = repeats_inner,
      min_per_group_inner = min_per_group_inner,
      min_per_group_outer = min_per_group_outer,
      n_boot = n_boot,
      seed = seed
    ),

    PD = nested_cv_one_comparison(
      dat = dat0,
      neg_labels = "PD",
      neg_name = "PD",
      pos_label = pos_label,
      thr_cols = thr_cols,
      group_var = group_var,
      score_var = score_var,
      v_outer = v_outer,
      repeats_outer = repeats_outer,
      v_inner = v_inner,
      repeats_inner = repeats_inner,
      min_per_group_inner = min_per_group_inner,
      min_per_group_outer = min_per_group_outer,
      n_boot = n_boot,
      seed = seed + 1000
    ),

    AD = nested_cv_one_comparison(
      dat = dat0,
      neg_labels = "AD",
      neg_name = "AD",
      pos_label = pos_label,
      thr_cols = thr_cols,
      group_var = group_var,
      score_var = score_var,
      v_outer = v_outer,
      repeats_outer = repeats_outer,
      v_inner = v_inner,
      repeats_inner = repeats_inner,
      min_per_group_inner = min_per_group_inner,
      min_per_group_outer = min_per_group_outer,
      n_boot = n_boot,
      seed = seed + 2000
    ),

    OTHERS = nested_cv_one_comparison(
      dat = dat0,
      neg_labels = c("CONT", "PD", "AD"),
      neg_name = "OTHERS",
      pos_label = pos_label,
      thr_cols = thr_cols,
      group_var = group_var,
      score_var = score_var,
      v_outer = v_outer,
      repeats_outer = repeats_outer,
      v_inner = v_inner,
      repeats_inner = repeats_inner,
      min_per_group_inner = min_per_group_inner,
      min_per_group_outer = min_per_group_outer,
      n_boot = n_boot,
      seed = seed + 3000
    )
  )
}


# =========================
# Run analysis
# =========================

res <- run_all_comparisons()


# =========================
# Export outputs
# =========================

walk(names(res), function(k) {
  write_csv(res[[k]]$fold_summary,       file.path(output_dir, paste0(prefix, "_FoldSummary_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$inner_eval,         file.path(output_dir, paste0(prefix, "_InnerValidationAUC_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$inner_summary,      file.path(output_dir, paste0(prefix, "_InnerSummary_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$inner_selected,     file.path(output_dir, paste0(prefix, "_InnerSelectedBestCriteria_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$selected_frequency, file.path(output_dir, paste0(prefix, "_SelectedCriteriaFrequency_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$oof,                file.path(output_dir, paste0(prefix, "_OOF_RowLevel_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$participant_oof,    file.path(output_dir, paste0(prefix, "_OOF_ParticipantLevel_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$boot_res,           file.path(output_dir, paste0(prefix, "_BootstrapMetrics_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$boot_ci,            file.path(output_dir, paste0(prefix, "_BootstrapCI_PSP_vs_", k, ".csv")))
  write_csv(res[[k]]$outer_distribution, file.path(output_dir, paste0(prefix, "_OuterDistribution_PSP_vs_", k, ".csv")))
})

summary_tbl <- bind_rows(
  res$CONT$pooled_metrics,
  res$PD$pooled_metrics,
  res$AD$pooled_metrics,
  res$OTHERS$pooled_metrics
)

outer_distribution_tbl <- bind_rows(
  res$CONT$outer_distribution,
  res$PD$outer_distribution,
  res$AD$outer_distribution,
  res$OTHERS$outer_distribution
)

selected_frequency_tbl <- bind_rows(
  res$CONT$selected_frequency,
  res$PD$selected_frequency,
  res$AD$selected_frequency,
  res$OTHERS$selected_frequency
)

write_csv(summary_tbl, file.path(output_dir, paste0(prefix, "_PooledOOF_BootstrapSummary_AllComparisons.csv")))
write_csv(outer_distribution_tbl, file.path(output_dir, paste0(prefix, "_OuterDistribution_AllComparisons.csv")))
write_csv(selected_frequency_tbl, file.path(output_dir, paste0(prefix, "_SelectedCriteriaFrequency_AllComparisons.csv")))

saveRDS(res, file.path(output_dir, paste0(prefix, "_FullResultsObject.rds")))

print(summary_tbl)
print(outer_distribution_tbl)
print(selected_frequency_tbl)

message("Analysis completed successfully. Outputs saved to: ", output_dir)
