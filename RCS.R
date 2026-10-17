
rm(list = ls())
graphics.off()

library(tidyverse)
library(readxl)
library(mice)
library(Hmisc)

# Data preparation
dat_raw <- read_excel("~/Desktop/混合效应模型.xlsx")

names(dat_raw) <- names(dat_raw) %>%
  stringr::str_replace_all("\\s+", "_")

dat_long <- dat_raw %>%
  pivot_longer(
    cols = matches("^(DQ|ZBMI|T)[1-4]$"),
    names_to = c(".value", "wave"),
    names_pattern = "^(DQ|ZBMI|T)([1-4])$"
  ) %>%
  mutate(
    wave = as.integer(wave),
    time_month = as.numeric(T)
  ) %>%
  select(-T)

dat_long <- dat_long %>%
  mutate(
    ID = factor(ID),
    DQ = as.numeric(DQ),
    ZBMI = as.numeric(ZBMI),
    MATERNAL_AGE = as.numeric(MATERNAL_AGE),
    WT500 = as.numeric(WT500),
    LEN5 = as.numeric(LEN5),
    SEX = factor(SEX),
    GESTATIONAL_WEEK = factor(GESTATIONAL_WEEK),
    NUMBER_OF_BIRTHS = factor(NUMBER_OF_BIRTHS),
    GRAVIDITY = factor(GRAVIDITY),
    DELIVERY_METHOD = factor(DELIVERY_METHOD),
    BF = factor(BF)
  )

lag_dat <- dat_long %>%
  arrange(ID, time_month) %>%
  group_by(ID) %>%
  mutate(
    time_now = time_month,
    time_next = lead(time_month),
    ZBMI_now = ZBMI,
    DQ_now = DQ,
    ZBMI_next = lead(ZBMI),
    DQ_next = lead(DQ),
    lag_interval = time_next - time_now
  ) %>%
  ungroup() %>%
  filter(!is.na(time_next))

cat("\n滞后配对检查：\n")
print(
  lag_dat %>%
    count(time_now, time_next, lag_interval)
)

time_now_mean <- mean(lag_dat$time_now, na.rm = TRUE)
time_now_sd <- sd(lag_dat$time_now, na.rm = TRUE)

lag_interval_mean <- mean(lag_dat$lag_interval, na.rm = TRUE)
lag_interval_sd <- sd(lag_dat$lag_interval, na.rm = TRUE)

lag_dat <- lag_dat %>%
  mutate(
    time_now_z = as.numeric((time_now - time_now_mean) / time_now_sd),
    lag_interval_z = as.numeric((lag_interval - lag_interval_mean) / lag_interval_sd)
  )

imp_vars <- c(
  "ID",
  "DQ_next",
  "ZBMI_next",
  "ZBMI_now",
  "DQ_now",
  "time_now_z",
  "lag_interval_z",
  "MATERNAL_AGE",
  "SEX",
  "GESTATIONAL_WEEK",
  "NUMBER_OF_BIRTHS",
  "GRAVIDITY",
  "DELIVERY_METHOD",
  "BF",
  "WT500",
  "LEN5"
)

lag_imp_base <- lag_dat %>%
  select(any_of(imp_vars))

cat("\n插补前缺失情况：\n")
missing_overview <- tibble(
  variable = names(lag_imp_base),
  n_missing = colSums(is.na(lag_imp_base)),
  pct_missing = round(colMeans(is.na(lag_imp_base)) * 100, 2)
) %>%
  arrange(desc(n_missing))
print(missing_overview, n = Inf)

# Multiple imputation
ini <- mice(lag_imp_base, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

meth["ID"] <- ""
pred[, "ID"] <- 0
pred["ID", ] <- 0

for (v in c("time_now_z", "lag_interval_z")) {
  if (v %in% names(meth)) {
    meth[v] <- ""
  }
}

num_vars <- names(lag_imp_base)[sapply(lag_imp_base, is.numeric)]
num_vars_miss <- num_vars[colSums(is.na(lag_imp_base[num_vars])) > 0]
meth[num_vars_miss] <- "pmm"

set.seed(2025)
imp_lag <- mice(
  lag_imp_base,
  m = 20,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  seed = 2025,
  printFlag = TRUE
)

plot(imp_lag)

imp_lag_list <- complete(imp_lag, action = "all")

# Helper functions
add_sig <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ ""
  )
}

get_mode_value <- function(x) {
  x_no_na <- x[!is.na(x)]
  if (length(x_no_na) == 0) return(NA)
  tab <- table(x_no_na)
  names(tab)[which.max(tab)]
}

add_rcs_basis <- function(dat, exposure, knots_value, prefix) {
  rcs_basis <- Hmisc::rcspline.eval(
    dat[[exposure]],
    knots = knots_value,
    inclx = TRUE
  )
  rcs_basis <- as.data.frame(rcs_basis)
  names(rcs_basis) <- paste0(prefix, "_s", seq_len(ncol(rcs_basis)))
  bind_cols(dat, rcs_basis)
}

get_fixed_knots <- function(data_list, exposure, knots_n) {
  all_x <- unlist(lapply(data_list, function(dat) dat[[exposure]]))
  all_x <- all_x[!is.na(all_x)]
  Hmisc::rcspline.eval(
    all_x,
    nk = knots_n,
    knots.only = TRUE
  )
}

make_formula <- function(outcome, rhs_terms) {
  as.formula(
    paste(
      outcome,
      "~",
      paste(rhs_terms, collapse = " + ")
    )
  )
}

select_best_knots_explicit <- function(
    data_list,
    outcome,
    exposure,
    adjust_vars,
    prefix,
    knot_candidates = 3:5
) {
  aic_table <- purrr::map_dfr(
    knot_candidates,
    function(k) {
      knots_value <- get_fixed_knots(
        data_list = data_list,
        exposure = exposure,
        knots_n = k
      )
      data_list_rcs <- lapply(
        data_list,
        add_rcs_basis,
        exposure = exposure,
        knots_value = knots_value,
        prefix = prefix
      )
      n_basis <- length(paste0(prefix, "_s", seq_len(k - 1)))
      rcs_terms <- paste0(prefix, "_s", seq_len(n_basis))
      formula_full <- make_formula(outcome, c(rcs_terms, adjust_vars))
      aics <- purrr::map_dbl(
        data_list_rcs,
        function(dat) {
          mod <- lm(formula_full, data = dat)
          AIC(mod)
        }
      )
      tibble(
        knots = k,
        mean_AIC = mean(aics, na.rm = TRUE),
        sd_AIC = sd(aics, na.rm = TRUE)
      )
    }
  )
  print(aic_table)
  best_k <- aic_table$knots[which.min(aic_table$mean_AIC)]
  cat(sprintf("平均AIC最小，选择 %d 个节点。\n", best_k))
  best_k
}

run_rcs_mi_explicit <- function(
    data_list,
    outcome,
    exposure,
    adjust_vars,
    knots_n,
    prefix
) {
  knots_value <- get_fixed_knots(
    data_list = data_list,
    exposure = exposure,
    knots_n = knots_n
  )

  data_list_rcs <- lapply(
    data_list,
    add_rcs_basis,
    exposure = exposure,
    knots_value = knots_value,
    prefix = prefix
  )

  n_basis <- ncol(
    Hmisc::rcspline.eval(
      data_list[[1]][[exposure]],
      knots = knots_value,
      inclx = TRUE
    )
  )

  linear_term <- paste0(prefix, "_s1")
  nonlinear_terms <- paste0(prefix, "_s", 2:n_basis)
  full_rcs_terms <- c(linear_term, nonlinear_terms)

  formula_full <- make_formula(outcome, c(full_rcs_terms, adjust_vars))
  formula_no_exposure <- make_formula(outcome, adjust_vars)
  formula_linear <- make_formula(outcome, c(linear_term, adjust_vars))

  models_full <- lapply(
    data_list_rcs,
    function(dat) lm(formula_full, data = dat)
  )
  models_no_exposure <- lapply(
    data_list_rcs,
    function(dat) lm(formula_no_exposure, data = dat)
  )
  models_linear <- lapply(
    data_list_rcs,
    function(dat) lm(formula_linear, data = dat)
  )

  d1_overall <- mice::D1(
    as.mira(models_full),
    as.mira(models_no_exposure)
  )
  d1_nonlinear <- mice::D1(
    as.mira(models_full),
    as.mira(models_linear)
  )

  list(
    knots = knots_value,
    data_list_rcs = data_list_rcs,
    linear_term = linear_term,
    nonlinear_terms = nonlinear_terms,
    full_rcs_terms = full_rcs_terms,
    formula_full = formula_full,
    formula_no_exposure = formula_no_exposure,
    formula_linear = formula_linear,
    models_full = models_full,
    models_no_exposure = models_no_exposure,
    models_linear = models_linear,
    d1_overall = d1_overall,
    d1_nonlinear = d1_nonlinear
  )
}

pool_predictions_lm <- function(
    rcs_result,
    exposure,
    prefix,
    adjust_vars,
    n_grid = 200
) {
  data_list_rcs <- rcs_result$data_list_rcs
  models <- rcs_result$models_full
  knots_value <- rcs_result$knots

  all_x <- unlist(lapply(data_list_rcs, function(dat) dat[[exposure]]))
  x_grid <- seq(
    min(all_x, na.rm = TRUE),
    max(all_x, na.rm = TRUE),
    length.out = n_grid
  )

  stacked <- bind_rows(data_list_rcs)

  base_new <- tibble(!!exposure := x_grid)

  for (v in adjust_vars) {
    if (is.numeric(stacked[[v]])) {
      base_new[[v]] <- mean(stacked[[v]], na.rm = TRUE)
    } else {
      mode_v <- get_mode_value(stacked[[v]])
      base_new[[v]] <- mode_v
    }
  }

  base_new <- add_rcs_basis(
    dat = base_new,
    exposure = exposure,
    knots_value = knots_value,
    prefix = prefix
  )

  pred_all <- purrr::map2_dfr(
    models,
    seq_along(models),
    function(mod, i) {
      nd <- base_new
      train_dat <- data_list_rcs[[i]]

      for (v in adjust_vars) {
        if (is.factor(train_dat[[v]])) {
          nd[[v]] <- factor(nd[[v]], levels = levels(train_dat[[v]]))
        }
      }

      pr <- predict(mod, newdata = nd, se.fit = TRUE)
      tibble(
        x = x_grid,
        yhat = as.numeric(pr$fit),
        se = as.numeric(pr$se.fit),
        .imp = i
      )
    }
  )

  m <- length(models)

  pred_pooled <- pred_all %>%
    group_by(x) %>%
    summarise(
      yhat = mean(yhat, na.rm = TRUE),
      within_var = mean(se^2, na.rm = TRUE),
      between_var = var(yhat, na.rm = TRUE),
      between_var = ifelse(is.na(between_var), 0, between_var),
      total_var = within_var + (1 + 1 / m) * between_var,
      se_total = sqrt(total_var),
      lower = yhat - 1.96 * se_total,
      upper = yhat + 1.96 * se_total,
      .groups = "drop"
    ) %>%
    rename(!!exposure := x)

  pred_pooled
}

plot_rcs_base <- function(
    pred_df,
    exposure,
    yhat_col = "yhat",
    lower_col = "lower",
    upper_col = "upper",
    xlab,
    ylab,
    main,
    p_overall,
    p_nonlinear,
    ref_vline = NULL,
    fixed_ylim = NULL
) {
  x <- pred_df[[exposure]]
  yhat <- pred_df[[yhat_col]]
  lower <- pred_df[[lower_col]]
  upper <- pred_df[[upper_col]]

  max_id <- which.max(yhat)
  peak_x <- x[max_id]
  peak_y <- yhat[max_id]

  par(
    bty = "n",
    las = 1,
    mar = c(5.5, 7, 4, 2)
  )

  x_range <- range(x, na.rm = TRUE)
  x_ticks <- pretty(x_range, n = 6)

  if (is.null(fixed_ylim)) {
    y_range <- range(c(lower, upper), na.rm = TRUE)
    y_ticks <- pretty(y_range, n = 6)
  } else {
    y_range <- fixed_ylim
    y_ticks <- seq(fixed_ylim[1], fixed_ylim[2], length.out = 5)
    y_ticks <- pretty(y_range, n = 5)
  }

  plot(
    x,
    yhat,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    cex.main = 0.95,
    xlim = x_range,
    ylim = y_range
  )

  abline(
    v = x_ticks,
    h = y_ticks,
    lty = 2,
    col = "grey88",
    lwd = 1
  )

  polygon(
    c(x, rev(x)),
    c(lower, rev(upper)),
    col = rgb(142, 136, 201, maxColorValue = 255, alpha = 70),
    border = NA
  )

  lines(
    x,
    yhat,
    col = "#8E88C9",
    lwd = 2.5
  )

  if (!is.null(ref_vline)) {
    abline(
      v = ref_vline,
      lty = 2,
      col = "grey45",
      lwd = 1.5
    )
  }

  abline(
    v = peak_x,
    lty = 1,
    col = rgb(142, 136, 201, maxColorValue = 255, alpha = 180),
    lwd = 1.2
  )

  points(
    x = peak_x,
    y = peak_y,
    pch = 19,
    col = "#8E88C9",
    cex = 1.2
  )

  axis(
    1,
    at = x_ticks,
    lwd = 0,
    lwd.ticks = 1,
    cex.axis = 1.2
  )

  axis(
    2,
    at = y_ticks,
    lwd = 0,
    lwd.ticks = 1,
    cex.axis = 1.2
  )

  mtext(
    xlab,
    side = 1,
    line = 3.2,
    cex = 1.6,
    font = 2
  )

  mtext(
    ylab,
    side = 2,
    line = 4.5,
    cex = 1.5,
    font = 2,
    las = 0
  )

  legend(
    "topleft",
    legend = c(
      paste0(
        "P-overall ",
        ifelse(
          p_overall < 0.001,
          "< 0.001",
          paste0("= ", sprintf("%.3f", p_overall))
        )
      ),
      paste0(
        "P-nonlinear ",
        ifelse(
          p_nonlinear < 0.001,
          "< 0.001",
          paste0("= ", sprintf("%.3f", p_nonlinear))
        )
      )
    ),
    bty = "o",
    box.col = "grey70",
    bg = "white",
    cex = 1.1,
    inset = 0.02
  )

  invisible(
    tibble(
      peak_x = peak_x,
      peak_y = peak_y
    )
  )
}

cat("\n==============================\n")
cat("分析1：ZBMI_now 预测 DQ_next\n")
cat("==============================\n")

adjust_vars_zbmi_to_dq <- c(
  "DQ_now",
  "time_now_z",
  "lag_interval_z",
  "MATERNAL_AGE",
  "SEX",
  "GESTATIONAL_WEEK",
  "NUMBER_OF_BIRTHS",
  "GRAVIDITY",
  "DELIVERY_METHOD",
  "BF",
  "WT500",
  "LEN5"
)

cat("自动选择最佳节点数：ZBMI_now -> DQ_next\n")
best_knots_zbmi_to_dq <- select_best_knots_explicit(
  data_list = imp_lag_list,
  outcome = "DQ_next",
  exposure = "ZBMI_now",
  adjust_vars = adjust_vars_zbmi_to_dq,
  prefix = "ZBMI"
)

rcs_zbmi_to_dq <- run_rcs_mi_explicit(
  data_list = imp_lag_list,
  outcome = "DQ_next",
  exposure = "ZBMI_now",
  adjust_vars = adjust_vars_zbmi_to_dq,
  knots_n = best_knots_zbmi_to_dq,
  prefix = "ZBMI"
)

cat("\nZBMI_now -> DQ_next：使用的RCS节点\n")
print(rcs_zbmi_to_dq$knots)

cat("\nZBMI_now -> DQ_next：完整模型公式\n")
print(rcs_zbmi_to_dq$formula_full)

cat("\nZBMI_now -> DQ_next：P-overall，RCS总体效应\n")
print(rcs_zbmi_to_dq$d1_overall)

cat("\nZBMI_now -> DQ_next：P-nonlinear，RCS非线性效应\n")
print(rcs_zbmi_to_dq$d1_nonlinear)

p_overall_zbmi_to_dq <- rcs_zbmi_to_dq$d1_overall$p.value[1]
p_nonlinear_zbmi_to_dq <- rcs_zbmi_to_dq$d1_nonlinear$p.value[1]

p_rcs_zbmi_to_dq_df <- pool_predictions_lm(
  rcs_result = rcs_zbmi_to_dq,
  exposure = "ZBMI_now",
  prefix = "ZBMI",
  adjust_vars = adjust_vars_zbmi_to_dq,
  n_grid = 200
)

cat("\nZBMI_now -> DQ_next：预测曲线数据前6行\n")
print(head(p_rcs_zbmi_to_dq_df))

peak_zbmi_to_dq <- p_rcs_zbmi_to_dq_df %>%
  slice(which.max(yhat)) %>%
  transmute(
    peak_zbmi = ZBMI_now,
    peak_dq = yhat
  )

cat("\nZBMI_now -> DQ_next：预测峰值，仅作图示辅助\n")
print(peak_zbmi_to_dq)

plot_rcs_base(
  pred_df = p_rcs_zbmi_to_dq_df,
  exposure = "ZBMI_now",
  xlab = "Current ZBMI",
  ylab = "Predicted DQ at next survey",
  main = "Restricted cubic spline analysis of ZBMI and subsequent DQ",
  p_overall = p_overall_zbmi_to_dq,
  p_nonlinear = p_nonlinear_zbmi_to_dq,
  ref_vline = 0,
  fixed_ylim = c(94, 102)
)

cat("\n==============================\n")
cat("分析2：DQ_now 预测 ZBMI_next\n")
cat("==============================\n")

adjust_vars_dq_to_zbmi <- c(
  "ZBMI_now",
  "time_now_z",
  "lag_interval_z",
  "MATERNAL_AGE",
  "SEX",
  "GESTATIONAL_WEEK",
  "NUMBER_OF_BIRTHS",
  "GRAVIDITY",
  "DELIVERY_METHOD",
  "BF",
  "WT500",
  "LEN5"
)

cat("自动选择最佳节点数：DQ_now -> ZBMI_next\n")
best_knots_dq_to_zbmi <- select_best_knots_explicit(
  data_list = imp_lag_list,
  outcome = "ZBMI_next",
  exposure = "DQ_now",
  adjust_vars = adjust_vars_dq_to_zbmi,
  prefix = "DQ"
)

rcs_dq_to_zbmi <- run_rcs_mi_explicit(
  data_list = imp_lag_list,
  outcome = "ZBMI_next",
  exposure = "DQ_now",
  adjust_vars = adjust_vars_dq_to_zbmi,
  knots_n = best_knots_dq_to_zbmi,
  prefix = "DQ"
)

cat("\nDQ_now -> ZBMI_next：使用的RCS节点\n")
print(rcs_dq_to_zbmi$knots)

cat("\nDQ_now -> ZBMI_next：完整模型公式\n")
print(rcs_dq_to_zbmi$formula_full)

cat("\nDQ_now -> ZBMI_next：P-overall，RCS总体效应\n")
print(rcs_dq_to_zbmi$d1_overall)

cat("\nDQ_now -> ZBMI_next：P-nonlinear，RCS非线性效应\n")
print(rcs_dq_to_zbmi$d1_nonlinear)

p_overall_dq_to_zbmi <- rcs_dq_to_zbmi$d1_overall$p.value[1]
p_nonlinear_dq_to_zbmi <- rcs_dq_to_zbmi$d1_nonlinear$p.value[1]

p_rcs_dq_to_zbmi_df <- pool_predictions_lm(
  rcs_result = rcs_dq_to_zbmi,
  exposure = "DQ_now",
  prefix = "DQ",
  adjust_vars = adjust_vars_dq_to_zbmi,
  n_grid = 200
)

cat("\nDQ_now -> ZBMI_next：预测曲线数据前6行\n")
print(head(p_rcs_dq_to_zbmi_df))

peak_dq_to_zbmi <- p_rcs_dq_to_zbmi_df %>%
  slice(which.max(yhat)) %>%
  transmute(
    peak_dq = DQ_now,
    peak_zbmi = yhat
  )

cat("\nDQ_now -> ZBMI_next：预测峰值，仅作图示辅助\n")
print(peak_dq_to_zbmi)

plot_rcs_base(
  pred_df = p_rcs_dq_to_zbmi_df,
  exposure = "DQ_now",
  xlab = "Current DQ",
  ylab = "Predicted ZBMI at next survey",
  main = "Restricted cubic spline analysis of DQ and subsequent ZBMI",
  p_overall = p_overall_dq_to_zbmi,
  p_nonlinear = p_nonlinear_dq_to_zbmi,
  ref_vline = 100,
  fixed_ylim = NULL
)

# Summary
rcs_mi_summary <- tibble(
  analysis = c(
    "ZBMI_now -> DQ_next",
    "ZBMI_now -> DQ_next",
    "DQ_now -> ZBMI_next",
    "DQ_now -> ZBMI_next"
  ),
  test = c(
    "P-overall",
    "P-nonlinear",
    "P-overall",
    "P-nonlinear"
  ),
  statistic = c(
    rcs_zbmi_to_dq$d1_overall$statistic[1],
    rcs_zbmi_to_dq$d1_nonlinear$statistic[1],
    rcs_dq_to_zbmi$d1_overall$statistic[1],
    rcs_dq_to_zbmi$d1_nonlinear$statistic[1]
  ),
  df1 = c(
    rcs_zbmi_to_dq$d1_overall$df1[1],
    rcs_zbmi_to_dq$d1_nonlinear$df1[1],
    rcs_dq_to_zbmi$d1_overall$df1[1],
    rcs_dq_to_zbmi$d1_nonlinear$df1[1]
  ),
  df2 = c(
    rcs_zbmi_to_dq$d1_overall$df2[1],
    rcs_zbmi_to_dq$d1_nonlinear$df2[1],
    rcs_dq_to_zbmi$d1_overall$df2[1],
    rcs_dq_to_zbmi$d1_nonlinear$df2[1]
  ),
  p_value = c(
    p_overall_zbmi_to_dq,
    p_nonlinear_zbmi_to_dq,
    p_overall_dq_to_zbmi,
    p_nonlinear_dq_to_zbmi
  )
) %>%
  mutate(
    sig = add_sig(p_value)
  )

cat("\n==============================\n")
cat("多重插补后 RCS 滞后分析汇总\n")
cat("==============================\n")
print(rcs_mi_summary)
