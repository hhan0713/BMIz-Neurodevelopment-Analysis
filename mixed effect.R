# 清空环境
rm(list = ls())
graphics.off()

# Packages and data
library(tidyverse)
library(readxl)
library(lme4)
library(mice)
library(broom.mixed)
library(ggeffects)

dat_raw <- read_excel("~/Desktop/混合效应模型.xlsx")
names(dat_raw) <- stringr::str_replace_all(names(dat_raw), "\\s+", "_")


# Variable preparation
aux_domain_vars <- names(dat_raw)[
  stringr::str_detect(names(dat_raw), "^(G|F|A|L|S)[1-4]$")
]

longitudinal_vars <- names(dat_raw)[
  stringr::str_detect(names(dat_raw), "^(DQ|ZBMI|T)[1-4]$")
]

covariate_vars <- c(
  "SEX", "MATERNAL_AGE", "GESTATIONAL_WEEK", "NUMBER_OF_BIRTHS",
  "GRAVIDITY", "DELIVERY_METHOD", "BF", "WT500", "LEN5"
)

needed_vars <- intersect(
  c(
    "ID", longitudinal_vars, aux_domain_vars,
    "GroupZBMI", "GroupCNBS", covariate_vars
  ),
  names(dat_raw)
)

dat_imp_base <- dat_raw %>%
  select(all_of(needed_vars)) %>%
  mutate(
    ID = factor(ID),
    across(matches("^(DQ|ZBMI|T|G|F|A|L|S)[1-4]$"), as.numeric),
    MATERNAL_AGE = as.numeric(MATERNAL_AGE),
    WT500 = as.numeric(WT500),
    LEN5 = as.numeric(LEN5),
    across(
      c(
        SEX, GESTATIONAL_WEEK, NUMBER_OF_BIRTHS,
        GRAVIDITY, DELIVERY_METHOD, BF
      ),
      factor
    ),
    GroupZBMI = factor(
      GroupZBMI,
      levels = c(1, 2, 3),
      labels = c(
        "Low-BMI group",
        "Normal-BMI group",
        "High-BMI group"
      )
    ),
    GroupCNBS = factor(
      GroupCNBS,
      levels = c(1, 2, 3),
      labels = c(
        "Low-level declining group",
        "Medium-level rising group",
        "High-level rising group"
      )
    ),
    GroupZBMI = relevel(GroupZBMI, ref = "Normal-BMI group"),
    GroupCNBS = relevel(GroupCNBS, ref = "High-level rising group")
  )

# Missing-data assessment and multiple imputation
missing_overview <- tibble(
  variable = names(dat_imp_base),
  n_missing = colSums(is.na(dat_imp_base)),
  pct_missing = round(colMeans(is.na(dat_imp_base)) * 100, 2)
) %>%
  arrange(desc(pct_missing))

print(missing_overview)
md.pattern(dat_imp_base)

ini <- mice(dat_imp_base, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

meth["ID"] <- ""
pred[, "ID"] <- 0
pred["ID", ] <- 0

for (v in c("GroupZBMI", "GroupCNBS")) {
  if (v %in% names(dat_imp_base)) {
    meth[v] <- ""
    pred[v, ] <- 0
  }
}

numeric_vars <- names(dat_imp_base)[vapply(dat_imp_base, is.numeric, logical(1))]
numeric_vars_missing <- numeric_vars[
  colSums(is.na(dat_imp_base[numeric_vars])) > 0
]
meth[numeric_vars_missing] <- "pmm"

set.seed(2025)
imp <- mice(
  dat_imp_base,
  m = 20,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  seed = 2025,
  printFlag = TRUE
)

plot(imp)

# Long-format data and standardized time
dat_long_observed_for_time <- dat_imp_base %>%
  pivot_longer(
    cols = matches("^(DQ|ZBMI|T)[1-4]$"),
    names_to = c(".value", "wave"),
    names_pattern = "^(DQ|ZBMI|T)([1-4])$"
  ) %>%
  mutate(
    wave = as.integer(wave),
    time_month = as.numeric(T)
  )

time_mean <- mean(
  dat_long_observed_for_time$time_month,
  na.rm = TRUE
)

time_sd <- sd(
  dat_long_observed_for_time$time_month,
  na.rm = TRUE
)

make_long_data <- function(data_wide) {
  data_wide %>%
    pivot_longer(
      cols = matches("^(DQ|ZBMI|T)[1-4]$"),
      names_to = c(".value", "wave"),
      names_pattern = "^(DQ|ZBMI|T)([1-4])$"
    ) %>%
    mutate(
      wave = as.integer(wave),
      time_month = as.numeric(T),
      ID = factor(ID),
      DQ = as.numeric(DQ),
      ZBMI = as.numeric(ZBMI),
      MATERNAL_AGE = as.numeric(MATERNAL_AGE),
      WT500 = as.numeric(WT500),
      LEN5 = as.numeric(LEN5),
      across(
        c(
          SEX, GESTATIONAL_WEEK, NUMBER_OF_BIRTHS,
          GRAVIDITY, DELIVERY_METHOD, BF
        ),
        factor
      ),
      GroupZBMI = factor(
        GroupZBMI,
        levels = c(
          "Low-BMI group",
          "Normal-BMI group",
          "High-BMI group"
        )
      ),
      GroupCNBS = factor(
        GroupCNBS,
        levels = c(
          "Low-level declining group",
          "Medium-level rising group",
          "High-level rising group"
        )
      ),
      GroupZBMI = relevel(GroupZBMI, ref = "Normal-BMI group"),
      GroupCNBS = relevel(GroupCNBS, ref = "High-level rising group"),
      time_z = as.numeric((time_month - time_mean) / time_sd),
      time_z2 = time_z^2,
      zbmi_low = as.integer(GroupZBMI == "Low-BMI group"),
      zbmi_high = as.integer(GroupZBMI == "High-BMI group"),
      zbmi_low_time = zbmi_low * time_z,
      zbmi_high_time = zbmi_high * time_z,
      zbmi_low_time2 = zbmi_low * time_z2,
      zbmi_high_time2 = zbmi_high * time_z2,
      cnbs_low = as.integer(GroupCNBS == "Low-level declining group"),
      cnbs_medium = as.integer(GroupCNBS == "Medium-level rising group"),
      cnbs_low_time = cnbs_low * time_z,
      cnbs_medium_time = cnbs_medium * time_z,
      cnbs_low_time2 = cnbs_low * time_z2,
      cnbs_medium_time2 = cnbs_medium * time_z2
    ) %>%
    select(-T)
}

# Mixed-model fitting function
fit_mixed_model <- function(data_wide, model_formula) {
  dat_long <- make_long_data(data_wide)

  lmer(
    model_formula,
    data = dat_long,
    REML = FALSE,
    na.action = na.omit,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )
}

imp_list <- complete(imp, action = "all")

# DQ outcome: BMI trajectory group as predictor
dq_full_formula <- DQ ~
  GroupZBMI * (time_z + I(time_z^2)) +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

models_dq_full <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = dq_full_formula
)

pool_dq <- pool(as.mira(models_dq_full))
dq_fixed_effects <- summary(pool_dq, conf.int = TRUE)
print(dq_fixed_effects)

dq_test_full_formula <- DQ ~
  zbmi_low + zbmi_high + time_z + time_z2 +
  zbmi_low_time + zbmi_high_time +
  zbmi_low_time2 + zbmi_high_time2 +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

dq_test_no_group_main_formula <- DQ ~
  time_z + time_z2 +
  zbmi_low_time + zbmi_high_time +
  zbmi_low_time2 + zbmi_high_time2 +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

dq_test_no_interaction_formula <- DQ ~
  zbmi_low + zbmi_high + time_z + time_z2 +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

models_dq_test_full <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = dq_test_full_formula
)

models_dq_test_no_group_main <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = dq_test_no_group_main_formula
)

models_dq_test_no_interaction <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = dq_test_no_interaction_formula
)

dq_group_main_test <- D1(
  as.mira(models_dq_test_full),
  as.mira(models_dq_test_no_group_main)
)

dq_interaction_test <- D1(
  as.mira(models_dq_test_full),
  as.mira(models_dq_test_no_interaction)
)

print(dq_group_main_test)
print(dq_interaction_test)
print(sapply(models_dq_full, isSingular, tol = 1e-4))

# ZBMI outcome: DQ trajectory group as predictor
zbmi_full_formula <- ZBMI ~
  GroupCNBS * (time_z + I(time_z^2)) +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

models_zbmi_full <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = zbmi_full_formula
)

pool_zbmi <- pool(as.mira(models_zbmi_full))
zbmi_fixed_effects <- summary(pool_zbmi, conf.int = TRUE)
print(zbmi_fixed_effects)

zbmi_test_full_formula <- ZBMI ~
  cnbs_low + cnbs_medium + time_z + time_z2 +
  cnbs_low_time + cnbs_medium_time +
  cnbs_low_time2 + cnbs_medium_time2 +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

zbmi_test_no_group_main_formula <- ZBMI ~
  time_z + time_z2 +
  cnbs_low_time + cnbs_medium_time +
  cnbs_low_time2 + cnbs_medium_time2 +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

zbmi_test_no_interaction_formula <- ZBMI ~
  cnbs_low + cnbs_medium + time_z + time_z2 +
  SEX + MATERNAL_AGE + GESTATIONAL_WEEK + NUMBER_OF_BIRTHS +
  GRAVIDITY + DELIVERY_METHOD + BF + WT500 + LEN5 +
  (1 + time_z | ID)

models_zbmi_test_full <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = zbmi_test_full_formula
)

models_zbmi_test_no_group_main <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = zbmi_test_no_group_main_formula
)

models_zbmi_test_no_interaction <- lapply(
  imp_list,
  fit_mixed_model,
  model_formula = zbmi_test_no_interaction_formula
)

zbmi_group_main_test <- D1(
  as.mira(models_zbmi_test_full),
  as.mira(models_zbmi_test_no_group_main)
)

zbmi_interaction_test <- D1(
  as.mira(models_zbmi_test_full),
  as.mira(models_zbmi_test_no_interaction)
)

print(zbmi_group_main_test)
print(zbmi_interaction_test)
print(sapply(models_zbmi_full, isSingular, tol = 1e-4))

# Model diagnostics
diagnose_model <- function(model, outcome_label) {
  par(mfrow = c(1, 2))

  qqnorm(
    resid(model),
    main = paste0(outcome_label, ": Normal Q-Q plot")
  )
  qqline(resid(model), col = "red", lwd = 2)

  plot(
    fitted(model),
    resid(model),
    xlab = "Fitted values",
    ylab = "Residuals",
    main = paste0(outcome_label, ": Residuals vs fitted")
  )
  abline(h = 0, lty = 2, col = "red", lwd = 2)

  par(mfrow = c(1, 1))
}

diagnose_model(models_dq_full[[1]], "DQ model")
diagnose_model(models_zbmi_full[[1]], "ZBMI model")

# Pooled marginal predictions
pool_predictions <- function(models, terms, time_mean, time_sd) {
  pred_all <- purrr::map2_dfr(
    models,
    seq_along(models),
    function(model, imputation_number) {
      out <- as.data.frame(
        ggpredict(model, terms = terms)
      )

      if (!"std.error" %in% names(out)) {
        out$std.error <- (
          out$conf.high - out$conf.low
        ) / (2 * 1.96)
      }

      out %>%
        mutate(.imp = imputation_number)
    }
  )

  m <- n_distinct(pred_all$.imp)

  pred_all %>%
    group_by(x, group) %>%
    summarise(
      predicted = mean(predicted, na.rm = TRUE),
      within_var = mean(std.error^2, na.rm = TRUE),
      between_var = var(predicted, na.rm = TRUE),
      between_var = replace_na(between_var, 0),
      total_var = within_var + (1 + 1 / m) * between_var,
      se_total = sqrt(total_var),
      df = if_else(
        between_var <= 1e-12,
        Inf,
        (m - 1) *
          (
            1 +
              within_var /
                ((1 + 1 / m) * between_var)
          )^2
      ),
      conf.low = predicted - qt(0.975, df) * se_total,
      conf.high = predicted + qt(0.975, df) * se_total,
      .groups = "drop"
    ) %>%
    mutate(
      time_month = x * time_sd + time_mean,
      time_month_plot = round(time_month)
    )
}

target_months <- c(12, 18, 24, 36)
target_time_z <- (target_months - time_mean) / time_sd

term_time <- paste0(
  "time_z [",
  paste(round(target_time_z, 6), collapse = ","),
  "]"
)

pred_dq_by_zbmi <- pool_predictions(
  models = models_dq_full,
  terms = c(term_time, "GroupZBMI"),
  time_mean = time_mean,
  time_sd = time_sd
) %>%
  mutate(
    GroupZBMI = factor(
      group,
      levels = c(
        "Low-BMI group",
        "Normal-BMI group",
        "High-BMI group"
      )
    )
  )

pred_zbmi_by_cnbs <- pool_predictions(
  models = models_zbmi_full,
  terms = c(term_time, "GroupCNBS"),
  time_mean = time_mean,
  time_sd = time_sd
) %>%
  mutate(
    GroupCNBS = factor(
      group,
      levels = c(
        "Low-level declining group",
        "Medium-level rising group",
        "High-level rising group"
      )
    )
  )

# Figures
zbmi_colors <- c(
  "Low-BMI group" = "#00A887",
  "Normal-BMI group" = "#8E88C9",
  "High-BMI group" = "#F26B1D"
)

p_dq_by_zbmi <- ggplot(
  pred_dq_by_zbmi,
  aes(
    x = time_month_plot,
    y = predicted,
    color = GroupZBMI,
    fill = GroupZBMI
  )
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.06,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = zbmi_colors) +
  scale_fill_manual(values = zbmi_colors) +
  scale_x_continuous(breaks = target_months) +
  scale_y_continuous(
    breaks = seq(89, 102, by = 2),
    limits = c(89, 102),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Month",
    y = "Developmental quotient (DQ)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.45,
      linetype = "dashed"
    ),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "grey25", size = 12),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(1.2, "cm"),
    plot.margin = margin(12, 15, 10, 12)
  )

cnbs_colors <- c(
  "Low-level declining group" = "#00A887",
  "Medium-level rising group" = "#8E88C9",
  "High-level rising group" = "#F26B1D"
)

p_zbmi_by_cnbs <- ggplot(
  pred_zbmi_by_cnbs,
  aes(
    x = time_month_plot,
    y = predicted,
    color = GroupCNBS,
    fill = GroupCNBS
  )
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.06,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = cnbs_colors) +
  scale_fill_manual(values = cnbs_colors) +
  scale_x_continuous(breaks = target_months) +
  scale_y_continuous(
    breaks = seq(-1.0, 1.5, by = 0.5),
    minor_breaks = NULL
  ) +
  coord_cartesian(ylim = c(-1.0, 1.5)) +
  labs(
    x = "Month",
    y = "ZBMI",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.45,
      linetype = "dashed"
    ),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "grey25", size = 12),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(1.2, "cm"),
    plot.margin = margin(12, 15, 10, 12)
  )

print(p_dq_by_zbmi)
print(p_zbmi_by_cnbs)

# Export results
write.csv(
  dq_fixed_effects,
  "pooled_results_DQ_by_GroupZBMI.csv",
  row.names = FALSE
)

write.csv(
  zbmi_fixed_effects,
  "pooled_results_ZBMI_by_GroupCNBS.csv",
  row.names = FALSE
)

write.csv(
  pred_dq_by_zbmi,
  "pooled_predictions_DQ_by_GroupZBMI.csv",
  row.names = FALSE
)

write.csv(
  pred_zbmi_by_cnbs,
  "pooled_predictions_ZBMI_by_GroupCNBS.csv",
  row.names = FALSE
)

ggsave(
  "figure_DQ_by_GroupZBMI_MI.png",
  p_dq_by_zbmi,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "figure_ZBMI_by_GroupCNBS_MI.png",
  p_zbmi_by_cnbs,
  width = 7,
  height = 5,
  dpi = 300
)
