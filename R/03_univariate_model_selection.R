# ==============================================================================
# 03_univariate_model_selection.R
# Univariate model of envelope extent: specification, inference and figures
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   Asks which forest configuration metrics accompany a larger occupied habitat
#   envelope (OHE), under a candidate set defined a priori. The goal is
#   inference, not prediction.
#
#   1. Six specifications (OLS, linear mixed model and Gamma GLM, each additive
#      or with the a priori interaction) are compared on a common response
#      scale by AICc, by BIC and by exact leave-one-out cross-validation, with
#      leave-one-genus-out beside it. Both criteria are reported because their
#      relative behavior depends on unobserved heterogeneity, which 67 units
#      from nine taxa carry (Brewer et al. 2016).
#   2. The full a priori model is fitted once, without selection. Its
#      coefficients are the inferential result, because inference on
#      coefficients is biased when it follows a selection step (Tredennick et
#      al. 2021; Yates et al. 2023).
#   3. The same model is refitted with the Moran eigenvectors retained by
#      script 02, because an information criterion computed on a non-spatial
#      regression is sensitive to spatial autocorrelation (Diniz-Filho et al.
#      2008).
#   4. The exhaustive search over the candidate subsets is reported as a
#      description of selection uncertainty, not as a second set of estimates.
#
#   The interaction between proximity and forest cover is the only one admitted,
#   because the fragmentation threshold hypothesis is a statement about that
#   interaction rather than about either term alone.
#
#   Three terms of the a priori set, patch density (PD), forest cover (PLAND)
#   and the product of proximity and forest cover, are divided by the window
#   area, which is the response, or built on a term that is. Script 05 shows
#   that their coefficients cannot be interpreted, and the manuscript reports
#   the model of script 05 (Table 4). This script is kept because the a priori
#   analysis is part of the record (Table 3 and the supplementary tables below).
#
# Input
#   Dados/Processados/Data_Raw_WithCoords.csv              written by script 01
#   Outputs/Manuscrito/PartII/TableII_2b_MEM_Vectors.csv   written by script 02
#
# Output   (in Outputs/Manuscrito/; docs/TABLE_MAP.md lists every file)
#   Table_S_PartIII_*.csv   Table 3 and Tables S15 to S25, S29 and S34
#   Figure_PIII_*.tiff      Figure 7 and Figures S2 to S7 and S10
#
# Run time
#   Less than a minute.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)        # project-relative paths
  library(tidyverse)   # data wrangling and plotting
  library(scales)      # axis formatting
  library(patchwork)   # multi-panel figures
  library(MuMIn)       # multimodel inference (dredge, model.avg)
  library(DHARMa)      # simulation-based residual diagnostics
  library(car)         # vif(), ncvTest(), durbinWatsonTest()
  library(lme4)        # linear mixed models, one of the six specifications
  library(performance) # r2() for mixed models
})

# Seed, set again before each stochastic call.
GLOBAL_SEED <- 123

N_UNITS_EXPECTED <- 67   # the design has 67 sampling units; Section 3 checks it

# Minimum number of units a genus must have for a genus-level association to be
# estimated in the heatmaps (Sections 14 and 15). Below it, the cell is shown
# as "n.e." (not estimable).
MIN_N_GENUS <- 4

# Diverging palette of the heatmaps (Okabe-Ito colors): blue for negative,
# near-white at zero, vermillion for positive associations.
heat_low  <- "#0072B2"   # blue       (negative association)
heat_mid  <- "#F7F7F7"   # near-white (no association)
heat_high <- "#D55E00"   # vermillion (positive association)

# Single color of the importance and partial-effect figures.
col_main <- "#0072B2"

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Predictor labels, as in Table 3 of the manuscript.
predictor_labels <- c(
  "PD_Forest"                   = "Patch density (PD)",
  "Mean_FRAC_Forest"            = "Shape complexity (FRAC_MN)",
  "Mean_ENN_Forest_m"           = "Isolation (ENN_MN)",
  "PROX_MN_Forest"              = "Forest proximity (PROX_MN)",
  "PLAND_Forest"                = "Forest cover (PLAND)",
  "PLAND_Forest:PROX_MN_Forest" = "Proximity \u00d7 cover",
  "PROX_MN_Forest:PLAND_Forest" = "Proximity \u00d7 cover"
)

# The same labels with units, in two lines, for the panels of Figure S4.
predictor_labels_units <- c(
  "PD_Forest"         = "Patch density\n(PD, patches/100 ha)",
  "Mean_FRAC_Forest"  = "Shape complexity\n(FRAC_MN)",
  "Mean_ENN_Forest_m" = "Isolation\n(ENN_MN, m)",
  "PROX_MN_Forest"    = "Forest proximity\n(PROX_MN)",
  "PLAND_Forest"      = "Forest cover\n(PLAND, %)"
)

# Remove the standardization suffix from every term of an interaction:
# "PLAND_Forest_z:PROX_MN_Forest_z" becomes "PLAND_Forest:PROX_MN_Forest".
strip_z <- function(x) gsub("_z(?=:|$)", "", x, perl = TRUE)

# Look up a label. When the lookup misses, fall back to the variable name with a
# warning, so that no axis reads "NA".
label_predictors <- function(x, lookup = predictor_labels) {
  out <- unname(lookup[x])
  gap <- is.na(out)
  if (any(gap)) {
    warning("No label defined for: ", paste(unique(x[gap]), collapse = ", "),
            ". Falling back to the variable name.", call. = FALSE)
    out[gap] <- x[gap]
  }
  out
}


# 2. Helper functions ----------------------------------------------------------

# First finite numeric element of x, or NA_real_ when x is NULL, empty or not
# numeric. Tames the different return shapes of the car and DHARMa test objects
# before rounding.
safe_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  v <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(v)) return(NA_real_)
  v
}

# A) Export a ggplot as a TIFF (500 dpi, LZW compression). The path and time of
# every figure written are reported. A failed write, most often a file left open
# in an image viewer, is reported as a warning and does not stop the run.
save_publication_plot <- function(plot_obj, file_name,
                                  width_mm = 140, height_mm = 110) {
  if (is.null(plot_obj)) {
    message(sprintf("   [skip] %s: plot object is NULL.", file_name))
    return(invisible(FALSE))
  }
  full_path <- file.path(output_path, file_name)
  ok <- tryCatch({
    ggsave(filename = file_name, plot = plot_obj, path = output_path,
           width = width_mm, height = height_mm, units = "mm",
           dpi = 500, device = "tiff", compression = "lzw", bg = "white")
    TRUE
  }, error = function(e) {
    message(sprintf(paste0(
      "   [WARN] Could not write %s\n",
      "          Reason: %s\n",
      "          If the file is open in an image viewer (Photos, Preview, etc.), ",
      "close it and re-run: a locked file cannot be overwritten."),
      normalizePath(full_path, mustWork = FALSE), conditionMessage(e)))
    FALSE
  })
  if (isTRUE(ok) && file.exists(full_path)) {
    message(sprintf("   Figure written: %s  (%s)",
                    normalizePath(full_path, mustWork = FALSE),
                    format(file.info(full_path)$mtime, "%Y-%m-%d %H:%M:%S")))
  }
  invisible(ok)
}

# B) Residual diagnostics for OLS, GLM and mixed models. Returns one row per
# test, with a flag indicating whether the assumption is met at alpha = 0.05.
run_residual_diagnostics <- function(model, label = "model", seed = GLOBAL_SEED) {

  # Shapiro-Wilk, the score test of non-constant variance and Durbin-Watson test
  # the Gaussian assumptions of normal residuals with constant variance. A Gamma
  # GLM makes neither assumption, so the three are reported as not applicable
  # outside the Gaussian families, and the simulation-based DHARMa
  # diagnostics, which are defined for any family, carry the verdict.
  is_gaussian <- inherits(model, "lm") && !inherits(model, "glm") ||
                 inherits(model, "lmerMod")

  # Shapiro-Wilk on raw residuals.
  sw <- if (is_gaussian) {
    tryCatch(shapiro.test(residuals(model)),
             error = function(e) list(statistic = NA_real_, p.value = NA_real_))
  } else list(statistic = NA_real_, p.value = NA_real_)

  # Breusch-Pagan / non-constant variance via car::ncvTest.
  bp <- if (is_gaussian) {
    tryCatch(car::ncvTest(model),
             error = function(e) list(ChiSquare = NA_real_, p = NA_real_))
  } else list(ChiSquare = NA_real_, p = NA_real_)

  # Durbin-Watson for residual serial dependence.
  dw <- if (is_gaussian) {
    tryCatch(car::durbinWatsonTest(model, reps = 999),
             error = function(e) list(dw = NA_real_, p = NA_real_))
  } else list(dw = NA_real_, p = NA_real_)

  # DHARMa simulation-based diagnostics.
  set.seed(seed)
  sim  <- DHARMa::simulateResiduals(fittedModel = model, n = 1000, seed = seed)
  disp <- DHARMa::testDispersion(sim, plot = FALSE)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  outl <- DHARMa::testOutliers(sim, plot = FALSE)

  diag_df <- data.frame(
    Model     = label,
    Test      = c("Shapiro-Wilk (normality)",
                  "Breusch-Pagan (homoscedasticity)",
                  "Durbin-Watson (independence)",
                  "DHARMa dispersion",
                  "DHARMa uniformity",
                  "DHARMa outliers"),
    Statistic = round(c(safe_num(sw$statistic),
                        safe_num(bp$ChiSquare),
                        safe_num(dw$dw),
                        safe_num(disp$statistic),
                        safe_num(unif$statistic),
                        safe_num(outl$statistic)), 4),
    P_value   = round(c(safe_num(sw$p.value),
                        safe_num(bp$p),
                        safe_num(dw$p),
                        safe_num(disp$p.value),
                        safe_num(unif$p.value),
                        safe_num(outl$p.value)), 4),
    Assumption_met = NA,
    stringsAsFactors = FALSE
  )

  # A test that does not apply to the family is reported as such, never as a
  # failed assumption.
  diag_df$Assumption_met <- ifelse(is.na(diag_df$P_value), "not applicable",
                                   ifelse(diag_df$P_value > 0.05, "TRUE", "FALSE"))
  diag_df
}

# C) Coefficients with 95% confidence intervals, for lm, glm, merMod and
# model.avg objects.
extract_coefs <- function(model, label = "final") {
  if (inherits(model, "averaging")) {
    # Conditional (subset) averages: each parameter is averaged only across the
    # models in which it occurs (Burnham & Anderson 2002).
    cm <- summary(model)$coefmat.subset
    ci <- confint(model, full = FALSE)
    df <- data.frame(
      Predictor = rownames(cm),
      Estimate  = cm[, "Estimate"],
      SE        = cm[, "Std. Error"],
      CI_lower  = ci[, 1],
      CI_upper  = ci[, 2],
      Z_value   = cm[, "z value"],
      P_value   = cm[, "Pr(>|z|)"],
      Source    = label,
      stringsAsFactors = FALSE
    )
  } else if (inherits(model, "merMod")) {
    # lme4 returns no p-value, and confint() on a merMod adds rows for the
    # random-effect standard deviations. Only the fixed effects are extracted,
    # and the p-value is the two-sided Wald test with the normal approximation.
    coefs <- summary(model)$coefficients
    ci    <- suppressMessages(
      confint(model, parm = "beta_", method = "Wald", oldNames = FALSE))
    ci    <- ci[rownames(coefs), , drop = FALSE]
    tval  <- coefs[, 3]
    df <- data.frame(
      Predictor = rownames(coefs),
      Estimate  = coefs[, 1],
      SE        = coefs[, 2],
      CI_lower  = ci[, 1],
      CI_upper  = ci[, 2],
      Z_value   = tval,
      P_value   = 2 * stats::pnorm(abs(tval), lower.tail = FALSE),
      Source    = label,
      stringsAsFactors = FALSE
    )
  } else {
    coefs <- summary(model)$coefficients
    ci    <- suppressMessages(confint(model))
    df <- data.frame(
      Predictor = rownames(coefs),
      Estimate  = coefs[, 1],
      SE        = coefs[, 2],
      CI_lower  = ci[, 1],
      CI_upper  = ci[, 2],
      Z_value   = coefs[, 3],
      P_value   = coefs[, 4],
      Source    = label,
      stringsAsFactors = FALSE
    )
  }
  rownames(df) <- NULL
  df
}

# D) Explained variance: deviance-based pseudo-R2 for a GLM (Cohen et al. 2003),
# adjusted R2 for an OLS, and marginal R2 for a mixed model.
explained_variance <- function(model) {
  if (inherits(model, "merMod")) {
    # Marginal R2 (Nakagawa & Schielzeth 2013): the share explained by the
    # fixed effects alone, comparable to the adjusted R2 of the OLS.
    r2 <- tryCatch(performance::r2(model)$R2_marginal,
                   error = function(e) NA_real_)
    as.numeric(r2)
  } else if (inherits(model, "glm")) {
    1 - (model$deviance / model$null.deviance)
  } else if (inherits(model, "lm")) {
    summary(model)$adj.r.squared
  } else {
    NA_real_
  }
}


# 3. Data ----------------------------------------------------------------------
message("\n[Data] Loading the dataset and preparing the predictors...")

# The five forest predictors of the a priori model.
predictors_part3 <- c("PD_Forest", "Mean_FRAC_Forest",
                      "Mean_ENN_Forest_m", "PROX_MN_Forest",
                      "PLAND_Forest")

# The a priori interaction. The fragmentation threshold hypothesis states that
# the spatial arrangement of habitat matters where habitat is scarce and ceases
# to matter where it is abundant (Andren 1994; Pardini et al. 2010; Villard &
# Metzger 2014). That is a statement about an interaction between connectivity
# and habitat amount, so the product of proximity and forest cover is admitted
# to the candidate set. It is the only interaction considered: with 67 units, a
# full factorial of the configuration metrics would leave too few observations
# per parameter. Predictors are standardized before the product is formed,
# which keeps the main-effect coefficients interpretable at the mean and removes
# the non-essential collinearity that products introduce (Aiken & West 1991;
# Schielzeth 2010).

# Same function as in scripts 00, 01 and 02; keep the four copies identical.
canonical_taxon <- function(x) {
  key <- gsub("[^a-z]", "", tolower(as.character(x)))
  out <- dplyr::case_when(
    key %in% c("aguariba", "aguriba", "alouattaguariba")           ~ "A_guariba",
    key %in% c("barach", "baracnoides", "barachnoides",
               "brachytelesarachnoides")                           ~ "B_arachnoides",
    key %in% c("bcrinitus", "btorquatus", "bradypuscrinitus")       ~ "B_crinitus",
    key %in% c("bhypox", "bhypoxantus", "bhypoxanthus",
               "brachyteleshypoxanthus")                           ~ "B_hypoxanthus",
    key %in% c("lwiedii", "leoparduswiedii")                        ~ "L_wiedii",
    key %in% c("mazama", "mazamasp")                                ~ "Mazama",
    key %in% c("mtrida", "mtridactyla", "myrmecophagatridactyla")   ~ "M_tridactyla",
    key %in% c("pconcolor", "concolor", "pumaconcolor")             ~ "P_concolor",
    key %in% c("tpecari", "tayassupecari")                          ~ "T_pecari",
    TRUE ~ NA_character_
  )
  if (any(is.na(out) & !is.na(x))) {
    stop("Unrecognized taxon name(s): ",
         paste(unique(x[is.na(out) & !is.na(x)]), collapse = ", "))
  }
  out
}

REQUIRED_INPUT <- "Data_Raw_WithCoords.csv"
primary_file   <- file.path(data_path, REQUIRED_INPUT)
if (!file.exists(primary_file)) {
  stop(sprintf(paste0("%s not found in %s. It is written by script 01. Run ",
                      "scripts 00 and 01 before this one."),
               REQUIRED_INPUT, data_path))
}
message(sprintf("   Data file: %s", basename(primary_file)))

# AOH_ha is the area of the occupied habitat envelope (OHE), in hectares.
data_raw <- read_csv(primary_file, show_col_types = FALSE) %>%
  rename(AOH_ha = AOH_unit_area_ha)

# Canonical taxon labels, so that SPECIES matches scripts 00, 01 and 02.
if ("SPECIES" %in% names(data_raw)) {
  data_raw$SPECIES <- canonical_taxon(data_raw$SPECIES)
}

# Categorical columns as factors.
data_raw <- data_raw %>%
  mutate(across(any_of(c("ORDER", "GENUS", "SPECIES",
                         "DIET", "LOCOMOTION", "UA_ID")), as.factor))

missing_preds <- setdiff(predictors_part3, names(data_raw))
if (length(missing_preds) > 0) {
  stop(sprintf(paste0("The following predictor(s) are missing from %s: %s. ",
                      "Rerun script 01."),
               REQUIRED_INPUT, paste(missing_preds, collapse = ", ")))
}

# Working data. An area of zero would be replaced by 0.001 ha so that the
# logarithm is defined; no unit has an area of zero, so this is a guard only.
# The same guard appears in the matrix sensitivity model (Section 10).
data_p3 <- data_raw %>%
  select(UA_ID, GENUS, SPECIES, AOH_ha, all_of(predictors_part3)) %>%
  drop_na(all_of(predictors_part3)) %>%
  mutate(AOH_ha   = ifelse(AOH_ha == 0, 0.001, AOH_ha),
         log_AOH  = log(AOH_ha))

# z-scores, so that the coefficients are directly comparable.
data_p3 <- data_p3 %>%
  mutate(across(all_of(predictors_part3),
                ~as.numeric(scale(.x)),
                .names = "{.col}_z"))

# UA_ID numbers the units within each taxon, so a sampling unit is identified by
# the pair (SPECIES, UA_ID).
n_obs <- nrow(data_p3)
n_su  <- nrow(distinct(data_p3, SPECIES, UA_ID))
message(sprintf("   Sample size: %d sampling units.", n_su))

# Guard: drop_na() above must remove nothing. A smaller sample means that the
# input was not written by script 01, or that script 00 was not run.
if (n_obs != N_UNITS_EXPECTED) {
  lost <- data_raw %>%
    select(SPECIES, UA_ID, all_of(predictors_part3)) %>%
    filter(!complete.cases(.)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value))
  if (nrow(lost) > 0) print(as.data.frame(lost), row.names = FALSE)
  stop(sprintf(paste0("Script 03 is running on %d units instead of %d. ",
                      "Run scripts 00 and 01, and make sure ",
                      "Data_Raw_WithCoords.csv is the file being read."),
               n_obs, N_UNITS_EXPECTED))
}

# Descriptive statistics (Table S15).
desc_stats <- data_p3 %>%
  select(AOH_ha, log_AOH, all_of(predictors_part3)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(N      = sum(!is.na(Value)),
            Mean   = round(mean(Value, na.rm = TRUE), 4),
            SD     = round(sd(Value, na.rm = TRUE),   4),
            Median = round(median(Value, na.rm = TRUE), 4),
            Min    = round(min(Value, na.rm = TRUE),   4),
            Max    = round(max(Value, na.rm = TRUE),   4),
            .groups = "drop") %>%
  # Exported labels use OHE, the term of the manuscript; internal column names
  # are unchanged.
  mutate(Variable = ifelse(Variable == "AOH_ha",  "OHE_ha",
                    ifelse(Variable == "log_AOH", "log_OHE", Variable)))
write.csv(desc_stats,
          file.path(output_path, "Table_S_PartIII_DescriptiveStatistics.csv"),
          row.names = FALSE)
message("   Descriptive statistics exported.")


# 4. Figure S2: the response before and after the log transformation ----------
message("\n[Figure S2] The response before and after the log transformation...")

plot_resp_raw <- ggplot(data_p3, aes(x = AOH_ha)) +
  geom_histogram(bins = 30, fill = "grey60", color = "black", linewidth = 0.3) +
  scale_x_continuous(labels = label_number(scale = 1e-3, big.mark = ",")) +
  labs(x = "OHE (thousand ha)", y = "Number of sampling units") +
  theme_classic(base_size = 9)

plot_resp_log <- ggplot(data_p3, aes(x = log_AOH)) +
  geom_histogram(bins = 30, fill = "#009E73", color = "black", linewidth = 0.3) +
  labs(x = "log(OHE)", y = "Number of sampling units") +
  theme_classic(base_size = 9)

plot_response_inspection <- plot_resp_raw + plot_resp_log +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 9))
save_publication_plot(plot_response_inspection,
                      "Figure_PIII_01_Response_Transformation.tiff",
                      width_mm = 140, height_mm = 65)


# 5. Choosing the model: six specifications, three criteria --------------------
# The family is chosen by an explicit comparison of the six specifications, not
# by a chain of rules that fits an OLS, tests its assumptions and falls back to
# another family when a test fails. Such a chain never compares the families on
# equal terms and lets one marginal diagnostic settle a structural question.
#
# A common response scale. A model fitted to log(OHE) and a model fitted to OHE
# have likelihoods defined for different response variables, so their
# information criteria are not comparable as they stand. The likelihood of a
# log-scale model is placed on the response scale by subtracting the Jacobian of
# the transformation, sum(log(y_i)). Without that correction the log-scale
# models appear about 1400 AICc units better, an artifact of the scale.
#
# The three criteria, in the order in which they are applied:
#   1. Residual diagnostics. A specification whose simulated residuals fail is
#      not considered, whatever its AICc. This is a gate, not a tie-breaker.
#   2. AICc weight on the common response scale.
#   3. Out-of-sample error in hectares, by exact leave-one-out
#      cross-validation, which settles ties that the AICc leaves open.
# Behind all three stands the nature of the response: the OHE is continuous,
# strictly positive and right-skewed, and its dispersion grows with its mean.
message("\n[Model choice] Comparing six specifications on a common response scale...")

options(na.action = "na.fail")   # required by MuMIn::dredge later on

ADD_RHS <- "PD_Forest_z + Mean_FRAC_Forest_z + Mean_ENN_Forest_m_z + PROX_MN_Forest_z + PLAND_Forest_z"
INT_RHS <- paste(ADD_RHS, "+ PROX_MN_Forest_z:PLAND_Forest_z")
JAC     <- sum(log(data_p3$AOH_ha))
n_obs   <- nrow(data_p3)
message(sprintf("   n = %d | Jacobian sum(log(AOH)) = %.2f", n_obs, JAC))

aicc_manual <- function(ll, k, n) -2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1)

# BIC is reported alongside AICc, not in its place. Brewer, Butler and Cooksley
# (2016) show that the relative performance of AIC, AICc and BIC depends on the
# unobserved heterogeneity in the data, and these data carry heterogeneity that
# no predictor captures: 67 units from 9 taxa. Where the two criteria agree, the
# choice is robust to that dependence; where they disagree, the disagreement is
# reported rather than resolved by preference.
bic_manual <- function(ll, k, n) -2 * ll + k * log(n)

fit_spec <- function(family_name, rhs) {
  if (family_name == "OLS")   return(lm(as.formula(paste("log_AOH ~", rhs)), data = data_p3))
  if (family_name == "LMM")   return(lme4::lmer(as.formula(paste("log_AOH ~", rhs, "+ (1 | GENUS)")),
                                                data = data_p3, REML = FALSE))
  if (family_name == "Gamma") return(glm(as.formula(paste("AOH_ha ~", rhs)),
                                         data = data_p3, family = Gamma(link = "log")))
}

specs <- expand.grid(
  family_name = c("OLS", "LMM", "Gamma"),
  structure   = c("additive", "interaction"),
  stringsAsFactors = FALSE
)

models <- list(); rows <- list()
for (i in seq_len(nrow(specs))) {
  fam <- specs$family_name[i]; est <- specs$structure[i]
  rhs <- if (est == "additive") ADD_RHS else INT_RHS
  m   <- fit_spec(fam, rhs)
  spec_label <- sprintf("%s, %s", fam, est)
  models[[spec_label]] <- m

  # number of estimated parameters, including the variance components
  k <- switch(fam,
              OLS   = length(coef(m)) + 1,
              LMM   = length(lme4::fixef(m)) + 2,
              Gamma = length(coef(m)) + 1)
  # likelihood on the response scale
  ll <- if (fam == "Gamma") as.numeric(logLik(m)) else as.numeric(logLik(m)) - JAC

  # residual diagnostics by simulation, which apply to every family alike
  set.seed(GLOBAL_SEED)
  sim <- tryCatch(DHARMa::simulateResiduals(m, n = 1000, plot = FALSE),
                  error = function(e) NULL)
  if (is.null(sim)) {
    p_unif <- NA_real_; p_disp <- NA_real_; p_out <- NA_real_
  } else {
    p_unif <- DHARMa::testUniformity(sim, plot = FALSE)$p.value
    p_disp <- DHARMa::testDispersion(sim, plot = FALSE)$p.value
    p_out  <- DHARMa::testOutliers(sim, plot = FALSE)$p.value
  }
  rows[[spec_label]] <- data.frame(
    Model = spec_label, Family = fam, Structure = est, k = k,
    logLik_response = round(ll, 3),
    AICc = round(aicc_manual(ll, k, n_obs), 3),
    BIC  = round(bic_manual(ll, k, n_obs), 3),
    DHARMa_uniformity = round(p_unif, 4),
    DHARMa_dispersion = round(p_disp, 4),
    DHARMa_outliers   = round(p_out, 4),
    stringsAsFactors = FALSE
  )
}
tab_comp <- do.call(rbind, rows)
tab_comp$Delta_AICc <- round(tab_comp$AICc - min(tab_comp$AICc), 3)
tab_comp$Delta_BIC  <- round(tab_comp$BIC  - min(tab_comp$BIC),  3)
tab_comp$Weight     <- round(exp(-0.5 * tab_comp$Delta_AICc) /
                               sum(exp(-0.5 * tab_comp$Delta_AICc)), 4)
# a specification passes the gate when no simulation-based test rejects
tab_comp$Residuals_OK <- with(tab_comp,
  (is.na(DHARMa_uniformity) | DHARMa_uniformity > 0.05) &
  (is.na(DHARMa_dispersion) | DHARMa_dispersion > 0.05) &
  (is.na(DHARMa_outliers)   | DHARMa_outliers   > 0.05))


# 5.1 Out-of-sample error, in hectares -----------------------------------------
# Computed only for the families whose predictions are defined on the hectare
# scale without conditioning on a random effect, which is what a comparison of
# predictive accuracy across families requires.
#
# Leave-one-out is used rather than random k-fold: Yates, Aandahl, Richards and
# Brook (2023) recommend exact or approximate leave-one-out to minimize bias, or
# k-fold with a bias correction when k < 10. Leave-one-out is exact here because
# the models are cheap and n is 67.
#
# A second scheme, leave-one-genus-out, is reported beside it. The 67 units come
# from 9 taxa and units of the same taxon share species-level traits, so a random
# split places relatives in the training and the test set at once. Comparing the
# two schemes shows how much that matters.
predict_ha <- function(family_name, rhs, tr, te) {
  if (family_name == "Gamma") {
    predict(glm(as.formula(paste("AOH_ha ~", rhs)), data = tr,
                family = Gamma(link = "log")), newdata = te, type = "response")
  } else {
    mm <- lm(as.formula(paste("log_AOH ~", rhs)), data = tr)
    # log-normal back-transformation, so that the two families are compared on
    # the hectare scale rather than on scales that are not the same quantity
    exp(predict(mm, newdata = te) + summary(mm)$sigma^2 / 2)
  }
}

# Both functions return the root mean squared error (RMSE) on the hectare scale.
# The error distribution is dominated by the largest units, so the RMSE differs
# from the mean absolute error here; report these values as RMSE.
cv_rmse_loo <- function(family_name, rhs) {
  err <- numeric(0)
  for (i in seq_len(n_obs)) {
    tr <- data_p3[-i, ]; te <- data_p3[i, , drop = FALSE]
    p  <- tryCatch(predict_ha(family_name, rhs, tr, te), error = function(e) NA_real_)
    err <- c(err, (te$AOH_ha - p)^2)
  }
  sqrt(mean(err, na.rm = TRUE))
}

cv_rmse_lgo <- function(family_name, rhs) {
  err <- numeric(0)
  for (g in unique(data_p3$GENUS)) {
    tr <- data_p3[data_p3$GENUS != g, ]; te <- data_p3[data_p3$GENUS == g, ]
    p  <- tryCatch(predict_ha(family_name, rhs, tr, te),
                   error = function(e) rep(NA_real_, nrow(te)))
    err <- c(err, (te$AOH_ha - p)^2)
  }
  sqrt(mean(err, na.rm = TRUE))
}

tab_comp$CV_RMSE_LOO_ha <- NA_real_
tab_comp$CV_RMSE_LGO_ha <- NA_real_
for (spec_label in tab_comp$Model) {
  fam <- tab_comp$Family[tab_comp$Model == spec_label]
  est <- tab_comp$Structure[tab_comp$Model == spec_label]
  if (fam == "LMM") next
  rhs <- if (est == "additive") ADD_RHS else INT_RHS
  tab_comp$CV_RMSE_LOO_ha[tab_comp$Model == spec_label] <- round(cv_rmse_loo(fam, rhs), 1)
  tab_comp$CV_RMSE_LGO_ha[tab_comp$Model == spec_label] <- round(cv_rmse_lgo(fam, rhs), 1)
}
# CV_RMSE_ha (leave-one-out) is the error used by the tie-break of Section 5.2.
tab_comp$CV_RMSE_ha <- tab_comp$CV_RMSE_LOO_ha

tab_comp <- tab_comp[order(tab_comp$AICc), ]
write.csv(tab_comp, file.path(output_path, "Table_S_PartIII_ModelComparison.csv"),
          row.names = FALSE)
message("\n   Comparison of the six specifications:")
print(tab_comp[, c("Model", "k", "AICc", "Delta_AICc", "Weight", "Delta_BIC",
                   "Residuals_OK", "CV_RMSE_LOO_ha", "CV_RMSE_LGO_ha")], row.names = FALSE)

if (tab_comp$Model[which.min(tab_comp$AICc)] != tab_comp$Model[which.min(tab_comp$BIC)]) {
  message("   NOTE. AICc and BIC point to different specifications:")
  message(sprintf("         AICc favors %s; BIC favors %s.",
                  tab_comp$Model[which.min(tab_comp$AICc)],
                  tab_comp$Model[which.min(tab_comp$BIC)]))
  message("         Report the disagreement; it reflects the unobserved heterogeneity")
  message("         among taxa (Brewer et al. 2016) and is not resolved by preference.")
} else {
  message("   AICc and BIC agree on the same specification.")
}


# 5.2 The decision -------------------------------------------------------------
# Applied in the stated order. Among the specifications whose residuals pass,
# those within two AICc units of the best are treated as indistinguishable by
# that criterion (Burnham & Anderson 2002); the tie is then settled by the
# out-of-sample error and, where none is available, by parsimony.
eligible <- tab_comp[tab_comp$Residuals_OK, , drop = FALSE]
if (nrow(eligible) == 0) {
  stop("No specification passed the residual diagnostics. Inspect Table_S_PartIII_ModelComparison.csv.")
}
tied <- eligible[eligible$Delta_AICc <= 2, , drop = FALSE]
if (all(is.na(tied$CV_RMSE_ha))) {
  chosen    <- tied$Model[which.min(tied$k)]
  criterion <- "fewest parameters among the models tied on AICc"
} else {
  chosen    <- tied$Model[which.min(tied$CV_RMSE_ha)]
  criterion <- "lowest cross-validated error among the models tied on AICc"
}
global_model       <- models[[chosen]]
model_family_label <- chosen

message(sprintf("\n[Decision] Retained specification: %s", chosen))
message(sprintf("            Criterion: %s.", criterion))
message(sprintf("            Delta AICc = %.2f, weight = %.3f, k = %d.",
                tab_comp$Delta_AICc[tab_comp$Model == chosen],
                tab_comp$Weight[tab_comp$Model == chosen],
                tab_comp$k[tab_comp$Model == chosen]))
message("            The structure of the response is the standing constraint:")
message("            the OHE is continuous, strictly positive and right-skewed,")
message("            and its dispersion grows with its mean.")

formula_global <- formula(global_model)

# Collinearity within the a priori set (Table S18). The VIF is computed on the
# linear predictor, the quantity it refers to in any of the three families.
# Values above 5 for the terms that make up the product are expected and are not
# an artifact of the product itself, since the predictors were standardized
# before it was formed (Aiken & West 1991; Schielzeth 2010).
vif_partIII <- car::vif(lm(as.formula(paste("log_AOH ~", INT_RHS)), data = data_p3))
write.csv(data.frame(Variable = names(vif_partIII), VIF = round(vif_partIII, 3)),
          file.path(output_path, "Table_S_PartIII_VIF.csv"), row.names = FALSE)
message(sprintf("   Maximum VIF among the terms of the retained model: %.2f", max(vif_partIII)))

# Diagnostics of the retained model (Table S17).
diag_final <- run_residual_diagnostics(global_model, label = model_family_label)
write.csv(diag_final,
          file.path(output_path, "Table_S_PartIII_Diagnostics.csv"), row.names = FALSE)


# 5.3 Taxonomic dependence, reported rather than decided by --------------------
# The 67 units come from 9 taxa, so units of the same taxon are not independent.
# The mixed specification, with a random intercept for genus, is in the
# comparison above and is not favored by it; the intraclass correlation is
# reported for completeness, not as a threshold that decides anything.
lmm_ref <- models[["LMM, interaction"]]
vc      <- as.data.frame(lme4::VarCorr(lmm_ref))
icc_su  <- vc$vcov[vc$grp == "GENUS"] / sum(vc$vcov)
write.csv(
  data.frame(Test  = c("ICC (genus), from the mixed specification",
                       "Delta AICc, retained model minus mixed model",
                       "Weight of the mixed specification"),
             Value = round(c(icc_su,
                             tab_comp$AICc[tab_comp$Model == chosen] -
                               tab_comp$AICc[tab_comp$Model == "LMM, interaction"],
                             tab_comp$Weight[tab_comp$Model == "LMM, interaction"]), 4)),
  file.path(output_path, "Table_S_PartIII_DependenceDiagnostic.csv"), row.names = FALSE)
message(sprintf("   ICC of the genus in the mixed specification: %.3f", icc_su))

# Summary of the retained model (one row).
glob_summary <- data.frame(
  Family            = model_family_label,
  N_obs             = n_obs,
  N_predictors      = length(all.vars(formula_global)) - 1,
  R2_or_DevExpl     = round(explained_variance(global_model), 4),
  AICc              = round(tab_comp$AICc[tab_comp$Model == chosen], 3),
  stringsAsFactors  = FALSE
)
write.csv(glob_summary,
          file.path(output_path, "Table_S_PartIII_GlobalModelSummary.csv"),
          row.names = FALSE)
message(sprintf("   Explained variance / deviance of the retained model: %.3f",
                glob_summary$R2_or_DevExpl))


# 6. The confirmatory model: one specification, fitted once (Table 3) ----------
# This is the inferential result of the script. Yates et al. (2023) show that
# inference on parameter estimates is biased when it is preceded by model
# selection, and that valid inference requires either a carefully specified
# single model or post-selection adjustments. Tredennick et al. (2021) make the
# same distinction by goal: what suits prediction does not suit inference.
#
# The goal here is inference, and the candidate set was defined a priori. The
# full a priori model is therefore fitted once, without selection, and its
# coefficients carry the confidence intervals and p-values. Section 8 explores
# the subsets as a description of selection uncertainty only. The family was
# chosen in Section 5 by residual behavior and out-of-sample error, properties
# of the specification rather than of any coefficient.
message("\n[Confirmatory] Fitting the full a priori model, without selection...")

confirm_model <- global_model
cf_conf <- summary(confirm_model)$coefficients
ci_conf <- confint.default(confirm_model)

confirm_table <- data.frame(
  Predictor = rownames(cf_conf),
  Estimate  = round(cf_conf[, 1], 4),
  SE        = round(cf_conf[, 2], 4),
  CI_lower  = round(ci_conf[, 1], 4),
  CI_upper  = round(ci_conf[, 2], 4),
  Statistic = round(cf_conf[, 3], 4),
  P_value   = signif(cf_conf[, 4], 4),
  Source    = "Full a priori model, no selection",
  row.names = NULL, stringsAsFactors = FALSE
)
write.csv(confirm_table,
          file.path(output_path, "Table_S_PartIII_ConfirmatoryModel.csv"),
          row.names = FALSE)
message("   Coefficients of the confirmatory model:")
print(confirm_table[, c("Predictor", "Estimate", "SE", "CI_lower", "CI_upper", "P_value")],
      row.names = FALSE)


# 7. Spatial sensitivity of the confirmatory model (Table S34) -----------------
# Diniz-Filho, Rangel and Bini (2008) show that an information criterion computed
# on a non-spatial regression is sensitive to spatial autocorrelation and yields
# unstable, overfitted minimum adequate models. Script 02 detected spatial
# structure in these data and retained Moran eigenvectors. The confirmatory
# model is refitted with those eigenvectors as covariates, and the coefficients
# are compared.
message("\n[Spatial] Refitting the confirmatory model with the retained Moran eigenvectors...")

mem_file <- file.path(output_path, "PartII", "TableII_2b_MEM_Vectors.csv")

if (file.exists(mem_file)) {
  mem_df <- read_csv(mem_file, show_col_types = FALSE) %>%
    mutate(UA_ID = as.character(UA_ID), SPECIES = as.character(SPECIES))
  mem_cols <- setdiff(names(mem_df), c("SPECIES", "UA_ID"))

  d_sp <- data_p3 %>%
    mutate(UA_ID = as.character(UA_ID), SPECIES = as.character(SPECIES)) %>%
    left_join(mem_df, by = c("SPECIES", "UA_ID"))

  if (anyNA(d_sp[, mem_cols])) {
    message("   Eigenvectors did not join to every unit; the check is skipped.")
  } else {
    rhs_sp  <- paste(INT_RHS, "+", paste(mem_cols, collapse = " + "))
    m_sp    <- glm(as.formula(paste("AOH_ha ~", rhs_sp)), data = d_sp,
                   family = Gamma(link = "log"))
    cf_sp   <- summary(m_sp)$coefficients
    shared  <- intersect(rownames(cf_conf), rownames(cf_sp))

    spatial_table <- data.frame(
      Predictor          = shared,
      Beta_non_spatial   = round(cf_conf[shared, 1], 4),
      Beta_with_MEMs     = round(cf_sp[shared, 1], 4),
      Absolute_change    = round(abs(cf_sp[shared, 1] - cf_conf[shared, 1]), 4),
      P_non_spatial      = signif(cf_conf[shared, 4], 4),
      P_with_MEMs        = signif(cf_sp[shared, 4], 4),
      stringsAsFactors   = FALSE, row.names = NULL
    )
    write.csv(spatial_table,
              file.path(output_path, "Table_S_PartIII_SpatialSensitivity.csv"),
              row.names = FALSE)
    print(spatial_table, row.names = FALSE)

    flips <- with(spatial_table,
                  sign(Beta_non_spatial) != sign(Beta_with_MEMs) &
                    Predictor != "(Intercept)")
    if (any(flips)) {
      message("   WARNING. At least one coefficient changes sign when the spatial")
      message("   structure is included. Report the spatial specification.")
    } else {
      message("   No coefficient changes sign when the spatial structure is included.")
    }
  }
} else {
  message("   TableII_2b_MEM_Vectors.csv not found. Run script 02 first; without it")
  message("   the spatial sensitivity of Diniz-Filho et al. (2008) cannot be checked.")
}


# 8. Exploratory multimodel analysis (Tables S18b, S19 and S20) ----------------
# The candidate subsets of the retained specification are ranked by AICc with
# MuMIn::dredge(). dredge() keeps an interaction only in models that contain both
# of its main effects, and no other rule is added, so forest cover also enters
# models without forest proximity. The candidate set therefore has 8 x 5 = 40
# models: the 8 subsets of PD, FRAC and ENN, crossed with none, PROX, PLAND,
# PROX + PLAND, and PROX + PLAND + their product. Proximity and cover each enter
# 24 of the 40 models, so their sums of Akaike weights (Table S20) are taken
# over the same number of models (Burnham & Anderson 2002).
#
# The competitive set is defined by Delta AICc < 2 (Burnham & Anderson 2002). On
# these data it holds two models, which differ only by Mean_ENN_Forest_m. When
# two or more models are competitive, their coefficients are averaged.
#
# This section describes selection uncertainty across the a priori candidate
# set. It is not the inferential result: its coefficients come after a search
# over 40 subsets, so their standard errors and p-values are conditional on a
# selection step and are biased as inference (Yates et al. 2023). The estimates
# that carry inference are those of Section 6.
#
# Both averaging schemes are reported, because they answer different questions.
# The subset, or conditional, average takes each parameter over the models that
# contain it and describes the effect where it appears (Grueber et al. 2011).
# The full average enters a zero for the models that omit the parameter and is
# the more conservative summary of its support across the whole set.
message("\n[Exploratory] AICc-based multimodel analysis of the candidate set...")

set.seed(GLOBAL_SEED)
dredge_table <- MuMIn::dredge(global_model, rank = "AICc", trace = FALSE)
message(sprintf("   Candidate subsets ranked: %d.", nrow(dredge_table)))

# Full dredge table (Table S19).
dredge_df <- as.data.frame(dredge_table) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))
write.csv(dredge_df,
          file.path(output_path, "Table_S_PartIII_Dredge_FullTable.csv"),
          row.names = FALSE)

# Competitive set (Delta AICc < 2; Burnham & Anderson 2002).
competitive <- MuMIn::get.models(dredge_table, subset = delta < 2)
n_competitive <- length(competitive)
message(sprintf("   Competitive set (Delta AICc < 2): %d model(s).",
                n_competitive))

# A single competitive model is used as it is; two or more are averaged, and
# extract_coefs() reports the conditional (subset) average.
if (n_competitive == 1) {
  final_model <- competitive[[1]]
  averaged    <- FALSE
  inf_label   <- "Single most-parsimonious model retained"
} else {
  final_model <- MuMIn::model.avg(competitive, fit = TRUE)
  averaged    <- TRUE
  inf_label   <- sprintf(
    "Conditional (subset) average across %d competitive models (Delta AICc < 2)",
    n_competitive
  )
}
message(sprintf("   Inference: %s.", inf_label))

# Variable importance: the sum of Akaike weights across the full dredge
# (Table S20).
sw_vec <- MuMIn::sw(dredge_table)
df_importance <- data.frame(
  Predictor   = names(sw_vec),
  Importance  = round(as.numeric(sw_vec), 4),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Importance))
write.csv(df_importance,
          file.path(output_path, "Table_S_PartIII_VariableImportance.csv"),
          row.names = FALSE)

# Coefficient table of the final model (single or averaged).
coef_df <- extract_coefs(final_model, label = inf_label) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

# Full (zero-method) average beside the conditional one (Table S18b), with the
# revised variance estimator of Burnham and Anderson (2004), which MuMIn
# implements as revised.var and which widens the interval to acknowledge
# selection uncertainty.
if (averaged) {
  avg_rev <- MuMIn::model.avg(competitive, revised.var = TRUE)
  ct_full <- MuMIn::coefTable(avg_rev, full = TRUE)
  ct_cond <- MuMIn::coefTable(avg_rev, full = FALSE)
  ci_full <- confint(avg_rev, full = TRUE)
  sw_all  <- MuMIn::sw(dredge_table)

  both <- data.frame(
    Predictor      = rownames(ct_full),
    Beta_full      = round(ct_full[, 1], 4),
    SE_full        = round(ct_full[, 2], 4),
    CI_lower_full  = round(ci_full[, 1], 4),
    CI_upper_full  = round(ci_full[, 2], 4),
    Beta_conditional = round(ct_cond[match(rownames(ct_full), rownames(ct_cond)), 1], 4),
    SE_conditional   = round(ct_cond[match(rownames(ct_full), rownames(ct_cond)), 2], 4),
    stringsAsFactors = FALSE, row.names = NULL
  )
  both$Sum_of_weights <- round(sw_all[match(both$Predictor, names(sw_all))], 4)
  write.csv(both,
            file.path(output_path, "Table_S_PartIII_Averaging_FullVsConditional.csv"),
            row.names = FALSE)
  message("\n   Full and conditional averages, with the sum of Akaike weights:")
  print(both, row.names = FALSE)
  message("   These describe selection uncertainty. They are not the inferential")
  message("   estimates; those are in Table_S_PartIII_ConfirmatoryModel.csv.")
}
write.csv(coef_df,
          file.path(output_path, "Table_S_PartIII_FinalCoefficients.csv"),
          row.names = FALSE)

# Compact table of the competitive models.
comp_df <- as.data.frame(subset(dredge_table, delta < 2))
metric_cols <- c("df", "logLik", "AICc", "delta", "weight")
var_cols    <- setdiff(names(comp_df), metric_cols)
formulas <- vapply(seq_len(nrow(comp_df)), function(i) {
  vp <- var_cols[!is.na(comp_df[i, var_cols])]
  vp <- vp[vp != "(Intercept)"]
  if (length(vp) == 0) "~ Null Model ~" else paste(vp, collapse = " + ")
}, character(1))
comp_clean <- data.frame(
  Rank               = seq_len(nrow(comp_df)),
  Selected_Variables = formulas,
  df                 = comp_df$df,
  LogLik             = round(comp_df$logLik, 3),
  AICc               = round(comp_df$AICc,   3),
  Delta_AICc         = round(comp_df$delta,  3),
  Akaike_Weight      = round(comp_df$weight, 3)
)
write.csv(comp_clean,
          file.path(output_path, "Table_S_PartIII_CompetitiveModels.csv"),
          row.names = FALSE)


# 9. Collinearity stability check: PD_Forest and Mean_FRAC_Forest (Table S21) --
# PD_Forest and Mean_FRAC_Forest are moderately correlated and enter the
# retained model with opposite signs (PD negative, FRAC positive). Two correlated
# predictors with opposite-sign effects can reinforce each other in a joint model
# (mutual suppression). The retained model is therefore refitted twice, once
# without PD_Forest and once without Mean_FRAC_Forest, and the coefficients are
# compared with the full model. If a remaining predictor keeps its sign and a
# comparable magnitude, the joint reading is defensible.
#
# The pairwise correlation is reported with Pearson, the measure of linear
# dependence that matches the VIF (Dormann et al. 2013), and with Spearman as a
# cross-check that is robust to skew and outliers (Zuur et al. 2010).
message("\n[Stability] Collinearity stability check (PD vs FRAC)...")

r_pd_frac          <- cor(data_p3$PD_Forest, data_p3$Mean_FRAC_Forest,
                          method = "pearson", use = "complete.obs")
rho_pd_frac        <- cor(data_p3$PD_Forest, data_p3$Mean_FRAC_Forest,
                          method = "spearman", use = "complete.obs")
message(sprintf("   r(PD_Forest, Mean_FRAC_Forest): Pearson = %.3f, Spearman = %.3f",
                r_pd_frac, rho_pd_frac))

# The models are fitted in the retained family, so that the coefficients are on
# the same scale as Table 3.
fit_retained_family <- function(rhs) {
  if (inherits(global_model, "glm")) {
    glm(as.formula(paste("AOH_ha ~", rhs)), data = data_p3, family = Gamma(link = "log"))
  } else {
    lm(as.formula(paste("log_AOH ~", rhs)), data = data_p3)
  }
}
lm_full_4    <- fit_retained_family(INT_RHS)
lm_drop_pd   <- fit_retained_family(sub("PD_Forest_z \\+ ", "", INT_RHS))
lm_drop_frac <- fit_retained_family(sub("Mean_FRAC_Forest_z \\+ ", "", INT_RHS))

# Coefficient of one term in a model, NA if absent. `term` comes first so that
# the function can be iterated with vapply() while the model is passed by name.
get_beta <- function(term, model) {
  cf <- coef(model)
  if (term %in% names(cf)) as.numeric(cf[term]) else NA_real_
}

stability_terms <- c("PD_Forest_z", "Mean_FRAC_Forest_z",
                     "Mean_ENN_Forest_m_z", "PROX_MN_Forest_z")

stability_df <- data.frame(
  Predictor      = stability_terms,
  Beta_full      = round(vapply(stability_terms, get_beta,
                                numeric(1), model = lm_full_4), 4),
  Beta_drop_PD   = round(vapply(stability_terms, get_beta,
                                numeric(1), model = lm_drop_pd), 4),
  Beta_drop_FRAC = round(vapply(stability_terms, get_beta,
                                numeric(1), model = lm_drop_frac), 4),
  row.names = NULL
)
# The two pairwise correlations are appended as the last rows of the table.
stability_df <- rbind(
  stability_df,
  data.frame(Predictor = "Pearson r(PD_Forest, Mean_FRAC_Forest)",
             Beta_full = round(r_pd_frac, 4),
             Beta_drop_PD = NA, Beta_drop_FRAC = NA),
  data.frame(Predictor = "Spearman rho(PD_Forest, Mean_FRAC_Forest)",
             Beta_full = round(rho_pd_frac, 4),
             Beta_drop_PD = NA, Beta_drop_FRAC = NA)
)
write.csv(stability_df,
          file.path(output_path, "Table_S_PartIII_StabilityCheck.csv"),
          row.names = FALSE)
message("   Stability check exported (full vs PD-dropped vs FRAC-dropped).")


# 10. Sensitivity model: proximity of the farming matrix (Table S23) -----------
# The proximity of the farming matrix is not part of the a priori model, which
# explains the response by forest configuration only. This sensitivity model
# adds it to the a priori model, in the retained family, runs the same AICc
# search and exports the coefficients of the best model. The model also
# contains PLAND and PD, so its matrix coefficient inherits the coupling that
# script 05 measures; script 05 refits the matrix term without them (Table S40).
message("\n[Sensitivity] Matrix-proximity sensitivity model...")

matrix_pred <- "PROX_MN_Agropecuaria"
if (matrix_pred %in% names(data_raw)) {

  data_sens <- data_raw %>%
    select(UA_ID, GENUS, AOH_ha, all_of(predictors_part3), all_of(matrix_pred)) %>%
    drop_na(all_of(c(predictors_part3, matrix_pred))) %>%
    mutate(AOH_ha  = ifelse(AOH_ha == 0, 0.001, AOH_ha),
           log_AOH = log(AOH_ha)) %>%
    mutate(across(all_of(c(predictors_part3, matrix_pred)),
                  ~as.numeric(scale(.x)), .names = "{.col}_z"))

  # A sensitivity check is informative only when the two models share a sample.
  if (nrow(data_sens) != n_obs) {
    stop(sprintf(paste0("The sensitivity model would use %d observations while ",
                        "the retained model uses %d. The two must share a ",
                        "sample to be comparable. Check PROX_MN_Agropecuaria in ",
                        "Data_Raw_WithCoords.csv."),
                 nrow(data_sens), n_obs))
  }
  message(sprintf("   Sensitivity model fitted on %d observations, as the retained model.",
                  nrow(data_sens)))

  # Same family as the retained model, so that the coefficients are comparable
  # with Table 3.
  rhs_sens <- paste(INT_RHS, "+ PROX_MN_Agropecuaria_z")
  lm_sens <- if (inherits(global_model, "glm")) {
    glm(as.formula(paste("AOH_ha ~", rhs_sens)), data = data_sens,
        family = Gamma(link = "log"))
  } else {
    lm(as.formula(paste("log_AOH ~", rhs_sens)), data = data_sens)
  }

  # VIF of the full sensitivity model; the matrix term may add collinearity.
  vif_sens <- car::vif(lm_sens)
  message(sprintf("   Maximum VIF in the sensitivity model: %.2f",
                  max(vif_sens)))

  # The same AICc search on the larger model.
  set.seed(GLOBAL_SEED)
  dredge_sens <- MuMIn::dredge(lm_sens, rank = "AICc", trace = FALSE)
  best_sens   <- MuMIn::get.models(dredge_sens, subset = 1)[[1]]

  coef_sens <- extract_coefs(best_sens,
                             label = "Sensitivity: best 5-predictor subset") %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))
  # VIF of each term in the full sensitivity model, as the last column.
  coef_sens$VIF_full_model <- round(
    vif_sens[match(stringr::str_remove(coef_sens$Predictor, "_z$") %>%
                     paste0("_z"), names(vif_sens))], 3)
  write.csv(coef_sens,
            file.path(output_path, "Table_S_PartIII_Sensitivity_MatrixModel.csv"),
            row.names = FALSE)
  message("   Matrix-proximity sensitivity model exported.")

} else {
  message(sprintf("   '%s' not found in the input data; sensitivity model skipped.",
                  matrix_pred))
}


# 11. Shape complexity against mean patch area (Table S22) ---------------------
# Checks whether shape complexity (Mean_FRAC_Forest) tracks the mean forest
# patch area (Mean_AREA_Forest_ha), which bears on how a positive FRAC
# coefficient can be read: if more complex shapes belong to smaller patches,
# the FRAC effect cannot be attributed to patch extent. Pearson and Spearman
# are both reported, because the marginal distributions can be skewed.
# Descriptive only; this section does not change the model.
message("\n[Descriptive] FRAC vs forest patch area correlation...")

area_pred <- "Mean_AREA_Forest_ha"
if (area_pred %in% names(data_raw)) {
  area_check <- data_raw %>%
    select(Mean_FRAC_Forest, all_of(area_pred)) %>%
    drop_na()
  r_frac_area_pearson  <- cor(area_check$Mean_FRAC_Forest,
                              area_check[[area_pred]], method = "pearson")
  r_frac_area_spearman <- cor(area_check$Mean_FRAC_Forest,
                              area_check[[area_pred]], method = "spearman")
  frac_area_df <- data.frame(
    Metric = c("Pearson r (FRAC vs patch area)",
               "Spearman rho (FRAC vs patch area)",
               "N"),
    Value  = round(c(r_frac_area_pearson, r_frac_area_spearman,
                     nrow(area_check)), 4)
  )
  write.csv(frac_area_df,
            file.path(output_path, "Table_S_PartIII_FRAC_Area_Correlation.csv"),
            row.names = FALSE)
  message(sprintf("   r(FRAC, patch area): Pearson = %.3f, Spearman = %.3f (n = %d).",
                  r_frac_area_pearson, r_frac_area_spearman, nrow(area_check)))
} else {
  message(sprintf("   '%s' not in input; FRAC-area check skipped.", area_pred))
}


# 12. Figure S10: coefficients of the a priori set -----------------------------
# The panel is split in two. Three of the six terms are divided by the window
# area, which is the response, or built on a term that is, and the permutation
# null of script 05 reproduces all three (Table S28b). They are drawn hollow and
# in gray and carry no significance mark. The three terms free of the response
# keep the mark.
message("\n[Figure S10] Coefficient plot of the a priori set...")

# The terms the response normalizes. PD is NP divided by the unit area, PLAND is
# the class area divided by the unit area, and the product inherits PLAND.
coupled_terms <- c("PD_Forest_z", "PLAND_Forest_z",
                   "PLAND_Forest_z:PROX_MN_Forest_z",
                   "PROX_MN_Forest_z:PLAND_Forest_z")

GRP_FREE    <- "Free of the response: interpretable"
GRP_COUPLED <- "Normalized by the response: not interpretable"

forest_df <- coef_df %>%
  filter(Predictor != "(Intercept)") %>%
  mutate(
    Predictor_raw = strip_z(Predictor),
    Label         = label_predictors(Predictor_raw),
    Group         = ifelse(Predictor %in% coupled_terms, GRP_COUPLED, GRP_FREE),
    Significance  = ifelse(Predictor %in% coupled_terms, "",
                           case_when(P_value < 0.001 ~ "***",
                                     P_value < 0.01  ~ "**",
                                     P_value < 0.05  ~ "*",
                                     TRUE            ~ ""))
  ) %>%
  arrange(Group, Estimate) %>%
  mutate(Label = factor(Label, levels = Label),
         Group = factor(Group, levels = c(GRP_FREE, GRP_COUPLED)))

plot_forest <- ggplot(forest_df,
                      aes(x = Estimate, y = Label,
                          colour = Group, fill = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                width = 0.18, linewidth = 0.6, orientation = "y") +
  geom_point(size = 2.8, shape = 21, stroke = 0.8) +
  geom_text(aes(label = Significance,
                x = ifelse(Estimate >= 0, CI_upper + 0.04, CI_lower - 0.04)),
            colour = "black", size = 3.4, fontface = "bold",
            hjust = ifelse(forest_df$Estimate >= 0, 0, 1),
            show.legend = FALSE) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_colour_manual(values = setNames(c("#0072B2", "#8A8A8A"),
                                        c(GRP_FREE, GRP_COUPLED)), guide = "none") +
  scale_fill_manual(values = setNames(c("#0072B2", "white"),
                                      c(GRP_FREE, GRP_COUPLED)), guide = "none") +
  labs(x = expression("Standardized coefficient (" * beta * ") with 95% CI"),
       y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.y        = element_text(colour = "black"),
        strip.text         = element_text(face = "bold", hjust = 0,
                                          margin = margin(3, 3, 3, 3)),
        strip.background   = element_rect(fill = "grey95", colour = NA),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing      = unit(6, "pt"))

save_publication_plot(plot_forest,
                      "Figure_PIII_02_Coefficient_ForestPlot.tiff",
                      width_mm = 140, height_mm = 90)


# 13. Figure S3: variable importance -------------------------------------------
message("\n[Figure S3] Variable importance...")

imp_plot_df <- df_importance %>%
  mutate(
    Predictor_clean = strip_z(Predictor),
    Label = label_predictors(Predictor_clean)
  ) %>%
  arrange(Importance) %>%
  mutate(Label = factor(Label, levels = Label))

plot_importance <- ggplot(imp_plot_df, aes(x = Importance, y = Label)) +
  geom_col(fill = col_main, width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", Importance)),
            hjust = -0.15, size = 2.8) +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Sum of Akaike weights", y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(colour = "black"))

save_publication_plot(plot_importance,
                      "Figure_PIII_03_Variable_Importance.tiff",
                      width_mm = 140, height_mm = 70)


# Axis labels shared by the two heatmaps --------------------------------------

# Metric labels in up to three lines, so that they fit under the columns
# without rotation.
metric_labels_axis <- c(
  "PD_Forest"         = "Patch\ndensity\n(PD)",
  "Mean_FRAC_Forest"  = "Shape\ncomplexity\n(FRAC_MN)",
  "Mean_ENN_Forest_m" = "Isolation\n(ENN_MN)",
  "PROX_MN_Forest"    = "Forest\nproximity\n(PROX_MN)",
  "PLAND_Forest"      = "Forest\ncover\n(PLAND)"
)

# Every predictor of the heatmaps must have a label.
stopifnot(all(predictors_part3 %in% names(metric_labels_axis)))

# Row labels as plotmath expressions: genus names in italics, the pooled row in
# bold, and the number of units of each row.
group_axis_labels <- function(groups, n_lookup, pooled) {
  txt <- ifelse(groups == pooled,
                sprintf('bold("%s")~"(n = %d)"', groups, n_lookup[groups]),
                sprintf('italic("%s")~"(n = %d)"', groups, n_lookup[groups]))
  parse(text = txt)
}


# 14. Figure S7: Spearman correlations by genus (Table S25) --------------------
# For every genus, and for the pooled data ("Overall"), each cell shows the
# Spearman correlation between a forest metric and log(OHE). Spearman is used
# because the number of units per genus is small and the association need not
# be linear within a genus. The figure is descriptive, not inferential: genera
# with fewer than MIN_N_GENUS units are shown as "n.e." (not estimable), and no
# p-value is reported per genus.
message("\n[Figure S7] Metric-association heatmap (by genus and overall)...")

# Long table of raw metric values with the response, one row per observation.
assoc_long <- data_p3 %>%
  select(GENUS, log_AOH, all_of(predictors_part3)) %>%
  pivot_longer(cols = all_of(predictors_part3),
               names_to = "Metric", values_to = "MetricValue")

# Per-genus Spearman correlations (only where n >= MIN_N_GENUS).
assoc_by_genus <- assoc_long %>%
  group_by(GENUS, Metric) %>%
  summarise(
    n   = sum(!is.na(MetricValue) & !is.na(log_AOH)),
    rho = if (n >= MIN_N_GENUS &&
              sd(MetricValue, na.rm = TRUE) > 0) {
      suppressWarnings(cor(MetricValue, log_AOH,
                           method = "spearman", use = "complete.obs"))
    } else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(Group = as.character(GENUS)) %>%
  select(Group, Metric, n, rho)

# Overall (pooled) Spearman correlations.
assoc_overall <- assoc_long %>%
  group_by(Metric) %>%
  summarise(
    n   = sum(!is.na(MetricValue) & !is.na(log_AOH)),
    rho = suppressWarnings(cor(MetricValue, log_AOH,
                               method = "spearman", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(Group = "Overall") %>%
  select(Group, Metric, n, rho)

assoc_df <- bind_rows(assoc_by_genus, assoc_overall)

write.csv(assoc_df,
          file.path(output_path, "Table_S_PartIII_MetricAssociation_ByGenus.csv"),
          row.names = FALSE)

# Number of units per row.
assoc_group_n <- assoc_df %>%
  group_by(Group) %>%
  summarise(n = max(n), .groups = "drop")
n_lookup_assoc <- setNames(assoc_group_n$n, assoc_group_n$Group)

# Rows from top to bottom: "Overall", then the genera in alphabetical order.
# The first factor level is drawn at the bottom of the y axis.
assoc_genera <- sort(setdiff(unique(assoc_df$Group), "Overall"))

assoc_plot_df <- assoc_df %>%
  mutate(
    Metric_lab = factor(metric_labels_axis[Metric],
                        levels = metric_labels_axis[predictors_part3]),
    Group = factor(Group, levels = c(rev(assoc_genera), "Overall")),
    cell_label = ifelse(is.na(rho), "n.e.", sprintf("%.2f", rho)),
    # Dark text on light or empty cells, white text on saturated cells.
    txt_col    = ifelse(is.na(rho) | abs(rho) < 0.55, "grey15", "white")
  )

plot_assoc <- ggplot(assoc_plot_df,
                     aes(x = Metric_lab, y = Group, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label, colour = txt_col), size = 2.6) +
  scale_colour_identity() +
  scale_fill_gradient2(low = heat_low, mid = heat_mid, high = heat_high,
                       midpoint = 0, limits = c(-1, 1),
                       na.value = "grey88",
                       name = expression("Spearman " * rho),
                       breaks = seq(-1, 1, 0.5)) +
  scale_y_discrete(labels = function(x) group_axis_labels(x, n_lookup_assoc, "Overall")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x  = element_text(colour = "black", vjust = 1, lineheight = 0.9),
        axis.text.y  = element_text(colour = "black"),
        panel.grid   = element_blank(),
        legend.position = "right")

save_publication_plot(plot_assoc,
                      "Figure_PIII_07_Metric_Association_Heatmap.tiff",
                      width_mm = 140, height_mm = 115)


# 15. Figure 7: standardized effect by genus (Table S24) -----------------------
# What is the association of each forest metric with log(OHE), for all species
# together and for each genus separately? The estimate is the standardized
# slope (beta) of a univariate OLS, log(OHE) ~ metric_z, fitted within each
# genus and for all species pooled. The metric keeps its global z-score, so
# beta is the change in log(OHE) per standard deviation of the metric, and all
# rows are on the same scale. The color encodes sign and magnitude (blue for
# negative, vermillion for positive).
#
# Why univariate, and not one model per genus: each genus has between four and
# twelve units, so a multivariable model or an AICc selection inside a genus
# would be overparameterized and unstable. Where a genus has fewer than
# MIN_N_GENUS units, or no variance in the metric, the cell is shown as "n.e."
#
# This figure is a descriptive companion to the models: it shows whether the
# pooled association is shared across genera or driven by a few of them. A
# univariate slope can differ in sign or magnitude from a multivariable
# coefficient, especially for PD_Forest and Mean_FRAC_Forest, which are
# correlated.
message("\n[Figure 7] Standardized-effect heatmap (by genus and all species)...")

# Significance marker from a p-value (descriptive at the genus level).
sig_star <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Univariate standardized slope of log(OHE) on one metric, within a subset.
# metric_z is the global z-score, so every beta is on a common "per global SD"
# scale and the genus rows and the pooled row are directly comparable.
univariate_slope <- function(df, metric) {
  z_var <- paste0(metric, "_z")
  x  <- df[[z_var]]
  y  <- df$log_AOH
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok]; y <- y[ok]
  n  <- length(y)
  if (n < MIN_N_GENUS || sd(x) == 0) {
    return(data.frame(n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  fit <- lm(y ~ x)
  cf  <- summary(fit)$coefficients
  if (!"x" %in% rownames(cf)) {
    return(data.frame(n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  data.frame(n    = n,
             beta = unname(cf["x", "Estimate"]),
             se   = unname(cf["x", "Std. Error"]),
             p    = unname(cf["x", "Pr(>|t|)"]))
}

genera_vec <- sort(as.character(unique(data_p3$GENUS)))
groups_vec <- c(genera_vec, "All species")

effect_rows <- list()
for (g in groups_vec) {
  sub <- if (g == "All species") data_p3 else dplyr::filter(data_p3, as.character(GENUS) == g)
  for (m in predictors_part3) {
    s <- univariate_slope(sub, m)
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(
      Group  = g,
      Metric = m,
      n      = s$n,
      beta   = round(s$beta, 4),
      se     = round(s$se, 4),
      p      = round(s$p, 4),
      stringsAsFactors = FALSE
    )
  }
}
effect_df <- dplyr::bind_rows(effect_rows)

write.csv(effect_df,
          file.path(output_path, "Table_S_PartIII_EffectByGenus.csv"),
          row.names = FALSE)


# 15.1 Genus slopes with the genus as a moderator (Table S29) ------------------
# The block above fits one model inside each genus, so a genus with four units
# estimates its own residual variance from four units. This block fits one
# model per metric over all 67 units, with the genus as a moderator, so every
# genus slope is read off a single fitted surface and shares one residual
# variance. Genera with few units borrow precision from the rest, and the two
# tables bracket the genus-level result.
#
# The family is Gamma with a log link, as elsewhere in this script, so the
# slopes are on the same multiplicative scale as the models. Each genus slope is
# the linear contrast that adds the reference slope to the interaction term of
# that genus. Its standard error comes from the full covariance matrix of the
# fit, because the two terms are correlated, and its interval uses the t
# quantile on the residual degrees of freedom.
slope_by_genus_moderator <- function(dat, metric) {
  x <- as.numeric(scale(dat[[metric]]))
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(NULL)
  d_fit <- data.frame(y = dat[["AOH_ha"]], x = x,
                      GENUS = factor(as.character(dat$GENUS)))
  d_fit <- d_fit[stats::complete.cases(d_fit), ]
  fit <- stats::glm(y ~ x * GENUS, family = stats::Gamma(link = "log"),
                    data = d_fit)
  bt <- stats::coef(fit)
  V  <- stats::vcov(fit)
  tq <- stats::qt(0.975, stats::df.residual(fit))
  levs <- levels(d_fit$GENUS)
  out <- lapply(levs, function(g) {
    cv <- rep(0, length(bt)); names(cv) <- names(bt)
    cv["x"] <- 1
    inter <- paste0("x:GENUS", g)
    if (inter %in% names(bt)) cv[inter] <- 1
    b  <- sum(cv * bt)
    se <- sqrt(as.numeric(t(cv) %*% V %*% cv))
    data.frame(
      Metric   = metric,
      Genus    = g,
      N_units  = sum(d_fit$GENUS == g),
      Slope    = round(b, 4),
      CI_lower = round(b - tq * se, 4),
      CI_upper = round(b + tq * se, 4),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

moderator_df <- dplyr::bind_rows(
  lapply(predictors_part3, function(m) slope_by_genus_moderator(data_p3, m))
)
moderator_df$Significant <- ifelse(
  moderator_df$CI_lower > 0 | moderator_df$CI_upper < 0, "yes", "no")

n_excluding_zero <- sum(moderator_df$Significant == "yes")
message(sprintf(
  "Moderator slopes: %d of %d intervals exclude zero.",
  n_excluding_zero, nrow(moderator_df)))

write.csv(moderator_df,
          file.path(output_path,
                    "Table_S_PartIII_SlopeByGenus_Moderator.csv"),
          row.names = FALSE)

# Rows from top to bottom: "All species", then the genera in alphabetical
# order. The first factor level is drawn at the bottom of the y axis.
group_levels <- c(rev(genera_vec), "All species")

# Number of units per row.
group_n   <- effect_df %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(n = max(n), .groups = "drop")
n_lookup  <- setNames(group_n$n, group_n$Group)

# Symmetric color limit for the diverging scale.
beta_lim <- max(abs(effect_df$beta), na.rm = TRUE)
beta_lim <- ifelse(is.finite(beta_lim), ceiling(beta_lim * 10) / 10, 1)

effect_plot_df <- effect_df %>%
  dplyr::mutate(
    Metric_lab = factor(metric_labels_axis[Metric],
                        levels = metric_labels_axis[predictors_part3]),
    Group      = factor(Group, levels = group_levels),
    cell_label = ifelse(is.na(beta), "n.e.",
                        paste0(sprintf("%.2f", beta),
                               vapply(p, sig_star, character(1)))),
    # Dark text on light cells, white text on strongly colored cells.
    txt_col    = ifelse(!is.na(beta) & abs(beta) > 0.6 * beta_lim,
                        "white", "grey15")
  )

plot_effect_heatmap <- ggplot(effect_plot_df,
                              aes(x = Metric_lab, y = Group, fill = beta)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label, colour = txt_col), size = 2.6) +
  scale_colour_identity() +
  scale_fill_gradient2(low = heat_low, mid = heat_mid, high = heat_high,
                       midpoint = 0, limits = c(-beta_lim, beta_lim),
                       na.value = "grey88",
                       name = expression("Standardized " * beta),
                       breaks = scales::pretty_breaks(n = 5)) +
  scale_y_discrete(labels = function(x) group_axis_labels(x, n_lookup, "All species")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x  = element_text(colour = "black", vjust = 1, lineheight = 0.9),
        axis.text.y  = element_text(colour = "black"),
        panel.grid   = element_blank(),
        legend.position = "right")

save_publication_plot(plot_effect_heatmap,
                      "Figure_PIII_08_Effect_Heatmap.tiff",
                      width_mm = 140, height_mm = 115)


# 16. Figure S4: partial effects -----------------------------------------------
# Predicted response for each predictor while the others are held at their mean
# (0 on the z-score scale), back-transformed to hectares.
message("\n[Figure S4] Partial effect curves...")

# Model used for prediction. When the competitive set holds a single model, it
# is used directly. When the competitive models are averaged, the curves come
# from the competitive model with the highest weight, because
# MuMIn::predict.averaging is sensitive to mismatches between the dredge term
# names and the new data.
predict_model <- if (averaged) competitive[[1]] else final_model
is_glm_pred   <- inherits(predict_model, "glm") &&
  family(predict_model)$family == "Gamma"

# Predictions on a grid for each predictor. Both branches predict on the link
# (log) scale, build the 95% confidence interval there, and only then
# back-transform with exp(), which avoids the asymmetry of intervals built on
# the response scale.
make_partial_df <- function(varname, model, data, raw_data) {
  z_var   <- paste0(varname, "_z")
  grid_z  <- seq(min(data[[z_var]], na.rm = TRUE),
                 max(data[[z_var]], na.rm = TRUE),
                 length.out = 100)
  newdata <- data %>%
    summarise(across(ends_with("_z"), ~0)) %>%
    slice(rep(1, length(grid_z)))
  newdata[[z_var]] <- grid_z

  # OLS on log(OHE): predict() returns log(OHE). Gamma GLM with a log link:
  # predict(type = "link") returns log(mu). Mixed model: predict.merMod has no
  # se.fit argument, so the population-level curve (random intercept at zero)
  # is built from the fixed-effect design matrix and its covariance.
  if (inherits(model, "merMod")) {
    for (g in names(lme4::getME(model, "flist"))) {
      newdata[[g]] <- data[[g]][1]          # placeholder level, unused below
    }
    fe_terms <- stats::delete.response(
      stats::terms(lme4::nobars(stats::formula(model))))
    X        <- stats::model.matrix(fe_terms, data = newdata)
    beta     <- lme4::fixef(model)
    X        <- X[, names(beta), drop = FALSE]
    V        <- as.matrix(stats::vcov(model))
    link_fit <- as.vector(X %*% beta)
    se_fit   <- sqrt(pmax(rowSums((X %*% V) * X), 0))
  } else if (is_glm_pred) {
    pred     <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
    link_fit <- pred$fit
    se_fit   <- pred$se.fit
  } else {
    pred     <- predict(model, newdata = newdata, se.fit = TRUE)
    link_fit <- pred$fit
    se_fit   <- pred$se.fit
  }
  link_lo   <- link_fit - 1.96 * se_fit
  link_up   <- link_fit + 1.96 * se_fit

  data.frame(
    Predictor = varname,
    x_real    = grid_z * sd(raw_data[[varname]], na.rm = TRUE) +
      mean(raw_data[[varname]], na.rm = TRUE),
    fit       = exp(link_fit),
    ci_low    = exp(link_lo),
    ci_up     = exp(link_up)
  )
}

# Panels in the order of predictors_part3.
partial_levels <- unname(predictor_labels_units[predictors_part3])

partial_df <- bind_rows(lapply(predictors_part3, make_partial_df,
                               model    = predict_model,
                               data     = data_p3,
                               raw_data = data_p3)) %>%
  mutate(Predictor_label = factor(label_predictors(Predictor, predictor_labels_units),
                                  levels = partial_levels))

raw_points_df <- data_p3 %>%
  select(AOH_ha, all_of(predictors_part3)) %>%
  pivot_longer(cols = all_of(predictors_part3),
               names_to = "Predictor", values_to = "x_real") %>%
  mutate(Predictor_label = factor(label_predictors(Predictor, predictor_labels_units),
                                  levels = partial_levels))

plot_partial <- ggplot() +
  geom_point(data = raw_points_df,
             aes(x = x_real, y = AOH_ha),
             alpha = 0.35, size = 1.3, color = "grey45") +
  geom_ribbon(data = partial_df,
              aes(x = x_real, ymin = ci_low, ymax = ci_up),
              fill = col_main, alpha = 0.20) +
  geom_line(data = partial_df,
            aes(x = x_real, y = fit),
            color = col_main, linewidth = 0.9) +
  facet_wrap(~ Predictor_label, scales = "free_x", ncol = 3) +
  scale_x_continuous(labels = label_number(big.mark = ",")) +
  scale_y_log10(labels = label_comma()) +
  labs(x = NULL,
       y = "OHE (ha, log scale)") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", lineheight = 0.9))

save_publication_plot(plot_partial,
                      "Figure_PIII_04_Partial_Effects.tiff",
                      width_mm = 190, height_mm = 125)


# 17. Figure S5: observed against predicted ------------------------------------
# Checks that the model captures the gradient of the OHE, not only its central
# tendency.
message("\n[Figure S5] Observed vs predicted...")

obs_pred_df <- data_p3 %>%
  mutate(predicted = if (is_glm_pred) {
    predict(predict_model, type = "response")
  } else if (inherits(predict_model, "merMod")) {
    # re.form = NA gives the population-level prediction, consistent with the
    # partial-effect curves above.
    exp(predict(predict_model, re.form = NA))
  } else {
    exp(predict(predict_model, type = "response"))
  })

r_squared_obs_pred <- cor(obs_pred_df$AOH_ha, obs_pred_df$predicted)^2

# The caption of Figure S5 quotes the value annotated in the panel and its
# counterpart on the logarithmic scale; both are exported.
r_squared_obs_pred_log <- cor(log(obs_pred_df$AOH_ha),
                              log(obs_pred_df$predicted))^2
write.csv(
  data.frame(
    Quantity = c("Squared correlation, observed against predicted, hectare scale",
                 "Squared correlation, observed against predicted, logarithmic scale",
                 "Model used for the figure",
                 "Largest sampling unit (ha)"),
    Value    = c(round(r_squared_obs_pred, 4),
                 round(r_squared_obs_pred_log, 4),
                 NA, round(max(obs_pred_df$AOH_ha), 2)),
    Note     = c("annotated in the panel", "",
                 model_family_label, ""),
    stringsAsFactors = FALSE
  ),
  file.path(output_path, "Table_S_PartIII_ObsPred_R2.csv"),
  row.names = FALSE, na = "")
message(sprintf("   R2 obs-pred: %.4f on the hectare scale, %.4f on the log scale",
                r_squared_obs_pred, r_squared_obs_pred_log))

plot_obs_pred <- ggplot(obs_pred_df, aes(x = predicted, y = AOH_ha)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.5) +
  geom_point(alpha = 0.55, size = 1.8, color = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE,
              color = "#D55E00", fill = "#FED9B7", linewidth = 0.7) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  annotate("text", x = min(obs_pred_df$predicted), y = max(obs_pred_df$AOH_ha),
           label = sprintf("R^2 == %.3f", r_squared_obs_pred),
           parse = TRUE, hjust = 0, vjust = 1, size = 3.2) +
  labs(x = "Predicted OHE (ha)", y = "Observed OHE (ha)") +
  theme_classic(base_size = 9) +
  theme(plot.margin = margin(5.5, 14, 5.5, 5.5))

save_publication_plot(plot_obs_pred,
                      "Figure_PIII_05_Observed_vs_Predicted.tiff",
                      width_mm = 140, height_mm = 125)


# 18. Figure S6: residual diagnostics ------------------------------------------
# Four panels for the retained model. They complement the simulation-based
# tests of Section 5 (Table S17).
message("\n[Figure S6] Residual diagnostic panel...")

tiff(file.path(output_path, "Figure_PIII_06_Residual_Diagnostics.tiff"),
     width = 140, height = 140, units = "mm", res = 500,
     compression = "lzw", pointsize = 9)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1),
    family = "sans", cex.lab = 0.95, cex.main = 1.0)
# A mixed model is drawn by hand, because plot.merMod does not accept the which=
# argument of plot.lm. Panels 1 to 3 mirror plot.lm; the fourth shows the
# residuals across the levels of the grouping factor, because leverage is not
# defined in the same way for a mixed model.
if (inherits(global_model, "merMod")) {
  res_g   <- stats::residuals(global_model)
  fit_g   <- stats::fitted(global_model)
  sres_g  <- res_g / stats::sd(res_g)
  grp_g   <- lme4::getME(global_model, "flist")[[1]]

  plot(fit_g, res_g, xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs fitted", pch = 19, col = "#00000088")
  abline(h = 0, lty = 2, col = "grey40")
  lines(stats::lowess(fit_g, res_g), col = "#D55E00", lwd = 2)

  stats::qqnorm(sres_g, main = "Normal Q-Q", pch = 19, col = "#00000088")
  stats::qqline(sres_g, lty = 2, col = "grey40")

  plot(fit_g, sqrt(abs(sres_g)), xlab = "Fitted values",
       ylab = expression(sqrt(abs("standardized residuals"))),
       main = "Scale-location", pch = 19, col = "#00000088")
  lines(stats::lowess(fit_g, sqrt(abs(sres_g))), col = "#D55E00", lwd = 2)

  boxplot(res_g ~ grp_g, xlab = "", ylab = "Residuals",
          main = "Residuals by grouping level", las = 2,
          cex.axis = 0.7, col = "grey90")
  abline(h = 0, lty = 2, col = "grey40")
} else {
  # sub.caption = "" suppresses the call that plot.lm writes in the outer margin.
  plot(global_model, which = c(1, 2, 3, 5), sub.caption = "")
}
invisible(dev.off())
message("   Residual diagnostic panel exported.")


# 19. Console summary ----------------------------------------------------------
message("\n============================================================")
message("  Script 03. Summary of results")
message("============================================================")
message(sprintf("  Sampling units                : %d", n_su))
message(sprintf("  Retained specification        : %s", model_family_label))
message(sprintf("  Global-model R^2 / DevExpl    : %.3f",
                glob_summary$R2_or_DevExpl))
if (exists("icc_su")) {
  message(sprintf("  ICC of the genus              : %.3f", icc_su))
}
message(sprintf("  Competitive set (Delta < 2)   : %d model(s)", n_competitive))
message(sprintf("  Multimodel summary            : %s", inf_label))
message("\n  Predictor importance (Akaike weights):")
for (i in seq_len(nrow(df_importance))) {
  message(sprintf("    - %-22s : %.3f",
                  df_importance$Predictor[i], df_importance$Importance[i]))
}
message("============================================================\n")


# 20. Files written ------------------------------------------------------------
message("[Export] Files written by script 03:")
message("  Tables:")
message("    - Table_S_PartIII_Averaging_FullVsConditional.csv")
message("    - Table_S_PartIII_CompetitiveModels.csv")
message("    - Table_S_PartIII_ConfirmatoryModel.csv")
message("    - Table_S_PartIII_DependenceDiagnostic.csv")
message("    - Table_S_PartIII_DescriptiveStatistics.csv")
message("    - Table_S_PartIII_Diagnostics.csv")
message("    - Table_S_PartIII_Dredge_FullTable.csv")
message("    - Table_S_PartIII_EffectByGenus.csv")
message("    - Table_S_PartIII_FRAC_Area_Correlation.csv")
message("    - Table_S_PartIII_FinalCoefficients.csv")
message("    - Table_S_PartIII_GlobalModelSummary.csv")
message("    - Table_S_PartIII_MetricAssociation_ByGenus.csv")
message("    - Table_S_PartIII_ModelComparison.csv")
message("    - Table_S_PartIII_ObsPred_R2.csv")
message("    - Table_S_PartIII_Sensitivity_MatrixModel.csv")
message("    - Table_S_PartIII_SlopeByGenus_Moderator.csv")
message("    - Table_S_PartIII_SpatialSensitivity.csv")
message("    - Table_S_PartIII_StabilityCheck.csv")
message("    - Table_S_PartIII_VIF.csv")
message("    - Table_S_PartIII_VariableImportance.csv")
message("  Figures (TIFF, 500 dpi):")
message("    - Figure_PIII_01_Response_Transformation.tiff")
message("    - Figure_PIII_02_Coefficient_ForestPlot.tiff")
message("    - Figure_PIII_03_Variable_Importance.tiff")
message("    - Figure_PIII_04_Partial_Effects.tiff")
message("    - Figure_PIII_05_Observed_vs_Predicted.tiff")
message("    - Figure_PIII_06_Residual_Diagnostics.tiff")
message("    - Figure_PIII_07_Metric_Association_Heatmap.tiff")
message("    - Figure_PIII_08_Effect_Heatmap.tiff")
message("    - session_info_part03.txt")


# 21. Session information ------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part03.txt"))


# References -------------------------------------------------------------------
# Aiken, L.S., West, S.G. (1991). Multiple Regression: Testing and Interpreting
#   Interactions. Sage, Newbury Park.
# Andren, H. (1994). Effects of habitat fragmentation on birds and mammals in
#   landscapes with different proportions of suitable habitat: a review.
#   Oikos, 71, 355-366.
# Bolker, B.M., Brooks, M.E., Clark, C.J., Geange, S.W., Poulsen, J.R., Stevens,
#   M.H.H., White, J.S.S. (2009). Generalized linear mixed models: a practical
#   guide for ecology and evolution. Trends in Ecology & Evolution, 24, 127-135.
# Borcard, D., Gillet, F., Legendre, P. (2018). Numerical Ecology with R, 2nd
#   ed. Springer.
# Brewer, M.J., Butler, A., Cooksley, S.L. (2016). The relative performance of
#   AIC, AICc and BIC in the presence of unobserved heterogeneity. Methods in
#   Ecology and Evolution, 7, 679-692.
# Burnham, K.P., Anderson, D.R. (2002). Model Selection and Multimodel
#   Inference: A Practical Information-Theoretic Approach, 2nd ed. Springer.
# Burnham, K.P., Anderson, D.R. (2004). Multimodel inference: understanding AIC
#   and BIC in model selection. Sociological Methods & Research, 33, 261-304.
# Cohen, J., Cohen, P., West, S.G., Aiken, L.S. (2003). Applied Multiple
#   Regression/Correlation Analysis for the Behavioral Sciences, 3rd ed.
#   Lawrence Erlbaum, Mahwah.
# Diniz-Filho, J.A.F., Rangel, T.F.L.V.B., Bini, L.M. (2008). Model selection
#   and information theory in geographical ecology. Global Ecology and
#   Biogeography, 17, 479-488.
# Dormann, C.F. et al. (2013). Collinearity: a review of methods to deal
#   with it and a simulation study evaluating their performance. Ecography,
#   36, 27-46.
# Grueber, C.E., Nakagawa, S., Laws, R.J., Jamieson, I.G. (2011). Multimodel
#   inference in ecology and evolution: challenges and solutions. Journal of
#   Evolutionary Biology, 24, 699-711.
# Nakagawa, S., Schielzeth, H. (2013). A general and simple method for
#   obtaining R2 from generalized linear mixed-effects models. Methods in
#   Ecology and Evolution, 4, 133-142.
# Pardini, R., Bueno, A.A., Gardner, T.A., Prado, P.I., Metzger, J.P. (2010).
#   Beyond the fragmentation threshold hypothesis: regime shifts in
#   biodiversity across fragmented landscapes. PLoS ONE, 5, e13666.
# Schielzeth, H. (2010). Simple means to improve the interpretability of
#   regression coefficients. Methods in Ecology and Evolution, 1, 103-113.
# Tredennick, A.T., Hooker, G., Ellner, S.P., Adler, P.B. (2021). A practical
#   guide to selecting models for exploration, inference, and prediction in
#   ecology. Ecology, 102, e03336.
# Villard, M.-A., Metzger, J.P. (2014). Beyond the fragmentation debate: a
#   conceptual model to predict when habitat configuration really matters.
#   Journal of Applied Ecology, 51, 309-318.
# Yates, L.A., Aandahl, Z., Richards, S.A., Brook, B.W. (2023). Cross
#   validation for model selection: a review with examples from ecology.
#   Ecological Monographs, 93, e1557.
# Zuur, A.F., Ieno, E.N., Walker, N., Saveliev, A.A., Smith, G.M. (2009). Mixed
#   Effects Models and Extensions in Ecology with R. Springer.
# Zuur, A.F., Ieno, E.N., Elphick, C.S. (2010). A protocol for data exploration
#   to avoid common statistical problems. Methods in Ecology and Evolution,
#   1, 3-14.
