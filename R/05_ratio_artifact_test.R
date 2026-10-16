# ==============================================================================
# 05_ratio_artifact_test.R
# Ratio artifact test: which coefficients survive the geometry of the design
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   Forest cover (PLAND) and patch density (PD) are divided by the area of the
#   sampling window, which is also the response: the area of the occupied
#   habitat envelope (OHE). A predictor built this way can correlate with the
#   response through its denominator alone. This is the ratio problem described
#   by Pearson (1897) and Kronmal (1993).
#
#   The script measures that effect with a permutation null. The forest class
#   area (CA) and the number of forest patches (NP) are shuffled together across
#   units, which keeps their mutual association, while the area of each unit
#   stays at its observed value. PLAND and PD are recomputed from the shuffled
#   numerators and the observed denominators, and the model is refitted 9,999
#   times. Under this null there is no relationship between the amount of forest
#   and the extent of the envelope, so the coefficients it produces measure the
#   construction alone. The same null is applied to the AICc difference between
#   the coupled model and the model without window-normalized predictors.
#
#   Six alternative specifications remove the normalized metrics singly and
#   jointly. Further checks ask whether the proximity index, which is not
#   divided by the window area but is built from the areas of neighboring
#   patches, can be separated from the amount of forest, and whether the
#   reported associations depend on overlapping units of different taxa. The
#   script ends by exporting the model the manuscript reports (Table 4 and
#   Figure 6).
#
# Input
#   Dados/Processados/Data_Raw_WithCoords.csv
#
# Output   (in Outputs/Manuscrito/PartV/)
#   TableS_PartV_Identity_and_PartialCorrelations.csv   (Table S28e)
#   TableS_PartV_Specifications.csv                     (Table S28a)
#   TableS_PartV_Matrix_Uncoupled.csv                   (Table S40)
#   TableS_PartV_Spatial_Thinning.csv                   (Table S41)
#   TableS_PartV_Interaction_Gain.csv                   (Table S28f)
#   TableS_PartV_NullSimulation.csv                     (Table S28b)
#   Figure_PV_01_Null_Distribution.tiff                 (Figure S9)
#   TableS_PartV_AICc_Under_Null.csv                    (Table S28c)
#   TableS_PartV_Proximity_Amount_Correlations.csv      (Table S30)
#   TableS_PartV_Proximity_Conditioned.csv              (Table S31)
#   TableS_PartV_Decision.csv                           (Table S28d)
#   TableS_PartV_Uncoupled_Model.csv                    (Table 4)
#   Figure_PV_02_Reported_Model_Coefficients.tiff       (Figure 6)
#   TableS_PartV_CV_And_BIC.csv                         (Table S35)
#
# Run time
#   One to two minutes on the machine used for the study. The null of the
#   coefficients runs N_SIM = 9,999 permutations; the null of the AICc
#   difference runs 2,000, because each permutation refits two models.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------
# Partial correlations are computed from residuals (function below), so the
# ppcor package, which attaches MASS and masks dplyr::select(), is not needed.
# dplyr verbs that other packages also export are called as dplyr::verb().

suppressPackageStartupMessages({
  library(here)        # project-relative paths
  library(tidyverse)   # data wrangling, ggplot2, purrr
  library(car)         # vif()
  library(MuMIn)       # AICc()
})

# Partial correlation of x and y given the columns of z: the Pearson
# correlation between the residuals of x ~ z and of y ~ z. With z = NULL it is
# the ordinary correlation.
partial_cor <- function(x, y, z = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (is.null(z)) return(cor(x, y))
  z <- as.data.frame(z)
  if (ncol(z) == 0) return(cor(x, y))
  ok <- stats::complete.cases(x, y, z)
  zz <- z[ok, , drop = FALSE]
  rx <- stats::resid(stats::lm(x[ok] ~ ., data = zz))
  ry <- stats::resid(stats::lm(y[ok] ~ ., data = zz))
  cor(rx, ry)
}

GLOBAL_SEED      <- 123
N_UNITS_EXPECTED <- 67
N_SIM            <- 9999   # permutations of the coefficient null
N_SIM_AIC_CAP    <- 2000   # permutations of the AICc null; each one refits two
                           # models, so the count is capped below N_SIM
THIN_RADIUS_M    <- 1000   # centroids closer than this join one cluster
N_THIN_REPS      <- 200    # thinned designs drawn from those clusters
set.seed(GLOBAL_SEED)

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartV")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


# 2. Data [1/12] ---------------------------------------------------------------
message("\n[1/12] Loading data...")

# AOH_ha is the area of the occupied habitat envelope (OHE), in hectares.
d <- read_csv(file.path(data_path, "Data_Raw_WithCoords.csv"),
              show_col_types = FALSE) %>%
  rename(AOH_ha = AOH_unit_area_ha)

# PROX_MN_Agropecuaria and PD_Agropecuaria enter only the farming-matrix models
# of Section 5. PD_Agropecuaria, like PD_Forest, is 100 * NP / window area.
needed <- c("AOH_ha", "CA_Forest", "NP_Forest", "PLAND_Forest", "PD_Forest",
            "Mean_FRAC_Forest", "Mean_ENN_Forest_m", "PROX_MN_Forest",
            "Mean_AREA_Forest_ha", "PROX_MN_Agropecuaria", "PD_Agropecuaria",
            "coord_x", "coord_y", "GENUS")
missing <- setdiff(needed, names(d))
if (length(missing) > 0) stop("Missing column(s): ", paste(missing, collapse = ", "))

d <- d %>%
  dplyr::select(GENUS, SPECIES, UA_ID, all_of(setdiff(needed, "GENUS"))) %>%
  tidyr::drop_na() %>%
  mutate(log_AOH = log(AOH_ha), GENUS = as.factor(GENUS))

if (nrow(d) != N_UNITS_EXPECTED) {
  stop(sprintf("Script 05 loaded %d units instead of %d. Rerun script 01.",
               nrow(d), N_UNITS_EXPECTED))
}
message(sprintf("   %d sampling units.", nrow(d)))

z <- function(x) as.numeric(scale(x))


# 3. The algebraic coupling [2/12] ---------------------------------------------
# Two checks. First, that the identity linking PD to the response,
# log(OHE) = log(NP) + log(100) - log(PD), holds on all 67 units, so that the
# concern is not hypothetical. Second, through partial correlations, how much of
# the association between forest proximity and the response is carried by the
# coupled variables.
message("\n[2/12] Verifying the algebraic identity and partial correlations...")

identity_error <- max(abs(log(d$NP_Forest) + log(100) - log(d$PD_Forest) - d$log_AOH))

pc <- partial_cor

coupling <- data.frame(
  Quantity = c(
    "max |log(NP) + log(100) - log(PD) - log(OHE)|",
    "Pearson r, PLAND_Forest with log(OHE)",
    "Spearman rho, PLAND_Forest with log(OHE)",
    "Pearson r, PD_Forest with log(OHE)",
    "Pearson r, log(CA_Forest) with log(OHE)",
    "Pearson r, PROX_MN_Forest with log(OHE)",
    "Partial r, PROX_MN with log(OHE), controlling PLAND",
    "Partial r, PROX_MN with log(OHE), controlling PLAND and PD",
    "Partial r, PROX_MN with log(OHE), controlling Mean_AREA_Forest_ha"
  ),
  Value = round(c(
    identity_error,
    cor(d$PLAND_Forest, d$log_AOH),
    cor(d$PLAND_Forest, d$log_AOH, method = "spearman"),
    cor(d$PD_Forest, d$log_AOH),
    cor(log(d$CA_Forest), d$log_AOH),
    cor(d$PROX_MN_Forest, d$log_AOH),
    pc(d$PROX_MN_Forest, d$log_AOH, d[, "PLAND_Forest"]),
    pc(d$PROX_MN_Forest, d$log_AOH, d[, c("PLAND_Forest", "PD_Forest")]),
    pc(d$PROX_MN_Forest, d$log_AOH, d[, "Mean_AREA_Forest_ha"])
  ), 4),
  Note = c(
    "if near zero, PD is an exact function of the response given NP",
    "marginal association, no conditioning",
    "distribution-free cross-check of the line above",
    "marginal association, no conditioning",
    "the response bounds this quantity from above",
    "PROX is not normalized by the window area",
    "conditioning on a ratio that contains the response",
    "conditioning on two quantities that contain the response",
    "conditioning on a quantity free of the response"
  )
)

write.csv(coupling,
          file.path(output_path, "TableS_PartV_Identity_and_PartialCorrelations.csv"),
          row.names = FALSE)
print(coupling, row.names = FALSE)

if (identity_error > 0.05) {
  message("   NOTE: the PD identity does not hold to within 0.05 on all units.")
  message("   Check the units of PD before interpreting the null simulation.")
} else {
  message(sprintf("   The PD identity holds to within %.4f on every unit.",
                  identity_error))
}


# 4. Competing specifications [3/12] -------------------------------------------
# Every specification is a Gamma GLM with a log link, the family retained by
# script 03, fitted on the same 67 units. Which one fits best is not the point:
# a predictor built from the response will always help to explain the response.
# The point is whether the coefficient of PROX_MN_Forest and the interaction
# keep their sign, magnitude and interval as the coupled predictors are removed.
#
#   S1  retained model                       reference
#   S2  without PLAND                        drops one coupled term
#   S3  without PD                           drops the other
#   S4  neither PLAND nor PD                 no window-normalized predictor
#   S5  PLAND replaced by log(CA_Forest)     absolute amount, not a ratio
#   S6  amount as mean forest patch area     free of the window
message("\n[3/12] Fitting the competing specifications...")

dz <- d %>%
  mutate(PLAND_z   = z(PLAND_Forest),
         PD_z      = z(PD_Forest),
         FRAC_z    = z(Mean_FRAC_Forest),
         ENN_z     = z(Mean_ENN_Forest_m),
         PROX_z    = z(PROX_MN_Forest),
         logCA_z   = z(log(CA_Forest)),
         AREAMN_z  = z(Mean_AREA_Forest_ha))

specs <- list(
  S1 = list(f = AOH_ha ~ ENN_z + FRAC_z + PD_z + PLAND_z + PROX_z + PROX_z:PLAND_z,
            label = "S1. Retained model (reference)",
            coupled = "PLAND and PD"),
  S2 = list(f = AOH_ha ~ ENN_z + FRAC_z + PD_z + PROX_z,
            label = "S2. Without PLAND",
            coupled = "PD only"),
  S3 = list(f = AOH_ha ~ ENN_z + FRAC_z + PLAND_z + PROX_z + PROX_z:PLAND_z,
            label = "S3. Without PD",
            coupled = "PLAND only"),
  S4 = list(f = AOH_ha ~ ENN_z + FRAC_z + PROX_z,
            label = "S4. No area-normalized predictor",
            coupled = "none"),
  S5 = list(f = AOH_ha ~ ENN_z + FRAC_z + PD_z + logCA_z + PROX_z + PROX_z:logCA_z,
            label = "S5. Amount as log(CA_Forest), absolute",
            coupled = "PD only; CA is bounded by the response"),
  S6 = list(f = AOH_ha ~ ENN_z + FRAC_z + AREAMN_z + PROX_z + PROX_z:AREAMN_z,
            label = "S6. Amount as mean forest patch area",
            coupled = "none")
)

fit_spec <- function(sp) {
  m  <- glm(sp$f, data = dz, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- suppressMessages(confint.default(m))

  int_row <- grep(":", rownames(cf), value = TRUE)
  amount  <- intersect(c("PLAND_z", "logCA_z", "AREAMN_z"), rownames(cf))

  # Largest VIF, computed on the linear model of log(OHE).
  vf <- tryCatch({
    lmfit <- lm(update(sp$f, log(AOH_ha) ~ .), data = dz)
    v <- car::vif(lmfit)
    max(v)
  }, error = function(e) NA_real_)

  data.frame(
    Specification   = sp$label,
    Coupled_terms   = sp$coupled,
    AICc            = round(MuMIn::AICc(m), 2),
    Max_VIF         = round(vf, 2),
    PROX_beta       = round(cf["PROX_z", "Estimate"], 4),
    PROX_CI_low     = round(ci["PROX_z", 1], 4),
    PROX_CI_high    = round(ci["PROX_z", 2], 4),
    PROX_p          = signif(cf["PROX_z", "Pr(>|t|)"], 4),
    Amount_term     = if (length(amount)) amount[1] else NA_character_,
    Amount_beta     = if (length(amount)) round(cf[amount[1], "Estimate"], 4) else NA_real_,
    Amount_p        = if (length(amount)) signif(cf[amount[1], "Pr(>|t|)"], 4) else NA_real_,
    Interaction_beta = if (length(int_row)) round(cf[int_row, "Estimate"], 4) else NA_real_,
    Interaction_CI_low  = if (length(int_row)) round(ci[int_row, 1], 4) else NA_real_,
    Interaction_CI_high = if (length(int_row)) round(ci[int_row, 2], 4) else NA_real_,
    Interaction_p    = if (length(int_row)) signif(cf[int_row, "Pr(>|t|)"], 4) else NA_real_,
    stringsAsFactors = FALSE
  )
}

spec_table <- purrr::map_dfr(specs, fit_spec)

write.csv(spec_table,
          file.path(output_path, "TableS_PartV_Specifications.csv"),
          row.names = FALSE)
print(as.data.frame(spec_table), row.names = FALSE)

# Summary of the specifications: does the proximity effect survive, and in
# which specifications does the interaction exclude zero?
prox_sig <- spec_table$PROX_CI_low * spec_table$PROX_CI_high > 0
message(sprintf("   PROX_MN_Forest ranges from %.3f to %.3f across the six specifications;",
                min(spec_table$PROX_beta), max(spec_table$PROX_beta)))
message(sprintf("   its 95%% CI excludes zero in %d of %d.",
                sum(prox_sig), length(prox_sig)))
int_rows <- spec_table %>% dplyr::filter(!is.na(Interaction_beta))
int_sig  <- int_rows$Interaction_CI_low * int_rows$Interaction_CI_high > 0
message(sprintf("   The interaction is present in %d specifications; its 95%% CI",
                nrow(int_rows)))
message(sprintf("   excludes zero in %d of them (%s).", sum(int_sig),
                if (any(int_sig)) paste(sub("\\..*$", "", int_rows$Specification[int_sig]),
                                        collapse = ", ") else "none"))


# 5. The farming matrix without window-normalized predictors [4/12] -----------
# Script 03 reports a sensitivity model that adds the proximity index of the
# farming matrix to the a priori set. That model also contains PLAND and PD, so
# its matrix coefficient inherits the coupling measured here. The proximity
# index of the farming matrix is not divided by the window area, so it can be
# added to S4, the specification without window-normalized predictors (M1).
# M1 is the matrix estimate free of the coupling. Matrix patch density (M2) is
# fitted beside it for contrast only: like PD_Forest, it is 100 * NP / window
# area, so M2 is not free of the coupling.
message("\n[4/12] The farming matrix, without window-normalized predictors...")

dz <- dz %>%
  mutate(PROXA_z = z(PROX_MN_Agropecuaria),
         PDA_z   = z(PD_Agropecuaria))

matrix_specs <- list(
  M1 = list(f = AOH_ha ~ ENN_z + FRAC_z + PROX_z + PROXA_z,
            label = "S4 plus farming matrix proximity",
            focal = "PROXA_z"),
  M2 = list(f = AOH_ha ~ ENN_z + FRAC_z + PROX_z + PDA_z,
            label = "S4 plus farming matrix patch density",
            focal = "PDA_z")
)

fit_matrix <- function(sp) {
  m  <- glm(sp$f, data = dz, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- suppressMessages(confint.default(m))
  vf <- tryCatch(max(car::vif(lm(update(sp$f, log(AOH_ha) ~ .), data = dz))),
                 error = function(e) NA_real_)
  keep <- setdiff(rownames(cf), "(Intercept)")
  data.frame(
    Specification = sp$label,
    Predictor     = keep,
    Estimate      = round(cf[keep, "Estimate"], 4),
    SE            = round(cf[keep, "Std. Error"], 4),
    CI_lower      = round(ci[keep, 1], 4),
    CI_upper      = round(ci[keep, 2], 4),
    P_value       = signif(cf[keep, 4], 4),
    OHE_ratio_per_SD = round(exp(cf[keep, "Estimate"]), 3),
    Deviance_explained = round(1 - m$deviance / m$null.deviance, 4),
    Max_VIF       = round(vf, 2),
    row.names     = NULL, stringsAsFactors = FALSE
  )
}

matrix_table <- purrr::map_dfr(matrix_specs, fit_matrix)
write.csv(matrix_table,
          file.path(output_path, "TableS_PartV_Matrix_Uncoupled.csv"),
          row.names = FALSE)
print(matrix_table, row.names = FALSE)
message("   Read the matrix coefficient from M1, not from the sensitivity model")
message("   of script 03, which also contains PLAND and PD.")


# 6. Spatial coincidence among units of different taxa [5/12] ------------------
# The envelopes of different taxa overlap, and some are nested, so the same
# landscape can enter the model more than once. The random intercept for genus
# does not account for this, because the overlapping units belong to different
# genera. The Moran eigenvectors account for it only in part, because proximity
# between centroids is not the same as overlap between polygons.
#
# Thinning test: units whose centroids lie closer than THIN_RADIUS_M join one
# cluster (single linkage), one unit is drawn at random from each cluster, the
# reported model (S4) is refitted on the thinned design, and the draw is
# repeated N_THIN_REPS times. This is a robustness check of the reported
# associations, not a second set of estimates.
message("\n[5/12] Thinning test for spatial coincidence among units of different taxa...")

xy   <- as.matrix(dz[, c("coord_x", "coord_y")])
dmat <- as.matrix(stats::dist(xy))
diag(dmat) <- Inf

# Single-linkage clustering, written out so that the rule is visible: two units
# join the same cluster when their centroids are closer than THIN_RADIUS_M.
parent <- seq_len(nrow(dz))
find_root <- function(a) { while (parent[a] != a) a <- parent[a]; a }
for (i in seq_len(nrow(dz) - 1)) {
  for (j in (i + 1):nrow(dz)) {
    if (dmat[i, j] < THIN_RADIUS_M) {
      ri <- find_root(i); rj <- find_root(j)
      if (ri != rj) parent[ri] <- rj
    }
  }
}
cluster_id <- vapply(seq_len(nrow(dz)), find_root, integer(1))
n_clusters <- length(unique(cluster_id))
n_in_group <- sum(table(cluster_id)[as.character(cluster_id)] > 1)

cross_taxon <- sum(vapply(unique(cluster_id), function(cl) {
  idx <- which(cluster_id == cl)
  length(idx) > 1 && length(unique(dz$SPECIES[idx])) > 1
}, logical(1)))

message(sprintf(paste0("   %d units fall in %d cluster(s) of coincident centroids ",
                       "at %d m; %d cluster(s) mix taxa. Thinned design: %d units."),
                n_in_group, sum(table(cluster_id) > 1), THIN_RADIUS_M,
                cross_taxon, n_clusters))

set.seed(GLOBAL_SEED)
thin_draws <- matrix(NA_real_, nrow = N_THIN_REPS, ncol = 6,
                     dimnames = list(NULL, c("PROX_beta", "PROX_p",
                                             "FRAC_beta", "FRAC_p",
                                             "ENN_beta", "ENN_p")))
for (r in seq_len(N_THIN_REPS)) {
  keep <- vapply(unique(cluster_id), function(cl) {
    idx <- which(cluster_id == cl)
    if (length(idx) == 1) idx else sample(idx, 1)
  }, integer(1))
  sub <- dz[sort(keep), , drop = FALSE]
  # Standardize again inside the thinned design, so that the coefficients stay
  # on the scale of one standard deviation of the data actually fitted.
  sub <- sub %>% mutate(ENN_z = z(Mean_ENN_Forest_m), FRAC_z = z(Mean_FRAC_Forest),
                        PROX_z = z(PROX_MN_Forest))
  fit <- tryCatch(glm(specs$S4$f, data = sub, family = Gamma(link = "log")),
                  error = function(e) NULL)
  if (is.null(fit)) next
  cf <- summary(fit)$coefficients
  thin_draws[r, ] <- c(cf["PROX_z", 1], cf["PROX_z", 4],
                       cf["FRAC_z", 1], cf["FRAC_z", 4],
                       cf["ENN_z",  1], cf["ENN_z",  4])
}

thin_summary <- data.frame(
  Predictor = c("Forest proximity (PROX_MN)", "Shape complexity (FRAC_MN)",
                "Isolation (ENN_MN)"),
  Beta_full_design = round(c(
    coef(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))[["PROX_z"]],
    coef(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))[["FRAC_z"]],
    coef(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))[["ENN_z"]]), 4),
  Beta_thinned_median = round(c(median(thin_draws[, "PROX_beta"], na.rm = TRUE),
                                median(thin_draws[, "FRAC_beta"], na.rm = TRUE),
                                median(thin_draws[, "ENN_beta"],  na.rm = TRUE)), 4),
  Beta_thinned_q2.5 = round(c(quantile(thin_draws[, "PROX_beta"], 0.025, na.rm = TRUE),
                              quantile(thin_draws[, "FRAC_beta"], 0.025, na.rm = TRUE),
                              quantile(thin_draws[, "ENN_beta"],  0.025, na.rm = TRUE)), 4),
  Beta_thinned_q97.5 = round(c(quantile(thin_draws[, "PROX_beta"], 0.975, na.rm = TRUE),
                               quantile(thin_draws[, "FRAC_beta"], 0.975, na.rm = TRUE),
                               quantile(thin_draws[, "ENN_beta"],  0.975, na.rm = TRUE)), 4),
  Pct_replicates_p_below_0.05 = round(100 * c(
    mean(thin_draws[, "PROX_p"] < 0.05, na.rm = TRUE),
    mean(thin_draws[, "FRAC_p"] < 0.05, na.rm = TRUE),
    mean(thin_draws[, "ENN_p"]  < 0.05, na.rm = TRUE)), 1),
  N_units_per_replicate = n_clusters,
  N_replicates = N_THIN_REPS,
  N_units_in_multi_unit_clusters = n_in_group,
  N_clusters_mixing_taxa = cross_taxon,
  Thinning_radius_m = THIN_RADIUS_M,
  row.names = NULL
)
write.csv(thin_summary,
          file.path(output_path, "TableS_PartV_Spatial_Thinning.csv"),
          row.names = FALSE)
print(thin_summary, row.names = FALSE)


# 7. Does the interaction pay for itself? [6/12] -------------------------------
# Script 03 compares additive and interaction specifications with PLAND in both
# arms, and PLAND carries the response in its denominator. The question here is
# whether the interaction still improves the fit when the amount term it
# interacts with is free of the window. Four models, stated in advance and all
# reported: two with PLAND and two with the mean forest patch area, the only
# measure of amount here that the response neither contains nor bounds. Each
# pair differs only in the product term. The gap between the two Delta AICc
# values is the part of the advantage of the interaction that the
# normalization produces.
message("\n[6/12] Interaction gain with and without a window-normalized amount term...")

pairs_int <- list(
  list(key = "Coupled amount (PLAND)",
       add = AOH_ha ~ ENN_z + FRAC_z + PD_z + PLAND_z + PROX_z,
       int = AOH_ha ~ ENN_z + FRAC_z + PD_z + PLAND_z + PROX_z + PROX_z:PLAND_z),
  list(key = "Amount free of the response (mean patch area)",
       add = AOH_ha ~ ENN_z + FRAC_z + AREAMN_z + PROX_z,
       int = AOH_ha ~ ENN_z + FRAC_z + AREAMN_z + PROX_z + PROX_z:AREAMN_z)
)

# AICc of a Gamma GLM; k counts the coefficients plus the dispersion parameter.
aicc_gamma <- function(f) {
  m <- glm(f, data = dz, family = Gamma(link = "log"))
  k <- length(coef(m)) + 1L
  n <- nobs(m)
  as.numeric(-2 * logLik(m) + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}

interaction_gain <- purrr::map_dfr(pairs_int, function(p) {
  a <- aicc_gamma(p$add)
  i <- aicc_gamma(p$int)
  data.frame(
    Amount_measure   = p$key,
    AICc_additive    = round(a, 2),
    AICc_interaction = round(i, 2),
    Delta_AICc       = round(i - a, 2),
    Verdict          = if (i - a <= -2) {
      "the interaction earns its parameter"
    } else if (i - a < 0) {
      "indistinguishable from the additive model; parsimony favors the additive one"
    } else {
      "the interaction costs more than it returns"
    },
    stringsAsFactors = FALSE
  )
})

write.csv(interaction_gain,
          file.path(output_path, "TableS_PartV_Interaction_Gain.csv"),
          row.names = FALSE)

print(as.data.frame(interaction_gain), row.names = FALSE)
message("   The first row is the comparison made in script 03. The second makes")
message("   the same comparison with an amount measure free of the response.")


# 8. Permutation null of the coefficients [7/12] -------------------------------
# The null keeps everything about the design except the association being
# tested. The response keeps its observed values. The pair (CA_Forest,
# NP_Forest) is permuted jointly across units, so the association between the
# two is kept and only their association with the window area is removed. PLAND
# and PD are then recomputed from the permuted numerators and the observed
# denominator:
#
#   PLAND* = 100 * CA* / OHE        PD* = NP* / (OHE / 100)
#
# Under this null the amount of forest is unrelated to the extent of the window,
# but PLAND* and PD* still carry the response in their denominators. Any
# coefficient they obtain is therefore an artifact of the ratio, and the
# distribution of those coefficients is the yardstick.
#
# Reading the result. If the observed coefficient falls inside the null
# distribution, the data cannot distinguish the effect from the artifact, and
# its sign should not be interpreted. If it falls outside, the association is
# stronger than the construction alone produces.
message(sprintf("\n[7/12] Null simulation, %d permutations...", N_SIM))

obs_fit  <- glm(specs$S1$f, data = dz, family = Gamma(link = "log"))
obs_cf   <- coef(obs_fit)
obs_pland <- obs_cf[["PLAND_z"]]
obs_pd    <- obs_cf[["PD_z"]]
obs_int   <- obs_cf[[grep(":", names(obs_cf), value = TRUE)]]

null_mat <- matrix(NA_real_, nrow = N_SIM, ncol = 3,
                   dimnames = list(NULL, c("PLAND", "PD", "Interaction")))

set.seed(GLOBAL_SEED)
pb_step <- max(1, floor(N_SIM / 10))

for (i in seq_len(N_SIM)) {
  idx <- sample.int(nrow(dz))
  ca_star <- dz$CA_Forest[idx]
  np_star <- dz$NP_Forest[idx]

  sim <- dz %>%
    mutate(PLAND_z = z(100 * ca_star / AOH_ha),
           PD_z    = z(np_star / (AOH_ha / 100)))

  fit <- tryCatch(
    glm(specs$S1$f, data = sim, family = Gamma(link = "log")),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(fit)) next

  cf <- coef(fit)
  ir <- grep(":", names(cf), value = TRUE)
  null_mat[i, ] <- c(cf[["PLAND_z"]], cf[["PD_z"]],
                     if (length(ir)) cf[[ir]] else NA_real_)

  if (i %% pb_step == 0) message(sprintf("   %d%%", round(100 * i / N_SIM)))
}

# Two quantities are computed, because they answer different questions.
#
# emp_p is the two-sided empirical p-value: the proportion of null draws whose
# absolute value is at least as large as the observed one. It compares
# magnitudes against a null that is not centered on zero (the ratio
# construction alone produces a negative coefficient), so it can return a small
# value for an observation that lies inside the body of an asymmetric null.
#
# pct_below is the percentile of the observed value within the null: the
# proportion of null draws at or below it.
#
# Decision rule, applied to every term: a coefficient is indistinguishable from
# the ratio artifact when it lies inside the central 95% interval of the null,
# between its 2.5th and 97.5th percentiles. The interval is used rather than
# emp_p because it does not assume a null symmetric about zero. Terms on which
# the two readings disagree are flagged in the column Rules_agree.
emp_p <- function(obs, null_vec) {
  nv <- null_vec[is.finite(null_vec)]
  (sum(abs(nv) >= abs(obs)) + 1) / (length(nv) + 1)
}

pct_below <- function(obs, null_vec) {
  nv <- null_vec[is.finite(null_vec)]
  100 * mean(nv <= obs)
}

null_summary <- data.frame(
  Term = c("PLAND_Forest", "PD_Forest", "PLAND x PROX interaction"),
  Observed = round(c(obs_pland, obs_pd, obs_int), 4),
  Null_mean = round(apply(null_mat, 2, mean, na.rm = TRUE), 4),
  Null_q2.5 = round(apply(null_mat, 2, quantile, 0.025, na.rm = TRUE), 4),
  Null_q97.5 = round(apply(null_mat, 2, quantile, 0.975, na.rm = TRUE), 4),
  Percentile_of_observed = round(c(pct_below(obs_pland, null_mat[, "PLAND"]),
                                   pct_below(obs_pd,    null_mat[, "PD"]),
                                   pct_below(obs_int,   null_mat[, "Interaction"])), 2),
  Empirical_p_two_sided = round(c(emp_p(obs_pland, null_mat[, "PLAND"]),
                                  emp_p(obs_pd,    null_mat[, "PD"]),
                                  emp_p(obs_int,   null_mat[, "Interaction"])), 4)
)
null_summary$Verdict <- ifelse(
  null_summary$Observed >= null_summary$Null_q2.5 &
    null_summary$Observed <= null_summary$Null_q97.5,
  "INSIDE the central 95% interval of the null: indistinguishable from the ratio artifact; do not interpret the sign",
  "OUTSIDE the central 95% interval of the null: stronger than the construction alone produces; may be interpreted with the caveat"
)
null_summary$Rules_agree <- ifelse(
  grepl("^INSIDE", null_summary$Verdict) & null_summary$Empirical_p_two_sided < 0.05,
  "no: inside the interval but two-sided p < 0.05; treat as the borderline case",
  "yes")
# Permutations whose fit ended with an error or a warning (usually a failure to
# converge) were skipped; this column counts the ones that entered the null.
null_summary$Permutations_used <- as.integer(colSums(is.finite(null_mat)))

write.csv(null_summary,
          file.path(output_path, "TableS_PartV_NullSimulation.csv"),
          row.names = FALSE)
print(null_summary, row.names = FALSE)

# Figure S9. Panel labels are set with a named vector rather than recode(),
# because car, attached for vif(), masks dplyr::recode().
term_labels <- c(PLAND       = "Forest cover (PLAND)",
                 PD          = "Patch density (PD)",
                 Interaction = "Proximity \u00d7 forest cover (PROX_MN \u00d7 PLAND)")

null_long <- as.data.frame(null_mat) %>%
  tidyr::pivot_longer(everything(), names_to = "Term", values_to = "Coefficient") %>%
  dplyr::filter(is.finite(Coefficient)) %>%
  dplyr::mutate(Term = factor(unname(term_labels[Term]),
                              levels = unname(term_labels)))

obs_long <- data.frame(
  Term = factor(unname(term_labels), levels = unname(term_labels)),
  Observed = c(obs_pland, obs_pd, obs_int)
)

p_null <- ggplot(null_long, aes(x = Coefficient)) +
  geom_histogram(bins = 60, fill = "grey70", colour = NA) +
  geom_vline(data = obs_long, aes(xintercept = Observed),
             colour = "#D55E00", linewidth = 0.8) +
  facet_wrap(~ Term, scales = "free", ncol = 1) +
  labs(x = "Standardized coefficient under the permutation null",
       y = "Number of permutations") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey25", colour = NA),
        strip.text = element_text(colour = "white", face = "bold"))

ggsave("Figure_PV_01_Null_Distribution.tiff", p_null, path = output_path,
       width = 140, height = 150, units = "mm", dpi = 500,
       device = "tiff", compression = "lzw", bg = "white")
message("   Figure_PV_01_Null_Distribution.tiff exported (Figure S9).")


# 9. AICc advantage of the coupled model under the null [8/12] -----------------
# The retained model S1 beats S4, the specification without window-normalized
# predictors, by a wide margin on AICc. That margin cannot by itself justify
# keeping S1: a predictor built from the response improves the fit to the
# response whether or not it carries ecological information.
#
# The same AICc difference is therefore recomputed under the permutation null.
# The numerators are shuffled and the denominator is held fixed, so S1 keeps its
# construction and loses any real information. If the null routinely
# reproduces the observed advantage, the advantage measures the construction,
# not the ecology, and the simpler specification is the one to report.
message("\n[8/12] AICc advantage of the coupled model, under the null...")

aic_obs_S1 <- MuMIn::AICc(glm(specs$S1$f, data = dz, family = Gamma(link = "log")))
aic_obs_S4 <- MuMIn::AICc(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))
delta_obs  <- aic_obs_S4 - aic_obs_S1

N_SIM_AIC <- min(N_SIM, N_SIM_AIC_CAP)
delta_null <- rep(NA_real_, N_SIM_AIC)

set.seed(GLOBAL_SEED)
for (i in seq_len(N_SIM_AIC)) {
  idx <- sample.int(nrow(dz))
  sim <- dz %>%
    mutate(PLAND_z = z(100 * dz$CA_Forest[idx] / AOH_ha),
           PD_z    = z(dz$NP_Forest[idx] / (AOH_ha / 100)))
  a1 <- tryCatch(MuMIn::AICc(glm(specs$S1$f, data = sim, family = Gamma(link = "log"))),
                 error = function(e) NA_real_, warning = function(w) NA_real_)
  a4 <- tryCatch(MuMIn::AICc(glm(specs$S4$f, data = sim, family = Gamma(link = "log"))),
                 error = function(e) NA_real_, warning = function(w) NA_real_)
  delta_null[i] <- a4 - a1
}

dn <- delta_null[is.finite(delta_null)]
aic_table <- data.frame(
  Quantity = c("AICc of the retained model (S1)",
               "AICc of the specification free of coupling (S4)",
               "Observed advantage of S1, in AICc units",
               "Mean advantage of S1 under the null",
               "Null 95% interval, lower",
               "Null 95% interval, upper",
               "Proportion of permutations reaching the observed advantage",
               "Permutations run",
               "Permutations in which both models fitted without error or warning"),
  Value = round(c(aic_obs_S1, aic_obs_S4, delta_obs, mean(dn),
                  quantile(dn, 0.025), quantile(dn, 0.975),
                  mean(dn >= delta_obs), N_SIM_AIC, length(dn)), 4)
)
write.csv(aic_table, file.path(output_path, "TableS_PartV_AICc_Under_Null.csv"),
          row.names = FALSE)
print(aic_table, row.names = FALSE)

if (mean(dn >= delta_obs) > 0.05) {
  message("   The null reproduces the observed AICc advantage. The information")
  message("   criterion measures the construction of the ratio, not ecology,")
  message("   and cannot justify keeping the coupled model.")
} else {
  message("   The observed AICc advantage exceeds what the construction alone")
  message("   produces. The coupled model carries information beyond the ratio,")
  message("   which does not license interpreting the coefficient of the ratio.")
}


# 10. Is forest proximity habitat amount measured again? [9/12] ----------------
# The proximity index is not normalized by the response, so it is not exposed
# to the ratio artifact. It is, however, built from the areas of the patches
# within a search radius divided by the square of their distances, so a
# landscape with more forest tends to score higher on it by construction. This
# is confounding with habitat amount, a different problem from the ratio
# artifact, and the concern raised by Fahrig (2013): an apparent configuration
# effect that is habitat amount in disguise.
#
# Two steps. First, the association between proximity and four measures of
# forest amount. Second, the coefficient of proximity with the admissible
# amount measures entered as covariates. The absolute class area is bounded by
# the response, so conditioning on it would be circular; it appears in the
# correlations but not among the covariates. The mean patch area and the number
# of patches are not divided by the response and are used as covariates.
message("\n[9/12] Is forest proximity confounded with habitat amount?")

amount_cor <- data.frame(
  Amount_measure = c("PLAND_Forest", "CA_Forest", "Mean_AREA_Forest_ha", "NP_Forest"),
  Pearson  = round(c(cor(d$PROX_MN_Forest, d$PLAND_Forest),
                     cor(d$PROX_MN_Forest, d$CA_Forest),
                     cor(d$PROX_MN_Forest, d$Mean_AREA_Forest_ha),
                     cor(d$PROX_MN_Forest, d$NP_Forest)), 4),
  Spearman = round(c(cor(d$PROX_MN_Forest, d$PLAND_Forest, method = "spearman"),
                     cor(d$PROX_MN_Forest, d$CA_Forest, method = "spearman"),
                     cor(d$PROX_MN_Forest, d$Mean_AREA_Forest_ha, method = "spearman"),
                     cor(d$PROX_MN_Forest, d$NP_Forest, method = "spearman")), 4),
  Contains_the_response = c("yes, in its denominator", "bounded by it", "no", "no")
)

dz$logNP_z <- z(log(dz$NP_Forest))

prox_models <- list(
  "No amount covariate"                 = AOH_ha ~ ENN_z + FRAC_z + PROX_z,
  "Plus mean forest patch area"         = AOH_ha ~ ENN_z + FRAC_z + PROX_z + AREAMN_z,
  "Plus log number of forest patches"   = AOH_ha ~ ENN_z + FRAC_z + PROX_z + logNP_z,
  "Plus both"                           = AOH_ha ~ ENN_z + FRAC_z + PROX_z + AREAMN_z + logNP_z
)

prox_table <- purrr::imap_dfr(prox_models, function(fm, lab) {
  m  <- glm(fm, data = dz, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- confint.default(m)
  data.frame(
    Specification = lab,
    PROX_beta     = round(cf["PROX_z", "Estimate"], 4),
    CI_lower      = round(ci["PROX_z", 1], 4),
    CI_upper      = round(ci["PROX_z", 2], 4),
    P_value       = signif(cf["PROX_z", 4], 4),
    OHE_ratio_per_SD = round(exp(cf["PROX_z", "Estimate"]), 3),
    stringsAsFactors = FALSE
  )
})

write.csv(amount_cor, file.path(output_path, "TableS_PartV_Proximity_Amount_Correlations.csv"),
          row.names = FALSE)
write.csv(prox_table, file.path(output_path, "TableS_PartV_Proximity_Conditioned.csv"),
          row.names = FALSE)
print(amount_cor, row.names = FALSE)
print(prox_table, row.names = FALSE)

surv <- all(prox_table$CI_lower * prox_table$CI_upper > 0)
shrink <- 1 - min(prox_table$PROX_beta) / prox_table$PROX_beta[1]
if (surv) {
  message(sprintf(paste0("   The 95%% CI of the proximity coefficient excludes zero under every ",
                         "amount covariate; the coefficient shrinks by at most %.0f%%."),
                  100 * shrink))
} else {
  message("   The 95% CI of the proximity coefficient includes zero under at least")
  message("   one amount covariate, so proximity cannot be separated from habitat amount.")
}


# 11. Decision table [10/12] ---------------------------------------------------
message("\n[10/12] Building the decision table...")

s6 <- spec_table %>% dplyr::filter(grepl("^S6", Specification))
s4 <- spec_table %>% dplyr::filter(grepl("^S4", Specification))

decision <- data.frame(
  Question = c(
    "Q1. Is the negative coefficient of forest cover more than the ratio artifact?",
    "Q2. Does the effect of forest proximity survive without any coupled predictor?",
    "Q3. Does the interaction survive with an amount measure free of the response?",
    "Q4. Is the coefficient of patch density more than the ratio artifact?",
    "Q5. Can AICc arbitrate between the coupled and the uncoupled specification?"
  ),
  Evidence = c(
    sprintf(paste0("Observed %.3f; null central 95%% interval [%.3f, %.3f]; ",
                   "observed at the %.1fth percentile of the null; ",
                   "two-sided empirical p = %.4f"),
            null_summary$Observed[1], null_summary$Null_q2.5[1],
            null_summary$Null_q97.5[1], null_summary$Percentile_of_observed[1],
            null_summary$Empirical_p_two_sided[1]),
    sprintf("S4 (no coupled predictor): PROX beta = %.3f [%.3f, %.3f], p = %s",
            s4$PROX_beta, s4$PROX_CI_low, s4$PROX_CI_high, format(s4$PROX_p)),
    sprintf("S6 (amount as mean patch area): interaction beta = %.3f [%.3f, %.3f], p = %s",
            s6$Interaction_beta, s6$Interaction_CI_low,
            s6$Interaction_CI_high, format(s6$Interaction_p)),
    sprintf(paste0("Observed %.3f; null central 95%% interval [%.3f, %.3f]; ",
                   "observed at the %.1fth percentile of the null; ",
                   "two-sided empirical p = %.4f"),
            null_summary$Observed[2], null_summary$Null_q2.5[2],
            null_summary$Null_q97.5[2], null_summary$Percentile_of_observed[2],
            null_summary$Empirical_p_two_sided[2]),
    sprintf("Observed advantage of the coupled model %.2f AICc units; the null reaches it in %.1f%% of permutations",
            delta_obs, 100 * mean(dn >= delta_obs))
  ),
  Verdict = c(
    null_summary$Verdict[1],
    ifelse(s4$PROX_CI_low * s4$PROX_CI_high > 0,
           "YES: the proximity effect does not require a coupled predictor",
           "NO: the proximity effect depends on the coupled predictors; report this"),
    ifelse(!is.na(s6$Interaction_p) &&
             s6$Interaction_CI_low * s6$Interaction_CI_high > 0,
           "YES: the threshold is reproduced without any area-normalized predictor",
           "NO: the threshold is not reproduced; report it as conditional on PLAND"),
    null_summary$Verdict[2],
    ifelse(mean(dn >= delta_obs) > 0.05,
           paste("NO: the construction of the ratio reproduces the observed AICc",
                 "advantage, so the criterion measures access to the response and",
                 "not ecological information. Report the uncoupled specification."),
           paste("The observed advantage exceeds the construction, which still does",
                 "not license interpreting the coefficient of the ratio."))
  ),
  stringsAsFactors = FALSE
)

write.csv(decision, file.path(output_path, "TableS_PartV_Decision.csv"),
          row.names = FALSE)

message("\n============================================================")
for (i in seq_len(nrow(decision))) {
  message("  ", decision$Question[i])
  message("     Evidence: ", decision$Evidence[i])
  message("     Verdict:  ", decision$Verdict[i])
  message("")
}
message("============================================================")


# 12. The reported model: Table 4 and Figure 6 [11/12] -------------------------
# The manuscript reports S4, the specification without window-normalized
# predictors, so its full coefficient table is exported here.
message("\n[11/12] Coefficients and figure of the reported model (S4)...")

m_free <- glm(specs$S4$f, data = dz, family = Gamma(link = "log"))
cf_free <- summary(m_free)$coefficients
ci_free <- confint.default(m_free)
free_table <- data.frame(
  Predictor = rownames(cf_free),
  Estimate  = round(cf_free[, "Estimate"], 4),
  SE        = round(cf_free[, "Std. Error"], 4),
  CI_lower  = round(ci_free[, 1], 4),
  CI_upper  = round(ci_free[, 2], 4),
  P_value   = signif(cf_free[, 4], 4),
  OHE_ratio_per_SD = round(exp(cf_free[, "Estimate"]), 3),
  row.names = NULL
)
write.csv(free_table, file.path(output_path, "TableS_PartV_Uncoupled_Model.csv"),
          row.names = FALSE)
print(free_table, row.names = FALSE)
message(sprintf("   Maximum VIF in this specification: %.2f",
                spec_table$Max_VIF[spec_table$Specification == specs$S4$label]))

# Figure 6: standardized coefficients of the reported model with 95% CI.
LAB <- c(PROX_z = "Forest proximity (PROX_MN)",
         FRAC_z = "Shape complexity (FRAC_MN)",
         ENN_z  = "Isolation (ENN_MN)")

coef_plot_df <- free_table %>%
  dplyr::filter(Predictor != "(Intercept)") %>%
  dplyr::mutate(Label = unname(LAB[Predictor]),
                Label = factor(Label, levels = rev(unname(LAB))),
                Supported = ifelse(CI_lower * CI_upper > 0,
                                   "95% CI excludes zero", "95% CI includes zero"))

p_coef <- ggplot(coef_plot_df, aes(x = Estimate, y = Label, colour = Supported)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.14,
                linewidth = 0.7, orientation = "y") +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c("95% CI excludes zero" = "#0072B2",
                                 "95% CI includes zero" = "#999999")) +
  labs(x = expression("Standardized coefficient (" * beta * ") with 95% CI"),
       y = NULL, colour = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top",
        axis.text.y     = element_text(colour = "black"))

ggsave("Figure_PV_02_Reported_Model_Coefficients.tiff", p_coef, path = output_path,
       width = 140, height = 75, units = "mm", dpi = 500,
       device = "tiff", compression = "lzw", bg = "white")
message("   Figure_PV_02_Reported_Model_Coefficients.tiff exported (Figure 6).")


# 13. Out-of-sample error and BIC [12/12] --------------------------------------
# Yates et al. (2023) recommend exact leave-one-out (LOO) cross-validation over
# random k-fold with k < 10, and warn that a random split places units of the
# same taxon in the training and the test set at once; leave-one-genus-out
# (LGO) is reported as well. Brewer, Butler and Cooksley (2016) show that the
# relative behavior of AICc and BIC depends on unobserved heterogeneity, so both
# criteria are reported.
#
# The coupled specification S1 will appear to predict better, because two of
# its predictors are built from the response. This is the artifact that Table
# S28c documents for AICc: out-of-sample error cannot arbitrate between the two
# specifications either.
message("\n[12/12] Out-of-sample error and BIC of the reported model...")

loo_rmse <- function(fml, dat) {
  err <- numeric(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    m <- tryCatch(glm(fml, data = dat[-i, ], family = Gamma(link = "log")),
                  error = function(e) NULL)
    err[i] <- if (is.null(m)) NA_real_ else
      (dat$AOH_ha[i] - predict(m, newdata = dat[i, , drop = FALSE],
                               type = "response"))^2
  }
  sqrt(mean(err, na.rm = TRUE))
}
lgo_rmse <- function(fml, dat) {
  err <- numeric(0)
  for (g in unique(dat$GENUS)) {
    tr <- dat[dat$GENUS != g, ]; te <- dat[dat$GENUS == g, ]
    m  <- tryCatch(glm(fml, data = tr, family = Gamma(link = "log")),
                   error = function(e) NULL)
    p  <- if (is.null(m)) rep(NA_real_, nrow(te)) else
      predict(m, newdata = te, type = "response")
    err <- c(err, (te$AOH_ha - p)^2)
  }
  sqrt(mean(err, na.rm = TRUE))
}

fit_free <- glm(specs$S4$f, data = dz, family = Gamma(link = "log"))
fit_coup <- glm(specs$S1$f, data = dz, family = Gamma(link = "log"))

cv_table <- data.frame(
  Specification = c("Reported model, free of coupling (S4)",
                    "Coupled model (S1)"),
  k        = c(length(coef(fit_free)) + 1, length(coef(fit_coup)) + 1),
  AICc     = round(c(MuMIn::AICc(fit_free), MuMIn::AICc(fit_coup)), 2),
  BIC      = round(c(stats::BIC(fit_free),  stats::BIC(fit_coup)), 2),
  LOO_RMSE_ha = round(c(loo_rmse(specs$S4$f, dz), loo_rmse(specs$S1$f, dz)), 1),
  LGO_RMSE_ha = round(c(lgo_rmse(specs$S4$f, dz), lgo_rmse(specs$S1$f, dz)), 1),
  Note = c("no predictor built from the response",
           "two predictors built from the response; every criterion favors it for that reason"),
  stringsAsFactors = FALSE
)
write.csv(cv_table, file.path(output_path, "TableS_PartV_CV_And_BIC.csv"),
          row.names = FALSE)
print(cv_table, row.names = FALSE)
message("   Read this table with Table S28c: neither the information criterion")
message("   nor the out-of-sample error can arbitrate between the two,")
message("   because both reward access to the response.")


# 14. Session information ------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part05.txt"))
