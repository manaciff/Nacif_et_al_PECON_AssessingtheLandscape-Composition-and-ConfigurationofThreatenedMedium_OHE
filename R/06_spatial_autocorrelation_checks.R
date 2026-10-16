# ==============================================================================
# 06_spatial_autocorrelation_checks.R
# Spatial autocorrelation checks of the variation partitioning and of the
# reported model
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   Script 02 models positive spatial autocorrelation only, with the
#   distance-based Moran eigenvector maps (dbMEMs) of positive eigenvalue. This
#   script checks that choice and the residuals of the models the manuscript
#   reports.
#   1. Tests the dbMEMs of negative eigenvalue together, by permutation, in the
#      same way script 02 tests the positive ones. This is a diagnostic: no
#      eigenvector is selected from it. Beside each p-value it writes the
#      Sidak-corrected value, 1 - (1 - p)^2, that a joint screen of the two
#      signs would use (Borcard et al. 2011, Section 7.4.4).
#   2. Tests Moran's I, two-sided, of log(envelope area) and of the residuals
#      of the variation partitioning models and of the reported model (Table
#      4), with the spatial weights of the dbMEMs. Regression residuals are not
#      exchangeable, so each series gets the test that suits it: a permutation
#      test for log(OHE); the exact moments of Moran's I for least-squares
#      residuals (spdep::lm.morantest) for the partition models; and a
#      permutation test on the simulated quantile residuals of DHARMa, which
#      are independent and uniform under a correct model, for the Gamma GLM.
#   3. Refits the reported model with the eigenvectors retained by script 02,
#      as script 03 does for the a priori model (Table S34).
#   4. Tests the same simulated quantile residuals of the reported model,
#      without and with the eigenvectors, for dispersion, uniformity and
#      outliers, as script 03 does for the a priori model (Table S17).
#
# Input
#   Dados/Processados/Data_Raw_WithCoords.csv                  written by script 01
#   Outputs/Manuscrito/PartII/TableII_1_VIF_VP.csv             written by script 02
#   Outputs/Manuscrito/PartII/TableII_2b_MEM_Vectors.csv       written by script 02
#   Outputs/Manuscrito/PartV/TableS_PartV_Uncoupled_Model.csv  written by script 05
#
# Output   (in Outputs/Manuscrito/PartII/)
#   TableII_6_MEM_Global_Tests.csv            (Table S38b)
#   TableII_6b_Reported_Model_With_MEMs.csv   (Table S34b)
#   TableII_6c_Moran_I_Residuals.csv          (Tables S34b and S38b)
#   TableII_6d_DHARMa_Reported_Model.csv      (Table S34b)
#   session_info_part06.txt
#
# Run time
#   Less than a minute.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)         # project-relative paths
  library(tidyverse)    # data wrangling
  library(vegan)        # rda(), anova.cca(), RsquareAdj()
  library(adespatial)   # dbmem()
  library(spdep)        # moran.mc(), lm.morantest()
  library(DHARMa)       # simulateResiduals()
  library(car)          # vif()
})

GLOBAL_SEED      <- 123
N_PERM           <- 9999
N_UNITS_EXPECTED <- 67

# Taxon names as script 01 writes them.
TAXA <- c("A_guariba", "B_arachnoides", "B_crinitus", "B_hypoxanthus",
          "L_wiedii", "M_tridactyla", "Mazama", "P_concolor", "T_pecari")

# Landscape predictors divided by the area of the sampling unit (script 02,
# Section 9.1).
RATIO_METRICS <- c("PD_Forest", "PD_Agropecuaria", "PLAND_Herbaceous")

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartII")
vif_file    <- file.path(output_path, "TableII_1_VIF_VP.csv")
mem_file    <- file.path(output_path, "TableII_2b_MEM_Vectors.csv")
table4_file <- here("Outputs", "Manuscrito", "PartV", "TableS_PartV_Uncoupled_Model.csv")


# 2. Data and inputs [1/5] -----------------------------------------------------
message("\n[1/5] Loading the data and the outputs of scripts 02 and 05...")

for (f in c(vif_file, mem_file, table4_file)) {
  if (!file.exists(f)) stop("Input not found: ", f, "\nRun scripts 02 and 05 first.")
}

d <- read_csv(file.path(data_path, "Data_Raw_WithCoords.csv"),
              show_col_types = FALSE) %>%
  mutate(UA_ID = as.integer(UA_ID))

if (nrow(d) != N_UNITS_EXPECTED) {
  stop(sprintf("Data_Raw_WithCoords.csv holds %d rows, expected %d. Rerun script 01.",
               nrow(d), N_UNITS_EXPECTED))
}
if (!all(d$SPECIES %in% TAXA)) {
  stop("Unrecognized taxon name(s): ",
       paste(setdiff(unique(d$SPECIES), TAXA), collapse = ", "),
       ". Rerun script 01, which writes the canonical names.")
}
if (anyNA(d$coord_x) || anyNA(d$coord_y)) {
  stop("Some units lack centroid coordinates. Rerun script 01.")
}

pred_vp <- read_csv(vif_file, show_col_types = FALSE)$Variable

mem_sel  <- read_csv(mem_file, show_col_types = FALSE) %>%
  mutate(UA_ID = as.integer(UA_ID))
mem_cols <- setdiff(names(mem_sel), c("SPECIES", "UA_ID"))

# The key (SPECIES, UA_ID) makes the join independent of row order.
d <- d %>% left_join(mem_sel, by = c("SPECIES", "UA_ID"))
if (nrow(d) != N_UNITS_EXPECTED || anyNA(d[, mem_cols])) {
  stop("The eigenvectors of script 02 did not join one to one to the 67 units.")
}

y  <- log(d$AOH_unit_area_ha)
xy <- as.matrix(d[, c("coord_x", "coord_y")])


# 3. Eigenvectors and spatial weights [2/5] ------------------------------------
# The construction of script 02. The Moran's I tests use the spatial weights of
# the dbMEMs, so the tests and the eigenvectors describe the same neighborhood.
message("\n[2/5] Building the positive and negative dbMEMs...")

mem_pos <- dbmem(xy, MEM.autocor = "positive", silent = TRUE)
mem_neg <- dbmem(xy, MEM.autocor = "negative", silent = TRUE)
lw      <- attr(mem_pos, "listw")

# Each eigenvector retained by script 02 must be the dbMEM of the same name
# rebuilt here. The sign of an eigenvector is arbitrary.
for (m in mem_cols) {
  if (!m %in% names(mem_pos) || abs(cor(d[[m]], mem_pos[[m]])) < 0.9999) {
    stop("Eigenvector ", m, " of TableII_2b_MEM_Vectors.csv does not match ",
         "the dbMEM rebuilt here. Rerun script 02.")
  }
}
message(sprintf("   %d positive and %d negative dbMEMs; retained by script 02: %s.",
                ncol(mem_pos), ncol(mem_neg), paste(mem_cols, collapse = ", ")))


# 4. Global tests of the positive and of the negative eigenvectors [3/5] -------
# Each set is tested as script 02 tests the positive one: an RDA of log(OHE) on
# all the eigenvectors of the set, tested by permutation.
message("\n[3/5] Global tests of the positive and of the negative eigenvectors...")

global_test <- function(M, label) {
  fit <- rda(y ~ ., data = as.data.frame(M))
  set.seed(GLOBAL_SEED)
  a <- anova(fit, permutations = N_PERM)
  p <- a$`Pr(>F)`[1]
  data.frame(Eigenvectors       = label,
             N_eigenvectors     = ncol(M),
             R2                 = round(RsquareAdj(fit)$r.squared, 4),
             Adj_R2             = round(RsquareAdj(fit)$adj.r.squared, 4),
             F                  = round(a$F[1], 4),
             p                  = p,
             p_Sidak_two_signs  = round(1 - (1 - p)^2, 4),
             Permutations       = N_PERM)
}

global_tab <- bind_rows(
  global_test(mem_pos, "Positive, modeled in script 02"),
  global_test(mem_neg, "Negative, diagnostic only")
)
write.csv(global_tab, file.path(output_path, "TableII_6_MEM_Global_Tests.csv"),
          row.names = FALSE)
print(global_tab, row.names = FALSE)


# 5. Moran's I and DHARMa tests of the residuals [4/5] -------------------------
# Two-sided tests, so that positive and negative autocorrelation are both
# detectable.
message("\n[4/5] Moran's I and DHARMa tests of the residuals...")

z <- function(x) as.numeric(scale(x))

X1  <- as.data.frame(lapply(d[, pred_vp], z))
X1c <- X1[, setdiff(pred_vp, RATIO_METRICS), drop = FALSE]
X2  <- as.data.frame(d[, mem_cols])

# The reported model (Table 4) and the same model with the eigenvectors.
dz <- d %>%
  mutate(ENN_z  = z(Mean_ENN_Forest_m),
         FRAC_z = z(Mean_FRAC_Forest),
         PROX_z = z(PROX_MN_Forest))
f_t4     <- AOH_unit_area_ha ~ ENN_z + FRAC_z + PROX_z
f_t4_mem <- as.formula(paste("AOH_unit_area_ha ~ ENN_z + FRAC_z + PROX_z +",
                             paste(mem_cols, collapse = " + ")))
m_t4     <- glm(f_t4,     data = dz, family = Gamma(link = "log"))
m_t4_mem <- glm(f_t4_mem, data = dz, family = Gamma(link = "log"))

# Permutation test of a series that is exchangeable under no autocorrelation.
moran_perm <- function(x, label, test) {
  set.seed(GLOBAL_SEED)
  mc <- moran.mc(x, lw, nsim = N_PERM, alternative = "two.sided")
  data.frame(Series = label, Test = test,
             Moran_I = round(unname(mc$statistic), 4),
             Expected_I = round(-1 / (length(x) - 1), 4),
             p_two_sided = round(mc$p.value, 4))
}

# Exact moments of Moran's I for least-squares residuals.
moran_lm <- function(fit, label) {
  mt <- lm.morantest(fit, lw, alternative = "two.sided")
  data.frame(Series = label, Test = "lm.morantest, exact moments",
             Moran_I = round(unname(mt$estimate[1]), 4),
             Expected_I = round(unname(mt$estimate[2]), 4),
             p_two_sided = round(mt$p.value, 4))
}

# Simulated quantile residuals of a GLM, as in the diagnostics of script 03. The
# same simulations give the Moran's I tests below and the dispersion, uniformity
# and outlier tests of the reported model.
sim_t4     <- simulateResiduals(fittedModel = m_t4,     n = 1000, seed = GLOBAL_SEED)
sim_t4_mem <- simulateResiduals(fittedModel = m_t4_mem, n = 1000, seed = GLOBAL_SEED)

moran_tab <- bind_rows(
  moran_perm(y, "log(OHE)", "permutation"),
  moran_lm(lm(y ~ ., data = X1),
           "Residuals, log(OHE) on the landscape set (10 predictors)"),
  moran_lm(lm(y ~ ., data = cbind(X1, X2)),
           "Residuals, log(OHE) on the landscape set and the MEMs (full partition)"),
  moran_lm(lm(y ~ ., data = cbind(X1c, X2)),
           "Residuals, log(OHE) on the landscape set without the three ratio metrics and the MEMs (Table S38)"),
  moran_perm(residuals(sim_t4),
             "Quantile residuals, reported model (Table 4)",
             "permutation of DHARMa residuals"),
  moran_perm(residuals(sim_t4_mem),
             "Quantile residuals, reported model with the MEMs (Table S34b)",
             "permutation of DHARMa residuals")
) %>%
  mutate(Permutations = ifelse(grepl("permutation", Test), N_PERM, NA))
write.csv(moran_tab, file.path(output_path, "TableII_6c_Moran_I_Residuals.csv"),
          row.names = FALSE)
print(moran_tab[, c("Series", "Moran_I", "p_two_sided")], row.names = FALSE)

# Dispersion, uniformity and outliers of the same simulated residuals: the tests
# of Table S17, in the same order and with the same defaults.
dharma_tests <- function(sim, label) {
  disp <- testDispersion(sim, plot = FALSE)
  unif <- testUniformity(sim, plot = FALSE)
  outl <- testOutliers(sim, plot = FALSE)
  data.frame(Model       = label,
             Test        = c("DHARMa dispersion", "DHARMa uniformity", "DHARMa outliers"),
             Statistic   = round(c(disp$statistic[[1]], unif$statistic[[1]],
                                   outl$statistic[[1]]), 4),
             P_value     = round(c(disp$p.value, unif$p.value, outl$p.value), 4),
             Simulations = 1000)
}

dharma_tab <- bind_rows(
  dharma_tests(sim_t4,     "Reported model (Table 4)"),
  dharma_tests(sim_t4_mem, "Reported model with the MEMs (Table S34b)")
)
write.csv(dharma_tab, file.path(output_path, "TableII_6d_DHARMa_Reported_Model.csv"),
          row.names = FALSE)
print(dharma_tab, row.names = FALSE)


# 6. The reported model with the eigenvectors [5/5] ----------------------------
# Script 03 refits the a priori model with the eigenvectors (Table S34). The
# reported model is refitted here in the same way; Table 4 itself is unchanged.
message("\n[5/5] Refitting the reported model (Table 4) with the eigenvectors...")

cf0 <- summary(m_t4)$coefficients
t4  <- read_csv(table4_file, show_col_types = FALSE)
i4  <- match(t4$Predictor, rownames(cf0))
if (anyNA(i4) || max(abs(round(cf0[i4, "Estimate"], 4) - t4$Estimate)) > 1e-4) {
  stop("The reported model fitted here does not reproduce ",
       "TableS_PartV_Uncoupled_Model.csv. Rerun script 05.")
}

cf1  <- summary(m_t4_mem)$coefficients
ci1  <- confint.default(m_t4_mem)
vif1 <- car::vif(m_t4_mem)
rows <- rownames(cf1)
i0   <- match(rows, rownames(cf0))

t4_mem_tab <- data.frame(
  Predictor                  = rows,
  Estimate_Table4            = round(cf0[i0, "Estimate"], 4),
  Estimate_with_MEMs         = round(cf1[, "Estimate"], 4),
  Absolute_change            = round(abs(cf1[, "Estimate"] - cf0[i0, "Estimate"]), 4),
  SE_with_MEMs               = round(cf1[, "Std. Error"], 4),
  CI_lower_with_MEMs         = round(ci1[, 1], 4),
  CI_upper_with_MEMs         = round(ci1[, 2], 4),
  P_Table4                   = signif(cf0[i0, 4], 4),
  P_with_MEMs                = signif(cf1[, 4], 4),
  OHE_ratio_per_SD_with_MEMs = round(exp(cf1[, "Estimate"]), 3),
  VIF_with_MEMs              = round(unname(vif1[match(rows, names(vif1))]), 3),
  row.names = NULL
)
write.csv(t4_mem_tab,
          file.path(output_path, "TableII_6b_Reported_Model_With_MEMs.csv"),
          row.names = FALSE)
print(t4_mem_tab[, c("Predictor", "Estimate_Table4", "Estimate_with_MEMs",
                     "CI_lower_with_MEMs", "CI_upper_with_MEMs", "P_with_MEMs")],
      row.names = FALSE)

flips <- with(t4_mem_tab, !is.na(Estimate_Table4) & Predictor != "(Intercept)" &
                sign(Estimate_Table4) != sign(Estimate_with_MEMs))
if (any(flips)) {
  message("   WARNING. A coefficient of Table 4 changes sign with the eigenvectors.")
} else {
  message("   No coefficient of Table 4 changes sign with the eigenvectors.")
}


# 7. Outputs and session information --------------------------------------------
message("\nOutputs in: ", output_path)
message("    TableII_6_MEM_Global_Tests.csv            (Table S38b)")
message("    TableII_6b_Reported_Model_With_MEMs.csv   (Table S34b)")
message("    TableII_6c_Moran_I_Residuals.csv          (Tables S34b and S38b)")
message("    TableII_6d_DHARMa_Reported_Model.csv      (Table S34b)")

writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part06.txt"))
