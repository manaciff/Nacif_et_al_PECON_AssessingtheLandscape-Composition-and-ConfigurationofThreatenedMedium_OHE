# ==============================================================================
# 04_geometric_coupling_diagnostic.R
# Geometric coupling between the response and the landscape predictors
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   The response, the area of the occupied habitat envelope (OHE), is also the
#   area of the window in which the landscape metrics are computed. Patch
#   density (PD, patches per 100 ha) is therefore tied to the response by an
#   exact identity:
#
#       log(OHE) = log(NP) + log(100) - log(PD)
#
#   The script checks this identity on the 67 sampling units and fits a
#   reference mixed model beside three variants. M2 adds forest cover (PLAND)
#   as a control for habitat amount. M3 and M4 add a term that is tied to the
#   response by construction, the number of patches (NP) or the forest class
#   area (CA). These two circular variants show the circularity; they are not
#   tests of the associations.
#
# Input
#   Dados/Processados/Data_Raw_WithCoords.csv   written by script 01; the same
#                                               file read by scripts 03 and 05
#
# Output
#   Outputs/Manuscrito/TableS28_Geometric_Coupling.csv   (Table S26)
#
# Run time
#   A few seconds.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)    # project-relative paths
  library(readr)   # read_csv()
  library(dplyr)   # data manipulation
  library(lme4)    # linear mixed models
})

GLOBAL_SEED <- 123
set.seed(GLOBAL_SEED)

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


# 2. Data [1/5] ----------------------------------------------------------------
message("\n[1/5] Loading data...")

input_file <- file.path(data_path, "Data_Raw_WithCoords.csv")
if (!file.exists(input_file)) {
  stop("Data_Raw_WithCoords.csv not found in ", data_path,
       ". Run scripts 00 and 01 before this one.", call. = FALSE)
}

# AOH_ha is the area of the occupied habitat envelope (OHE), in hectares.
d <- read_csv(input_file, show_col_types = FALSE) %>%
  rename(AOH_ha = AOH_unit_area_ha) %>%
  mutate(log_AOH = log(AOH_ha),
         log_NP  = log(NP_Forest),
         log_CA  = log(pmax(CA_Forest, 1)),
         GENUS   = as.factor(GENUS))

# Guard against a wrong input file. With the corrected data and the response at
# full precision, mean log(OHE) = 10.3735 and sd = 0.9392.
message(sprintf("   n = %d | mean log(OHE) = %.4f | sd = %.4f",
                nrow(d), mean(d$log_AOH), sd(d$log_AOH)))
if (abs(mean(d$log_AOH) - 10.3735) > 0.01) {
  warning("mean log(OHE) does not match the corrected dataset. ",
          "Confirm that script 01 has been run and that ",
          "Data_Raw_WithCoords.csv is the file it wrote, before ",
          "interpreting this table.")
}

# Standardize every variable that enters a model, so that the coefficients are
# comparable.
z_vars <- c("Mean_FRAC_Forest", "PD_Forest", "PROX_MN_Forest",
            "PLAND_Forest", "log_NP", "log_CA")
d <- d %>% mutate(across(all_of(z_vars), ~as.numeric(scale(.x)), .names = "{.col}_z"))


# 3. The four mixed models [2/5] -----------------------------------------------
message("[2/5] Fitting the four mixed models...")

base3 <- c("Mean_FRAC_Forest_z", "PD_Forest_z", "PROX_MN_Forest_z")

fit_mixed <- function(preds) {
  f <- as.formula(paste("log_AOH ~", paste(preds, collapse = " + "), "+ (1 | GENUS)"))
  # Maximum likelihood (REML = FALSE), so that the models are comparable.
  lme4::lmer(f, data = d, REML = FALSE)
}

model_set <- list(
  list(name = "M1. Reference model",
       preds = base3,
       note  = "reference"),
  list(name = "M2. + PLAND (proportion of forest)",
       preds = c(base3, "PLAND_Forest_z"),
       note  = "habitat-amount control; PLAND = 100 x CA / OHE"),
  list(name = "M3. + log(NP) (number of patches)",
       preds = c(base3, "log_NP_z"),
       note  = "circular: log(OHE) = log(NP) + log(100) - log(PD)"),
  list(name = "M4. + log(CA) (absolute forest area)",
       preds = c(base3, "log_CA_z"),
       note  = "circular: CA = PLAND x OHE / 100")
)

rows <- list()

for (spec in model_set) {
  m  <- fit_mixed(spec$preds)
  cf <- summary(m)$coefficients          # Estimate, Std. Error, t value
  # lme4 does not return p-values. The two-sided Wald test with the normal
  # approximation is used, as in script 03.
  for (v in spec$preds) {
    est <- cf[v, 1]; se <- cf[v, 2]; tv <- cf[v, 3]
    rows[[length(rows) + 1]] <- data.frame(
      Model     = spec$name,
      Predictor = sub("_z$", "", v),
      Estimate  = round(est, 4),
      SE        = round(se, 4),
      P_value   = round(2 * pnorm(abs(tv), lower.tail = FALSE), 4),
      Note      = spec$note,
      stringsAsFactors = FALSE
    )
  }
  # Intraclass correlation: share of the variance left by the fixed effects
  # that lies between genera.
  vc        <- as.data.frame(lme4::VarCorr(m))
  var_genus <- vc$vcov[vc$grp == "GENUS"]
  var_resid <- vc$vcov[vc$grp == "Residual"]
  rows[[length(rows) + 1]] <- data.frame(
    Model = spec$name, Predictor = "ICC (genus)",
    Estimate = round(var_genus / (var_genus + var_resid), 4),
    SE = NA_real_, P_value = NA_real_, Note = spec$note,
    stringsAsFactors = FALSE
  )
  message(sprintf("   %-38s ICC = %.4f", spec$name, var_genus / (var_genus + var_resid)))
}


# 4. Partial correlations of PROX_MN with log(OHE) [3/5] -----------------------
message("[3/5] Partial correlations...")

# Correlation between the residuals of x and y, each regressed on the controls.
partial_cor <- function(x, y, ctrl) {
  if (length(ctrl) == 0) return(cor(d[[x]], d[[y]]))
  fx <- as.formula(paste(x, "~", paste(ctrl, collapse = " + ")))
  fy <- as.formula(paste(y, "~", paste(ctrl, collapse = " + ")))
  cor(residuals(lm(fx, data = d)), residuals(lm(fy, data = d)))
}

controls <- list(
  list(v = character(0),                lab = "none",            note = "no control"),
  list(v = "PLAND_Forest",              lab = "PLAND",           note = "habitat-amount control"),
  list(v = "log_NP",                    lab = "log(NP)",         note = "circular control"),
  list(v = c("PLAND_Forest", "log_NP"), lab = "PLAND + log(NP)", note = "circular control")
)

for (ctl in controls) {
  r <- partial_cor("PROX_MN_Forest", "log_AOH", ctl$v)
  rows[[length(rows) + 1]] <- data.frame(
    Model = "Partial correlation of PROX_MN with log(OHE)",
    Predictor = paste("controlling for", ctl$lab),
    Estimate = round(r, 4), SE = NA_real_, P_value = NA_real_,
    Note = ctl$note, stringsAsFactors = FALSE
  )
  message(sprintf("   controlling for %-16s r = %.4f", ctl$lab, r))
}


# 5. Numerical check of the identity [4/5] -------------------------------------
message("[4/5] Algebraic identity check...")

# The identity is exact for the window in which the metrics were computed. The
# small discrepancy arises because the metrics are counted on 30 m pixels,
# whereas the response is the area of the polygon.
max_discrepancy <- max(abs(log(d$NP_Forest) + log(100) - log(d$PD_Forest) - d$log_AOH))
rows[[length(rows) + 1]] <- data.frame(
  Model = "Algebraic identity check",
  Predictor = "max |log(NP)+log(100)-log(PD) - log(OHE)|",
  Estimate = round(max_discrepancy, 4), SE = NA_real_, P_value = NA_real_,
  Note = paste("PD = NP/(area/100) holds exactly on the pixels;",
               "the residual is pixel against polygon area"),
  stringsAsFactors = FALSE
)
message(sprintf("   maximum absolute discrepancy across %d units: %.4f", nrow(d), max_discrepancy))


# 6. Export [5/5] --------------------------------------------------------------
message("[5/5] Exporting...")

out_table <- do.call(rbind, rows)
write.csv(out_table, file.path(output_path, "TableS28_Geometric_Coupling.csv"),
          row.names = FALSE, na = "")

message("\n============================================================")
message("  Script 04. Geometric coupling: summary")
message("============================================================")
print(out_table, row.names = FALSE)
message("\n  Exported: TableS28_Geometric_Coupling.csv")
message("  M3 and M4 add a term tied to the response by construction. They are")
message("  reported to make the circularity visible, not to test the associations.")
message("============================================================\n")


# 7. Session information -------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part04.txt"))
