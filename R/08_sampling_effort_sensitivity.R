# ==============================================================================
# 08_sampling_effort_sensitivity.R
# Sampling effort: does the number of occurrence records inside an envelope
# account for the two reported associations?
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   The response is the area of a kernel utilization distribution at the 99%
#   isopleth, computed with a single smoothing parameter (h). At a fixed h, this
#   area grows with the number and the spread of the records behind it. Records
#   concentrate near roads, research stations and protected areas, where forest
#   also persists. Part of the association between forest configuration and
#   envelope extent could therefore reflect sampling effort rather than
#   landscape structure.
#
#   The script counts the retained occurrence records inside each sampling unit
#   and refits the reported model (Table 4 of the manuscript) with
#   log(records + 1), standardized, as an additional covariate.
#
# How to read the result
#   The record count is not a pure nuisance variable. Larger and more stable
#   populations are expected to maintain broader areas of use (Section 2.3 of
#   the manuscript), and they also produce more records, so the count partly
#   measures the same quantity as the response. Conditioning on it removes
#   sampling effort and part of the ecological signal. The coefficients of the
#   augmented model are therefore a lower bound on the associations, and those
#   of the reported model are an upper bound.
#
# Limits of the measure
#   Only the Rio de Janeiro subset of the occurrence compilation is archived
#   here. The kernels were estimated on the biome-wide compilation of Macedo et
#   al. (2019), so the count is the number of records inside each envelope, not
#   the complete input to the kernel. Ten units hold no record of the Rio de
#   Janeiro subset; their envelopes were produced by records from outside the
#   state. This is why log(records + 1) is used, and why the analysis is
#   repeated on the 57 units that hold at least one record.
#
# Input
#   Dados/Brutos/Coordenates_sp_Data.xlsx
#   Dados/Processados/Data_Raw_WithCoords.csv
#   08_Dados_Especies/Dados_geo_especies/Sp_data_singlepart/*.shp
#
# Output   (in Outputs/Manuscrito/PartVI/)
#   Records_Per_Unit.csv
#   TableS42_Sampling_Effort.csv               (Table S42)
#   TableS42b_Sampling_Effort_Robustness.csv   (Table S42b)
#   session_info_part08.txt
#
# Run time
#   Less than a minute. The jackknife refits the model 67 times.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)      # project-relative paths
  library(readxl)    # read_excel()
  library(dplyr)     # data manipulation
  library(sf)        # polygons, projection, point in polygon
  library(car)       # vif()
  library(DHARMa)    # simulation-based residual diagnostics
})

GLOBAL_SEED <- 123
set.seed(GLOBAL_SEED)

data_raw    <- here("Dados", "Brutos")
data_path   <- here("Dados", "Processados")
shape_path  <- here("08_Dados_Especies", "Dados_geo_especies", "Sp_data_singlepart")
output_path <- here("Outputs", "Manuscrito", "PartVI")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# The nine focal taxa as named in the occurrence file (Spp), in the analysis
# table (SPECIES) and in the polygon files (shp). The polygon file of Bradypus
# crinitus keeps the name B_torquatus, from before B. crinitus was split from
# B. torquatus.
TAXA <- tibble::tribble(
  ~Spp,                       ~SPECIES,          ~shp,
  "Alouatta guariba",         "A_guariba",       "A_guariba_utm",
  "Brachyteles arachnoides",  "B_arachnoides",   "B_arachnoides_utm",
  "Brachyteles hypoxanthus",  "B_hypoxanthus",   "B_hypoxantus_utm",
  "Bradypus crinitus",        "B_crinitus",      "B_torquatus_utm",
  "Leopardus wiedii",         "L_wiedii",        "L_wiedii_utm",
  "Mazama americana",         "Mazama",          "Mazama_utm",
  "Myrmecophaga tridactyla",  "M_tridactyla",    "M_tridactyla_utm",
  "Puma concolor",            "P_concolor",      "P_concolor_utm",
  "Tayassu pecari",           "T_pecari",        "T_pecari_utm"
)


# 2. Occurrence records [1/6] --------------------------------------------------
# Filtered as in Section 2.3.1 of the manuscript: collection date from 1990 to
# 2020, valid coordinates, and no duplication of taxon, coordinate pair and year.
message("\n[1/6] Occurrence records...")

occ_file <- file.path(data_raw, "Coordenates_sp_Data.xlsx")
if (!file.exists(occ_file)) {
  stop("Coordenates_sp_Data.xlsx not found in ", data_raw, call. = FALSE)
}

occ <- readxl::read_excel(occ_file) %>%
  filter(Spp %in% TAXA$Spp) %>%
  filter(!is.na(Ano), Ano >= 1990, Ano <= 2020) %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  left_join(TAXA, by = "Spp") %>%
  mutate(Latitude = round(as.numeric(Latitude), 6),
         Longitude = round(as.numeric(Longitude), 6)) %>%
  distinct(SPECIES, Latitude, Longitude, Ano, .keep_all = TRUE)

message(sprintf("   %d records retained across %d taxa",
                nrow(occ), dplyr::n_distinct(occ$SPECIES)))


# 3. Sampling-unit polygons and the match to UA_ID [2/6] -----------------------
# The polygon files carry no UA_ID. Each polygon is matched to a row of the
# analysis table by its planar area, and the match within each taxon must be
# one-to-one. The script stops if any match exceeds the area tolerance, rather
# than proceeding on a mismatched join. With the full-precision areas written
# by script 01, the relative error of every match is below 1e-5
# (column area_rel_err of Records_Per_Unit.csv).
message("[2/6] Sampling-unit polygons and the match to UA_ID...")

ref <- read.csv(file.path(data_path, "Data_Raw_WithCoords.csv"),
                stringsAsFactors = FALSE)
if (!"AOH_unit_area_ha" %in% names(ref)) {
  stop("Column AOH_unit_area_ha not found in Data_Raw_WithCoords.csv. ",
       "Run scripts 00 and 01 before this one.", call. = FALSE)
}
stopifnot(nrow(ref) == 67)

MAX_AREA_REL_ERR <- 0.02   # tolerance on the relative area error of a match

counts <- list()
for (i in seq_len(nrow(TAXA))) {
  tx  <- TAXA[i, ]
  shp <- sf::st_read(file.path(shape_path, paste0(tx$shp, ".shp")), quiet = TRUE)
  shp <- sf::st_zm(shp, drop = TRUE)          # the files are PolygonZ

  cand <- ref[ref$SPECIES == tx$SPECIES, ]
  if (nrow(cand) != nrow(shp)) {
    stop(sprintf("%s: %d polygons but %d rows in the analysis table.",
                 tx$SPECIES, nrow(shp), nrow(cand)), call. = FALSE)
  }

  area_shp <- as.numeric(sf::st_area(shp)) / 1e4      # hectares
  # Greedy assignment on relative area error, smallest error first.
  free <- rep(TRUE, nrow(cand)); assign_id <- integer(nrow(shp))
  err  <- outer(area_shp, cand$AOH_unit_area_ha,
                function(a, b) abs(a - b) / pmax(b, 1))
  for (k in order(apply(err, 1, min))) {
    j <- which(free)[which.min(err[k, free])]
    assign_id[k] <- j; free[j] <- FALSE
  }
  rel_err <- err[cbind(seq_len(nrow(shp)), assign_id)]
  if (max(rel_err) > MAX_AREA_REL_ERR) {
    stop(sprintf("%s: the polygon-to-UA_ID match is not safe. ",
                 tx$SPECIES),
         sprintf("Largest relative area error %.4f exceeds %.2f.",
                 max(rel_err), MAX_AREA_REL_ERR), call. = FALSE)
  }
  stopifnot(!anyDuplicated(assign_id))

  # Records of this taxon that fall inside each of its polygons.
  pts <- occ %>% filter(SPECIES == tx$SPECIES)
  n_in <- rep(0L, nrow(shp))
  if (nrow(pts)) {
    sp <- sf::st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326)
    sp <- sf::st_transform(sp, sf::st_crs(shp))
    n_in <- lengths(sf::st_intersects(shp, sp))
  }

  counts[[tx$SPECIES]] <- data.frame(
    SPECIES      = tx$SPECIES,
    UA_ID        = cand$UA_ID[assign_id],
    GENUS        = cand$GENUS[assign_id],
    A_ha_table   = cand$AOH_unit_area_ha[assign_id],
    A_ha_polygon = round(area_shp, 1),
    area_rel_err = round(rel_err, 5),
    n_records    = as.integer(n_in),
    stringsAsFactors = FALSE
  )
  message(sprintf("   %-16s %2d units | %3d records inside | max area error %.4f",
                  tx$SPECIES, nrow(shp), sum(n_in), max(rel_err)))
}

rec <- do.call(rbind, counts)
rec <- rec[order(rec$SPECIES, rec$UA_ID), ]
stopifnot(nrow(rec) == 67)
write.csv(rec, file.path(output_path, "Records_Per_Unit.csv"), row.names = FALSE)
message(sprintf("   %d of %d records fall inside the 67 envelopes; %d units hold none.",
                sum(rec$n_records), nrow(occ), sum(rec$n_records == 0)))


# 4. The reported model with and without sampling effort [3/6] -----------------
message("[3/6] Fitting the two models...")

z <- function(x) as.numeric(scale(x))

d <- ref %>%
  left_join(rec[, c("SPECIES", "UA_ID", "n_records")], by = c("SPECIES", "UA_ID")) %>%
  mutate(AOH_ha  = AOH_unit_area_ha,
         lognrec = log(n_records + 1),
         ENN_z   = z(Mean_ENN_Forest_m),
         FRAC_z  = z(Mean_FRAC_Forest),
         PROX_z  = z(PROX_MN_Forest),
         NREC_z  = z(lognrec))
stopifnot(!any(is.na(d$n_records)))

LAB <- c("(Intercept)" = "Intercept",
         ENN_z  = "Isolation (ENN_MN)",
         FRAC_z = "Shape complexity (FRAC_MN)",
         PROX_z = "Forest proximity (PROX_MN)",
         NREC_z = "Sampling effort, log(records + 1)")

f_reported <- AOH_ha ~ ENN_z + FRAC_z + PROX_z
f_effort   <- AOH_ha ~ ENN_z + FRAC_z + PROX_z + NREC_z

# Gamma GLM with a log link, as in Table 4. The VIF is computed on the linear
# model of log(OHE) with the same predictors.
tidy_glm <- function(f, label) {
  m  <- glm(f, data = d, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- confint.default(m)
  vf <- car::vif(lm(update(f, log(AOH_ha) ~ .), data = d))
  data.frame(
    Model              = label,
    Predictor          = unname(LAB[rownames(cf)]),
    Estimate           = round(cf[, 1], 4),
    SE                 = round(cf[, 2], 4),
    CI_lower           = round(ci[, 1], 4),
    CI_upper           = round(ci[, 2], 4),
    P_value            = signif(cf[, 4], 4),
    OHE_ratio_per_SD   = ifelse(rownames(cf) == "(Intercept)", NA,
                                round(exp(cf[, 1]), 3)),
    VIF                = ifelse(rownames(cf) == "(Intercept)", NA,
                                round(vf[match(rownames(cf), names(vf))], 3)),
    Deviance_explained = ifelse(rownames(cf) == "(Intercept)",
                                round(1 - m$deviance / m$null.deviance, 4), NA),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

M1_LABEL <- "M1. Reported model (Table 4)"
M2_LABEL <- "M2. Reported model plus sampling effort"
tab <- rbind(tidy_glm(f_reported, M1_LABEL),
             tidy_glm(f_effort,   M2_LABEL))
write.csv(tab, file.path(output_path, "TableS42_Sampling_Effort.csv"),
          row.names = FALSE, na = "")
print(tab, row.names = FALSE)


# 5. Residual diagnostics of the augmented model [4/6] -------------------------
message("[4/6] Residual diagnostics (DHARMa)...")
set.seed(GLOBAL_SEED)
m_eff <- glm(f_effort, data = d, family = Gamma(link = "log"))
sim   <- DHARMa::simulateResiduals(m_eff, n = 1000)
disp  <- DHARMa::testDispersion(sim, plot = FALSE)
unif  <- DHARMa::testUniformity(sim, plot = FALSE)
outl  <- DHARMa::testOutliers(sim, plot = FALSE)
message(sprintf("   dispersion p = %.4f | uniformity p = %.4f | outliers p = %.4f",
                disp$p.value, unif$p.value, outl$p.value))


# 6. Robustness [5/6] ----------------------------------------------------------
message("[5/6] Robustness: jackknife and two subsets...")

rb <- function(q, v, n = "") data.frame(Quantity = q, Value = as.character(v),
                                        Note = n, stringsAsFactors = FALSE)
lo <- log(d$AOH_ha)
rows <- list(
  rb("Records inside the 67 envelopes, of those retained", sum(rec$n_records),
     sprintf("%d records retained for 1990 to 2020, deduplicated by taxon, coordinate and year", nrow(occ))),
  rb("Sampling units with zero records inside", sum(rec$n_records == 0),
     "the kernels were estimated on the biome-wide compilation; these units hold no record of the Rio de Janeiro subset"),
  rb("Median records per unit", median(rec$n_records),
     sprintf("range %d to %d", min(rec$n_records), max(rec$n_records))),
  rb("Pearson r, log(records + 1) with log(OHE)", round(cor(d$lognrec, lo), 4),
     "compare with the two rows marked 'for comparison'"),
  rb("Spearman rho, records with OHE",
     round(cor(d$n_records, d$AOH_ha, method = "spearman"), 4),
     "distribution-free cross-check"),
  rb("Pearson r, log(records + 1) with PROX_MN", round(cor(d$lognrec, d$PROX_MN_Forest), 4), ""),
  rb("Pearson r, log(records + 1) with FRAC_MN", round(cor(d$lognrec, d$Mean_FRAC_Forest), 4), ""),
  rb("Pearson r, PROX_MN with log(OHE)", round(cor(d$PROX_MN_Forest, lo), 4), "for comparison"),
  rb("Pearson r, PLAND_Forest with log(OHE)", round(cor(d$PLAND_Forest, lo), 4), "for comparison"),
  rb("DHARMa dispersion, uniformity and outliers of the augmented model",
     sprintf("p = %.4f, %.4f, %.4f", disp$p.value, unif$p.value, outl$p.value),
     "1000 simulated residual sets; p > 0.05 indicates no detectable departure")
)

# Jackknife: drop one unit at a time, standardize again and refit.
jk <- matrix(NA_real_, nrow = nrow(d), ncol = 4,
             dimnames = list(NULL, c("ENN_z", "FRAC_z", "PROX_z", "NREC_z")))
for (i in seq_len(nrow(d))) {
  dd <- d[-i, ]
  dd$ENN_z  <- z(dd$Mean_ENN_Forest_m); dd$FRAC_z <- z(dd$Mean_FRAC_Forest)
  dd$PROX_z <- z(dd$PROX_MN_Forest);    dd$NREC_z <- z(dd$lognrec)
  jk[i, ] <- coef(glm(f_effort, data = dd, family = Gamma(link = "log")))[-1]
}
for (v in colnames(jk)) {
  rows[[length(rows) + 1]] <- rb(
    sprintf("Jackknife of %s, one unit removed at a time", LAB[[v]]),
    sprintf("%.4f to %.4f", min(jk[, v]), max(jk[, v])),
    sprintf("median %.4f", median(jk[, v])))
}

# Refit the augmented model on a subset of units, standardizing again.
refit_subset <- function(keep, label) {
  dd <- d[keep, ]
  dd$ENN_z  <- z(dd$Mean_ENN_Forest_m); dd$FRAC_z <- z(dd$Mean_FRAC_Forest)
  dd$PROX_z <- z(dd$PROX_MN_Forest);    dd$NREC_z <- z(dd$lognrec)
  m <- glm(f_effort, data = dd, family = Gamma(link = "log"))
  ci <- confint.default(m)
  for (v in c("PROX_z", "FRAC_z")) {
    rows[[length(rows) + 1]] <<- rb(
      sprintf("%s %s", LAB[[v]], label),
      sprintf("%.4f [%.4f, %.4f]", coef(m)[v], ci[v, 1], ci[v, 2]),
      sprintf("n = %d", nrow(dd)))
  }
}
refit_subset(d$n_records >= 1, "with units of zero records removed")
refit_subset(rank(-d$n_records, ties.method = "first") > 4,
             "with the four best-sampled units removed")

rb_tab <- do.call(rbind, rows)
write.csv(rb_tab, file.path(output_path, "TableS42b_Sampling_Effort_Robustness.csv"),
          row.names = FALSE, na = "")
print(rb_tab, row.names = FALSE)


# 7. Summary [6/6] -------------------------------------------------------------
message("[6/6] Done.")
message("\n============================================================")
message("  Script 08. Sampling effort: summary")
message("============================================================")
message("  Coefficient without and with the record count as a covariate,")
message("  and the 95% CI of the second:")
for (v in c("PROX_z", "FRAC_z")) {
  a <- tab[tab$Model == M1_LABEL & tab$Predictor == LAB[[v]], ]
  b <- tab[tab$Model == M2_LABEL & tab$Predictor == LAB[[v]], ]
  message(sprintf("  %-28s %.3f -> %.3f  [%.3f, %.3f]",
                  LAB[[v]], a$Estimate, b$Estimate, b$CI_lower, b$CI_upper))
}
message("  Read the two estimates as an upper and a lower bound, not as an")
message("  estimate and a correction: the record count partly measures the same")
message("  quantity as the response.")
message("============================================================\n")


# 8. Session information -------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part08.txt"))
