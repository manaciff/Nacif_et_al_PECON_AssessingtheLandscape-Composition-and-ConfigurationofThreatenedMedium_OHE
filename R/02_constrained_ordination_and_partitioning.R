# ==============================================================================
# 02_constrained_ordination_and_partitioning.R
# Constrained ordination and variation partitioning
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   1. Screens the candidate predictors for collinearity (Pearson matrix and an
#      iterative variance inflation factor, VIF, with a threshold of 5; Zuur et
#      al. 2010).
#   2. Partitions the variation in log(envelope area) between a landscape set
#      and a spatial set of distance-based Moran eigenvector maps (MEMs). MEMs
#      are selected only if the global spatial model is significant, with the
#      double stopping criterion of Blanchet et al. (2008). The partition is
#      repeated without the predictors divided by the area of the sampling
#      window, which is the response.
#   3. Runs a redundancy analysis (RDA) of land use and cover composition, the
#      Hellinger-transformed class area of forest, farming, and herbaceous and
#      shrub vegetation (Legendre & Gallagher 2001), with marginal permutation
#      tests of each predictor.
#
# Input
#   Dados/Processados/Data_Raw_WithCoords.csv   written by script 01
#
# Output   (in Outputs/Manuscrito/PartII/)
#   FigS1_Collinearity_Pearson.tiff                     (Figure S1)
#   TableII_1_VIF_RDA.csv, TableII_1_VIF_VP.csv         (Tables S9 and S10)
#   TableII_1c_Correlation_PD_FRAC.csv
#   TableII_2_MEMs_Selected.csv                         (Table S38b)
#   TableII_2a_MEM_Global_Test.csv                      (Table S38b)
#   TableII_2b_MEM_Vectors.csv                          (read by scripts 03 and 06)
#   TableII_3_VP_Fractions.csv, TableII_3_VP_Significance.csv
#   TableII_3c_VP_Fractions_Without_Ratio_Metrics.csv   (Table S38)
#   TableII_3d_VP_Sensitivity_Summary.csv               (both partitions)
#   Figure_VP_Venn.tiff                                 (Figure S8)
#   TableII_4_RDA_Anova_Marginal.csv                    (Table 2)
#   TableII_4_RDA_Anova_Terms_Sequential.csv            (Table S39b)
#   TableII_4b_RDA_Transformation_Comparison.csv        (Table S39)
#   TableII_5_Response_Size_Dependence.csv              (Table S39)
#   TableII_4_RDA_Anova_Axes.csv, TableII_4_RDA_GlobalSummary.csv
#   TableII_4c_RDA_Biplot_Scores.csv, TableII_4d_RDA_Axis1_Forest.csv
#   TableII_0_Sample_Sizes.csv                          (Table S27)
#   Figure_RDA_Biplot.tiff                              (Figure 5)
#
# Guards
#   The script checks its input and repairs nothing. It stops if any unit is
#   incomplete over the analysis pool, if any unit lacks coordinates, or if
#   edge density is zero in a unit where the class is present.
#
# Figures
#   Titles belong in the captions, not inside the images.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)          # project-relative paths
  library(tidyverse)     # data wrangling and ggplot2
  library(vegan)         # varpart(), rda(), anova.cca(), RsquareAdj()
  library(adespatial)    # dbmem(), forward.sel(), give.thresh()
  library(car)           # vif()
  library(ggrepel)       # non-overlapping labels in the RDA biplot
})

GLOBAL_SEED <- 123
N_PERM      <- 9999

# The full sampling design. Every analysis of this script must use 67 units.
N_UNITS_EXPECTED <- 67

# Same function as in scripts 00, 01 and 03; keep the four copies identical.
# It stops on an unknown taxon name instead of returning NA.
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

# Short codes of the predictors, used in Figure S1 and in the RDA biplot (as in
# Table 2 of the manuscript). F = forest, A = farming matrix, H = herbaceous and
# shrub vegetation.
PRED_CODE <- c(
  PLAND_Forest         = "PLAND-F",
  PLAND_Agropecuaria   = "PLAND-A",
  PLAND_Herbaceous     = "PLAND-H",
  PD_Forest            = "PD-F",
  ED_Forest            = "ED-F",
  Mean_FRAC_Forest     = "FRAC_MN-F",
  Mean_ENN_Forest_m    = "ENN_MN-F",
  PROX_MN_Forest       = "PROX_MN-F",
  Mean_AREA_Forest_ha  = "AREA_MN-F",
  PD_Agropecuaria      = "PD-A",
  ED_Agropecuaria      = "ED-A",
  PROX_MN_Agropecuaria = "PROX_MN-A",
  Mean_Elevation_m     = "Mean elevation",
  Mean_Pop_Density     = "Human population density"
)
pred_code <- function(v) {
  out <- unname(PRED_CODE[v])
  ifelse(is.na(out), v, out)
}

# Okabe-Ito colors and a fixed symbol for each genus, the same as in script 01.
colors_genus <- c(
  "Puma"         = "#D55E00",
  "Alouatta"     = "#009E73",
  "Mazama"       = "#CC79A7",
  "Brachyteles"  = "#E69F00",
  "Tayassu"      = "#0072B2",
  "Myrmecophaga" = "#F0E442",
  "Leopardus"    = "#56B4E9",
  "Bradypus"     = "#999999"
)
shapes_genus <- c(Alouatta = 1, Brachyteles = 2, Bradypus = 3, Leopardus = 4,
                  Myrmecophaga = 5, Mazama = 6, Puma = 7, Tayassu = 8)

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartII")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


# 2. Data and input checks [1/9] -----------------------------------------------
message("\n[1/9] Loading Data_Raw_WithCoords.csv (written by script 01)...")

coords_file <- file.path(data_path, "Data_Raw_WithCoords.csv")
if (!file.exists(coords_file)) {
  stop("Data_Raw_WithCoords.csv not found in: ", data_path, "\n",
       "Run script 01 first to generate this file.")
}

data_raw <- read_csv(coords_file, show_col_types = FALSE)
message(sprintf("   Loaded: %d rows.", nrow(data_raw)))

data_raw <- data_raw %>%
  mutate(across(c(ORDER, GENUS, SPECIES, DIET, LOCOMOTION), as.factor),
         UA_ID = as.integer(UA_ID))

# The input is checked, not repaired: script 01 fills the metric gaps and
# recodes structural zeros. If the input does not hold all 67 complete units,
# the script stops and names the cells at fault, instead of quietly analyzing
# a subset.
data_raw$SPECIES <- canonical_taxon(data_raw$SPECIES)
data_raw$UA_ID   <- as.integer(data_raw$UA_ID)

if (nrow(data_raw) != N_UNITS_EXPECTED) {
  stop(sprintf(paste0("Data_Raw_WithCoords.csv holds %d rows, expected %d. ",
                      "Rerun script 01."),
               nrow(data_raw), N_UNITS_EXPECTED))
}

ANALYSIS_POOL <- c(
  "PLAND_Forest", "PLAND_Agropecuaria", "PLAND_Herbaceous",
  "PD_Forest", "ED_Forest", "Mean_FRAC_Forest", "Mean_ENN_Forest_m",
  "PROX_MN_Forest", "Mean_AREA_Forest_ha",
  "PD_Agropecuaria", "ED_Agropecuaria", "PROX_MN_Agropecuaria",
  "Mean_Elevation_m", "Mean_Pop_Density"
)
pool_present <- intersect(ANALYSIS_POOL, names(data_raw))
missing_pool <- setdiff(ANALYSIS_POOL, names(data_raw))
if (length(missing_pool) > 0) {
  stop("Predictor(s) absent from the input: ", paste(missing_pool, collapse = ", "))
}

incomplete <- !complete.cases(data_raw[, pool_present])
if (any(incomplete)) {
  offenders <- data_raw %>%
    filter(incomplete) %>%
    select(SPECIES, UA_ID, all_of(pool_present)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value))
  print(as.data.frame(offenders), row.names = FALSE)
  stop(sprintf(paste0("%d unit(s) are incomplete over the analysis pool and ",
                      "would be dropped by complete.cases(). Run scripts 00 ",
                      "and 01 again."),
               sum(incomplete)))
}

if (!all(c("coord_x", "coord_y") %in% names(data_raw)) ||
    any(is.na(data_raw$coord_x)) || any(is.na(data_raw$coord_y))) {
  stop("Some units lack centroid coordinates. Script 01 stops on this ",
       "condition; rerun it rather than filtering here.")
}

# A false zero is harder to see than a missing value, so edge density is
# checked explicitly against the class area.
for (cls in c("Forest", "Agropecuaria", "Herbaceous")) {
  ca_col <- paste0("CA_", cls); ed_col <- paste0("ED_", cls)
  if (all(c(ca_col, ed_col) %in% names(data_raw))) {
    suspect <- data_raw[[ca_col]] > 0 & data_raw[[ed_col]] == 0
    if (any(suspect, na.rm = TRUE)) {
      print(as.data.frame(data_raw[which(suspect), c("SPECIES", "UA_ID",
                                                     ca_col, ed_col)]),
            row.names = FALSE)
      stop(sprintf(paste0("ED_%s is zero in unit(s) where the class is present. ",
                          "Rerun script 01, which recodes a structural zero only ",
                          "where the class is absent."), cls))
    }
  }
}

message(sprintf("   Input verified: %d units, complete over %d predictors, all with coordinates.",
                nrow(data_raw), length(pool_present)))
message(sprintf("   %d observations across %d species and %d genera.",
                nrow(data_raw),
                n_distinct(data_raw$SPECIES),
                n_distinct(data_raw$GENUS)))


# 3. Pearson correlation matrix of the candidate predictors [2/9] --------------
# Figure S1 shows every pairwise correlation among the candidates before any
# variable is removed.
message("\n[2/9] Building the Pearson correlation matrix of the candidate predictors...")

all_candidates <- c(
  "PLAND_Forest", "PLAND_Agropecuaria", "PLAND_Herbaceous",
  "PD_Forest",    "ED_Forest",          "Mean_FRAC_Forest",
  "Mean_ENN_Forest_m", "PROX_MN_Forest", "Mean_AREA_Forest_ha",
  "PD_Agropecuaria", "ED_Agropecuaria",  "PROX_MN_Agropecuaria",
  "Mean_Elevation_m", "Mean_Pop_Density"
)
all_candidates <- intersect(all_candidates, names(data_raw))

cor_mat_candidates <- cor(data_raw[, all_candidates],
                          use = "complete.obs", method = "pearson")

# Lower triangle without the diagonal, in the order of all_candidates. Blue for
# negative and vermillion for positive correlations, as in the heatmaps of
# script 03; coefficients with |r| >= 0.7 are written in white.
cor_long <- as.data.frame(as.table(cor_mat_candidates)) %>%
  setNames(c("Row", "Col", "r")) %>%
  mutate(i = match(Row, all_candidates), j = match(Col, all_candidates)) %>%
  filter(i > j) %>%
  mutate(Row = factor(pred_code(as.character(Row)),
                      levels = rev(pred_code(all_candidates[-1]))),
         Col = factor(pred_code(as.character(Col)),
                      levels = pred_code(all_candidates[-length(all_candidates)])))

p_s1 <- ggplot(cor_long, aes(x = Col, y = Row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) >= 0.7), size = 2.1) +
  scale_fill_gradient2(low = "#0072B2", mid = "#F7F7F7", high = "#D55E00",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.5), name = "Pearson r") +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid  = element_blank())

ggsave("FigS1_Collinearity_Pearson.tiff", p_s1, path = output_path,
       width = 140, height = 112, units = "mm", dpi = 500,
       device = "tiff", compression = "lzw", bg = "white")
message("   FigS1_Collinearity_Pearson.tiff exported (Figure S1).")


# 4. Analytical sample size ----------------------------------------------------
# The number of units entering each analysis is written to disk in Section 10,
# so that the manuscript quotes one source.
data_sp <- data_raw
n_sp    <- nrow(data_sp)

stopifnot(n_sp == N_UNITS_EXPECTED,
          all(!is.na(data_sp$coord_x)), all(!is.na(data_sp$coord_y)))
message(sprintf("\n   %d sampling units retained for spatial analyses, all with coordinates.",
                n_sp))


# 5. Predictor pools [3/9] -----------------------------------------------------
# VP pool:  composition (PLAND) + forest configuration + matrix configuration
#           + environment. The response is log(envelope area). PLAND is divided
#           by the window area, the response; Section 9.1 repeats the partition
#           without such predictors.
# RDA pool: forest configuration + matrix configuration + environment. PLAND is
#           excluded because the RDA response is the land cover composition
#           itself.
#
# Not offered to either pool:
#   - LPI of forest and of the farming matrix (redundant with PLAND of the
#     same class)
#   - PLAND of the non-vegetated class (excluded at the kernel stage)
#   - CA and NP (they scale with the area of the sampling unit)
#   - Water (median PLAND = 0.22%)
message("\n[3/9] Defining predictor pools...")

# Composition: the three focal land cover classes
pool_composition <- c(
  "PLAND_Forest",
  "PLAND_Agropecuaria",
  "PLAND_Herbaceous"
)

# Forest configuration
pool_config_forest <- c(
  "PD_Forest",           # patch density
  "ED_Forest",           # edge density
  "Mean_FRAC_Forest",    # mean fractal dimension (shape complexity)
  "Mean_ENN_Forest_m",   # mean Euclidean nearest-neighbor distance
  "PROX_MN_Forest",      # mean proximity index
  "Mean_AREA_Forest_ha"  # mean patch area
)

# Farming matrix configuration (Agropecuaria = farming)
pool_config_matrix <- c(
  "PD_Agropecuaria",
  "ED_Agropecuaria",
  "PROX_MN_Agropecuaria"
)

# Environmental covariates
pool_environment <- c(
  "Mean_Elevation_m",
  "Mean_Pop_Density"
)

pool_VP  <- c(pool_composition, pool_config_forest, pool_config_matrix, pool_environment)
pool_RDA <- c(pool_config_forest, pool_config_matrix, pool_environment)

message(sprintf("   VP pool: %d predictors | RDA pool: %d predictors.",
                length(pool_VP), length(pool_RDA)))

missing_pred <- setdiff(unique(c(pool_VP, pool_RDA)), names(data_sp))
if (length(missing_pred) > 0) {
  stop("Predictors missing from data: ", paste(missing_pred, collapse = ", "))
}


# 6. Collinearity: iterative VIF reduction [4/9] -------------------------------
# The predictor with the highest VIF is dropped and the model is refitted until
# every remaining VIF is <= 5 (Zuur et al. 2010).
message("\n[4/9] Iterative VIF reduction...")

vif_iterative <- function(df, vars, response_expr,
                          threshold = 5, max_iter = 30) {
  retained <- vars
  for (it in seq_len(max_iter)) {
    f <- as.formula(
      paste0(response_expr, " ~ ", paste(retained, collapse = " + "))
    )
    m <- tryCatch(lm(f, data = df), error = function(e) NULL)
    if (is.null(m)) {
      stop("VIF model failed to fit. Check for perfectly collinear inputs.")
    }
    v <- tryCatch(car::vif(m), error = function(e) NULL)
    if (is.null(v) || any(!is.finite(v))) {
      cor_mat  <- cor(df[, retained], use = "complete.obs")
      drop_var <- retained[which.max(rowSums(abs(cor_mat)) - 1)]
      message(sprintf("   Iter %02d: perfect collinearity; dropping %s", it, drop_var))
    } else {
      if (max(v) <= threshold) break
      drop_var <- names(v)[which.max(v)]
      message(sprintf("   Iter %02d: dropping %s (VIF = %.2f)", it, drop_var, max(v)))
    }
    retained <- setdiff(retained, drop_var)
    if (length(retained) < 2) break
  }
  # Recompute the VIFs of the retained set, so that the exported factors
  # describe the set that was kept even when the loop stopped early.
  v_final <- tryCatch(
    car::vif(lm(as.formula(paste(response_expr, "~",
                                 paste(retained, collapse = " + "))), data = df)),
    error = function(e) NULL)
  list(retained = retained, vif_final = if (is.null(v_final)) v else v_final)
}

# VP pool
message("\n   VP pool (response: log(AOH_unit_area_ha))")
vif_VP        <- vif_iterative(data_sp, pool_VP,
                               response_expr = "log(AOH_unit_area_ha)")
predictors_VP <- vif_VP$retained

# RDA pool. The VIF depends only on the predictors, so any response serves; the
# first principal coordinate of the composition (Bray-Curtis) is used.
message("\n   RDA pool (response proxy: PCoA1 of compositional dissimilarity)")
comp_mat_vif <- data_sp %>%
  transmute(Forest     = sqrt(replace_na(CA_Forest,       0)),
            Farming    = sqrt(replace_na(CA_Agropecuaria, 0)),
            ShrubHerb  = sqrt(replace_na(CA_Herbaceous,   0))) %>%
  as.matrix()
bray_dist  <- vegan::vegdist(comp_mat_vif, method = "bray")
pcoa1      <- cmdscale(bray_dist, k = 1)[, 1]
data_sp$comp_pcoa1 <- pcoa1

vif_RDA        <- vif_iterative(data_sp, pool_RDA,
                                response_expr = "comp_pcoa1")
predictors_RDA <- vif_RDA$retained

write.csv(
  data.frame(Variable = names(vif_VP$vif_final),
             VIF      = round(as.numeric(vif_VP$vif_final), 3),
             Pool     = "Variation Partitioning"),
  file.path(output_path, "TableII_1_VIF_VP.csv"), row.names = FALSE
)
write.csv(
  data.frame(Variable = names(vif_RDA$vif_final),
             VIF      = round(as.numeric(vif_RDA$vif_final), 3),
             Pool     = "RDA constrained ordination"),
  file.path(output_path, "TableII_1_VIF_RDA.csv"), row.names = FALSE
)

message(sprintf("\n   Predictors retained, VP (n = %d): %s",
                length(predictors_VP), paste(predictors_VP, collapse = ", ")))
message(sprintf("   Predictors retained, RDA (n = %d): %s",
                length(predictors_RDA), paste(predictors_RDA, collapse = ", ")))

# PD_Forest and Mean_FRAC_Forest are both retained because their VIFs stay below
# 5, although their pairwise correlation is close to the 0.7 screening threshold
# (Dormann et al. 2013). The VIF is the multivariable diagnostic; a pairwise
# correlation only flags a pair for inspection. The pair is examined again in
# script 03 (stability check).
cor_pd_frac <- data_sp %>%
  select(PD_Forest, Mean_FRAC_Forest) %>%
  filter(complete.cases(.))

pearson_pd_frac  <- cor(cor_pd_frac$PD_Forest, cor_pd_frac$Mean_FRAC_Forest,
                        method = "pearson")
spearman_pd_frac <- cor(cor_pd_frac$PD_Forest, cor_pd_frac$Mean_FRAC_Forest,
                        method = "spearman")

message(sprintf("\n   PD_Forest x Mean_FRAC_Forest on n = %d: Pearson r = %.3f | Spearman rho = %.3f",
                nrow(cor_pd_frac), pearson_pd_frac, spearman_pd_frac))

write.csv(
  data.frame(Pair    = "PD_Forest x Mean_FRAC_Forest",
             N       = nrow(cor_pd_frac),
             Pearson = round(pearson_pd_frac, 4),
             Spearman = round(spearman_pd_frac, 4)),
  file.path(output_path, "TableII_1c_Correlation_PD_FRAC.csv"), row.names = FALSE
)


# 7. Standardization -----------------------------------------------------------
data_z          <- data_sp
vars_to_scale   <- unique(c(predictors_VP, predictors_RDA))
data_z[, vars_to_scale] <- scale(data_sp[, vars_to_scale])


# 8. Spatial eigenvectors (MEMs) [5/9] -----------------------------------------
# dbmem() builds distance-based Moran eigenvector maps from the centroids of the
# sampling units. The global model regresses log(envelope area) on all the
# eigenvectors and is tested by permutation. This screens for spatial
# autocorrelation among sampling units: geographically close units holding
# similar values of the response.
#
# It does not test the dependence among units of the same genus. That is a
# different form of non-independence, examined in script 03 by comparing models
# with and without a random intercept for genus.
#
# Only positive autocorrelation is modeled (MEM.autocor = "positive"). It is the
# signature of spatially contagious processes, the structure expected among
# envelopes that share regional landscapes (Borcard et al. 2011, Section
# 7.4.4). Script 06 screens the eigenvectors of negative autocorrelation as a
# diagnostic; they are not used here.
#
# Two-step selection (Blanchet et al. 2008): individual eigenvectors are
# selected by forward selection only if the global test is significant. The
# selection then stops when the next eigenvector is not significant at
# alpha = 0.05 or would bring the adjusted R-squared above that of the global
# model (argument adjR2thresh). If the global test is not significant, no
# eigenvector is retained and the spatial fraction is empty.
message("\n[5/9] Building MEMs with adespatial::dbmem and forward selection...")

y_log      <- log(data_z$AOH_unit_area_ha)
coords_xy  <- as.matrix(data_z[, c("coord_x", "coord_y")])

set.seed(GLOBAL_SEED)
mem_full <- dbmem(coords_xy, MEM.autocor = "positive", silent = TRUE)
message(sprintf("   %d positive MEMs generated by dbmem.", ncol(mem_full)))

mem_global     <- rda(y_log ~ ., data = as.data.frame(mem_full))
set.seed(GLOBAL_SEED)
mem_global_anova <- anova(mem_global, permutations = N_PERM)
mem_global_p     <- mem_global_anova$`Pr(>F)`[1]
mem_global_R2a   <- RsquareAdj(mem_global)$adj.r.squared

if (is.na(mem_global_p) || mem_global_p > 0.05) {
  message(sprintf(
    "   Global spatial model not significant (p = %.4f). No MEM retained.",
    mem_global_p
  ))
  mems_selected <- NULL
  mem_mat       <- matrix(numeric(0), nrow = n_sp, ncol = 0)
} else {
  message(sprintf(
    "   Global spatial model significant (p = %.4f, R2a = %.3f).",
    mem_global_p, mem_global_R2a
  ))
  # forward.sel() prints why it stopped; the message is kept for
  # TableII_2a_MEM_Global_Test.csv.
  set.seed(GLOBAL_SEED)
  fs_log <- capture.output(
    mem_sel <- forward.sel(y_log, as.matrix(mem_full),
                           alpha       = 0.05,
                           adjR2thresh = mem_global_R2a,
                           nperm       = N_PERM))
  mem_stop <- grep("Procedure stopped", fs_log, value = TRUE)
  mem_stop <- if (length(mem_stop)) trimws(mem_stop[1]) else
    "all candidate eigenvectors entered"
  message("   ", mem_stop)
  mems_selected <- mem_sel$variables
  mem_mat       <- as.matrix(mem_full[, mems_selected, drop = FALSE])
  message(sprintf("   MEMs retained: %d (%s)",
                  length(mems_selected), paste(mems_selected, collapse = ", ")))
  write.csv(mem_sel,
            file.path(output_path, "TableII_2_MEMs_Selected.csv"),
            row.names = FALSE)

  # The selected eigenvectors are exported for the spatial sensitivity check of
  # script 03: model selection on data with spatial autocorrelation should be
  # checked with the spatial structure as a covariate (Diniz-Filho, Rangel and
  # Bini 2008). The key (SPECIES, UA_ID) makes the join independent of row order.
  write.csv(
    data.frame(SPECIES = as.character(data_sp$SPECIES),
               UA_ID   = data_sp$UA_ID,
               mem_mat, check.names = FALSE),
    file.path(output_path, "TableII_2b_MEM_Vectors.csv"), row.names = FALSE)
  message("   TableII_2b_MEM_Vectors.csv exported for the spatial sensitivity check of script 03.")
}

# The global test and the outcome of the selection, for Table S38b. The
# truncation distance is the one dbmem() uses by default: the longest edge of
# the minimum spanning tree of the centroids.
mem_global_tab <- data.frame(
  Item = c("Sampling units",
           "Positive dbMEMs tested (MEM.autocor = positive)",
           "Truncation distance (m)",
           "Global model R2",
           "Global model adjusted R2",
           "Global model F",
           "Global model p",
           "Permutations",
           "MEMs retained",
           "Forward selection"),
  Value = c(n_sp,
            ncol(mem_full),
            round(give.thresh(dist(coords_xy)), 1),
            round(RsquareAdj(mem_global)$r.squared, 4),
            round(mem_global_R2a, 4),
            round(mem_global_anova$F[1], 4),
            mem_global_p,
            N_PERM,
            if (is.null(mems_selected)) "none" else paste(mems_selected, collapse = ", "),
            if (is.null(mems_selected)) "not run: global test not significant" else mem_stop)
)
write.csv(mem_global_tab,
          file.path(output_path, "TableII_2a_MEM_Global_Test.csv"),
          row.names = FALSE)
message("   TableII_2a_MEM_Global_Test.csv exported (Table S38b).")


# 9. Variation partitioning [6/9] ----------------------------------------------
# Landscape (X1) against space (X2), on adjusted R-squared. The fractions of
# varpart() are read by their row names, not by position, and two identities
# that any two-set partition must satisfy are checked before anything is
# written:
#
#     [a] + [b] = adjusted R2 of X1 alone
#     [b] + [c] = adjusted R2 of X2 alone
#
# Permutation tests use the rda(Y, X, Z) interface.
message("\n[6/9] Variation partitioning: landscape (X1) vs space (X2)...")

# Guard: rda() rejects missing values, and a unit with an NA would be dropped.
# Script 01 guarantees complete data, so this filter must remove nothing.
X1_all <- as.matrix(data_z[, predictors_VP])
vp_ok  <- complete.cases(X1_all)
n_excl <- sum(!vp_ok)
if (n_excl > 0) {
  print(as.data.frame(data_z[!vp_ok, c("SPECIES", "UA_ID")]), row.names = FALSE)
  stop(sprintf(paste0("%d sampling unit(s) carry an NA among the retained VP ",
                      "predictors and would be dropped. Run scripts 00 and 01 ",
                      "before this one."), n_excl))
}
stopifnot(sum(vp_ok) == N_UNITS_EXPECTED)
message(sprintf("   Variation partitioning fitted on %d sampling units.", sum(vp_ok)))
y_log_vp <- y_log[vp_ok]
X1       <- X1_all[vp_ok, , drop = FALSE]
X2       <- if (ncol(mem_mat) > 0) mem_mat[vp_ok, , drop = FALSE] else mem_mat

# varpart() needs at least two tables. If no MEM was retained, only the
# landscape R-squared is computed, and vp_res = NULL tells Section 9.2 to
# skip the Venn diagram.
if (ncol(X2) == 0) {

  message("   No spatial component retained: computing the landscape R-squared only.")

  set.seed(GLOBAL_SEED)
  rda_full  <- rda(y_log_vp, X1)
  test_full <- anova(rda_full, permutations = N_PERM)
  r2a_x1    <- RsquareAdj(rda_full)$adj.r.squared

  vp_table <- data.frame(
    Fraction       = c("[a] Landscape (X1)", "[d] Residuals (unexplained)"),
    Adj_R2         = round(c(r2a_x1, 1 - r2a_x1), 4),
    Pct_of_total   = round(100 * c(r2a_x1, 1 - r2a_x1), 2),
    Interpretation = c("Total landscape effect (no spatial component retained)",
                       "Variance not explained by landscape")
  )

  vp_tests <- data.frame(
    Test  = "Full landscape model (X1)",
    F_obs = round(test_full$F[1], 4),
    p_val = test_full$`Pr(>F)`[1]
  )

  vp_res <- NULL

} else {

  vp_res <- varpart(y_log_vp, X1, X2)
  ind    <- vp_res$part$indfract

  # Row names written by vegan: "[a] = X1|X2", "[b] = X2|X1", "[c]", "[d] =
  # Residuals". Row two is the pure spatial fraction and row three the shared
  # one, which is why the fractions are selected by name.
  rn         <- rownames(ind)
  row_pure1  <- grep("X1\\|X2", rn)
  row_pure2  <- grep("X2\\|X1", rn)
  row_shared <- setdiff(grep("^\\[[abc]\\]", rn), c(row_pure1, row_pure2))
  if (length(row_pure1) != 1 || length(row_pure2) != 1 || length(row_shared) != 1)
    stop("varpart() returned fraction labels this script does not recognize: ",
         paste(rn, collapse = " | "))

  frac_a      <- ind[row_pure1,  "Adj.R.squared"]   # landscape, space removed
  frac_c      <- ind[row_pure2,  "Adj.R.squared"]   # space, landscape removed
  frac_shared <- ind[row_shared, "Adj.R.squared"]   # jointly explained
  frac_d      <- 1 - (frac_a + frac_shared + frac_c)

  r2a_X1 <- RsquareAdj(rda(y_log_vp, X1))$adj.r.squared
  r2a_X2 <- RsquareAdj(rda(y_log_vp, X2))$adj.r.squared
  tol    <- 1e-6
  if (abs((frac_a + frac_shared) - r2a_X1) > tol ||
      abs((frac_shared + frac_c) - r2a_X2) > tol) {
    stop(sprintf(paste0("Variation partitioning failed its identity check.\n",
                        "  [a] + shared = %.6f, adjusted R2 of X1 alone = %.6f\n",
                        "  shared + [c] = %.6f, adjusted R2 of X2 alone = %.6f\n",
                        "Inspect the row names of varpart()$part$indfract."),
                 frac_a + frac_shared, r2a_X1, frac_shared + frac_c, r2a_X2))
  }
  message(sprintf(paste0("   Partition identity check passed: [a] + shared = ",
                         "%.4f = adjR2(X1); shared + [c] = %.4f = adjR2(X2)."),
                  r2a_X1, r2a_X2))

  vp_table <- data.frame(
    Fraction       = c("[a] Pure landscape (X1 | X2)",
                       "[b] Shared landscape and space",
                       "[c] Pure space (X2 | X1)",
                       "[d] Residuals (unexplained)"),
    Adj_R2         = round(c(frac_a, frac_shared, frac_c, frac_d), 4),
    Pct_of_total   = round(100 * c(frac_a, frac_shared, frac_c, frac_d), 2),
    Interpretation = c("Landscape effect independent of geography",
                       "Landscape variation that is spatially structured",
                       "Spatial structure without landscape correspondence",
                       "Variance not explained by either set")
  )

  # Permutation tests
  set.seed(GLOBAL_SEED)
  rda_full            <- rda(y_log_vp, cbind(X1, X2))
  test_full           <- anova(rda_full, permutations = N_PERM)

  set.seed(GLOBAL_SEED)
  rda_pure_landscape  <- rda(y_log_vp, X1, X2)   # [a]: X1 | X2
  test_pure_landscape <- anova(rda_pure_landscape, permutations = N_PERM)

  set.seed(GLOBAL_SEED)
  rda_pure_space      <- rda(y_log_vp, X2, X1)   # [c]: X2 | X1
  test_pure_space     <- anova(rda_pure_space, permutations = N_PERM)

  vp_tests <- data.frame(
    Test  = c("Full model (X1 + X2)",
              "Pure landscape (X1 | X2)",
              "Pure space (X2 | X1)"),
    F_obs = round(c(test_full$F[1],
                    test_pure_landscape$F[1],
                    test_pure_space$F[1]), 4),
    p_val = c(test_full$`Pr(>F)`[1],
              test_pure_landscape$`Pr(>F)`[1],
              test_pure_space$`Pr(>F)`[1])
  )
}

write.csv(vp_table,
          file.path(output_path, "TableII_3_VP_Fractions.csv"),
          row.names = FALSE)
message("   VP fractions:")
print(vp_table[, c("Fraction", "Adj_R2", "Pct_of_total")], row.names = FALSE)

write.csv(vp_tests,
          file.path(output_path, "TableII_3_VP_Significance.csv"),
          row.names = FALSE)

message("   Permutation test results:")
print(vp_tests, row.names = FALSE)


# 9.1 The partition without the predictors divided by the window area ----------
# Three predictors of the retained landscape set are divided by the area of the
# sampling unit, which is the response of this partition:
#
#     PD_Forest        = 100 * NP_Forest       / unit area   (patches per 100 ha)
#     PD_Agropecuaria  = 100 * NP_Agropecuaria / unit area   (patches per 100 ha)
#     PLAND_Herbaceous = 100 * CA_Herbaceous   / unit area   (percent)
#
# Each carries the response in its denominator, as script 05 shows for
# PD_Forest, so part of the variation they explain is arithmetic. The partition
# is repeated without the three, keeping only predictors whose formulas do not
# divide by the window: the shape, isolation, proximity and mean patch area of
# forest, the proximity of the farming matrix, and the two environmental
# covariates. This is the rule that defines the reported model (Table 4 of the
# manuscript, script 05).
#
# The predictors to drop are found by name, among the class-level metrics that
# landscapemetrics divides by the total landscape area (PLAND, PD, ED and LPI).
# The script stops if the set found differs from the three above, because the
# manuscript would then no longer describe this partition.
if (!is.null(vp_res)) {

  ratio_pattern  <- "^(PLAND|PD|ED|LPI)_"
  ratio_expected <- c("PD_Forest", "PD_Agropecuaria", "PLAND_Herbaceous")
  ratio_metrics  <- grep(ratio_pattern, colnames(X1), value = TRUE)
  if (!setequal(ratio_metrics, ratio_expected)) {
    found <- if (length(ratio_metrics) > 0) paste(ratio_metrics, collapse = ", ") else "none"
    stop(sprintf(paste0("Section 9.1 expected %s to be the landscape predictors ",
                        "divided by the unit area, but the set holds: %s. ",
                        "Update this block and the manuscript text together."),
                 paste(ratio_expected, collapse = ", "), found))
  }

  X1_clean <- X1[, setdiff(colnames(X1), ratio_metrics), drop = FALSE]
  vp_clean <- varpart(y_log_vp, X1_clean, X2)
  indc     <- vp_clean$part$indfract
  rnc      <- rownames(indc)

  # Fractions by row name, as in Section 9.
  row_a <- grep("X1\\|X2", rnc)
  row_c <- grep("X2\\|X1", rnc)
  row_b <- setdiff(grep("^\\[[abc]\\]", rnc), c(row_a, row_c))
  if (length(row_a) != 1 || length(row_c) != 1 || length(row_b) != 1)
    stop("varpart() returned fraction labels this script does not recognize: ",
         paste(rnc, collapse = " | "))
  ac <- indc[row_a, "Adj.R.squared"]   # landscape, space removed
  cc <- indc[row_c, "Adj.R.squared"]   # space, landscape removed
  bc <- indc[row_b, "Adj.R.squared"]   # jointly explained

  # The two identities of Section 9, on the reduced landscape set.
  r2a_X1c <- RsquareAdj(rda(y_log_vp, X1_clean))$adj.r.squared
  r2a_X2c <- RsquareAdj(rda(y_log_vp, X2))$adj.r.squared
  if (abs((ac + bc) - r2a_X1c) > 1e-6 || abs((bc + cc) - r2a_X2c) > 1e-6) {
    stop(sprintf(paste0("Section 9.1 failed its identity check.\n",
                        "  [a] + shared = %.6f, adjusted R2 of X1 alone = %.6f\n",
                        "  shared + [c] = %.6f, adjusted R2 of X2 alone = %.6f"),
                 ac + bc, r2a_X1c, bc + cc, r2a_X2c))
  }

  set.seed(GLOBAL_SEED)
  tc_full  <- anova(rda(y_log_vp, cbind(X1_clean, X2)), permutations = N_PERM)
  set.seed(GLOBAL_SEED)
  tc_land  <- anova(rda(y_log_vp, X1_clean, X2),        permutations = N_PERM)
  set.seed(GLOBAL_SEED)
  tc_space <- anova(rda(y_log_vp, X2, X1_clean),        permutations = N_PERM)

  vp_clean_tab <- data.frame(
    Fraction = c("[a] Pure landscape (X1 | X2)",
                 "[b] Shared landscape and space",
                 "[c] Pure space (X2 | X1)",
                 "[d] Residuals (unexplained)",
                 "Full model (X1 + X2)"),
    Adj_R2   = round(c(ac, bc, cc, 1 - (ac + bc + cc), ac + bc + cc), 4),
    F_obs    = c(round(tc_land$F[1], 4), NA, round(tc_space$F[1], 4), NA,
                 round(tc_full$F[1], 4)),
    p_val    = c(tc_land$`Pr(>F)`[1], NA, tc_space$`Pr(>F)`[1], NA,
                 tc_full$`Pr(>F)`[1]),
    Note     = paste0("landscape set without the predictors divided by the area ",
                      "of the sampling unit (", paste(ratio_metrics, collapse = ", "),
                      "); remaining: ", paste(colnames(X1_clean), collapse = ", "))
  )
  write.csv(vp_clean_tab,
            file.path(output_path,
                      "TableII_3c_VP_Fractions_Without_Ratio_Metrics.csv"),
            row.names = FALSE)
  message(sprintf("   Partition without the %d predictors divided by the unit area (%s):",
                  length(ratio_metrics), paste(ratio_metrics, collapse = ", ")))
  print(vp_clean_tab[, c("Fraction", "Adj_R2", "F_obs", "p_val")], row.names = FALSE)

  # The two partitions side by side, from the unrounded fractions.
  vp_compare <- data.frame(
    Fraction              = c("Landscape predictors (n)", vp_clean_tab$Fraction),
    Full_landscape_set    = c(ncol(X1),
                              round(c(frac_a, frac_shared, frac_c, frac_d,
                                      frac_a + frac_shared + frac_c), 4)),
    Without_ratio_metrics = c(ncol(X1_clean),
                              round(c(ac, bc, cc, 1 - (ac + bc + cc),
                                      ac + bc + cc), 4))
  )
  write.csv(vp_compare,
            file.path(output_path, "TableII_3d_VP_Sensitivity_Summary.csv"),
            row.names = FALSE)
}


# 9.2 Venn diagram (Figure S8) [7/9] -------------------------------------------
message("\n[7/9] Building the variation partitioning figure...")

if (is.null(vp_res)) {
  message("   No spatial eigenvector retained: the Venn diagram is not drawn.")
} else {
  tiff(filename = file.path(output_path, "Figure_VP_Venn.tiff"),
       width = 140, height = 95, units = "mm", res = 500,
       compression = "lzw", pointsize = 9)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  # Set names are drawn inside the circles with text(), so that they are never
  # clipped by the frame. The circles of plot.varpart() are centered at x = 0
  # (landscape) and x = 1 (space).
  plot(vp_res,
       digits  = 2,
       bg      = c("#0072B2", "#E69F00"),
       alpha   = 80,
       Xnames  = c("", ""),
       cex     = 1.1)
  text(x = c(0, 1), y = c(0.45, 0.45), labels = c("Landscape", "Space"),
       cex = 1.2, font = 2)
  invisible(dev.off())
  message("   Figure_VP_Venn.tiff exported (Figure S8).")
}


# 10. Redundancy analysis of land cover composition [8/9] ----------------------
# Response: Hellinger-transformed class area of the three focal classes
#           (forest, farming, and herbaceous and shrub vegetation).
# Predictors: configuration and environmental variables retained after the VIF
#             screening (PLAND excluded, because it is the composition itself).
#
# The absolute class area is the product of the class proportion and the area
# of the sampling unit, which is the response of the univariate analyses. An
# ordination of the square root of the absolute class area would therefore
# order the units by the size of the window rather than by their composition.
# The Hellinger transformation, the square root of the relative class area,
# removes this size component and is the standard preparation for a linear
# ordination of composition data (Legendre & Gallagher 2001). The square root
# of the absolute class area is also analyzed, for comparison only (Table S39).
message("\n[8/9] Constrained ordination: RDA of land cover composition...")

LULC_CA <- data_z %>%
  transmute(Forest       = replace_na(CA_Forest,       0),
            Farming      = replace_na(CA_Agropecuaria, 0),
            `Shrub/Herb` = replace_na(CA_Herbaceous,   0)) %>%
  as.matrix()

# Hellinger: sqrt(x_ij / row total).
LULC_hell <- vegan::decostand(LULC_CA, method = "hellinger")

# Square root of the absolute class area, for the comparison exported below.
LULC_sqrt <- sqrt(LULC_CA)

# Guard: every unit must enter the ordination.
X_RDA_all <- as.matrix(data_z[, predictors_RDA])
rda_ok    <- complete.cases(X_RDA_all)
if (sum(!rda_ok) > 0) {
  excluded <- data_z %>% filter(!rda_ok) %>%
    transmute(unit_label = paste(SPECIES, UA_ID)) %>% pull(unit_label)
  stop(sprintf(paste0("%d unit(s) would be excluded from the RDA for missing ",
                      "predictors: %s. Run scripts 00 and 01 again."),
               sum(!rda_ok), paste(excluded, collapse = "; ")))
}
stopifnot(sum(rda_ok) == N_UNITS_EXPECTED)
message(sprintf("   RDA fitted on %d sampling units.", sum(rda_ok)))
n_rda     <- sum(rda_ok)
X_RDA     <- X_RDA_all[rda_ok, , drop = FALSE]
LULC_hell <- LULC_hell[rda_ok, , drop = FALSE]
LULC_sqrt <- LULC_sqrt[rda_ok, , drop = FALSE]

# Every object used after the ordination carries the same rows as the
# ordination itself.
data_rda <- data_z[rda_ok, , drop = FALSE]
stopifnot(nrow(data_rda) == n_rda, nrow(LULC_hell) == n_rda,
          nrow(LULC_sqrt) == n_rda, nrow(X_RDA) == n_rda)

# How much of each candidate response is the size of the window: the share of
# inertia on the first principal component, and its correlation with the
# square root of the unit area.
pc_share <- function(Y) {
  p <- prcomp(Y, center = TRUE, scale. = FALSE)
  c(pc1_pct = 100 * p$sdev[1]^2 / sum(p$sdev^2),
    r_with_sqrt_area = abs(cor(p$x[, 1], sqrt(data_rda$AOH_unit_area_ha))))
}
size_diag <- rbind(`sqrt(class area), superseded` = pc_share(LULC_sqrt),
                   `Hellinger, reported`          = pc_share(LULC_hell))
write.csv(data.frame(Response = rownames(size_diag), round(size_diag, 4)),
          file.path(output_path, "TableII_5_Response_Size_Dependence.csv"),
          row.names = FALSE)
message(sprintf(paste0("   Size dependence of the ordination response: ",
                       "sqrt(CA) PC1 = %.1f%% of inertia, r = %.3f with ",
                       "sqrt(unit area); Hellinger PC1 = %.1f%%, r = %.3f."),
                size_diag[1, 1], size_diag[1, 2],
                size_diag[2, 1], size_diag[2, 2]))

rda_landscape    <- rda(LULC_hell ~ ., data = as.data.frame(X_RDA))

set.seed(GLOBAL_SEED)
rda_global_anova <- anova(rda_landscape, permutations = N_PERM)
# by = "margin" tests each term adjusted for all the others; this is what the
# manuscript reports (Table 2). by = "terms" is a sequential test whose result
# depends on the order of the terms in the formula; it is exported for
# description only, because its variances sum to the constrained inertia.
set.seed(GLOBAL_SEED)
rda_margin_anova <- anova(rda_landscape, by = "margin", permutations = N_PERM)
set.seed(GLOBAL_SEED)
rda_terms_anova  <- anova(rda_landscape, by = "terms", permutations = N_PERM)
set.seed(GLOBAL_SEED)
rda_axes_anova   <- anova(rda_landscape, by = "axis",  permutations = N_PERM)

rda_R2  <- RsquareAdj(rda_landscape)$r.squared
rda_R2a <- RsquareAdj(rda_landscape)$adj.r.squared

message(sprintf("   RDA: R2 = %.3f, R2-adj = %.3f, p = %.4f",
                rda_R2, rda_R2a, rda_global_anova$`Pr(>F)`[1]))

# Table 2: marginal tests, ordered by F.
rda_margin_df <- as.data.frame(rda_margin_anova) %>%
  rownames_to_column("Predictor") %>%
  arrange(desc(.data$F))
write.csv(rda_margin_df,
          file.path(output_path, "TableII_4_RDA_Anova_Marginal.csv"),
          row.names = FALSE)

rda_terms_df <- as.data.frame(rda_terms_anova) %>%
  rownames_to_column("Predictor")
write.csv(rda_terms_df,
          file.path(output_path, "TableII_4_RDA_Anova_Terms_Sequential.csv"),
          row.names = FALSE)

# The same two tests on the square root of the absolute class area, to show the
# effect of the transformation on the ranking of the predictors (Table S39).
rda_superseded <- rda(LULC_sqrt ~ ., data = as.data.frame(X_RDA))
set.seed(GLOBAL_SEED)
sup_margin <- as.data.frame(anova(rda_superseded, by = "margin",
                                  permutations = N_PERM)) %>%
  rownames_to_column("Predictor")
set.seed(GLOBAL_SEED)
sup_terms  <- as.data.frame(anova(rda_superseded, by = "terms",
                                  permutations = N_PERM)) %>%
  rownames_to_column("Predictor")
write.csv(
  bind_rows(mutate(sup_margin, Test = "marginal", Response = "sqrt(class area)"),
            mutate(sup_terms,  Test = "sequential", Response = "sqrt(class area)"),
            mutate(rda_margin_df, Test = "marginal", Response = "Hellinger"),
            mutate(rda_terms_df,  Test = "sequential", Response = "Hellinger")),
  file.path(output_path, "TableII_4b_RDA_Transformation_Comparison.csv"),
  row.names = FALSE)

rda_axes_df <- as.data.frame(rda_axes_anova) %>%
  rownames_to_column("Axis")
write.csv(rda_axes_df,
          file.path(output_path, "TableII_4_RDA_Anova_Axes.csv"),
          row.names = FALSE)

# Global summary of the RDA and of the MEM screening.
rda_global_summary <- data.frame(
  Metric = c("R2 (constrained proportion)", "Adjusted R2",
             "Constrained inertia", "Residual inertia",
             "Global permutation p",
             "MEM global model p", "MEM global model adjusted R2"),
  Value  = c(round(rda_R2, 4), round(rda_R2a, 4),
             round(sum(rda_landscape$CCA$eig), 4),
             round(sum(rda_landscape$CA$eig), 4),
             rda_global_anova$`Pr(>F)`[1],
             mem_global_p,
             round(mem_global_R2a, 4))
)
write.csv(rda_global_summary,
          file.path(output_path, "TableII_4_RDA_GlobalSummary.csv"),
          row.names = FALSE)
message(sprintf("   RDA global exported: R2 = %.4f | adj R2 = %.4f | p = %.4f",
                rda_R2, rda_R2a, rda_global_anova$`Pr(>F)`[1]))

# Biplot scores quoted in Section 3.3 of the manuscript, exported both as the
# correlations stored by vegan (CCA$biplot) and as the scaling 2 scores drawn in
# Figure_RDA_Biplot.tiff.
bp_cor    <- rda_landscape$CCA$biplot[, 1:2, drop = FALSE]
bp_plot   <- scores(rda_landscape, display = "bp", choices = 1:2, scaling = 2)
rda_bp_df <- data.frame(
  Predictor        = rownames(bp_cor),
  RDA1_correlation = round(bp_cor[, 1], 4),
  RDA2_correlation = round(bp_cor[, 2], 4),
  RDA1_plotted     = round(bp_plot[rownames(bp_cor), 1], 4),
  RDA2_plotted     = round(bp_plot[rownames(bp_cor), 2], 4),
  row.names        = NULL
)
write.csv(rda_bp_df,
          file.path(output_path, "TableII_4c_RDA_Biplot_Scores.csv"),
          row.names = FALSE)

# Correlation of the RDA1 site scores with forest cover (Section 3.3).
rda1_wa <- scores(rda_landscape, display = "sites", choices = 1, scaling = 2)[, 1]
rda1_lc <- scores(rda_landscape, display = "lc",    choices = 1, scaling = 2)[, 1]
write.csv(
  data.frame(
    Scores              = c("RDA1 site scores (weighted averages)",
                            "RDA1 linear-combination scores"),
    r_with_PLAND_Forest = round(c(cor(rda1_wa, data_rda$PLAND_Forest),
                                  cor(rda1_lc, data_rda$PLAND_Forest)), 4),
    n                   = n_rda
  ),
  file.path(output_path, "TableII_4d_RDA_Axis1_Forest.csv"),
  row.names = FALSE)
message(sprintf("   RDA1 site scores vs PLAND_Forest: r = %.4f",
                cor(rda1_wa, data_rda$PLAND_Forest)))

# Sample size used by each analysis (Table S27).
write.csv(
  data.frame(
    Analysis = c("Rows loaded", "With centroid coordinates",
                 "MEM screening", "Variation partitioning", "RDA",
                 "Correlation PD_Forest x Mean_FRAC_Forest"),
    N        = c(nrow(data_raw), n_sp, length(y_log), sum(vp_ok), n_rda,
                 nrow(cor_pd_frac))
  ),
  file.path(output_path, "TableII_0_Sample_Sizes.csv"), row.names = FALSE
)

# If any analysis ran on a different sample, the script stops here.
sizes_II <- c(nrow(data_raw), n_sp, length(y_log), sum(vp_ok), n_rda,
              nrow(cor_pd_frac))
if (any(sizes_II != N_UNITS_EXPECTED)) {
  stop("Script 02 analyses did not all use ", N_UNITS_EXPECTED,
       " units: ", paste(sizes_II, collapse = ", "))
}
message(sprintf("   All analyses of script 02 ran on %d sampling units.", N_UNITS_EXPECTED))
message(sprintf("   MEM global: p = %.4f | adj R2 = %.4f",
                mem_global_p, mem_global_R2a))

# Biplot, scaling 2 (Figure 5). Solid arrows are predictors and dashed arrows
# are the land cover classes of the response.
rda_site_scores    <- as.data.frame(scores(rda_landscape, display = "sites",   scaling = 2))
rda_species_scores <- as.data.frame(scores(rda_landscape, display = "species", scaling = 2))
rda_biplot_arrows  <- as.data.frame(scores(rda_landscape, display = "bp",      scaling = 2))
rda_site_scores$Genus <- data_rda$GENUS
stopifnot(nrow(rda_site_scores) == nrow(data_rda))

rda_biplot_arrows$Label  <- pred_code(rownames(rda_biplot_arrows))
rda_species_scores$Label <- unname(c(Forest = "Forest", Farming = "Farming",
                                     `Shrub/Herb` = "Herbaceous and shrub")[
                                       rownames(rda_species_scores)])

axis_lab <- function(i) {
  pct <- round(100 * rda_landscape$CCA$eig[i] /
                 sum(c(rda_landscape$CCA$eig, rda_landscape$CA$eig)), 1)
  sprintf("RDA%d (%.1f%%)", i, pct)
}

genera_present <- unique(as.character(rda_site_scores$Genus))
colors_subset  <- colors_genus[names(colors_genus) %in% genera_present]

# Arrow labels. Each label starts beyond the tip of its arrow, at least 0.35
# units from the origin, and is joined to the tip by a thin line when moved.
# The site points, and points spaced along every arrow, enter the same layer
# with empty labels, so that the labels are repelled from the points and from
# the arrows as well (ggrepel does not draw empty labels).
radial_nudge <- function(x, y, mult = 1.08, min_r = 0) {
  r <- sqrt(x^2 + y^2)
  k <- ifelse(r > 0, pmax(r * mult, min_r) / r, 1)
  data.frame(nudge_x = x * k - x, nudge_y = y * k - y)
}
along_arrows <- function(x, y, at = seq(0.1, 0.95, by = 0.05)) {
  data.frame(x = as.vector(outer(at, x)), y = as.vector(outer(at, y)),
             label = "", face = "plain")
}
rda_labels <- rbind(
  data.frame(x = rda_biplot_arrows$RDA1, y = rda_biplot_arrows$RDA2,
             label = rda_biplot_arrows$Label, face = "bold"),
  data.frame(x = rda_species_scores$RDA1, y = rda_species_scores$RDA2,
             label = rda_species_scores$Label, face = "plain"),
  data.frame(x = rda_site_scores$RDA1, y = rda_site_scores$RDA2,
             label = "", face = "plain"),
  along_arrows(c(rda_biplot_arrows$RDA1, rda_species_scores$RDA1),
               c(rda_biplot_arrows$RDA2, rda_species_scores$RDA2))
)
rda_nudge <- radial_nudge(rda_labels$x, rda_labels$y, min_r = 0.35)
rda_nudge[rda_labels$label == "", ] <- 0

p_rda <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(data = rda_site_scores,
             aes(x = RDA1, y = RDA2, colour = Genus, shape = Genus),
             size = 2.2, stroke = 0.8, alpha = 0.9) +
  geom_segment(data = rda_biplot_arrows,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.18, "cm")),
               colour = "grey15", linewidth = 0.6) +
  geom_segment(data = rda_species_scores,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               colour = "black", linewidth = 0.5, linetype = "dashed") +
  ggrepel::geom_text_repel(data = rda_labels,
                           aes(x = x, y = y, label = label, fontface = face),
                           nudge_x = rda_nudge$nudge_x,
                           nudge_y = rda_nudge$nudge_y,
                           size = 2.6, colour = "black",
                           segment.colour = "grey60",
                           min.segment.length = 0.2, box.padding = 0.2,
                           max.overlaps = Inf, seed = GLOBAL_SEED) +
  scale_colour_manual(values = colors_subset) +
  scale_shape_manual(values  = shapes_genus[genera_present]) +
  coord_equal() +
  labs(x = axis_lab(1), y = axis_lab(2)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "right",
        legend.text     = element_text(face = "italic"))

ggsave(filename    = "Figure_RDA_Biplot.tiff",
       plot        = p_rda,
       path        = output_path,
       width       = 190, height = 160, units = "mm",
       dpi         = 500, device = "tiff", compression = "lzw", bg = "white")
message("   Figure_RDA_Biplot.tiff exported (Figure 5).")


# 11. Output inventory [9/9] ---------------------------------------------------
message("\n[9/9] All analyses of script 02 completed.")
message("\nOutputs in: ", output_path)
message("  Figures:")
message("    FigS1_Collinearity_Pearson.tiff      (Figure S1)")
message("    Figure_VP_Venn.tiff                  (Figure S8)")
message("    Figure_RDA_Biplot.tiff               (Figure 5)")
message("  Tables:")
message("    TableII_1_VIF_VP.csv")
message("    TableII_1_VIF_RDA.csv")
message("    TableII_1c_Correlation_PD_FRAC.csv")
message("    TableII_0_Sample_Sizes.csv           (n per analysis; all = 67)")
message("    TableII_2_MEMs_Selected.csv          (if MEMs were retained)")
message("    TableII_2a_MEM_Global_Test.csv       (global test and selection)")
message("    TableII_2b_MEM_Vectors.csv")
message("    TableII_3_VP_Fractions.csv           (identity-checked)")
message("    TableII_3_VP_Significance.csv")
message("    TableII_3c_VP_Fractions_Without_Ratio_Metrics.csv")
message("    TableII_3d_VP_Sensitivity_Summary.csv (the two partitions)")
message("    TableII_4_RDA_Anova_Marginal.csv     (Table 2)")
message("    TableII_4_RDA_Anova_Terms_Sequential.csv")
message("    TableII_4b_RDA_Transformation_Comparison.csv")
message("    TableII_4_RDA_Anova_Axes.csv")
message("    TableII_4_RDA_GlobalSummary.csv")
message("    TableII_4c_RDA_Biplot_Scores.csv")
message("    TableII_4d_RDA_Axis1_Forest.csv")
message("    TableII_5_Response_Size_Dependence.csv")


# 12. Session information ------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part02.txt"))
