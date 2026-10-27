# What each script does, and why

Manuscript: *Assessing the Landscape Composition and Configuration of Threatened Medium- and Large-Sized
Mammals using Habitat Envelopes as Sampling Windows*

One section per script. Each states the question the script answers, the decisions it encodes, the
checks it enforces and what it writes. `TABLE_MAP.md` maps every output file to the table or figure of
the manuscript that reports it.

---

## `00_recompute_class_metrics.R`

**Question.** Do the class-level metrics stored in the consolidated dataset reproduce from the per-unit
rasters, and can the empty cells of two sampling units be filled without changing the method?

**Why it exists.** Two units, *A. guariba* UA 6 and *B. arachnoides* UA 1, reached the consolidated
dataset with empty cells in edge density and largest patch index of every class, and in the number of
patches, patch density and proximity index of the farming class. Without this script, later steps would
drop both units.

**What it does.** Recomputes every class-level metric of the 67 units with `landscapemetrics`, with the
same reclassification (pasture and agriculture merged into farming before the metrics are computed) and
the same eight-neighbor rule used for the stored values. The proximity index is recovered from the
FRAGSTATS patch tables. The script then compares, metric by metric, the recomputed and the stored values.
Nothing is imputed: every filled value is measured again from the same raster.

**Run time.** More than 99% of each raster is empty padding, so every raster is cropped with
`terra::trim()` before the metrics are computed. The result of each unit is cached, so an interrupted run
resumes where it stopped.

**Checks.** Stops unless 67 rasters are found, and on a duplicated (taxon, unit) key.

**Writes.** `Dados/Processados/Tabela_Metricas_Recomputada.csv`, `Part0/TableS_Recompute_Agreement.csv`
(Table S32) and `Part0/TableS_Repaired_Cells.csv` (Table S33).

---

## `01_prepare_data_and_ordination.R`

**Question.** What is the composition of the landscape inside each occupied habitat envelope, and does it
differ among taxonomic and functional groups?

**What it does, in order.**

1. **Canonical taxon labels.** `canonical_taxon()` maps every spelling variant found in the project onto
   one label, and an unknown name stops the run. The variant that mattered was *Brachyteles hypoxanthus*,
   spelled `B_hypoxantus` in the table and correctly in the shapefile.
2. **Response at full precision.** If `AOH_area_full_precision.csv` exists (written by
   `arcpy_AOH_area_full_precision.py`), it replaces the rounded areas of `Data_Raw_FINAL.csv`.
3. **Metric repair.** Fills the empty cells of two units from the table of script 00.
4. **Conditional structural zeros.** FRAGSTATS omits the row of a class that is absent from a unit, and
   that omission became `NA`. A metric is set to zero only when the class area of its class is zero; an
   `NA` in a class that is present is a missing value, and the run stops.
5. **Centroids** of the sampling units, from the per-taxon shapefiles. Invalid geometries are repaired,
   and the point on the surface replaces a centroid that falls outside its polygon.
6. **Ordination and tests.** NMDS on Bray-Curtis dissimilarities of the square root of the class area of
   the three focal classes; `envfit` of forest cover and six configuration metrics of forest;
   `betadisper` with a permutation test before PERMANOVA, because the design is unbalanced and PERMANOVA
   responds to differences in dispersion as well as in location.

**Figures.** Figure 2 (composition) and Figure 3 (configuration) show one density ridge per genus. Each
density is evaluated only within the possible range of its metric (for example, 0 to 100% for the
percentage of landscape), so that no ridge extends into impossible values.

**Checks.** Stops unless the input has 67 rows, all 67 carry coordinates, and all 67 are complete over
the fourteen-predictor analysis pool used by scripts 02 and 03.

**Writes.** `Data_Raw_WithCoords.csv`, the input of scripts 02 to 08; Table 1; Figures 2 to 4; Tables S3
to S8 and S11 to S14; `TableS_NMDS_Summary.csv`, `TableS_Metric_Repair_Log.csv` and
`TableS_Sample_Sizes_PartI.csv`.

---

## `02_constrained_ordination_and_partitioning.R`

**Question.** Which landscape and environmental predictors structure the composition, and how much of the
variation in envelope extent is explained by the landscape once spatial structure is accounted for?

**What it does.**

- **Collinearity.** Figure S1 shows the Pearson correlations among the fourteen candidate predictors. An
  iterative variance inflation factor (VIF) reduction with a threshold of 5 (Zuur et al., 2010) then
  removes forest and farming edge density from the ordination pool, and the same two plus the percentage
  of forest and of farming from the partitioning pool (Tables S9 and S10).
- **Variation partitioning** of log(envelope extent) between the landscape set and a spatial set of
  distance-based Moran eigenvector maps (MEMs) of the unit centroids, modeling positive autocorrelation
  only. Individual eigenvectors are selected only if the global spatial model is significant, and the
  forward selection stops at the first eigenvector that is not significant or that would bring the
  adjusted R-squared above that of the global model (Blanchet et al., 2008). The global test and the
  reason the selection stopped are written to `TableII_2a_MEM_Global_Test.csv`. The fractions are read by
  the row names that `vegan` writes, and two arithmetic identities are checked before anything is
  exported.
- **Partition without the metrics normalized by the window.** Three predictors of the landscape set,
  PD_Forest, PD_Agropecuaria and PLAND_Herbaceous, are divided by the area of the sampling unit, which is
  the response of the partition. Section 9.1 repeats the partition without them (Table S38) and writes
  both partitions side by side to `TableII_3d_VP_Sensitivity_Summary.csv`. The script stops if the
  predictors divided by the unit area are not those three.
- **Redundancy analysis (RDA)** of the Hellinger-transformed class area of the three focal classes,
  constrained by configuration and environmental predictors. The percentage of landscape is excluded,
  because it is the composition itself. Term tests are marginal (Table 2); the sequential tests (Table
  S39b) and both tests on the square root of the absolute class area (Table S39) are exported for
  comparison.

**Checks.** The input is checked, not repaired. The script stops if any unit is incomplete over the
analysis pool, if any unit lacks coordinates, or if edge density is zero in a unit where the class is
present.

**Writes.** Tables 2, S9, S10, S38, S39 and S39b, and part of Table S38b; Figures 5, S1 and S8; the Moran
eigenvectors read by scripts 03 and 06 (`TableII_2b_MEM_Vectors.csv`).

---

## `03_univariate_model_selection.R`

**Question.** Which forest configuration metrics accompany a larger envelope, under the candidate set
defined a priori?

**What it does.**

- **Six specifications** on a common response scale: ordinary least squares, a linear mixed model with a
  random intercept for genus, and a Gamma generalized linear model with a log link, each with and without
  the a priori interaction. They are compared by AICc, by BIC, and by exact leave-one-out
  cross-validation, with leave-one-genus-out beside it (Brewer et al., 2016; Yates et al., 2023). The
  retained specification is the Gamma model with the interaction.
- **The full a priori model is fitted once, without selection** (Table 3), because inference on
  coefficients is biased when it follows a selection step (Tredennick et al., 2021; Yates et al., 2023).
- **Spatial sensitivity.** The same model is refitted with the Moran eigenvectors of script 02 (Table
  S34; Diniz-Filho et al., 2008).
- **Selection uncertainty.** `MuMIn::dredge()` ranks 40 candidate subsets: the product of forest proximity
  and forest cover enters only models that contain both, and no other rule is added, so proximity and cover
  each enter 24 of the 40 models. The conditional
  and full averages of the competitive models and the sums of Akaike weights are reported as a
  description of selection uncertainty, not as a second set of estimates (Tables S18b to S20).
- **Genus-level slopes, by two designs.** Section 15 fits a univariate regression inside each genus
  (Figure 7, Table S24). Section 15.1 fits one Gamma model per metric to all 67 units with genus as a
  moderator, so that small genera borrow precision from the rest (Table S29). With four to twelve units
  per genus, neither is confirmatory.

**The a priori interaction.** The fragmentation threshold hypothesis is a statement about an interaction
between connectivity and habitat amount, so the product of the proximity index and forest cover is the
only interaction admitted.

**Important.** Three terms of the a priori set (PD, PLAND and the product) are divided by the window
area, which is the response, or built on a term that is. Script 05 shows that their coefficients cannot
be interpreted, and the manuscript reports the model of script 05 (Table 4). Figure S10 separates the
three terms that can be read from the three that cannot.

---

## `04_geometric_coupling_diagnostic.R`

**Question.** Which predictors are algebraically tied to the response, and by how much?

**What it does.** Checks, on the 67 units, the identity that ties patch density to the response:
log(OHE) = log(NP) + log(100) - log(PD). It then fits four linear mixed models: a reference model, the
same model with forest cover (PLAND) as a control for habitat amount, and two models conditioned on the
number of patches or on the forest class area, which are tied to the response by construction and are
reported to make that circularity visible (Table S26).

---

## `05_ratio_artifact_test.R`

**Question.** Are the coefficients of forest cover and patch density, and the interaction built on forest
cover, more than the ratio construction produces on its own?

**This script selects the model the manuscript reports (Table 4, Figure 6).**

**The permutation null.** The forest class area and the number of forest patches are shuffled together
across sampling units, which keeps their mutual association, while the area of each unit stays at its
observed value. Forest cover and patch density are recomputed from the shuffled numerators and the
observed denominators, and the model is refitted 9,999 times. Under this null the amount of forest is
unrelated to the extent of the envelope, yet the ratio construction remains, so the coefficients it
produces measure the construction alone (Tables S28b and S28d, Figure S9).

**The same null is applied to the information criterion** (Table S28c), with 2,000 permutations,
because each one refits two models. A predictor built from the response improves the fit to the
response whether or not it carries ecological information, so the AICc advantage of the coupled model
cannot justify keeping it.

**Six specifications** remove the normalized metrics singly and jointly, or measure habitat amount by the
absolute class area or by the mean forest patch area (Table S28a). A separate comparison asks how much of
the advantage of the interaction remains when the amount term is not divided by the window (Table S28f).

**Further checks.** Whether the proximity index, which is built from the areas of neighboring patches,
can be separated from the amount of forest (Tables S30 and S31); the farming matrix without
window-normalized predictors (Table S40); a thinning test for overlapping units of different taxa (Table
S41); and out-of-sample error and BIC (Table S35).

---

## `06_spatial_autocorrelation_checks.R`

**Question.** Does the choice to model positive spatial autocorrelation only leave negative
autocorrelation, or any residual autocorrelation, unaccounted for?

**What it does.**

- **Negative autocorrelation.** Tests the distance-based eigenvectors of negative eigenvalue together, by
  permutation, in the same way script 02 tests the positive ones. The test is a diagnostic and selects no
  eigenvector. Beside each p-value it writes the Sidak-corrected value that a joint screen of the two signs
  would use (Borcard et al., 2011).
- **Moran's I, two-sided,** of log(envelope extent) and of the residuals of the partition models and of the
  reported model, with the spatial weights of the eigenvectors. Least-squares residuals are tested with the
  exact moments of `spdep::lm.morantest`; the Gamma model is tested on the simulated quantile residuals of
  `DHARMa`, which are independent and uniform under a correct model.
- **The reported model with the eigenvectors.** Refits the model of Table 4 with the eigenvectors retained
  by script 02, as script 03 does for the a priori model (Table S34). Table 4 itself is unchanged.
- **Residual tests of the reported model.** Tests the same simulated quantile residuals, without and with
  the eigenvectors, for dispersion, uniformity and outliers, as script 03 does for the a priori model
  (Table S17).

**Checks.** The script stops if the eigenvectors written by script 02 do not match the ones it rebuilds,
or if its fit of the reported model does not reproduce `PartV/TableS_PartV_Uncoupled_Model.csv`.

**Writes.** `PartII/TableII_6_MEM_Global_Tests.csv`, `PartII/TableII_6b_Reported_Model_With_MEMs.csv`,
`PartII/TableII_6c_Moran_I_Residuals.csv` and `PartII/TableII_6d_DHARMa_Reported_Model.csv` (Tables S34b
and S38b).

---

## `08_sampling_effort_sensitivity.R`

**Question.** Does the number of occurrence records inside an envelope account for the two reported
associations?

At a fixed bandwidth, the area of a kernel utilization distribution grows with the number and spread of
the records behind it, and records concentrate near roads, research stations and protected areas, where
forest also persists. The script counts the retained records inside each sampling unit and refits the
reported model (Table 4) with log(records + 1) as an additional covariate.

**How to read the result.** The record count is not a pure nuisance variable: a larger and more stable
population uses a broader area and also produces more records. Conditioning on the count removes
sampling effort and part of the ecological signal, so the augmented coefficients are a lower bound and
those of Table 4 an upper bound.

**Writes.** `PartVI/TableS42_Sampling_Effort.csv`, `PartVI/TableS42b_Sampling_Effort_Robustness.csv` and
`PartVI/Records_Per_Unit.csv`.

---

## `run_all.R` and `arcpy_AOH_area_full_precision.py`

`run_all.R` runs scripts 02 to 08 in order (and scripts 00 and 01 when `RUN_SCRIPT_00` and
`RUN_SCRIPT_01` are set to `TRUE`), counts the files each script writes and keeps a log in
`Outputs/Manuscrito/run_all_log.txt`.

`arcpy_AOH_area_full_precision.py` runs in ArcGIS Pro and exports the area of every sampling unit at full
precision; script 01 reads the file it writes.

---

## What has no script here

The study area map (Figure 1) and the weighted conservation value index (Figure 8, Tables S36 and S37)
were built in ArcGIS Pro 3.1. Section 5 of the supplementary material states the steps: the nine envelope
layers are intersected in the South America Albers Equal Area Conic projection, the sensitivity weights
are summed over the resulting polygons, and the result is intersected with the protected area network.
The input polygons are archived with the repository, so the index can be rebuilt.

---

## Checks common to the chain

Every script stops rather than analyzing a reduced sample. The checks are:

- the input has 67 rows;
- all 67 units carry centroid coordinates;
- all 67 units are complete over the analysis pool;
- no metric was set to zero in a class that is present in the unit;
- the sample size of every analysis is 67, checked and written to disk.
