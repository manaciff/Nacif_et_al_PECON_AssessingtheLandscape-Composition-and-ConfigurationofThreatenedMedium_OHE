# Assessing the Landscape Composition and Configuration of Threatened Medium- and Large-Sized Mammals using Habitat Envelopes as Sampling Windows

Code and data for the analysis of landscape composition and configuration inside the occupied habitat
envelopes of nine threatened medium- and large-sized mammal taxa in the Atlantic Forest of Rio de
Janeiro State, Brazil.

**Author:** Maria Eduarda Nacif
**Status:** manuscript in preparation. Last revision of the analysis chain: 15 October 2026.
**R version:** 4.5.2

---

## What the study asks

**Methodological question.** Landscape studies often measure habitat configuration inside a window whose
area is also the response variable: home ranges, utilization distributions, territories and occupied
polygons are all used this way. Any metric normalized by the area of that window then carries the
response in its own denominator. This repository measures the consequence and shows that an information
criterion cannot arbitrate between a model that contains such a metric and one that does not.

**Empirical question.** Which features of the remaining forest accompany a larger occupied habitat
envelope, once the metrics exposed to that dependence are set aside?

**The goal is inference, not prediction.** Exploration, inference and prediction require different model
selection procedures (Tredennick et al., 2021). The estimates reported here come from single models
specified a priori and fitted once. The exhaustive search over subsets is reported as a description of
selection uncertainty, not as a second set of estimates, because inference on coefficients is biased
when it follows a selection step (Yates et al., 2023).

## The response variable

The response is the **area of the occupied habitat envelope (OHE)**, in hectares: the 99% isopleth of a
kernel utilization distribution fitted to occurrence records, masked, clipped at the boundary of Rio de
Janeiro State and split into spatially separate parts. L. Macedo estimated the kernels with the protocol
of Macedo et al. (2019), who report a bandwidth of 0.02684 for their kernel maps; the repository starts
from the resulting polygons (`08_Dados_Especies/`) and does not reproduce the kernel step. Ten units
contain no retained record inside the state, because their kernels come from records outside it.

It is **not** the Area of Habitat (AOH) of Brooks et al. (2019), the suitable habitat within a species
range; in these data the forest inside the envelope is the closer analog of that quantity, and it
covers a median of 54.9% of the envelope. It is **not** the Area of Occupancy (AOO) of the IUCN Red List,
which is measured on a 2 x 2 km grid. No IUCN threshold is applied to any value here. Some column names in
the code keep the earlier label AOH (for example, `AOH_ha`); they always mean the OHE.

## Sampling design

67 sampling units from 9 taxa in 8 genera. Every analysis uses all 67 units; each script stops, and
names the cells at fault, if a step would reduce that number.

---

## Repository layout

```
R/                                 analysis scripts, numbered in execution order
run_all.R                          runs the chain in order (see "How to run")
Dados/Processados/                 tabular inputs, and the files the scripts write back
Dados/FRAGSTATS_RESULT/            FRAGSTATS patch and class tables, read by script 00
Dados/Rasters/UA_RASTER/           67 per-unit rasters, read by script 00 (not tracked; see below)
Dados/Shapefiles/                  hexagonal grid and the conservation value index layers
08_Dados_Especies/                 per-taxon envelope polygons, read by scripts 01 and 08
Outputs/                           archived tables and figures of the run of 16 October 2026, one subfolder
                                   per part; the scripts write new copies to Outputs/Manuscrito/
arcpy_AOH_area_full_precision.py   exports the envelope areas at full precision from ArcGIS Pro
```

Paths are resolved with `here::here()`. The `.here` file marks the root of the repository, so a clone
runs on its own.

Two inputs are not tracked, because of their size or origin:

- `Dados/Rasters/UA_RASTER/`, the 67 per-unit rasters (2.6 GB), read only by script 00.
  `Dados/Rasters/UA_RASTER/READ_THIS_FIRST.md` explains how to obtain them. Script 00 writes
  `Dados/Processados/Tabela_Metricas_Recomputada.csv`, which is tracked, so the other scripts run
  without the rasters.
- `Dados/Brutos/Coordenates_sp_Data.xlsx`, the occurrence records compiled by Macedo et al. (2019), read
  only by script 08.

## The scripts

| Script | What it does | Reads | Main outputs |
|---|---|---|---|
| `00_recompute_class_metrics.R` | Recomputes the class-level metrics from the rasters and fills the empty cells of two units | rasters, FRAGSTATS patch tables, `Data_Raw_FINAL.csv` | `Tabela_Metricas_Recomputada.csv`, Tables S32 and S33 |
| `01_prepare_data_and_ordination.R` | Builds the analytical dataset; NMDS, envfit, betadisper and PERMANOVA of composition | `Data_Raw_FINAL.csv`, script 00 table, shapefiles | `Data_Raw_WithCoords.csv`, Table 1, Figures 2 to 4 |
| `02_constrained_ordination_and_partitioning.R` | Collinearity screening, variation partitioning (landscape against space) and RDA of composition | `Data_Raw_WithCoords.csv` | Table 2, Figures 5, S1 and S8 |
| `03_univariate_model_selection.R` | Six specifications compared; the full a priori model; spatial sensitivity; selection uncertainty; genus-level slopes | `Data_Raw_WithCoords.csv`, Moran eigenvectors of script 02 | Table 3, Figure 7, Figures S2 to S7 and S10 |
| `04_geometric_coupling_diagnostic.R` | Checks the identity that ties patch density to the response | `Data_Raw_WithCoords.csv` | Table S26 |
| `05_ratio_artifact_test.R` | Permutation test of the ratio artifact; selects and reports the model of the manuscript | `Data_Raw_WithCoords.csv` | Table 4, Figure 6, Tables S28a to S28f, Figure S9 |
| `06_spatial_autocorrelation_checks.R` | Screens negative spatial autocorrelation, tests Moran's I of the residuals, refits the reported model with the Moran eigenvectors and tests its residuals with DHARMa | `Data_Raw_WithCoords.csv`, outputs of scripts 02 and 05 | Tables S34b and S38b |
| `08_sampling_effort_sensitivity.R` | Refits the reported model with the number of occurrence records as a covariate | `Data_Raw_WithCoords.csv`, occurrence records, shapefiles | Tables S42 and S42b |

`docs/SCRIPTS.md` describes each script in more detail, and `docs/TABLE_MAP.md` lists the file behind
every table and figure of the manuscript and of the supplementary material.

## How to run

1. Install the packages listed below.
2. Open R at the root of the repository (or open `ANALISES_TCC.Rproj` in the working folder of the
   study, which holds the repository).
3. Run the scripts in numerical order, or run `source("run_all.R")`. By default, `run_all.R` runs scripts
   02 to 08; set `RUN_SCRIPT_00` and `RUN_SCRIPT_01` to `TRUE` at its top to include scripts 00 and 01.
   Script 00 needs the rasters.

Script 02 must run before scripts 03 and 06, because it writes the Moran eigenvectors that they read;
script 06 also reads the reported model written by script 05.
Script 00 must run before script 01 at least once: without its table, two sampling units keep empty
cells. On the machine used for the study, scripts 02 to 08 ran in about two minutes, script 05 being
the slowest; script 00 takes a few minutes.

## Reproducibility

- The seed is fixed at 123 in every script.
- Permutation tests use 9,999 permutations, except the null of the AICc difference in script 05,
  which uses 2,000 because each permutation refits two models. Permutations whose fit ends with an error
  or a warning are skipped, and Tables S28b and S28c report how many entered each null.
- Cross-validation is exact leave-one-out, with leave-one-genus-out reported beside it.
- AICc and BIC are both reported.
- Term tests in the redundancy analysis are marginal; the sequential tests are exported for comparison.
- The redundancy analysis uses Hellinger-transformed class areas, so that the size of the sampling window
  does not enter the analysis of composition.
- `sessionInfo()` is printed at the end of every script and written to `session_info_part*.txt` in the
  output folder.
- Figures are written as TIFF files at 500 dpi with LZW compression, with Okabe-Ito colors, which remain
  distinguishable under the common forms of color vision deficiency.

### Packages

`here`, `terra`, `landscapemetrics`, `sf`, `tidyverse`, `readxl`, `vegan`, `adespatial`, `spdep`, `car`,
`MuMIn`, `DHARMa`, `lme4`, `performance`, `ggridges`, `ggrepel`, `patchwork` and `scales`.

`MASS` exports `select()` and `car` exports `recode()`, which mask the `dplyr` functions of the same
name. The scripts therefore call `dplyr::select()` where the masking could occur and do not use
`recode()`.

## Data sources

- **Land cover:** MapBiomas Collection 10.1, reference year 2024, 30 m resolution.
- **Occurrence records:** compiled by Macedo et al. (2019), filtered to records from 1990 to 2020.
- **Landscape metrics:** `landscapemetrics` for the class-level metrics; FRAGSTATS for the proximity
  index, which requires a search radius that `landscapemetrics` does not implement.

## Known limitations of the design

1. The sampling window is the response, so metrics normalized by it cannot be separated from that
   construction. Script 05 measures how far this reaches.
2. The occurrence records span 1990 to 2020, while the land cover map is from 2024.
3. A single kernel bandwidth, calibrated on *Panthera onca*, is applied to all taxa.
4. Genera contribute four to twelve units each, so genus-level estimates are exploratory.
5. The mean proximity index is not divided by the window area, but its search radius (28,230 m)
   exceeds the equivalent circular diameter of 54 of the 67 envelopes. Patches outside an envelope are
   not counted, so the index can rise with envelope size; the permutation null of script 05 does not
   test this dependence.

## What is not in this repository

The study area map (Figure 1) and the weighted conservation value index (Figure 8, Tables S36 and S37)
were built in ArcGIS Pro, not in R. Section 5 of the supplementary material states the geoprocessing
steps, and the input polygons are archived here, so both can be rebuilt.

## Changes

The commit history records the changes of each revision of the code.

## Citation

Please cite the manuscript once it is published. Until then, cite this repository.

## References

Brooks, T.M. et al. (2019) Measuring terrestrial Area of Habitat (AOH) and its utility for the IUCN Red
List. *Trends in Ecology & Evolution* 34, 977-986.
Macedo, L., Monjeau, A., Neves, A. (2019) Assessing the most irreplaceable protected areas for the
conservation of mammals in the Atlantic Forest: lessons for the governance of mosaics. *Sustainability*
11, 3029. https://doi.org/10.3390/su11113029
Tredennick, A.T., Hooker, G., Ellner, S.P., Adler, P.B. (2021) A practical guide to selecting models for
exploration, inference, and prediction in ecology. *Ecology* 102, e03336.
Yates, L.A., Aandahl, Z., Richards, S.A., Brook, B.W. (2023) Cross validation for model selection: a
review with examples from ecology. *Ecological Monographs* 93, e1557.
