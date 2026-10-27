# Where each table and figure comes from

The files written by the scripts are named after the part of the chain that produces them, while the
manuscript and the supplementary material number their tables and figures in reading order. The two
numbering schemes do not coincide, so this file maps one onto the other. Every number reported in the
two documents can be traced to a row of a file listed here.

Paths are relative to `Outputs/Manuscrito/`.

## Main text

| Item | Produced by | File |
|---|---|---|
| Figure 1, study area map | ArcGIS Pro, outside the R chain | not in the repository; see Supplementary Section 5 |
| Figure 2, composition ridgelines | `01` | `PartI/Figure2_Ridgeline_Composition.tiff` |
| Figure 3, configuration ridgelines | `01` | `PartI/Figure3_Ridgeline_Configuration.tiff` |
| Figure 4, NMDS | `01` | `PartI/Figure1_NMDS_Composition_Genus.tiff` |
| Figure 5, RDA biplot | `02` | `PartII/Figure_RDA_Biplot.tiff` |
| Figure 6, coefficients of the reported model | `05` | `PartV/Figure_PV_02_Reported_Model_Coefficients.tiff` |
| Figure 7, genus-level effect heatmap | `03` | `Figure_PIII_08_Effect_Heatmap.tiff` |
| Figure 8, conservation value index | ArcGIS Pro, outside the R chain | not in the repository |
| Table 1, global PERMANOVA | `01` | `PartI/TableS4_PERMANOVA_Global.csv` |
| Table 2, marginal tests of the RDA | `02` | `PartII/TableII_4_RDA_Anova_Marginal.csv` |
| Section 3.3, RDA biplot scores and the correlation of the first-axis site scores with forest cover | `02` | `PartII/TableII_4c_RDA_Biplot_Scores.csv`, `PartII/TableII_4d_RDA_Axis1_Forest.csv` |
| Section 3.3, variation partitioning | `02` | `PartII/TableII_3_VP_Fractions.csv`, `PartII/TableII_3_VP_Significance.csv`; the global spatial test is in `PartII/TableII_4_RDA_GlobalSummary.csv` and, with the outcome of the selection, in `PartII/TableII_2a_MEM_Global_Test.csv` |
| Table 3, full a priori model | `03` | `Table_S_PartIII_ConfirmatoryModel.csv`; the sums of Akaike weights of the last column are in `Table_S_PartIII_VariableImportance.csv` |
| Table 4, reported model of envelope extent | `05` | `PartV/TableS_PartV_Uncoupled_Model.csv` |

## Supplementary material

| Item | Produced by | File |
|---|---|---|
| S1, S2 | written by hand | no producing script |
| S3, S4 | `01` | `PartI/Table3_Medians_Composition_PLAND.csv`, `PartI/Table3b_Medians_MinorClasses.csv`, `PartI/Table1_Composition_By_Species.csv`, `PartI/Table1b_MinorClasses_And_LPI_By_Species.csv` |
| S5, S6 | `01` | `PartI/Table4_Medians_Configuration.csv`, `PartI/Table2_Configuration_By_Species.csv` |
| S7, envfit | `01` | `PartI/TableS2_Envfit_Configuration.csv` |
| S8, pairwise PERMANOVA by genus | `01` | `PartI/TableS5_PERMANOVA_Pairwise_Genus.csv` |
| S9, VIF of the RDA pool | `02` | `PartII/TableII_1_VIF_RDA.csv` |
| S10, VIF of the partitioning pool | `02` | `PartII/TableII_1_VIF_VP.csv` |
| S11, betadisper | `01` | `PartI/TableS3_Betadisper.csv` |
| S12, S13, S14, pairwise PERMANOVA | `01` | `PartI/TableS5_PERMANOVA_Pairwise_Species.csv`, `..._Diet.csv`, `..._Locomotion.csv` |
| S15, descriptive statistics | `03` | `Table_S_PartIII_DescriptiveStatistics.csv` |
| S16, comparison of the six specifications | `03` | `Table_S_PartIII_ModelComparison.csv` |
| S17, residual diagnostics | `03` | `Table_S_PartIII_Diagnostics.csv` |
| S18, VIF of the a priori model | `03` | `Table_S_PartIII_VIF.csv` |
| S18a, full a priori model | `03` | `Table_S_PartIII_ConfirmatoryModel.csv` |
| S18b, model averaging | `03` | `Table_S_PartIII_Averaging_FullVsConditional.csv` |
| S19, selection table | `03` | `Table_S_PartIII_Dredge_FullTable.csv` |
| S20, variable importance | `03` | `Table_S_PartIII_VariableImportance.csv` |
| S21, stability check | `03` | `Table_S_PartIII_StabilityCheck.csv` |
| S22, FRAC against patch area | `03` | `Table_S_PartIII_FRAC_Area_Correlation.csv` |
| S23, matrix sensitivity | `03` | `Table_S_PartIII_Sensitivity_MatrixModel.csv` |
| S24, genus-level regressions | `03` | `Table_S_PartIII_EffectByGenus.csv` |
| S25, genus-level Spearman correlations | `03` | `Table_S_PartIII_MetricAssociation_ByGenus.csv` |
| S26, mixed-model coupling controls | `04` | `TableS28_Geometric_Coupling.csv` |
| S27, sample sizes | `01`, `02` | `PartI/TableS_Sample_Sizes_PartI.csv`, `PartII/TableII_0_Sample_Sizes.csv` |
| S28a, six specifications | `05` | `PartV/TableS_PartV_Specifications.csv` |
| S28b, permutation null | `05` | `PartV/TableS_PartV_NullSimulation.csv` |
| S28c, AICc under the null | `05` | `PartV/TableS_PartV_AICc_Under_Null.csv` |
| S28d, decision table | `05` | `PartV/TableS_PartV_Decision.csv` |
| S28e, algebraic coupling | `05` | `PartV/TableS_PartV_Identity_and_PartialCorrelations.csv` |
| S28f, gain from the interaction | `05` | `PartV/TableS_PartV_Interaction_Gain.csv` |
| S29, genus as moderator | `03` | `Table_S_PartIII_SlopeByGenus_Moderator.csv` |
| S30, proximity against amount | `05` | `PartV/TableS_PartV_Proximity_Amount_Correlations.csv` |
| S31, proximity conditioned on amount | `05` | `PartV/TableS_PartV_Proximity_Conditioned.csv` |
| S32, agreement of the recomputation | `00` | `Part0/TableS_Recompute_Agreement.csv` |
| S33, repaired cells | `00` | `Part0/TableS_Repaired_Cells.csv` |
| S34, spatial sensitivity | `03` | `Table_S_PartIII_SpatialSensitivity.csv` |
| S34b, reported model refitted with the Moran eigenvectors | `06` | `PartII/TableII_6b_Reported_Model_With_MEMs.csv`; the Moran's I of its residuals, without and with the eigenvectors, are in `PartII/TableII_6c_Moran_I_Residuals.csv`, and the DHARMa tests of dispersion, uniformity and outliers in `PartII/TableII_6d_DHARMa_Reported_Model.csv` |
| S35, criteria and out-of-sample error | `05` | `PartV/TableS_PartV_CV_And_BIC.csv` |
| S36, S37, conservation value index | ArcGIS Pro, outside the R chain | not in the repository |
| S38, partition without the metrics normalized by the window | `02` | `PartII/TableII_3c_VP_Fractions_Without_Ratio_Metrics.csv`; the two partitions side by side are in `PartII/TableII_3d_VP_Sensitivity_Summary.csv` |
| S38b, screening of the spatial eigenvectors | `02`, `06` | `PartII/TableII_2a_MEM_Global_Test.csv`, `PartII/TableII_2_MEMs_Selected.csv`, `PartII/TableII_6_MEM_Global_Tests.csv`, `PartII/TableII_6c_Moran_I_Residuals.csv` |
| S39, effect of the ordination transformation | `02` | `PartII/TableII_4b_RDA_Transformation_Comparison.csv`, `PartII/TableII_5_Response_Size_Dependence.csv` |
| S39b, sequential term tests of the RDA | `02` | `PartII/TableII_4_RDA_Anova_Terms_Sequential.csv` |
| S40, farming matrix without window-normalized predictors | `05` | `PartV/TableS_PartV_Matrix_Uncoupled.csv`; the second specification in the file, with matrix patch density, contains a ratio of the window and is a contrast, not an uncoupled estimate |
| S41, thinning test for spatial coincidence | `05` | `PartV/TableS_PartV_Spatial_Thinning.csv` |
| S42, sampling effort bounds | `08` | `PartVI/TableS42_Sampling_Effort.csv` |
| S42b, sampling effort robustness | `08` | `PartVI/TableS42b_Sampling_Effort_Robustness.csv` |
| Figure S1, correlation matrix | `02` | `PartII/FigS1_Collinearity_Pearson.tiff` |
| Figure S2, response transformation | `03` | `Figure_PIII_01_Response_Transformation.tiff` |
| Figure S3, variable importance | `03` | `Figure_PIII_03_Variable_Importance.tiff` |
| Figure S4, partial effects | `03` | `Figure_PIII_04_Partial_Effects.tiff` |
| Figure S5, observed against predicted | `03` | `Figure_PIII_05_Observed_vs_Predicted.tiff`; the two values quoted in the caption are in `Table_S_PartIII_ObsPred_R2.csv` |
| Figure S6, residual diagnostics | `03` | `Figure_PIII_06_Residual_Diagnostics.tiff` |
| Figure S7, genus-level Spearman heatmap | `03` | `Figure_PIII_07_Metric_Association_Heatmap.tiff` |
| Figure S8, variation partitioning | `02` | `PartII/Figure_VP_Venn.tiff` |
| Figure S9, permutation null | `05` | `PartV/Figure_PV_01_Null_Distribution.tiff` |
| Figure S10, coefficients of the a priori set | `03` | `Figure_PIII_02_Coefficient_ForestPlot.tiff` |

## Codes used in the figures

Figure S1 and the RDA biplot (Figure 5) use the codes of Table 2: the suffix F marks a metric of forest,
A of the farming matrix, and H of herbaceous and shrub vegetation. Elevation and human population
density are written in full.

## Which model the figures of script 03 show

Section 5.2 of script `03` chooses the retained specification among the models within two AICc units,
by the smaller cross-validated error and, where none is available, by the smaller number of parameters.
On these data that is the Gamma model with a logarithmic link and the a priori interaction, so
`global_model` is a `glm`, and Figures S5 and S6 show that model. Figure S4 shows the competitive model
with the highest Akaike weight.
