# Gradients-3-GBT
This contains the full scripts for the targeted metabolite data processing for the GBT Fate experiments conducted on the Gradients 3 cruise

## Script 1-A
"skyline.standards.qc.R"
This takes several different files that are skyline integration output for various compounds, combines them, and does quality control.
This script is for the injections that were spiked with internal standards.
The user has to specify whether the LCMS files of interest were the reversed phase (CyanoAq) or HILIC chromatography and the ionization mode (HILICPos or HILICNeg).
For each the output will be a "Fraction_QC_skyline_IS_data.csv" which is the completed QC file, and "Fraction_frequency_of_standards_IS_data.csv" which is a counts table for the instances of each metabolite (only considering unlabled isotopologues) that we observed in each experiment.

## Script 1-B
"skyline.targeted.isotopes.qc.R"
This script is for the injections that were not spiked with internal standards.
This takes several different files that are skyline integration output for various compounds, combines them, and does quality control.
The user has to specify whether the LCMS files of interest were the reversed phase (CyanoAq) or HILIC chromatography and the ionization mode (HILICPos or HILICNeg).
For each the output will be a "Fraction_QC_skyline_isotope_data.csv" which is the completed QC file, "Fraction_QC_skyline_isotope_data_blankSubtract.csv" which has the QC data with the blank subtracted values, and "Fraction_frequency_of_isotope_data.csv" which is a counts table for the instances of each metabolite that we observed in each experiment.

## Script 2
"Fraction_skyline_isotopes_BMIS.rmd"
This script normalizes the targeted metabolomics data using the best-matched internal standard method (Boysen and Heal, Analytical Chemistry, 2018).
There is a version of this for each fraction (CyanoAq, HILICNeg, HILICPos).
Inputs: "Fraction_QC_skyline_IS_data.csv" generated from Script 1. Also needs the Sample list for each fraction.
Outputs: 
1) "Fraction_QC_skyline_IS_data_BMISed.csv" which has normalized values. Again, this is only done for unlabeled metabolites.
2) "Fraction__BMIS Model.csv" which just has each analyte and its corresopnding internal standard that is best for normalization
3) "Fraction_BMIS_multipliers.csv" which has the multiplier that should be used when normalizing to any of the internal standards that elute in that fraction/ionization for each sample (aka the peak are for that internal standard in that sample over the average peak area for that internal standard)
4) "Fraction_HILICNeg_IS Means.csv" which has the mean peak area for each internal standard.

## Script 3
"MID.calculations.R"
This script uses the quality controlled skyline output from 1-B to calculate the mass isotopologue distributions (MID) for each fraction.
Inputs
1) "Isotope_Possibilities_IngallsStandards_trim.csv" and "All_C-N_Isotope_possibilities_Ingalls_Standards.csv" which have the possible isotopologues we targeted
2) Outputs of Script 1B: "Fraction_QC_skyline_isotope_data_blankSubtract.csv" and "Fraction_frequency_of_isotope_data.csv"
Outputs:
1) "MID_mean_QC_Targeted_blkSubtract_Fraction.csv" and "MID_QC_Targeted_blkSubtract_Fraction.csv" which have the MID values for each isotopologue for all the samples (latter) and the mean of the triplicates (former).
2) "Fraction_goodCompounds_MID_lm_significance_test.csv" that has the results of an ANOVA and of a linear model of MID vs time to determine which isotopologues have significant trends in the two experiments. 
3) "MID_QC_Targeted_blkSubtract_Fraction.pdf" that has bar plots of the MID for the metabolites over time in both experiments.

## Script 4
"DMG.plot.r"
This plot generates Figure S4 from the manuscript - the peak area of the 13C-4, 15N-1 dimethylglycine isotopologue over time in both experiments
Input:
1) possible isotopologues: "Isotope_Possibilities_IngallsStandards_trim.csv", "All_C-N_Isotope_possibilities_Ingalls_Standards.csv",
2) Quality control data from script 1-B "HILICPos_QC_skyline_isotope_data_blankSubtract.csv" and "HILICPos_frequency_of_isotope_data.csv"
3) Sample log "Sample_Log_GBTFateExperiments.csv"
Outputs: "DMG_labeled_noRects.pdf"

## Script 5
"Quantify_GBT.R"
Self explanitory...
Inputs:
1) "HILICPos_QC_skyline_IS_data_BMISed.csv" internal standard spiked normalized data from Script 2
2) "MID_QC_Targeted_blkSubtractHILICPos.csv" MID calculations from Script 3
3) The sample log
Outputs:
1) "GBT_concentrations.csv" which has all the concentrations calculated for all isotopolouges. "GBT.nM" is the the nM concentration of that single GBT isotopologue in marine particles. The error is propagated so a standard deviation is included as "GBT.nM.sd"
2) "GBT_concentrations_aves.csv" has the average concentrations for all isotopologues considering the biological triplicates
3) "GBT_Concentration_Over_Time_withInset.pdf" Which is similar to the bottom panels of Figure 2 and the top panels of Figure S2 in the manuscript with GBT concentrations over time for the various isotopologues in the two experiments.
4) "Turnovertime.particulate.GBT.csv" is a simple estimate of turnover time from calculated uptake rates and initial particulate GBT concentrations.


## Script 6
"Quantify_knowns.R"
This script takes in the QC data from the internal standard spiked runs (from Script 1-A) and the unspiked runs (from script 1-B) 
as well as the BMIS data (Script 2) and the MID data (Script 3).
It also combines the best estimates for all the other compounds with that for GBT, so it reads in the output from Script 5
Using this, it calculates an estimated nM concentration for the known comopunds (nM particulate metabolite).
Outputs:
1) "IE_RF_Fraction.csv" which has the ionization efficieny and response factors for each of our metabolite standards
2) "Estimates_of_uMInVial_IS_samples_fraction.csv" which has the ESTIMATED concentrations as uM in the LCMS vial. Note that this is not the final estimate since compounds with internal standards and external standard curves will still be refined.
3) "Estimates_of_nMInSample_IS_samples_Fraction.csv" which just has ESTIMATED concentrations as nM in sample. Note that this is not the final estimate.
4) "nM.concentrations.Fraction.csv" this is the best and final estimate concentrations including the appropriate estimates for GBT and for compounds that we have internal standards for

## Script 7
"Combine.Fractions.R"
This takes in the nM data, ie.rf data, MID data, frequency data, from all the fractions and combines them into single files.
Outputs:
1) "nM.concentrations.allFractions.csv" nM concentration of all metabolites and isotopologues in the samples
2) "IE_RF_allFractions.csv" ionization efficiency and response factors for all of the standard metabolites
3) "MID_mean_blkSubtract_allFractions.csv" mean MIDs for the triplicates
4) "MID_blkSubtract_allFractions.csv" MIDs for all replicates
5) "SuppTab_MID_significant_trend.csv" Isotopologues that have significant trends over time in either expeirment
6) "SuppTab_Average_nM_and_MID.csv" Combined MID values and nM values for all isotopologues 
7) "SuppTab_list_of_metabolites_we_lookedAt.csv" A list of metabolites that had significant trends in MID in either experiment

## Script 8
"MID.line.plots.no.base.peak.R"
This code generates line plots of MID values over time in both experiments for a variety of metabolites. 
The plots do not include the fully unlabled (all 12C and 14C) isotopologue because that overwhelms the signal and makes the y-axis hard to read for the remaining isotopologues.
This code was used in generating Figure 5 in the manuscript.
Inputs:
1) possible isotopologues: "Isotope_Possibilities_IngallsStandards_trim.csv", "All_C-N_Isotope_possibilities_Ingalls_Standards.csv",
2) All MID data - "MID_blkSubtract_allFractions.csv" from Script 7
3) Natural abundance data- estimated: "Natural_abundance_of_Standards.csv" and measured: "Natural_abundance_HILICPos_from_G3.csv"
4) Sample log
Outputs:
1) "MID.Lines.Small.and.Nuc.pdf" This is Figure 5 in the manuscript
2) "MID.Lines.Big.Change.New.S.N.order.2.pdf" This is figure S5 in the manuscript

## Script 9
"Plot.Known.Conc.And.Rates.R"
This script plots the concentration (in nM, nM C, nM 13C, nm N, or nM 15N) of metabolites and specific isotopologues over time in both experiments.
This script was used in generating Figures 3, S1, S3, and S6 in the manuscript.
Inputs:
1) Sample Log
2) "nM.concentrations.allFractions.csv" from Script 7
Outputs:
1) "isotopologues.rates.csv" which is a simple calculation of 'instantanous' rates of uptake/production of various isotopologues
2) "Total.conc.nM.N.select.metabs.v.time.new.Order.pdf" which is figure S6 in the manuscript.
3) "Stacked_13C_15N_relative_concentration_OverTime_noRects_switchOrder.pdf" which is similar to Figure S3 in the manuscript (different color pallet)
4) "Stacked_13C_15N_concentration_OverTime_noRects_switchOrder.pdf" which is Figure 3 in the manuscript
5) "Osmolytes.v.time.bars.abs.and.rel.no.rect.pdf" which is Figure S1 in the manuscript

## Script 10
GBT.concentration.and.MID.plot.R
This script generates the GBT concentration and MID plots vs time
Inputs:
1) "GBT_concentrations.csv" generated from script 5 above
2) "Natural_abundance_of_Standards.csv" which is provided in the raw data folder
3) "MID_mean_QC_Targeted_blkSubtractHILICPos.csv" which is generated in script 3 above.
4) the sample log
5) "Natural_abundance_HILICPos_from_G3.csv" which is provided in the raw data folder.
Outputs:
1) GBT_Conc_MID.new.order.pdf is Supplemental Figure 2 in the manuscript - GBT concentrations vs time and MID vs time in the two fate experiments
2) GBT_Conc_v_time_newColors.pdf is the bottom panels of Figure 2 in the manuscript - GBT concentrations vs time for the two experiments 

## Script 11
Boxplots.of.FullLabeled.R
This makes the heatmaps of the MID of the 'fully labeled'/directly transformed from GBT isotopologues used in Figures 4 and 7 (pathway figures)
Inputs:
1) MID averages from Script 3
2) MID of all samples from Script 3
3) Sample log
Outputs:
1) "heatmaps.of.direct.gbt.transformation.metabolites.pdf"