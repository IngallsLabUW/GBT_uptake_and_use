# Gradients-3-GBT
Total Hydrolyzable Amino Acid (THAA) data processing for Glycine Betaine Fate Experiments from Gradients 3 cruise

## Script 1:
"THAA.quantification.Step1.R"

This script uses raw peak areas from skyline on the unlabled amino acids that were run at the same time as standard additions so that we can quantify the concentration in the various samples
The input is:
"StandardCurvePeakAreaData.csv" which is the necessary peak area data from skyline and
"Sample.Log.csv" which has necessary metadata to calculate the final concentrations in seawater rather than in the lcms vial

the output is:
"Quantified_THAA_inSeaWater.csv" which has the concentrations of the amino acids in the different samples


## Script 2:
"THAA.1.QC.and.MID.Calc.R"

This reads in the possible isotopologues (stored in: "RawData/THAA_Isotope_possibilities.csv") 
and the Skyline report (stored in: "RawData/GBT_Fate_THAA_Isotopes_skyline_report.csv")

The script does a quality control filtering of the raw Skyline output and exports 
1) "QC.isotopes.THAA.csv" which has all the data with the various QC flags
2) "QC.isotopes.THAA.blk.sub.csv" which is the same as the prior output with additional columns that have the peak areas less that in the associated blank samples
3) "THAA.isotopologue.frequency.csv" which compiles the number of observations of each THAA isotopologue in each experiment

The script then calculates mass isotopologue distributions (MID) for each THAA that was measured in multiple samples
The outputs of this section are
1) "MID_QC_blkSubtract.pdf" which is a document that has bar charts of the MID for each THAA
2) "MID_QC_blkSubtract.csv" whih has the MID data for each individual sample
3) "MID_mean_QC_blkSubtract.csv" which has the MID data averaged across the biological triplicates for each time point
4) "THAA_MID_t.test_significance_test.csv" which has the results for significance tests to determine if the MID was significantly differet at T0 and Tfinal
5) "THAA_MID_increasing_Isotopologues_new.scale.pdf" which plots select THAA isotopologues that had MID changes in the southern experiment (as in Figure 6 in the manuscript)

## Script 3:
"Plot.THAA.conc.T0.TLong.R"

This makes several plots that compare individual THAAs and total THAA concentrations for the initial and final time points in both experiments.
This includes "THAA_total_C.pdf" which is Supplemental Figure 7 in the manuscript.

Inputs:
1) "Quantified_THAA_inSeaWater.csv" generated in script 1
2) "AA_Molecular_Formulas.csv" which just has a parsed molecular formula for each amino acid




