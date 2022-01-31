# Gradients-3-GBT
This folder contains the code and raw data and example output for the targeted metabolite data and analyses done in the GBT Kinetics experiments conducted on the Gradients 3 cruise in 2019.

## Script 1
"Quantify_GBT_Aug8_2019.R"
This script reads in the output from skyline peak integration of TQS runs measuring 5C13, 1N15-GBT and C12, N14-GBT from the kinetics experiments and the associated standard addition curves that were run.
Using that data it quantifies the two forms of GBT in the samples and then calculates the uptake kinetics using a few methods. 

Inputs:
1) "GBT Kintetics Gradients 3 Sample List TQS for R code use.csv" Sample list
2) "QC_output_All_TQS_GBT_IS-GBT.csv" the quality control filtered skyline output
3) "BiologicNormalizationInformation.csv" the relevant normalization values for the different samples - dilution factor due to dilutions before running and water volume filtered, etc

Outputs:
1) "All.Standard.Addition.Models.csv" - different models for quantifying GBT. We ran several standard addition curves so the models use various combinations to calculate final GBT concentrations.
2) "GBT_Concentrations_and_Uptake_Rates_w_StdError.csv" - the best quantified GBT data
3) "Comparing.quantification.estimates.uptake.rates.pdf"
4) "GBT_Uptake_Rates_w_StdError.csv" - update rates calculated from the quantified GBT data
5) "GBT K1 Model.rdata", "GBT K1 Model no explicit dGBT.rdata", "GBT K1 Model no 2000.rdata" different uptake kinetics models for the first experiment
6) "GBT K1 Model coeff.csv" the model coefficients for the GBT uptake kinetics from the first experiment
7) "GBT K2 Model.rdata", "GBT K2 Model no explicit dGBT.rdata", "GBT K2 Model no 2000.rdata" different uptake kinetics models for the second experiment
8) "GBT K2 Model coeff.csv" the model coefficients for the GBT uptake kinetics from the second experiment

## Script 2
"Monte.Carlo.Uptake.Model.R"
This script takes the initial uptake rate data as a starting place and does a monte carlo error analysis to determine the uptake kintetics and associated uncertainty
Note that this takes a long time to run since it generates 1000 models
It also does the analysis two ways for each experiment - parameterizing dGBT and without parameterizing dGBT.

Inputs:
1) "GBT_Uptake_Rates_w_StdError.csv"  from Script 1

Output:
1) "GBT Uptake K1 all Monte Carlo Model Resuls.csv" and "GBT Uptake K2 all Monte Carlo Model Resuls.csv"
2) "GBT Uptake K1 all Monte Carlo Model Resuls no dGBT.csv" and "GBT Uptake K2 all Monte Carlo Model Resuls no dGBT.csv"
3) "monteCarlo_uptake_models_K1.pdf" and "monteCarlo_uptake_models_K2.pdf"  which shows outliers calculated a variety of ways
4) "K1-MonteCarlo_OutliersRemoved_summary_withError.csv" and "K2-MonteCarlo_OutliersRemoved_summary_withError.csv"
5) "K1-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv" and "K2-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv"

## Script 3
"TurnoverTime.R"
turnover time calculation based on the Wright-Hobbie linear transformation
Plot t/f (incubation time in hour divided by ration of added that was taken up) vs A (added concentration in nM). 
Fit least squares regression. Then Kt+S = -x intercept and Vmax = 1/slope and turnover time = y-intercept

Inputs:
1) "GBT_Concentrations_and_Uptake_Rates_w_StdError.csv" from Script 1 above
2) "K1-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv" and "K2-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv" from script 2 above

Outputs:
1) "GBT K1 linear transformation all Monte Carlo Model Resuls.csv" and "GBT K2 linear transformation all Monte Carlo Model Resuls.csv" writes out the model and the monte carlo error results for the Wright-Hobbie linear transformation for the two experiments.
2) "K1-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv" and "k2-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv" writes out the model and the monte carlo error results for the Wright-Hobbie linear transformation for the two experiments after removing outliers.


## Script 4
"Plot.Final.uptake.models.R"
This generates the plots and tables associated with the GBT uptake kinetics that are in the manuscript.

Inputs:
1) "GBT_Uptake_Rates_w_StdError.csv"  from script 1
2) "K1-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv" and "K2-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv" from script 2 above
3) "GBT K1 Model coeff.csv" and "GBT K2 Model coeff.csv" from script 1 above
4) "K1-MonteCarlo_OutliersRemoved_summary_withError.csv" and "K2-MonteCarlo_OutliersRemoved_summary_withError.csv" from script 2 above
5) "K1-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv" and "k2-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv" from script 3 above
6) "CTD Temperature (2).csv" which is provided in the raw data folder
7) "CTD Chloropigment.csv" which is provided in the raw data folder
8) "PCPN.Dat.R" which is provided in the raw data folder

Outputs:
1) "Uptake-Kinetics-with-rangeOfFit-K1andK2_insets.pdf" which is a plot like the one in the manuscript - figure 2 top pannels
2) "xtable1.tex" is a latex table similar to that in the manuscript. It has for each station and each kinetics experiment the oceanographic parameters: temperature, chl and the model parameters and error estimates.
3) "Experiment_Conditions_new.tex" is similar but has more oceanographic information including the flow cytometry based cell biomass for different plankton groups and total particulat carbon and nitrogen.
4) "Uptake_Kinetic_Parameters.tex" is a table with just the uptake kinetics parameters
5) "Uptake_Kinetic_2_vs_3_Parameters.tex" is a table showing the diference in the uptake kinetic parameters when modeled with two or three parameters. This is a supplemental table in the mansucript.
