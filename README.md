# Glycine betaine uptake and use in marine microbial communities

A repository for the code and analyses performed in Boysen et al. (in prep). Largely written in R

## Repo structure

This repository can be cloned or forked to replicate the analyses run in the untargeted portion of the GBT update and use manuscript.
It is a stable pipeline that accepts as input the mzXML files that can be downloaded from [Metabolights/MetabolomicsWorkbench (link soon!)]
and will produce several outputs. The pipeline is managed by the `output_control.Rmd` script in this repository, which has several sections
of code that must be changed to perform the full analysis. These sections are marked out in the Rmarkdown document (not yet). To run
the full analysis, both positive and negative mode HILIC mzXMLs and positive-mode RP mzXMLs should be downloaded into the mzXMLs_pos, mzXMLs_neg, and 
mzXMLs_cyano folders, respectively. Then, modify the script to point to the correct subset (pos/neg/cyano) and polarity (pos/neg). The script will
also need to be run twice, once for each experiment (north/south). The "output_folder_*" directories contain useful debugging and intermediate objects
from the xcms and custom code process. After the 
complete pipeline output is produced, figures and statistics for the manuscript can be produced by running the `qc_and_vis.Rmd` Rmarkdown.

Final repo structure after complete analysis:
  - mzXMLs_pos
  - mzXMLs_neg
  - mzXMLs_cyano
  - output_folder_cyano
  - output_folder_cyano_exp1
  - output_folder_pos
  - output_folder_pos_exp1
  - output_folder_neg
  - output_folder_neg_exp1
  - scripts
  - output_control.Rmd
  - qc_and_vis.Rmd

## Citation

TBD

## Contact

Please feel free to use this Github Issues page to raise questions or request clarification.
