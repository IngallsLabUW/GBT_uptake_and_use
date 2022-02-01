# Glycine betaine uptake and use in marine microbial communities

A repository for the code and analyses performed in Boysen et al. (in prep). Largely written in R

## Repo structure

This repository contains data and scripts for the targeted and untargeted portions of the GBT uptake and use manuscript.

### Raw data download

The raw data analyzed during this experiment can be obtained from Metabolomics Workbench at [Study ID ST002008](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR001273). These files contain both the raw data in .mzXML/.raw forms as well as the targeted and untargeted peaklists. 

Raw data:
  - ST002008_CYANO_QE_HF_Orbitrap_POS.zip (2.6G)
    - Contains .mzXML files for samples run on the CYANO HPLC column and QE HF Orbitrap in positive mode
  - ST002008_HILIC_QE_HF_Orbitrap_NEG.zip (620.6M)
    - Contains .mzXML files for samples run on the HILIC HPLC column and QE HF Orbitrap in negative mode
  - ST002008_HILIC_QE_HF_Orbitrap_POS.zip (2.2G)
    - Contains .mzXML files for samples run on the HILIC HPLC column and QE HF Orbitrap in positive mode
  - ST002008_HILIC_Waters_TQS_triple-quadrupole_POS.zip (15.6M)
    - Contains .raw files for samples run on the HILIC HPLC column and Waters TQS triple-quad

Peaklists:
  - ST002008_AN003272_Results.txt (4.3M)
    - Untargeted data from the "GBT Fate" experiments (both positive and negative mode, HILIC only)
  - ST002008_AN003273_Results.txt (104.3K)
    - Untargeted data from the "GBT Fate" experiments (negative mode, HILIC only)
  - ST002008_AN003274_Results.txt (4.5M)
    - Untargeted data from the "GBT Fate" experiments (positive mode, CYANO only)

### Targeted analysis

The targeted subdirectory contains two sub-folders, one for the GBT *uptake* experiments run on the TQS triple-quad and the other for the GBT *fate* analysis run on the QE HF Orbitrap.



### Untargeted analysis

The untargeted subdirectory contains a pipeline that accepts as input the mzXML files downloaded above and will produce several outputs.

The script expects to find the mzXML files in three separate directories of the untargeted folder, named "mzXMLs_cyano", "mzXMLs_neg", and "mzXMLs_pos", respectively. Zipped directories downloaded from Metabolomics Workbench can be extracted directly into the "untargeted" folder, which should automatically create the folders with the correct file names.

The pipeline is managed by the `output_control.Rmd` script in this subdirectory, which has several sections of code that must be changed to perform the full analysis. These sections are marked out in the Rmarkdown document (not yet). To run the full analysis, both positive and negative mode HILIC mzXMLs and positive-mode RP mzXMLs should be downloaded into the mzXMLs_pos, mzXMLs_neg, and mzXMLs_cyano folders, respectively. Then, modify the script to point to the correct subset (pos/neg/cyano), polarity (pos/neg), and experiment number (1 for northern experiment, 2 for the southern). The "output_folder_*" directories contain useful debugging and intermediate objects from the xcms and custom code process. After the complete pipeline output is produced, figures and statistics for the manuscript can be produced by running the `qc_and_vis.Rmd` Rmarkdown.

Final repo structure after complete analysis:
  - targeted
    - FateTidy
      - Cloned from github
    - UptakeKinetics
      - Cloned from github

  - untargeted
    - mzXMLs_pos
      - Created manually by extracting Metabolomics Workbench zipped file ST002008_HILIC_QE_HF_Orbitrap_POS.zip
    - mzXMLs_neg
      - Created manually by extracting Metabolomics Workbench zipped file ST002008_HILIC_QE_HF_Orbitrap_NEG.zip
    - mzXMLs_cyano
      - Created manually by extracting Metabolomics Workbench zipped file ST002008_CYANO_QE_HF_Orbitrap_POS.zip
    - output_folder_cyano
      - Created by running output_control.Rmd with setup chunk options polarity = "pos", identifier = "cyano", and exp.num = 2
    - output_folder_cyano_exp1
      - Created by running output_control.Rmd with setup chunk options polarity = "pos", identifier = "cyano", and exp.num = 1
    - output_folder_pos
      - Created by running output_control.Rmd with setup chunk options polarity = "pos", identifier = "pos", and exp.num = 2
    - output_folder_pos_exp1
      - Created by running output_control.Rmd with setup chunk options polarity = "pos", identifier = "pos", and exp.num = 1
    - output_folder_neg
      - Created by running output_control.Rmd with setup chunk options polarity = "neg", identifier = "neg", and exp.num = 2
    - output_folder_neg_exp1
      - Created by running output_control.Rmd with setup chunk options polarity = "neg", identifier = "neg", and exp.num = 2
    - scripts
      - Cloned from GitHub
    - output_control.Rmd
      - Cloned from GitHub
    - qc_and_vis.Rmd
      - Cloned from GitHub

## Citation

TBD

## Contact

Please feel free to use this Github Issues page to raise questions or request clarification.
