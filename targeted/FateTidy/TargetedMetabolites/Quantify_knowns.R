## quantify knowns

## Concentration in vial = Peak Area / ionization efficiency *  1 / RF-ratio
## Where RF-ratio (response factor ratio) is the difference between the peak area of a standard spiked into environmental matrix less the peak area in un-spiked matrix, divided by the peak area of the standard in water, 
## and the IE (ionization efficiency) is the ratio of the peak area of the standard in water over the concentration of the standard.

## libs ----------
library(tidyverse)
library(stats)
library(RCurl)
library(here)
library(readr)


## set column -------
# column.mode = "HILICPos"
# column.mode = "HILICNeg"
column.mode = "CyanoAq"

## load QC-ed, MID, and BMISed data -------------
if(column.mode == "HILICPos"){
  QC.dat = read_csv("FateTidy/output/TargetedMetabolites/HILICPos_QC_skyline_isotope_data_blankSubtract.csv")
  QC.dat.IS = read_csv("FateTidy/output/TargetedMetabolites/HILICPos_QC_skyline_IS_data.csv")
  BMIS.dat = read_csv("FateTidy/output/TargetedMetabolites/HILICPos_QC_skyline_IS_data_BMISed.csv")
  MID = read_csv("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_blkSubtractHILICPos.csv")
} else if(column.mode == "HILICNeg"){
  QC.dat = read_csv("FateTidy/output/TargetedMetabolites/HILICNeg_QC_skyline_isotope_data_blankSubtract.csv")
  QC.dat.IS = read_csv("FateTidy/output/TargetedMetabolites/HILICNeg_QC_skyline_IS_data.csv")
  BMIS.dat = read_csv("FateTidy/output/TargetedMetabolites/HILICNeg_QC_skyline_IS_data_BMISed.csv")
  MID = read_csv("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_blkSubtractHILICNeg.csv")
}else if(column.mode == "CyanoAq"){
  QC.dat = read_csv("FateTidy/output/TargetedMetabolites/CyanoAq_QC_skyline_isotope_data_blankSubtract.csv")
  QC.dat.IS = read_csv("FateTidy/output/TargetedMetabolites/CyanoAq_QC_skyline_IS_data.csv") %>%
    filter(!grepl("-pos",`Replicate Name`))
  BMIS.dat = read_csv("FateTidy/output/TargetedMetabolites/CyanoAq_QC_skyline_IS_data_BMISed.csv")
  MID = read_csv("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_blkSubtractCyanoAq.csv")
}

#Get list of better names -------------
urlfile <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"

stds.dat<-read_csv(url(urlfile)) %>%
  dplyr::rename(Compound.Name= Compound_Name)

if(column.mode == "HILICPos"){
stds.dat.to.use <- stds.dat %>%
  filter(Column == "HILIC",
         z ==1) %>%
  dplyr::rename(Identification = Compound.Name,
         BestMatch = Compound_Name_Figure) %>%
  dplyr::select(BestMatch, Identification, HILIC_Mix, Concentration_uM) %>% 
    unique() %>%
  mutate(Identification = ifelse(BestMatch=="Glucosylglycerol",BestMatch,Identification))
} else if(column.mode == "HILICNeg"){
  stds.dat.to.use <- stds.dat %>%
    filter(Column == "HILIC",
           z ==-1) %>%
    dplyr::rename(Identification = Compound.Name,
                  BestMatch = Compound_Name_Figure) %>%
    dplyr::select(BestMatch, Identification, HILIC_Mix, Concentration_uM) %>% 
    unique()
} else if(column.mode == "CyanoAq"){
  stds.dat.to.use <- stds.dat %>%
    filter(Column == "RP",
           z ==1) %>%
    dplyr::rename(Identification = Compound.Name,
                  BestMatch = Compound_Name_Figure) %>%
    dplyr::select(BestMatch, Identification, HILIC_Mix, Concentration_uM) %>% 
    unique()
}

## get out the Stds from the QC data --------
if(column.mode =="HILICPos"){
QC.dat.IS.stds = QC.dat.IS %>%
  filter(type=="Std",
         `Protein Name` != "Internal Standards") %>%
  left_join(., stds.dat.to.use, by = c("Protein Name"="Identification")) %>%
  filter(!is.na(Concentration_uM)) 
DMSP.etc = QC.dat.IS %>%
  filter(type=="Std",
         `Protein Name` != "Internal Standards") %>%
  left_join(., stds.dat.to.use, by = c("Protein Name"="BestMatch"))%>%
  filter(!is.na(Concentration_uM)) %>%
  filter(grepl("Glucosylglycero",Identification)|
           grepl("Dimethylsulfoniop",Identification)|
           grepl("Serine",Identification)|
           grepl("Asparagine",Identification)|
           grepl("Aspartic acid",Identification)|
           grepl("Hydroxylysine",Identification))
QC.dat.IS.stds <- full_join(QC.dat.IS.stds,DMSP.etc)%>%
  filter(!is.na(Concentration_uM)) %>%
  mutate(MixID = ifelse(grepl("Mix1", `Replicate Name`),"Mix1",
                        ifelse(grepl("Mix2", `Replicate Name`),"Mix2","noMix")),
         solventID = ifelse(grepl("InMatrix",`Replicate Name`),"Matrix","Water"))
} else if(column.mode =="HILICNeg"){
  QC.dat.IS.stds = QC.dat.IS %>%
    filter(type=="Std",
           `Protein Name` != "Internal Standards") %>%
    left_join(., stds.dat.to.use, by = c("Protein Name"="Identification")) %>%
    filter(!is.na(Concentration_uM)) %>%
    mutate(MixID = ifelse(grepl("Mix1", `Replicate Name`),"Mix1",
                          ifelse(grepl("Mix2", `Replicate Name`),"Mix2","noMix")),
           solventID = ifelse(grepl("InMatrix",`Replicate Name`),"Matrix","Water"))
} else if(column.mode =="CyanoAq"){
  QC.dat.IS.stds = QC.dat.IS %>%
    filter(type=="Std",
           `Protein Name` != "Internal Standards") %>%
    left_join(., stds.dat.to.use, by = c("Protein Name"="Identification")) %>%
    filter(!is.na(Concentration_uM)) %>%
    mutate(MixID = ifelse(grepl("Mix1", `Replicate Name`),"Mix1",
                          ifelse(grepl("Mix2", `Replicate Name`),"Mix2","noMix")),
           solventID = ifelse(grepl("InMatrix",`Replicate Name`),"Matrix","Water"),
           HILIC_Mix = "noMix")
}
## determine the IE: peak area of the standard in water to the concentration added -----------
IE = QC.dat.IS.stds %>%
  filter(HILIC_Mix == MixID,
         solventID == "Water") %>%
  mutate(IE = Area/Concentration_uM) %>%
  group_by(`Protein Name`, `Precursor Ion Name`,
           Fraction1, BestMatch) %>%
  summarise(IE.rsd = sd(IE)/mean(IE),
            IE = mean(IE))

## determine the RF ratio ----------------------------
## Where RF-ratio (response factor ratio) is the difference between the peak area of a standard spiked into environmental matrix less the peak area in un-spiked matrix, divided by the peak area of the standard in water, 
if (column.mode =="CyanoAq"){
  RF = QC.dat.IS.stds %>%
    mutate(Std.key = ifelse(grepl("H2OinMatrix",`Replicate Name`),
                            "H2OinMatrix",
                            ifelse(grepl("StdsInMatrix",`Replicate Name`),
                                   "StdsInMatrix","StdsInWater"
                                   ))) %>%
    filter(HILIC_Mix == MixID | MixID == "noMix") %>%
    group_by(Std.key, solventID,`Protein Name`, `Precursor Ion Name`, BestMatch) %>%
    summarise(Area = mean(Area, na.rm = T)) %>%
    mutate(Area = ifelse(is.na(Area),0,Area)) %>%
    ungroup() %>%
    dplyr::select( -solventID) %>%
    spread(Std.key, Area) %>%
    mutate(RF = ((StdsInMatrix - H2OinMatrix)/StdsInWater)) %>%
    dplyr::select(`Protein Name`, `Precursor Ion Name`, BestMatch, RF)
} else {
  RF = QC.dat.IS.stds %>%
    filter(HILIC_Mix == MixID | MixID == "noMix") %>%
    group_by(MixID, solventID,`Protein Name`, `Precursor Ion Name`, BestMatch) %>%
    summarise(Area = mean(Area, na.rm = T)) %>%
    mutate(Area = ifelse(is.na(Area),0,Area),
           ID = ifelse(MixID == "noMix",MixID, solventID)) %>%
    ungroup() %>%
    dplyr::select(-MixID, -solventID) %>%
    spread(ID, Area) %>%
    mutate(RF = ((Matrix - noMix)/Water)) %>%
    dplyr::select(`Protein Name`, `Precursor Ion Name`, BestMatch, RF)
}

## apply to get estimated concentrations in the vial for the base peaks -------------------
Quant.bp = BMIS.dat %>%
  filter(type != "Std") %>%
  full_join(., RF, by = c("MassFeature" = "Precursor Ion Name")) %>%
  full_join(., IE, by = c("MassFeature" = "Precursor Ion Name",
                          "Protein Name", "BestMatch")) %>%
  mutate(Concentration_uM.in.vial = Adjusted_Area /IE * (1/RF)) %>%
  mutate(Experiment = ifelse(grepl("Fate1",Run.Cmpd),"Fate1",
                             ifelse(grepl("Fate2",Run.Cmpd),"Fate2",NA)),
         Timepoint = ifelse(grepl("T0-",Run.Cmpd), 0,
                                  ifelse(grepl("T4-",Run.Cmpd), 4,
                                         ifelse(grepl("T5-",Run.Cmpd), 5,
                                                ifelse(grepl("T8-",Run.Cmpd), 8,
                                                       ifelse(grepl("T12-",Run.Cmpd), 12,
                                                              ifelse(grepl("T24-",Run.Cmpd), 24,
                                                                     ifelse(grepl("T36-",Run.Cmpd),36,
                                                                            ifelse(grepl("T48-",Run.Cmpd), 48,
                                                                                   ifelse(grepl("T50-",Run.Cmpd), 50,
                                                                                          ifelse(grepl("Tlong-",Run.Cmpd), 96,NA))))))))))) %>%
  dplyr::select(MassFeature, `Protein Name`, Experiment, Timepoint, type, replicate, Concentration_uM.in.vial) %>%
  filter(!is.na(`Protein Name`),
         !is.na(Experiment))

## write RF and IE --------------------
(IE.RF <- full_join(IE, RF) %>%
   unique())

write_csv(IE.RF, path = paste0("FateTidy/output/TargetedMetabolites/IE_RF_",column.mode,".csv"))
write_csv(Quant.bp, path = paste0("FateTidy/output/TargetedMetabolites/Estimates_of_uMInVial_IS_samples_",column.mode,".csv"))
## use the base peak values to get the estimated concentrations in the vial for all isotopologues ----------
conversion.in.base.peak.to.uM <- left_join(Quant.bp, MID,
                                           by = c("Protein Name", "MassFeature" = "Precursor Ion Name",
                                                  "Experiment", "Timepoint", "type", 
                                                  "replicate" = "rep")) %>%
  filter(!is.na(Concentration_uM.in.vial),
         !is.na(all.area),
         C.lab == 0,
         N.lab == 0) %>%
  mutate(Conc.All.uM.in.vial = Concentration_uM.in.vial/MID) %>%
  dplyr::select(Experiment, Timepoint, replicate, type, `Protein Name`,
                all.area, Conc.All.uM.in.vial) %>%
  unique()

All.uM.in.vial <- left_join(conversion.in.base.peak.to.uM, MID,
                            by = c("Protein Name",
                                   "Experiment", "Timepoint", "type", 
                                   "replicate" = "rep", "all.area")) %>%
  mutate(Conc.uM.in.vial = MID * Conc.All.uM.in.vial)
  
## ditch estimates for the compounds that we have internal standards for -----------
IS.cmpds <- QC.dat.IS %>%
  filter(`Protein Name`=="Internal Standards"|
         `Protein Name`=="Internal Standards_pos"|
           `Protein Name`=="Internal standards_pos"|
         `Protein Name`=="Internal Standards_neg")
unique(IS.cmpds$`Precursor Ion Name`)

Int.stds <- c("L-Alanine","L-Histidine","L-Proline",
              "L-Valine", "L-Isoleucine",
              "L-Methionine", "Glycine betaine","Adenine",
              "Cytosine","Guanine","Arsenobetaine",
              "Isethionic acid","Sucrose",'Trehalose',
              "L-Phenylalanine","L-Tryptophan")

All.uM.in.vial.no.intstdcmpds <- All.uM.in.vial %>%
  filter(!(`Protein Name` %in% Int.stds))

## do the internal standard correction to get the concentration in the vial for those compounds ------
## conc.in.vial.cmpd =  int.std.con.in.vial * cmpd.area / int.std.area

Int.std.areas = QC.dat.IS %>%
  filter(type =="Smp"|type=="Blk",
         # !grepl("-pos",`Replicate Name`),
         `Precursor Ion Name` %in% 
           unique(IS.cmpds$`Precursor Ion Name`)) %>%
  dplyr::select(`Protein Name`,
                `Precursor Ion Name`,
                Area,
                Timepoint,
                Experiment,
                rep, SampID) %>%
  unique() %>%
  spread(`Precursor Ion Name`, Area)

stds.dat.IS = stds.dat.to.use %>%
  filter((BestMatch %in% 
           unique(IS.cmpds$`Precursor Ion Name`))|
           (Identification %in% 
              unique(IS.cmpds$`Precursor Ion Name`))|
           grepl("2H",BestMatch)|
           grepl("13C",BestMatch)) %>%
  mutate(Conc.uM.int.std.in.vial = Concentration_uM, 
         `Protein Name` =  ifelse(grepl("Alanine", BestMatch),"L-Alanine",
                                  ifelse(grepl("Histidine", BestMatch),"L-Histidine",
                                         ifelse(grepl("Proline", BestMatch), "L-Proline",
                                                ifelse(grepl("Valine",BestMatch),"L-Valine",
                                                       ifelse(grepl("Isoleucine",BestMatch),"L-Isoleucine",
                                                              ifelse(grepl("Methionine", BestMatch),"L-Methionine",
                                                                     ifelse(grepl("Arsenobetaine",BestMatch),"Arsenobetaine",
                                                                            ifelse(grepl("Adenine",BestMatch),"Adenine",
                                                                                   ifelse(grepl("Cytosine",BestMatch),"Cytosine",
                                                                                          ifelse(grepl("Guanine",BestMatch),"Guanine",
        ifelse(grepl("etaine",BestMatch),"Glycine betaine",
               ifelse(grepl("Trehalose",BestMatch),"Trehalose",
                      ifelse(grepl("Sucrose",BestMatch),"Sucrose",
                      ifelse(grepl("Isethionic",BestMatch),"Isethionic acid",
                             ifelse(grepl("Phenylalanine",BestMatch),"L-Phenylalanine",
                                    ifelse(grepl("Tryptophan",BestMatch),"L-Tryptophan",
                                           NA))))))))))))))))) %>%
  dplyr::select(BestMatch,`Protein Name`, Conc.uM.int.std.in.vial )
           
Int.std.cmpds.bp.conc = QC.dat.IS %>%
  filter(type =="Smp"|type=="Blk",
         N.lab == 0,
         C.lab == 0,
         `Protein Name` %in% Int.stds) %>%
  full_join(stds.dat.IS) %>%
  full_join(., Int.std.areas %>% dplyr::select(-`Protein Name`)) %>%
  filter(!is.na(Area)) %>%
  mutate(Conc.uM.in.vial = ifelse(grepl("Alanine", BestMatch), Area * Conc.uM.int.std.in.vial / `DL-Alanine, D3`,
                                  ifelse(grepl("Histidine", BestMatch), Area * Conc.uM.int.std.in.vial / `DL-Histidine, 15N`,
                                         ifelse(grepl("Proline", BestMatch),Area * Conc.uM.int.std.in.vial / `DL-Proline, D7`,
                                                ifelse(grepl("Valine",BestMatch),Area * Conc.uM.int.std.in.vial / `DL-Valine, D8`,
             ifelse(grepl("Isoleucine",BestMatch),Area * Conc.uM.int.std.in.vial / `L-Isoleucine, 15N`,
                      ifelse(grepl("Methionine", BestMatch),Area * Conc.uM.int.std.in.vial / `L-Methionine, D3`,
                          ifelse(grepl("Arsenobetaine",BestMatch),Area * Conc.uM.int.std.in.vial / `Arsenobetaine, 13C2`,
          ifelse(grepl("Adenine",BestMatch),Area * Conc.uM.int.std.in.vial / `Adenine, 15N2`,
               ifelse(grepl("Cytosine",BestMatch),Area * Conc.uM.int.std.in.vial / `Cytosine, 13C2-15N3`,
                      ifelse(grepl("Guanine",BestMatch),Area * Conc.uM.int.std.in.vial / `Guanine, 13C-15N2`,
                           ifelse(grepl("Trehalose",BestMatch), Area * Conc.uM.int.std.in.vial / `Trehalose, 13C`,
                                  ifelse(grepl("Isethionic",BestMatch), Area * Conc.uM.int.std.in.vial / `Isethionic Acid, 13C2`,
                                         ifelse(grepl("Phenylalanine",BestMatch), Area * Conc.uM.int.std.in.vial / `L-Phenylalanine, D8`,
                                                ifelse(grepl("Tryptophan",BestMatch), Area * Conc.uM.int.std.in.vial / `L-Tryptophan, D3`,
                                                       NA))))))))))))))) %>%
  filter(!is.na(Conc.uM.in.vial)) %>%
  dplyr::select(`Protein Name`, `N.lab`, `C.lab`,
                Experiment, Timepoint, rep, Conc.uM.in.vial, type) %>%
  unique()
## apply that to the various isotopologues ---------
conversion.in.base.peak.to.uM.intstdcmpds <- left_join(Int.std.cmpds.bp.conc, MID,
                                           by = c("Protein Name", "C.lab", "N.lab",
                                                  "Experiment", "Timepoint", "type", 
                                                 "rep")) %>%
  filter(!is.na(Conc.uM.in.vial),
         !is.na(all.area),
         C.lab == 0,
         N.lab == 0) %>%
  mutate(Conc.All.uM.in.vial = Conc.uM.in.vial/MID) %>%
  dplyr::select(Experiment, Timepoint, rep, type, `Protein Name`,
                all.area, Conc.All.uM.in.vial) %>%
  unique()

All.uM.in.vial.intstdcmpds <- left_join(conversion.in.base.peak.to.uM.intstdcmpds, MID,
                            by = c("Protein Name",
                                   "Experiment", "Timepoint", "type", 
                                    "rep", "all.area")) %>%
  mutate(Conc.uM.in.vial = MID * Conc.All.uM.in.vial) %>%
  unique()


## combine all that data together ------
All.uM.in.vial.to.use <- full_join(All.uM.in.vial.no.intstdcmpds, All.uM.in.vial.intstdcmpds,
                                   by = c("Experiment", "Timepoint", "type", "Replicate Name",
                                          "replicate" = "rep",
                                          "Protein Name",  "Precursor Ion Name", "all.area", "Conc.All.uM.in.vial", 
                                          "Retention Time", "Area", "Background", "Height", "Fraction1",
                                          "N.lab", "C.lab", "runDate", "SampID", "ionization", 
                                          "Blk.flag", "ppm.flag", "size.flag", "RT.check", "all.flags",
                                          "rawArea", "iso.mz", "RT..min.", "Column", "Blk.Max.Area", 
                                          "Area.less.blk", "Area.filled.w.background",
                                          "MID", "Conc.uM.in.vial")) %>%
  filter(`Protein Name` != "Glycine betaine") %>%
  unique()

## convert to estimated concentration in the actual samples ----------
#  from GBT quantification:  total.GBT.nM.in.sample = total.GBT.nM.in.vial * (37.078/30) *(400e-6/`Vol Filtered (L)`),
Sample.log <- read_csv("FateTidy/RawData/Sample_Log_GBTFateExperiments.csv") %>%
  filter(grepl("Metab",`Sample Type`) |grepl("metab",`Sample Type`))
Sample.log.to.join <- Sample.log  %>%
  filter(!grepl("Blk",`Sample ID`))%>%
  dplyr::select(Experiment, Timepoint,`Time of day`, rep, `Vol Filtered (L)`, `Incubation time (hr)`,
                `Time since last sample (hr)`) %>%
  unique()

All.nM.inSamples = full_join(All.uM.in.vial.to.use, Sample.log.to.join,
                             by = c("Experiment", "Timepoint",
                                    "replicate" = "rep")) %>%
  mutate(nM.in.sample = Conc.uM.in.vial * 1000 * (37.078/30) *(400e-6/`Vol Filtered (L)`)) %>%  ## 30 ratio because we spiked in water or standards after aliquoting 30 um. 400e-6 because that was the original reconsitution volume
  unique()



## add gbt data ----
if(column.mode=="HILICPos"){
  GBT.quant <- read_csv("FateTidy/output/TargetedMetabolites/GBT_concentrations.csv") %>%
    dplyr::rename(nM.in.sample = GBT.nM) %>%
    unique()
  
  All.nM.inSamples <- full_join(All.nM.inSamples, GBT.quant,
                              by = c("Experiment", "Timepoint", "type", "Protein Name",
                                     "replicate" = "rep",
                                     "all.area", "Replicate Name", "Precursor Ion Name", 
                                     "Retention Time", "Area", "Background", "Height", "Fraction1", "N.lab", "C.lab", 
                                     "runDate", "SampID", "ionization", "Blk.flag", "ppm.flag", "size.flag", "RT.check",
                                     "all.flags", "rawArea", "iso.mz", "RT..min.", "Column", "Blk.Max.Area", "Area.less.blk", 
                                     "Area.filled.w.background", "MID", "Time of day",
                                     "Vol Filtered (L)", "Incubation time (hr)", "Time since last sample (hr)",
                                     'nM.in.sample'))
} 

tally.up = All.nM.inSamples %>%
  group_by(Experiment, `Precursor Ion Name`) %>%
  summarise(n = n())
## save out data ----------
write_csv(All.nM.inSamples, path = paste0("FateTidy/output/TargetedMetabolites/nM.concentrations.",column.mode,".csv"))


## take vial concentration to estimated in situ concentration for the base peak estimates from the -is samples (this work is just mostly for glycine) ----------
IS.BP.nM.inSamples = full_join(Quant.bp, Sample.log.to.join,
                             by = c("Experiment", "Timepoint",
                                    "replicate" = "rep")) %>%
  mutate(nM.in.sample = Concentration_uM.in.vial * 1000 * (37.078/30) *(400e-6/`Vol Filtered (L)`)) %>%  ## 30 ratio because we spiked in water or standards after aliquoting 30 um. 400e-6 because that was the original reconsitution volume
  unique()

write_csv(IS.BP.nM.inSamples, path = paste0("FateTidy/output/TargetedMetabolites/Estimates_of_nMInSample_IS_samples_",column.mode,".csv"))

