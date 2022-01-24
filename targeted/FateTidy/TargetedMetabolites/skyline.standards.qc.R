## QC for the -IS samples from GBT Fate QE runs
library(tidyverse)
library(RCurl)
library (readr)

## step 0: set some parameters:---------------
# column.mode = "HILICPos"
# column.mode = "HILICNeg"
column.mode = "CyanoAq"

if(column.mode=="CyanoAq"){
  RT.flex = 0.5
  RT.flex.iso = 0.25
  min.value = 1e3
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 10
}else if (column.mode=="HILICNeg"){
  RT.flex = 1.5
  RT.flex.iso = 0.3
  min.value = 1e3
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 10
} else {
  RT.flex = 1
  RT.flex.iso = 0.3
  min.value = 1e3
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 10
}

## Step 1: read in data -----

iso.poss <- read_csv("FateTidy/Isotope_Possibilities_IngallsStandards_trim.csv")
iso.poss <- iso.poss %>%
  dplyr::select(Compound.and.Iso.name, iso.mz, RT..min., Column, Fraction1, HILICMix) 

## get standards
if(column.mode == "HILICPos"){
  dat.raw = read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_IS_QE_Transition Results.csv",
                     na = "#N/A") %>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(!(`Protein Name`%in% c("L-Methionine","Carnitine","O-Acetylcarnitine")))
  
  dat.add.me <- read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_pos_QE_Transition Results_methionine_carnitine_acetylcarnitine.csv",
                         na = "#N/A")%>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(`Protein Name`%in% c("L-Methionine","Carnitine","O-Acetylcarnitine"),
           grepl("-IS",`Replicate Name`) | grepl("Poo", `Replicate Name`) | grepl("Std", `Replicate Name`))
  dat.add.me.too <- read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_IS_QE_Transition Results_some_more_Osmos.csv",
                             na = "#N/A")%>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(`Protein Name`%in% c("DMSP","Hydroxylysine","Gonyol","Glucosylglycerol"),
           grepl("-IS",`Replicate Name`) | grepl("Poo", `Replicate Name`) | grepl("Std", `Replicate Name`))
  
  add.me.three <- read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_pos_QE_4AminoAcids_Nov19.csv",
                           na = "#N/A") %>%
    mutate(Fraction1 = "HILICPos",
           `Protein` = `Precursor Ion Name`,
           `Protein Name` = `Precursor Ion Name`,
          `Precursor Ion Name` = paste0(`Precursor Ion Name`,"_13C-0 15N-0"))
  
  dat.raw <- full_join(dat.raw, dat.add.me) %>%
    full_join(.,dat.add.me.too)%>%
    full_join(.,add.me.three)
} else if (column.mode == "HILICNeg"){
  dat.raw = read_csv("FateTidy/RawData/Targeted HILIC/HILICNeg_GBTFate_IS_only_a_few_OSMOS_NOV02.csv",
                     na = "#N/A") %>%
    mutate(Fraction1 = "HILICNeg")
  stds.dat = read_csv("FateTidy/RawData/Targeted HILIC/HILICNeg_GBTFate_pos_only_a_few_OSMOS_NOV02.csv",
                      na = "#N/A") %>%
    filter(grepl("Std",`Replicate Name`)) %>%
    mutate(Fraction1 = "HILICNeg")
  
  dat.raw <- full_join(dat.raw, stds.dat)
} else if (column.mode=="CyanoAq") {
  dat.raw = read_csv("FateTidy/RawData/Targeted CyanoAq/CyanoAq_GBTFate_pos_and_IS_isopologues_Nov02.csv",
                     na = "#N/A") %>%
    mutate(Fraction1 = "CyanoAq")
}


joined.dat <- dat.raw %>%
  left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name", "Fraction1")) %>%
  unique()%>%
  mutate(Area = as.numeric(Area),
         N.lab =  ifelse(grepl("15N",`Precursor Ion Name`),
                         str_sub(`Precursor Ion Name`, start = -1, -1),
                         0),
         C.lab =  ifelse(grepl("15N",`Precursor Ion Name`),
                         str_sub(`Precursor Ion Name`, start = -7, -7),
                         str_sub(`Precursor Ion Name`, start = -1, -1)),
         name.to.break = `Replicate Name`) %>%
  separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
  mutate(Timepoint = str_sub(SampID, start = 11, -1),
         Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
         Timepoint = as.numeric(Timepoint),
         Experiment = str_sub(SampID, start = 4,8))

dat.no.IS <- joined.dat %>%
  filter(`Protein Name` != "Internal Standards",
         `Protein Name` != "Internal Standards_neg",
         `Protein Name` !="Internal standards_pos")%>%
  mutate(N.lab = as.numeric(N.lab),
         C.lab = as.numeric(C.lab))

dat.no.is.base.peak.only <- dat.no.IS %>%
  filter(C.lab == 0,
         N.lab ==0)

## get internal standards
IS <- dat.raw %>%
  filter(`Protein Name` == "Internal Standards"|
         `Protein Name` == "Internal Standards_neg"|
           `Protein Name` == "Internal standards_pos")%>%
  mutate(Area = as.numeric(Area),
         # N.lab = str_sub(`Precursor Ion Name`, start = -1, -1),
         # C.lab = str_sub(`Precursor Ion Name`, start = -7, -7),
         name.to.break = `Replicate Name`) %>%
  separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
  mutate(Timepoint = str_sub(SampID, start = 11, -1),
         Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
         Timepoint = as.numeric(Timepoint),
         Experiment = str_sub(SampID, start = 4,8))

## get master list with compound mix info -----
# AllStds_LabInfo <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"), header = T) 
urlfile="https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"

AllStds_LabInfo<-read_csv(url(urlfile))

## get RT range for comopunds (in Mix 1 and 2)------
## Step 2: Quality control by using the blanks--------
## get blank data for base peak
blank.dat <- dat.no.is.base.peak.only %>%
  filter(grepl("Blk", `Replicate Name`)) 
## get max value
max.blank <- blank.dat %>%
  group_by(`Protein Name`, `Precursor Ion Name`, 
           iso.mz, RT..min., Column, N.lab, C.lab) %>%
  summarise(Blk.Max.Area = max(Area, na.rm = T)) %>%
  mutate(Blk.Max.Area = ifelse(is.infinite(Blk.Max.Area),0,Blk.Max.Area))
## compare base peak to blank
blank.check <- full_join(dat.no.is.base.peak.only, max.blank) %>%
  mutate(Blk.ratio = Area/Blk.Max.Area,
         Blk.flag = ifelse(Blk.ratio < This.much.bigger.than.max.Blank, "comparable to blank","ok"))

## compare ppm to our defined threshold
ppm.check <- blank.check %>%
  mutate(ppm.flag = ifelse(abs(`Mass Error PPM`) > ppm.thresh, "bad ppm", "ok"))

## size.check
size.check <- ppm.check %>%
  mutate(size.flag = ifelse(Area > min.value, "ok", "too small"))

## compare the RT of the base peak to the expected (from standards)
std.rt <- dat.no.is.base.peak.only %>%
  mutate(correct.mix = ifelse(is.na(HILICMix),TRUE,
                              ifelse(HILICMix == "Mix1", grepl("Mix1",SampID),
                                     ifelse(HILICMix == "Mix2",  grepl("Mix2",SampID), NA)))) %>%
  filter(grepl("Std",`Replicate Name`),
         SampID != "H2OinMatrix",
         N.lab == 0,
         C.lab == 0,
         correct.mix == TRUE) %>%
  ungroup() %>%
  group_by(`Protein Name`, `Precursor Ion Name`, iso.mz, RT..min., Column, N.lab, C.lab) %>%
  summarise(minRT = min(`Retention Time`),
            maxRT = max(`Retention Time`)) %>%
  gather(RT.option, Value, -`Protein Name`, -`Precursor Ion Name`, -iso.mz, -RT..min., -Column, -N.lab, -C.lab) %>%
  summarise(minRT = min(Value),
            maxRT = max(Value))
rt.check.1 <- full_join(size.check, std.rt) %>%
  mutate(RT.check = ifelse((`Retention Time` > minRT - RT.flex) & (`Retention Time` < maxRT + RT.flex),
                           "ok","bad RT"))
## combine all flags 
qc.dat <- rt.check.1 %>%
  mutate(all.flags = ifelse(ppm.flag=="ok" & RT.check =="ok" & Blk.flag=="ok" &
                              size.flag=="ok",
                            "all ok","something is wrong")) %>%
  dplyr::select(`Protein Name`, `Replicate Name`, `Precursor Ion Name`,
                `Retention Time`, `Area`, `Background`, `Height`,
                `Fraction1`,`N.lab`,`C.lab`,`runDate`,`type`,
                `SampID`, `ionization`,`rep`, `Timepoint`, `Experiment`,`Blk.flag`,
                `ppm.flag`, size.flag, RT.check, all.flags) %>%
  mutate(rawArea = Area,
         Area = ifelse(all.flags == "all ok",Area,NA))

## combine with Int std dat 
qc.dat.full <- qc.dat %>%
  full_join(., IS  %>%
              dplyr::select(`Protein Name`, `Replicate Name`, `Precursor Ion Name`,
                            `Retention Time`, `Area`, `Background`, `Height`,
                            `Fraction1`,`runDate`,`type`,
                            `SampID`, `ionization`,`rep`, `Timepoint`, `Experiment`))

## write out QC data------
write_csv(qc.dat.full, path = paste0("FateTidy/output/TargetedMetabolites/",column.mode,"_QC_skyline_IS_data.csv"))

## get number of times a compound is ok 
cmpd.freq.tab <- qc.dat.full %>%
  filter(type == "Smp",
         !is.na(Area)) %>%
  ungroup() %>%
  group_by(`Protein Name`, `Precursor Ion Name`, Experiment) %>%
  summarise(n = n())

write_csv(cmpd.freq.tab, path = paste0("FateTidy/output/TargetedMetabolites/",column.mode,"_frequency_of_standards_IS_data.csv"))
