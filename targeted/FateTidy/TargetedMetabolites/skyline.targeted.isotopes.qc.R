## Processing Targeted Output
library(tidyverse)
library(RCurl)
## step 0: set some parameters:---------------
# column.mode = "HILICPos"
# column.mode = "Lipid"
column.mode = "HILICNeg"
# column.mode = "CyanoAq"

if(column.mode=="CyanoAq"){
  RT.flex = 0.5
  RT.flex.iso = 0.25
  min.value = 1e3
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 10
} else if (column.mode=="HILICPos"){
  RT.flex = 1
  RT.flex.iso = 0.3
  min.value = 1e3
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 10
}else if (column.mode=="HILICNeg"){
  RT.flex = 1.5
  RT.flex.iso = 0.3
  min.value = 1e3
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 10
} else if (column.mode=="Lipid"){
  RT.flex = 0.5
  RT.flex.iso = 0.5
  min.value = 1e4
  This.much.bigger.than.max.Blank = 2
  ppm.thresh = 15
}

## Step 1: read in data -----
if(column.mode == "Lipid"){
  iso.poss <- read_csv("FateTidy/All_Chl_C-N_Isotope_possibilities.csv") %>%
    mutate(HILICMix = NA)
} else {
  iso.poss <- read_csv("FateTidy/Isotope_Possibilities_IngallsStandards_trim.csv",
                       col_types = cols(
                         .default = col_character(),
                         iso.mz = col_double(),
                         RT..min. = col_double(),
                         z = col_double(),
                         m.z = col_double(),
                         C = col_double(),
                         this.N = col_double(),
                         possibilities = col_double(),
                         X = col_double(),
                         Fraction2 = col_character(),
                         Conc..uM = col_double(),
                         Date.added = col_double()
                       )) 
  iso.poss <- iso.poss %>%
    dplyr::select(Compound.and.Iso.name, iso.mz, RT..min., Column, Fraction1, HILICMix) 
}
if(column.mode == "HILICPos"){
  dat.raw = read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_pos_QE_Transition Results.csv",
                                na = "#N/A") %>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(!(`Protein Name`%in% c("L-Methionine","Carnitine","O-Acetylcarnitine")))
  
  ## select the few extra isotopologues integrated on 8/28/2020
  more.data <- read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_pos_QE_Transition Results_justAFewBonusIsotopes.csv",
                        na = "#N/A") %>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(`Precursor Ion Name` == "Glycerophosphocholine_13C-6 15N-0"|
             `Precursor Ion Name` == "Glycerophosphocholine_13C-7 15N-0"|
             `Precursor Ion Name` == "Glycerophosphocholine_13C-8 15N-0"|
             `Precursor Ion Name` == "Glycerophosphocholine_13C-6 15N-1"|
             `Precursor Ion Name` == "Glycerophosphocholine_13C-7 15N-1"|
             `Precursor Ion Name` == "(3-Carboxypropyl)trimethylammonium_13C-6 15N-1"
    )
  
  
  dat.add.me <- read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_pos_QE_Transition Results_methionine_carnitine_acetylcarnitine.csv",
                         na = "#N/A")%>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(`Protein Name`%in% c("L-Methionine","Carnitine","O-Acetylcarnitine"),
           grepl("-pos",`Replicate Name`) | grepl("Poo", `Replicate Name`) | grepl("Std", `Replicate Name`))
  
  dat.add.me.too <- read_csv("FateTidy/RawData/Targeted HILIC/HILICPos_GBTFate_pos_QE_Transition Results_some_more_Osmos.csv",
                             na = "#N/A")%>%
    mutate(Fraction1 = "HILICPos") %>%
    filter(`Protein Name`%in% c("DMSP","Hydroxylysine","Gonyol","Glucosylglycerol"),
           grepl("-pos",`Replicate Name`) | grepl("Poo", `Replicate Name`) | grepl("Std", `Replicate Name`))
  
  dat.raw <- full_join(dat.raw, more.data) %>%
    full_join(., dat.add.me) %>%
    full_join(.,dat.add.me.too)
  
  joined.dat <- dat.raw %>%
    left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name","Fraction1")) %>%
    unique()%>%
    mutate(Area = as.numeric(Area),
           name.to.break = `Replicate Name`,
           labels.to.break = `Precursor Ion Name`,
           N.lab = ifelse(grepl("15N",`Precursor Ion Name`),
                          str_sub(`Precursor Ion Name`, start = -1, -1),0),
           C.lab =  ifelse(grepl("15N",`Precursor Ion Name`),
                           str_sub(`Precursor Ion Name`, start = -7, -7),
                           str_sub(`Precursor Ion Name`, start = -1, -1))) %>%
    # dplyr::select(-Other) %>%
    separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
    mutate(Timepoint = str_sub(SampID, start = 11, -1),
           Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
           Timepoint = as.numeric(Timepoint),
           Experiment = str_sub(SampID, start = 4,8))
  
} else if (column.mode == "HILICNeg"){
  dat.raw = read_csv("FateTidy/RawData/Targeted HILIC/HILICNeg_GBTFate_pos_only_a_few_OSMOS_NOV02.csv",
                     na = "#N/A") %>%
    mutate(Fraction1 = "HILICNeg")
  joined.dat <- dat.raw %>%
    left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name","Fraction1")) %>%
    unique()%>%
    mutate(Area = as.numeric(Area),
           name.to.break = `Replicate Name`,
           labels.to.break = `Precursor Ion Name`,
           N.lab = ifelse(grepl("15N",`Precursor Ion Name`),
                          str_sub(`Precursor Ion Name`, start = -1, -1),0),
           C.lab =  ifelse(grepl("15N",`Precursor Ion Name`),
                           str_sub(`Precursor Ion Name`, start = -7, -7),
                           str_sub(`Precursor Ion Name`, start = -1, -1))) %>%
    # dplyr::select(-Other) %>%
    separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
    mutate(Timepoint = str_sub(SampID, start = 11, -1),
           Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
           Timepoint = as.numeric(Timepoint),
           Experiment = str_sub(SampID, start = 4,8))
  
} else if (column.mode=="CyanoAq") {
  dat.raw = read_csv("FateTidy/RawData/Targeted CyanoAq/CyanoAq_GBTFate_pos_isopologues_Nov03.csv",
                     na = "#N/A") %>%
    filter(grepl("-pos",`Replicate Name`)) %>%
    mutate(Fraction1 = "CyanoAq")
  poo.stds = read_csv("FateTidy/RawData/Targeted CyanoAq/CyanoAq_GBTFate_pos_and_IS_isopologues_Nov02.csv",
                     na = "#N/A") %>%
    filter(grepl("Poo",`Replicate Name`)|
             grepl("Std",`Replicate Name`)) %>%
    mutate(Fraction1 = "CyanoAq")
  dat.raw <- full_join(dat.raw, poo.stds)
  joined.dat <- dat.raw %>%
    left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name","Fraction1")) %>%
    unique()%>%
    mutate(Area = as.numeric(Area),
           name.to.break = `Replicate Name`,
           labels.to.break = `Precursor Ion Name`,
           N.lab = ifelse(grepl("15N",`Precursor Ion Name`),
                          str_sub(`Precursor Ion Name`, start = -1, -1),0),
           C.lab =  ifelse(grepl("15N",`Precursor Ion Name`),
                           str_sub(`Precursor Ion Name`, start = -7, -7),
                           str_sub(`Precursor Ion Name`, start = -1, -1))) %>%
    # dplyr::select(-Other) %>%
    separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
    mutate(Timepoint = str_sub(SampID, start = 11, -1),
           Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
           Timepoint = as.numeric(Timepoint),
           Experiment = str_sub(SampID, start = 4,8))
  
} else if (column.mode=="Lipid") {
  filename  = "FateTidy/RawData/Targeted Lipid/Lipid_Skyline_Chl_New_RT_Isotopologues_QE_Transition Results.csv"
  # filename  = "FateTidy/RawData/Targeted Lipid/Lipid_Skyline_Chl_Isotopologues_QE_Transition Results.csv"
  
  dat.raw = read_csv(filename,
                     na = "#N/A") %>%
    mutate(Fraction1 = "Lipid") %>%
  filter(grepl("2006",`Replicate Name`),
         `Protein Name` =="Chlorophyll") %>%
    mutate(name.to.break = ifelse(`Replicate Name` =="200617_Smp_GBTFate1MT0-pos_A_20200624191823",
                                     "200617_Smp_GBTFate1MT0-pos_A",
                                     ifelse(grepl("TruePooGBT_Full-pos",`Replicate Name`),
                                            str_replace(`Replicate Name`, "TruePooGBT_Full-pos","TruePooGBT_pos-Full"),
                                            ifelse(grepl("TruePooGBT_Half-pos",`Replicate Name`),
                                                   str_replace(`Replicate Name`, "TruePooGBT_Half-pos","TruePooGBT_pos-Half"),
                                                   ifelse(grepl("TruePooGBT-Full_DDApos", `Replicate Name`),
                                                          str_replace(`Replicate Name`, "TruePooGBT-Full_DDApos", "TruePooGBT-pos-FullDDA"),
                                     ifelse(grepl("pos1",`Replicate Name`),
                                            str_replace(`Replicate Name`,"pos1","pos-1"),
                                            ifelse(grepl("pos2",`Replicate Name`),
                                                   str_replace(`Replicate Name`,"pos2","pos-2"),
                                                   ifelse(grepl("pos3",`Replicate Name`),
                                                          str_replace(`Replicate Name`,"pos3","pos-3"),
                                                          ifelse(grepl("posNew",`Replicate Name`),
                                                                 str_replace(`Replicate Name`,"posNew","pos-New"),
                                            `Replicate Name`)))))))))
  
  joined.dat <- dat.raw %>%
    left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name")) %>%
    unique()%>%
    mutate(Area = as.numeric(Area),
           N.lab = str_sub(`Precursor Ion Name`, start = -1, -1),
           C.lab = str_sub(`Precursor Ion Name`, start = -7, -7),
           labels.to.break = `Precursor Ion Name`) %>%
    separate(labels.to.break, into = c("Other","C.lab","N.lab"), sep = "_") %>%
    separate(C.lab, into = c("C.lab","N.lab"), sep = " ") %>%
    separate(C.lab, into = c("Other","C.lab"), sep ="C-")%>%
    separate(N.lab, into = c("Other","N.lab"), sep ="N-") %>%
    dplyr::select(-Other) %>%
    separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
    mutate(Timepoint = str_sub(SampID, start = 11, -1),
           Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
           Timepoint = as.numeric(Timepoint),
           Experiment = str_sub(SampID, start = 4,8),
           # RT..min. = 22.7) ## old RT - might be right
           RT..min. = 25.5) %>% ## newer RT - might be right
  filter(
    C.lab < 5, ## these filters are for data generated 9/16/2020 when I just integrated a subset with the new RT
    N.lab < 2 ## these filters are for data generated 9/16/2020 when I just integrated a subset with the new RT
  )
}


dat.no.IS <- joined.dat %>%
  filter(`Protein Name` != "Internal Standards",
         `Protein Name` != "Internal Standards_neg",
         `Protein Name` != "Internal standards_pos")%>%
  mutate(N.lab = as.numeric(N.lab),
         C.lab = as.numeric(C.lab))
## Step 2: Quality control by using the blanks--------
## get blank data for base peak
blank.dat <- dat.no.IS %>%
  filter(grepl("Blk", `Replicate Name`)) 
## get max value
Choline.Arginine.blank.max <- blank.dat %>%
  filter(`Protein Name`=="Choline" | `Protein Name`=="L-Arginine") %>%
  filter(abs(`Retention Time` - RT..min.) < 5) %>%
  group_by(`Protein Name`, `Precursor Ion Name`, 
           iso.mz, RT..min., Column, N.lab, C.lab) %>%
  summarise(Blk.Max.Area = max(Area, na.rm = T)) %>%
  mutate(Blk.Max.Area = ifelse(is.infinite(Blk.Max.Area),0,Blk.Max.Area))
  
max.blank <- blank.dat %>%
  filter(`Protein Name`!="Choline" ,
         `Protein Name`!="L-Arginine") %>%
  filter(abs(`Retention Time` - RT..min.) < RT.flex) %>%
  group_by(`Protein Name`, `Precursor Ion Name`, 
           iso.mz, RT..min., Column, N.lab, C.lab) %>%
  summarise(Blk.Max.Area = max(Area, na.rm = T)) %>%
  mutate(Blk.Max.Area = ifelse(is.infinite(Blk.Max.Area),0,Blk.Max.Area)) %>%
  full_join(Choline.Arginine.blank.max)
## compare base peak to blank
blank.check <- full_join(dat.no.IS, max.blank) %>%
  mutate(Blk.ratio = Area/Blk.Max.Area,
         Blk.flag = ifelse(Blk.ratio < This.much.bigger.than.max.Blank, "comparable to blank","ok"))

## compare ppm to our defined threshold
ppm.check <- blank.check %>%
  mutate(ppm.flag = ifelse(abs(`Mass Error PPM`) > ppm.thresh, "bad ppm", "ok"))

## size.check
size.check <- ppm.check %>%
  mutate(size.flag = ifelse(Area > min.value, "ok", "too small"))


## compare the RT of the base peak to the expected (from standards)
std.rt <- dat.no.IS %>%
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
## compare the RT of the isotopes to that of the base peak
base.peak.RT <- rt.check.1 %>%
  ungroup() %>%
  filter(N.lab == 0,
         C.lab == 0) %>%
  select(`Replicate Name`, `Protein Name`, `Retention Time`, RT.check) %>%
  unique() %>%
  rename(Base.Peak.RT = `Retention Time`,
         Base.Peak.ok = RT.check)
rt.check.2 <- full_join(rt.check.1, base.peak.RT) %>%
  dplyr::select(-maxRT, -minRT) %>%
  full_join(., std.rt %>%
              ungroup() %>%
              dplyr::select( -N.lab, -`Precursor Ion Name`, -iso.mz)) %>%
  mutate(RT.check = ifelse(!is.na(RT.check),RT.check,
                           ifelse(is.na(`Retention Time`), "bad RT",
                           ifelse(Base.Peak.ok == "bad RT" &
                                  (`Retention Time` > minRT - RT.flex) & (`Retention Time` < maxRT + RT.flex),
                           "ok",
                           ifelse((`Retention Time` > Base.Peak.RT - RT.flex.iso) & (`Retention Time` < Base.Peak.RT + RT.flex.iso),
                           "ok","bad RT")))))
## combine all flags 
qc.dat <- rt.check.2 %>%
  mutate(all.flags = ifelse(ppm.flag=="ok" & RT.check =="ok" & Blk.flag=="ok" &
                              size.flag=="ok",
                            "all ok","something is wrong")) %>%
  dplyr::select(`Protein Name`, `Replicate Name`, `Precursor Ion Name`,
                `Retention Time`, `Area`, `Background`, `Height`,
                `Fraction1`,`N.lab`,`C.lab`,`runDate`,`type`,
                `SampID`, `ionization`,`rep`, `Timepoint`, `Experiment`,`Blk.flag`,
                `ppm.flag`,size.flag, RT.check, all.flags) %>%
  mutate(rawArea = Area,
         Area = ifelse(all.flags == "all ok",Area,NA))

## combine with Int std dat 
qc.dat.full <- qc.dat %>%
  full_join(., joined.dat %>%
              filter(`Protein Name` == "Internal Standards"|
                       `Protein Name`=="Internal Standards_neg")%>%
              mutate(N.lab = as.numeric(N.lab),
                     C.lab = as.numeric(C.lab)) %>%
              dplyr::select(`Protein Name`, `Replicate Name`, `Precursor Ion Name`,
                            `Retention Time`, `Area`, `Background`, `Height`,
                            `Fraction1`,`N.lab`,`C.lab`,`runDate`,`type`,
                            `SampID`, `ionization`,`rep`, `Timepoint`, `Experiment`))

## ditch compounds I know are bad
# (13C3-15N1 is not Alanine, wrong RT)
qc.dat.full.trim <- qc.dat.full %>%
  filter(`Precursor Ion Name` !="L-Alanine_13C-3 15N-1",
         `Precursor Ion Name` !="Proline betaine_13C-3 15N-0")

## Step 3: save qc dat------------
write_csv(qc.dat.full.trim, path = paste0("FateTidy/output/TargetedMetabolites/",column.mode,"_QC_skyline_isotope_data.csv"))
# write_csv(qc.dat.full.trim, path = paste0("Fate/output/skyline/",column.mode,"_newRT_Chl_QC_skyline_isotope_data.csv")) ## just for the one with new RT chl

## get number of times a compound is ok 
cmpd.freq.tab <- qc.dat.full.trim %>%
  filter(type == "Smp",
         !is.na(Area)) %>%
  ungroup() %>%
  group_by(`Protein Name`, C.lab, N.lab, `Precursor Ion Name`, Experiment) %>%
  summarise(n = n())

write_csv(cmpd.freq.tab, path = paste0("FateTidy/output/TargetedMetabolites/",column.mode,"_frequency_of_isotope_data.csv"))
# write_csv(cmpd.freq.tab, path = paste0("Fate/output/skyline/",column.mode,"_newRT_Chl_frequency_of_isotope_data.csv")) ## just for the one with new RT chl


## step 4: do blank subtraction --------------------
blk.sub.dat <- qc.dat.full.trim %>%
  full_join(max.blank)  %>%
  mutate(Area.less.blk = Area - Blk.Max.Area)

write_csv(blk.sub.dat, path = paste0("FateTidy/output/TargetedMetabolites/",column.mode,"_QC_skyline_isotope_data_blankSubtract.csv"))
# write_csv(blk.sub.dat, path = paste0("Fate/output/skyline/",column.mode,"_newRT_Chl_QC_skyline_isotope_data_blankSubtract.csv"))## just for the one with new RT chl

