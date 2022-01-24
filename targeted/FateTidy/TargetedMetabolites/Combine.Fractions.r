# Combine Fractions

library(tidyverse)

##concentrations: --------------
Cyano.nM = read_csv("FateTidy/output/TargetedMetabolites/nM.concentrations.CyanoAq.csv")
HILICNeg.nM = read_csv("FateTidy/output/TargetedMetabolites/nM.concentrations.HILICNeg.csv")
HILICPos.nM = read_csv("FateTidy/output/TargetedMetabolites/nM.concentrations.HILICPos.csv",
                       col_types = cols(
                         .default = col_double(),
                         Experiment = col_character(),
                         replicate = col_character(),
                         type = col_character(),
                         `Protein Name` = col_character(),
                         `Replicate Name` = col_character(),
                         `Precursor Ion Name` = col_character(),
                         Fraction1 = col_character(),
                         SampID = col_character(),
                         ionization = col_character(),
                         Blk.flag = col_character(),
                         ppm.flag = col_character(),
                         size.flag = col_character(),
                         RT.check = col_character(),
                         all.flags = col_character(),
                         Column = col_character(),
                         total.GBT.nM.in.sample = col_double(),
                         total.GBT.nM.in.sample.sd = col_double(),
                         GBT.nM.sd = col_double()
                       ))

all.nM = full_join(Cyano.nM, HILICNeg.nM) %>%
  full_join(HILICPos.nM)
write_csv(all.nM, "FateTidy/output/TargetedMetabolites/nM.concentrations.allFractions.csv")

## IERF------------
Cyano.IERF = read_csv("FateTidy/output/TargetedMetabolites//IE_RF_CyanoAq.csv")
HILICNeg.IERF = read_csv("FateTidy/output/TargetedMetabolites//IE_RF_HILICNeg.csv")
HILICPos.IERF = read_csv("FateTidy/output/TargetedMetabolites//IE_RF_HILICPos.csv")
all.IERF = full_join(Cyano.IERF, HILICNeg.IERF) %>%
  full_join(HILICPos.IERF)
write_csv(all.IERF, "FateTidy/output/TargetedMetabolites//IE_RF_allFractions.csv")


## MID mean ------
Cyano.MID = read_csv("FateTidy/output/TargetedMetabolites//MID_mean_QC_Targeted_blkSubtractCyanoAq.csv")
HILICNeg.MID = read_csv("FateTidy/output/TargetedMetabolites//MID_mean_QC_Targeted_blkSubtractHILICNeg.csv")
HILICPos.MID = read_csv("FateTidy/output/TargetedMetabolites//MID_mean_QC_Targeted_blkSubtractHILICPos.csv")
all.MID = full_join(Cyano.MID, HILICNeg.MID) %>%
  full_join(HILICPos.MID)
write_csv(all.MID, "FateTidy/output/TargetedMetabolites//MID_mean_blkSubtract_allFractions.csv")

## MID all ------
Cyano.MID = read_csv("FateTidy/output/TargetedMetabolites//MID_QC_Targeted_blkSubtractCyanoAq.csv")
HILICNeg.MID = read_csv("FateTidy/output/TargetedMetabolites//MID_QC_Targeted_blkSubtractHILICNeg.csv")
HILICPos.MID = read_csv("FateTidy/output/TargetedMetabolites//MID_QC_Targeted_blkSubtractHILICPos.csv")
all.MID = full_join(Cyano.MID, HILICNeg.MID) %>%
  full_join(HILICPos.MID)
write_csv(all.MID, "FateTidy/output/TargetedMetabolites//MID_blkSubtract_allFractions.csv")


## MID sig ----
Cyano.MID.sig <- read_csv("FateTidy/output/TargetedMetabolites//CyanoAq_goodCompounds_MID_lm_significance_test.csv")
HILICPos.MID.sig <- read_csv("FateTidy/output/TargetedMetabolites//HILICPos_goodCompounds_MID_lm_significance_test.csv")
HILICNeg.MID.sig <- read_csv("FateTidy/output/TargetedMetabolites//HILICNeg_goodCompounds_MID_lm_significance_test.csv")

all.MID.sig <- full_join(Cyano.MID.sig, HILICNeg.MID.sig ) %>%
  full_join(HILICPos.MID.sig)
library(stats)
all.MID.sig.to.write <- all.MID.sig %>%
  dplyr::select(`Precursor Ion Name`, Experiment, p.value, sig) %>%
  dplyr::rename(Isotopologue =`Precursor Ion Name` ) %>%
  unique() %>%
  arrange(Isotopologue, Experiment)# %>%
  # mutate(p.fdr = p.adjust(p = p.value, method = "fdr"))

write_csv(all.MID.sig.to.write,"FateTidy/output/TargetedMetabolites//SuppTab_MID_significant_trend.csv" )

## MID all means, all nM averaged, gly, asp, asn nM ------------
all.nM.ave <- all.nM %>%
  group_by(Experiment, Timepoint, type, `Protein Name`, `Precursor Ion Name`, N.lab, C.lab, ionization, Column,
           `Time of day`, `Incubation time (hr)`) %>%
  summarise(nM.in.Sample.ave = mean(nM.in.sample, na.rm = T),
            nM.in.sample.sd = sd(nM.in.sample, na.rm = T)) %>%
  filter(ionization!="IS")

aminoaciddata <- read_csv("FateTidy/output/TargetedMetabolites//Estimates_of_nMInSample_IS_samples_HILICPos.csv")
cmpds <- c("Glycine","Asparagine","Aspartic acid")
aminoaciddata <- aminoaciddata %>%
  filter(`Protein Name` %in% cmpds)%>%
  mutate(N.lab =0, 
         C.lab =0,
         Column = "HILIC",
         ionization = "pos") %>%
  dplyr::rename(`Precursor Ion Name`=MassFeature) %>%
  group_by(Experiment, Timepoint, type, `Protein Name`, `Precursor Ion Name`, N.lab, C.lab, ionization, Column,
           `Time of day`, `Incubation time (hr)`) %>%
  summarise(nM.in.Sample.ave = mean(nM.in.sample, na.rm = T),
            nM.in.sample.sd = sd(nM.in.sample, na.rm = T)) %>%
  filter(type !="Blk")


mean.MID.all <- read_csv("FateTidy/output/TargetedMetabolites//MID_mean_blkSubtract_allFractions.csv")%>%
  filter(ionization!="IS")

freq.dat.ca <- read_csv("FateTidy/output/TargetedMetabolites//CyanoAq_frequency_of_isotope_data.csv")
freq.dat.hn <- read_csv("FateTidy/output/TargetedMetabolites//HILICNeg_frequency_of_isotope_data.csv")
freq.dat.hp <- read_csv("FateTidy/output/TargetedMetabolites//HILICPos_frequency_of_isotope_data.csv")
all.freq <- full_join(freq.dat.ca, freq.dat.hn) %>% full_join(freq.dat.hp) 
has.bp <- all.freq %>% 
  filter(C.lab==0, N.lab==0, n>2) 
all.freq.trim <- all.freq %>% 
  filter(n>2,
         `Protein Name` %in% has.bp$`Protein Name`) 

big.tab.of.conc.and.MID.to.write <- full_join(mean.MID.all, all.nM.ave) %>%
  filter(`Precursor Ion Name` %in% all.freq.trim$`Precursor Ion Name`,
         `Protein Name` !="S-Farnesyl-L-cysteine methyl ester") %>%
  full_join(aminoaciddata) %>%
  unique() %>%
  dplyr::rename(Isotopologue = `Precursor Ion Name`,
                Metabolite = `Protein Name`) %>%
  dplyr::select(Experiment, Timepoint, `Incubation time (hr)`,Metabolite,
                Isotopologue, C.lab, N.lab, iso.mz, ionization, 
                nM.in.Sample.ave, nM.in.sample.sd,
                MID.mean, MID.sd) %>%
  arrange(Metabolite, iso.mz, Experiment)

write_csv(big.tab.of.conc.and.MID.to.write, "FateTidy/output/TargetedMetabolites//SuppTab_Average_nM_and_MID.csv")


## get list of metabs that we measured isotopes in -----
unique(all.freq.trim$`Protein Name`)
list.of.metabs.w.any.sig <- big.tab.of.conc.and.MID.to.write %>%
  dplyr::select(Metabolite, Isotopologue, Experiment) %>%
  unique() %>%
  left_join(all.MID.sig.to.write) %>%
  filter(sig ==T) %>%
  dplyr::select(Metabolite) %>%
  unique() %>%
  mutate(Isotopologue.Significant.Change = T) %>%
  full_join(big.tab.of.conc.and.MID.to.write %>%
              dplyr::select(Metabolite) %>%
              unique()) %>%
  mutate(Isotopologue.Significant.Change = ifelse(is.na(Isotopologue.Significant.Change),F,Isotopologue.Significant.Change))

write_csv(list.of.metabs.w.any.sig, "FateTidy/output/TargetedMetabolites//SuppTab_list_of_metabolites_we_lookedAt.csv")

# ## combine big MID table with untargeted -------
# Timepoint.key = big.tab.of.conc.and.MID.to.write %>%
#   dplyr::select(Experiment, Timepoint, `Incubation time (hr)`) %>%
#   unique() 
# untargeted_MID = read_csv("FateTidy/output/TargetedMetabolites//Combined.Untargeted.Significant.MF.MIDs.csv")
# combined.MID.tab = full_join(big.tab.of.conc.and.MID.to.write ,
#                              untargeted_MID %>% left_join(Timepoint.key)) %>%
#   mutate(Experiment = ifelse(Experiment == "Fate1","North",
#                              ifelse(Experiment =="Fate2","South",Experiment)))
# write_csv(combined.MID.tab, "FateTidy/output/TargetedMetabolites//SuppTab_Average_nM_and_MID_Targeted_and_Untargeted.csv")
# 


# ## combine big significant compound table with untargeted -------
# Targeted.MID.sig.to.write <- all.MID.sig %>%
#   dplyr::select(`Precursor Ion Name`, Experiment, sig) %>%
#   # dplyr::rename(Isotopologue =`Precursor Ion Name` ) %>%
#   # mutate(Isotopologue =`Precursor Ion Name`) %>% 
#   separate(`Precursor Ion Name`, into = c("Metabolite.Name","Isotopologue"), sep = "_13") %>%
#   mutate(Isotopologue = paste0("13",Isotopologue),
#          Isotopologue = ifelse(grepl("15N",Isotopologue),Isotopologue, paste0(Isotopologue, " 15N-0"))) %>%
#   unique() %>%
#   filter(sig,
#          Isotopologue!='13C-0 15N-0') %>%
#   ungroup() %>%
#   group_by(Metabolite.Name, Experiment) %>%
#   summarise(isotopes_w_trends = paste(Isotopologue,collapse = ";"))
# 
# untargeted_sig = read_csv("FateTidy/output/TargetedMetabolites//Combined.Untargeted.Significant.MF.Info.csv")
# combined.sig.tab = full_join(Targeted.MID.sig.to.write ,
#                              untargeted_sig %>%
#                                mutate(Experiment = ifelse(Experiment == 1,"Fate1","Fate2"),
#                                       Metabolite.Name = ifelse(Experiment == "Fate1",
#                                                                paste0("N-",feature),
#                                                                paste0("S-",feature)),
#                                       LC.column = ifelse(LC.column=="cyano","cyano","HILIC")) %>%
#                                dplyr::rename(untargeted_compound_ID = compound_name) %>%
#                                dplyr::select(-feature))%>%
#   mutate(Experiment = ifelse(Experiment == "Fate1","North",
#                              ifelse(Experiment =="Fate2","South",Experiment))) %>%
#   dplyr::select(Experiment, everything())
# write_csv(combined.sig.tab, "FateTidy/output/TargetedMetabolites//SuppTab_significant_Trends_Targeted_and_Untargeted.csv")


