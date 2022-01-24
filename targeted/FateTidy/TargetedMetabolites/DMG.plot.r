## DMG plot

library(tidyverse)
library(cowplot)
library(ggplot2)
library(wesanderson)
library(ggpubr)
library(RColorBrewer)

facet.exp.names <- c("Fate1" = "North",
                     "Fate2" = "South")

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
iso.poss.plus <- read_csv("FateTidy/All_C-N_Isotope_possibilities_Ingalls_Standards.csv",
                          col_types = cols(
                            .default = col_character(),
                            m.z = col_double(),
                            z = col_double(),
                            C = col_double(),
                            this.N = col_double(),
                            possibilities = col_double(),
                            iso.mz = col_double(),
                            X = col_double(),
                            RT..min. = col_double(),
                            Conc..uM = col_double(),
                            HILICMix = col_character(),
                            Date.added = col_double()
                          )) %>%
  filter(Column == "HILIC") %>%
  filter(Compound.and.Iso.name == "Glycerophosphocholine_13C-6 15N-0"|
           Compound.and.Iso.name == "Glycerophosphocholine_13C-7 15N-0"|
           Compound.and.Iso.name == "Glycerophosphocholine_13C-8 15N-0"|
           Compound.and.Iso.name == "Glycerophosphocholine_13C-6 15N-1"|
           Compound.and.Iso.name == "Glycerophosphocholine_13C-7 15N-1"|
           Compound.and.Iso.name == "(3-Carboxypropyl)trimethylammonium_13C-6 15N-1"|
           Compound.and.Iso.name == "Carnitine_13C-6 15N-1"|
           Compound.and.Iso.name == "Carnitine_13C-6 15N-0"|
           Compound.and.Iso.name == "Carnitine_13C-7 15N-1"
  )

iso.poss.use <- full_join(iso.poss,iso.poss.plus) %>%
  dplyr::select(Compound.and.Iso.name, iso.mz, RT..min., Column, Fraction1) %>%
  filter(Fraction1 == "HILICPos",
         grepl("Dimethylglycine",Compound.and.Iso.name))

qc.dat <- read_csv("FateTidy/output/TargetedMetabolites//HILICPos_QC_skyline_isotope_data_blankSubtract.csv") %>%
  filter(grepl("Dimethylglycine",`Precursor Ion Name`))
freq.dat <- read_csv("FateTidy/output/TargetedMetabolites/HILICPos_frequency_of_isotope_data.csv") %>%
  filter(grepl("Dimethylglycine",`Precursor Ion Name`))
freq.dat.3ormore <- freq.dat %>%
  filter(n >2)

Sample.log <- read_csv("FateTidy/RawData/Sample_Log_GBTFateExperiments.csv") %>%
  filter(grepl("Metab",`Sample Type`) |grepl("metab",`Sample Type`))
Sample.log.to.join <- Sample.log %>%
  filter(!grepl("Blk",`Sample ID`))%>%
  dplyr::select(Experiment, Timepoint,`Time of day`, rep, `Vol Filtered (L)`, `Incubation time (hr)`,
                `Time since last sample (hr)`) %>%
  unique()

DMG.restored <- qc.dat %>%
  mutate(Area = ifelse(RT.check=="bad RT" & Blk.flag=="ok", rawArea, Area)) %>%
  full_join(Sample.log.to.join) %>%
  left_join(iso.poss.use) %>%
  arrange(iso.mz) %>%
  filter(Experiment == "Fate1"|Experiment=="Fate2") %>%
  mutate(Experiment = factor(Experiment, levels = c("Fate2","Fate1")))
DMG.mean = DMG.restored %>%
  group_by(`Precursor Ion Name`, Experiment, SampID, iso.mz, Timepoint, `Incubation time (hr)`, `Time since last sample (hr)`,
           `Time of day`) %>%
  summarise(meanDMG = mean(Area, na.rm = T),
            sdDMG = sd(Area, na.rm = T))

rects <- tibble(ymin = -Inf, ymax = Inf, xmin = c(11,11, 11+24, 11+24, 11+48, 11+48, 11+24*3, 11+24*3),
                xmax = c(23, 23, 23+24, 23+24, 23+24*2, 23+24*2, 23+24*3, 23+24*3),
                Experiment = factor(rep(c("Fate1","Fate2"),4), levels = c("Fate2","Fate1")))


just.13c4.15n1 <- ggplot(DMG.restored %>% filter(`Precursor Ion Name` == "Dimethylglycine_13C-4 15N-1")) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_point(aes(x = `Incubation time (hr)`, y = Area)) +
  geom_line(data = DMG.mean%>% filter(`Precursor Ion Name` == "Dimethylglycine_13C-4 15N-1"),
            aes(x = `Incubation time (hr)`, y = meanDMG))+
  geom_errorbar(data = DMG.mean%>% filter(`Precursor Ion Name` == "Dimethylglycine_13C-4 15N-1"), 
                aes(x = `Incubation time (hr)`, ymin = meanDMG - sdDMG, ymax = meanDMG + sdDMG)) +
  facet_grid(`Precursor Ion Name`~Experiment, scales = "free",
             labeller = labeller( .cols = facet.exp.names)) +
  theme_bw()

just.13c4.15n1

# ggsave(filename = "Fate/output/Figures/DMG_labeled.pdf", width = 7, height = 4)

just.13c4.15n1.noRect <- ggplot(DMG.restored %>% filter(`Precursor Ion Name` == "Dimethylglycine_13C-4 15N-1")) +
  geom_point(aes(x = `Incubation time (hr)`, y = Area)) +
  geom_line(data = DMG.mean%>% filter(`Precursor Ion Name` == "Dimethylglycine_13C-4 15N-1"),
            aes(x = `Incubation time (hr)`, y = meanDMG))+
  geom_errorbar(data = DMG.mean%>% filter(`Precursor Ion Name` == "Dimethylglycine_13C-4 15N-1"), 
                aes(x = `Incubation time (hr)`, ymin = meanDMG - sdDMG, ymax = meanDMG + sdDMG))+
  facet_grid(`Precursor Ion Name`~Experiment, scales = "free",
             labeller = labeller( .cols = facet.exp.names)) +
  theme_bw()

just.13c4.15n1.noRect

ggsave(filename = "FateTidy/output/TargetedMetabolites/DMG_labeled_noRects.pdf", width = 7, height = 4)
