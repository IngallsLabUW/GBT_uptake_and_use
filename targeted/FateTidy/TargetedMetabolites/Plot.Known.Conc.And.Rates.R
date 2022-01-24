## quantify knowns plot

## libs ----------
library(tidyverse)
library(stats)
library(RCurl)
library(here)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

## get sample log -----
Sample.log <- read_csv("FateTidy/RawData/Sample_Log_GBTFateExperiments.csv") %>%
  filter(grepl("Metab",`Sample Type`) |grepl("metab",`Sample Type`))
Sample.log.to.join <- Sample.log %>%
  filter(!grepl("Blk",`Sample ID`))%>%
  dplyr::select(Experiment, Timepoint,`Time of day`, rep, `Vol Filtered (L)`, `Incubation time (hr)`,
                `Time since last sample (hr)`) %>%
  unique()
## get data ----
All.nM.inSamples <- read_csv("FateTidy/output/TargetedMetabolites/nM.concentrations.allFractions.csv",
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
                             )) %>%
  mutate(Experiment = factor(Experiment, levels = c("Fate2","Fate1")))
## pick compounds -----------------
cmpds <- c("Glycine betaine",
           "Sarcosine",
           "Creatine",
           "Trimethylamine N−oxide",
           "Choline",
           "Glycerophosphocholine",
           "Carnitine",
           "(3−Carboxypropyl)trimethylammonium",
           "L-Glutamic acid",
           "L-Glutamine")

isotopologues <- c("Glycine betaine_13C-5 15N-1",
                     "Sarcosine_13C-3 15N-1",
                     "Creatine_13C-3 15N-1",
                     "Carnitine_13C-5 15N-1",
                     "(3-Carboxypropyl)trimethylammonium_13C-5 15N-1",
                     "Choline_13C-5 15N-1",
                     "Glycerophosphocholine_13C-5 15N-1",
                     "Trimethylamine N-oxide_13C-3 15N-1",
                   "L-Glutamic acid_13C-0 15N-1",
                   "L-Glutamic acid_13C-1 15N-1",
                   "L-Glutamic acid_13C-1 15N-0",
                   "L-Glutamine_13C-0 15N-1",
                   "L-Glutamine_13C-1 15N-0")
## get std info -----------

Exp.labs <-  c("North", "South")
names(Exp.labs) <- c("Fate1", "Fate2")

rects <- tibble(ymin = -Inf, ymax = Inf, xmin = c(11,11, 11+24, 11+24, 11+48, 11+48, 11+24*3, 11+24*3),
                xmax = c(23, 23, 23+24, 23+24, 23+24*2, 23+24*2, 23+24*3, 23+24*3),
                Experiment = factor(rep(c("Fate1","Fate2"),4), levels = c("Fate2", 'Fate1')))


library (readr)

std.url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv"

stds.dat <- read.csv(url(std.url), header = T)

# stds.dat <- read.csv(text = getURL(std.url), header = T)
stds.dat.to.use <- stds.dat %>%
  dplyr::rename(Compound.Name = Compound_Name_Original,
                Emperical.Formula = Empirical_Formula)%>%
  mutate(Compound.Name = ifelse(Compound_Name == "beta-Alaninebetaine", "beta-Alaninebetaine", 
                                ifelse(Compound_Name == "L-Alanine", "L-Alanine", 
                                       ifelse(Compound_Name == "L-Arginine", "L-Arginine", 
                                              ifelse(Compound_Name == "L-Glutamic acid", "L-Glutamic acid",  
                                                     ifelse(Compound_Name == "(3-Carboxypropyl)trimethylammonium","(3-Carboxypropyl)trimethylammonium",
                                                     ifelse(Compound_Name == "L-Glutamine", "L-Glutamine", Compound.Name )))))))%>%
  filter(Column == "HILIC",
         !grepl("13C",Compound.Name),
         !grepl("15N",Compound.Name),
         !grepl("D",Emperical.Formula)) %>%
  full_join(.,stds.dat %>%
              dplyr::rename(Compound.Name = Compound_Name_Original,
                            Emperical.Formula = Empirical_Formula) %>%
              mutate(Compound.Name = ifelse(Compound_Name == "beta-Alaninebetaine", "beta-Alaninebetaine", 
                                            ifelse(Compound_Name == "L-Alanine", "L-Alanine", 
                                                   ifelse(Compound_Name == "L-Arginine", "L-Arginine", 
                                                          ifelse(Compound_Name == "L-Glutamic acid", "L-Glutamic acid",  
                                                                 ifelse(Compound_Name == "(3-Carboxypropyl)trimethylammonium","(3-Carboxypropyl)trimethylammonium",
                                                                        ifelse(Compound_Name == "L-Glutamine", "L-Glutamine", Compound.Name )))))))%>%
              filter(Column == "RP",
                     !grepl("13C",Compound.Name),
                     !grepl("15N",Compound.Name),
                     !grepl("D",Emperical.Formula),
                     Compound.Name=="L-Phenylalanine"|
                       Compound.Name=="L-Tryptophan"|
                       Compound.Name=="Methylthioadenosine"|
                       Compound.Name == "Butyrylcarnitine"|
                       Compound.Name=="	Xanthine") ) %>%
  dplyr::rename(Identification = Compound.Name,
                BestMatch = Compound_Name_Figure) %>%
  dplyr::select(BestMatch, Identification, Emperical.Formula, HILIC_Mix, Concentration_uM) %>% 
  unique() %>%
  mutate(Emp.to.split = Emperical.Formula) %>%
  separate(col = Emp.to.split, sep = "H", into = c("C","Other")) %>%
  separate(col = C, sep = "C", into = c("C","C.num")) %>%
  mutate(C.num = ifelse(C.num == "",1,ifelse(is.na(C.num),1,as.numeric(C.num)))) %>%
  mutate(N.temp = ifelse(grepl("N",Other),1,0)) %>%
  separate(col = Other, sep = "N", into = c("other","N.2")) %>% 
  separate(col = N.2, sep = "O", into = c("N.3","other")) %>%
  mutate(N.num  = ifelse(N.temp == 0, N.temp,
                         ifelse(is.na(other) & is.na(N.3),0,
                                ifelse(is.na(N.3), N.temp,
                                       ifelse(N.3=="",N.temp,
                                ifelse(is.na(other) & N.3 == "",N.temp,
                                       ifelse(other =="" & N.3 == "", N.temp, as.numeric(N.3)))))))) %>%
  dplyr::select(BestMatch, Identification, Emperical.Formula, C.num, N.num) %>%
  mutate(Identification = ifelse(BestMatch == "Glycine betaine", "Glycine betaine", Identification))


## get isotopologue concentration over time -------------

isotopo.conc.all <- All.nM.inSamples %>%
  # filter(`Precursor Ion Name` %in% isotopologues) %>%
  left_join(., stds.dat.to.use, by = c("Protein Name" = "Identification")) %>%
  mutate(nM.C.lab.in.isotopologue = nM.in.sample * C.lab,
         nM.N.lab.in.isotopologue = nM.in.sample * N.lab) 
isotopo.conc.aves.all <- isotopo.conc.all %>%
  group_by(Experiment, SampID, Timepoint, `Time of day`, `Time since last sample (hr)`, `Incubation time (hr)`,
           `Precursor Ion Name`, `Protein Name`, C.lab, N.lab, C.num, N.num) %>%
  summarise(nM.13C.in.isotopologue.sd = sd(nM.C.lab.in.isotopologue),
            nM.13C.in.isotopologue = mean(nM.C.lab.in.isotopologue),
            nM.15N.in.isotopologue.sd = sd(nM.N.lab.in.isotopologue),
            nM.15N.in.isotopologue = mean(nM.N.lab.in.isotopologue))

isotopo.conc <-isotopo.conc.all %>%
  filter(`Precursor Ion Name` %in% isotopologues) %>%
  mutate(`Precursor Ion Name` = factor(`Precursor Ion Name`, levels = isotopologues))
isotopo.conc.aves <- isotopo.conc.aves.all %>%
  filter(`Precursor Ion Name` %in% isotopologues) %>%
  mutate(`Precursor Ion Name` = factor(`Precursor Ion Name`, levels = isotopologues))

ggplot() +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar(data = isotopo.conc.aves,
            aes(x = `Incubation time (hr)`, ymin = nM.13C.in.isotopologue - nM.13C.in.isotopologue.sd,
                ymax = nM.13C.in.isotopologue + nM.13C.in.isotopologue.sd)) +
  geom_line(data = isotopo.conc.aves,
            aes(x = `Incubation time (hr)`, y =nM.13C.in.isotopologue )) +
  facet_grid(`Precursor Ion Name`~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))
# ggsave(filename = "FateTidy/output/skyline/nM.C.in.select.isotopologues.v.time.pdf",
#        height = 8, width = 6.5)

ggplot() +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar(data = isotopo.conc.aves,
                aes(x = `Incubation time (hr)`, ymin = nM.15N.in.isotopologue - nM.15N.in.isotopologue.sd,
                    ymax = nM.15N.in.isotopologue + nM.15N.in.isotopologue.sd)) +
  geom_line(data = isotopo.conc.aves,
            aes(x = `Incubation time (hr)`, y =nM.15N.in.isotopologue )) +
  facet_grid(`Precursor Ion Name`~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))
# ggsave(filename = "FateTidy/output/skyline/nM.N.in.select.isotopologues.v.time.pdf",
#        height = 8, width = 6.5)


## get rate ----------------
previou.time.point = isotopo.conc.aves.all %>%
  mutate(previous.time.point.ID = ifelse(Experiment=="Fate1" & Timepoint==4, "GBTFate1MT0",
                                         ifelse(Experiment=="Fate1" & Timepoint==8, "GBTFate1MT4",
                                                ifelse(Experiment=="Fate1" & Timepoint==12, "GBTFate1MT8", 
                                                       ifelse(Experiment=="Fate1" & Timepoint==36, "GBTFate1MT12", 
                                                              ifelse(Experiment=="Fate1" & Timepoint==50, "GBTFate1MT36",
                                                                     ifelse(Experiment=="Fate1" & Timepoint==96,"GBTFate1MT50",
                                                                            ifelse(Experiment=="Fate2" & Timepoint==5,  "GBTFate2MT0", 
                                                                                   ifelse(Experiment=="Fate2" & Timepoint==12,"GBTFate2MT5",
                                                                                          ifelse(Experiment=="Fate2" & Timepoint==24,"GBTFate2MT12",
                                                                                                 ifelse(Experiment=="Fate2" & Timepoint==48,"GBTFate2MT24",
                                                                                                        ifelse(Experiment=="Fate2" & Timepoint==96, "GBTFate2MT48",
                                                                                                               NA))))))))))))
previous.time.point.conc = isotopo.conc.aves %>%
  ungroup() %>%
  dplyr::rename(previous.con.nM.13C = nM.13C.in.isotopologue,
                previous.con.nM.15N = nM.15N.in.isotopologue,
                previous.time.point.ID = SampID) %>%
  dplyr::select(previous.con.nM.13C,previous.con.nM.15N,  previous.time.point.ID, `Precursor Ion Name`) %>%
  full_join(.,previou.time.point) %>%
  filter(!is.na(SampID)) %>%
  mutate(previous.con.nM.13C = ifelse(is.na(previous.con.nM.13C),0,previous.con.nM.13C),
         previous.con.nM.15N = ifelse(is.na(previous.con.nM.15N),0,previous.con.nM.15N))
Uptake.rate.of.full.labeled.all <- previous.time.point.conc %>%
  mutate(uptake.overall.nM.13C.per.h = ifelse(!grepl("Glycine betaine", `Protein Name`) & Timepoint == 0,
                                            0,
                                            nM.13C.in.isotopologue/`Incubation time (hr)`),
         uptake.since.last.point.nM.13C.per.h = (nM.13C.in.isotopologue-previous.con.nM.13C)/`Time since last sample (hr)`,
         uptake.since.last.point.nM.13C.per.h = ifelse(!grepl("Glycine betaine", `Protein Name`) & Timepoint == 0,
                                                     0,
                                                     ifelse(is.infinite(uptake.since.last.point.nM.13C.per.h),
                                                     uptake.overall.nM.13C.per.h,
                                                   uptake.since.last.point.nM.13C.per.h)),
         uptake.overall.nM.15N.per.h = ifelse(!grepl("Glycine betaine", `Protein Name`) & Timepoint == 0,
                                            0,
                                            nM.15N.in.isotopologue/`Incubation time (hr)`),
         uptake.since.last.point.nM.15N.per.h = (nM.15N.in.isotopologue-previous.con.nM.15N)/`Time since last sample (hr)`,
         uptake.since.last.point.nM.15N.per.h = ifelse(!grepl("Glycine betaine", `Protein Name`) & Timepoint == 0,
                                                     0,
                                                     ifelse(is.infinite(uptake.since.last.point.nM.15N.per.h),
                                                     uptake.overall.nM.15N.per.h,
                                                     uptake.since.last.point.nM.15N.per.h))) 

Uptake.rate.of.full.labeled <- Uptake.rate.of.full.labeled.all %>%
  filter(`Precursor Ion Name` %in% isotopologues) %>%
  mutate(`Precursor Ion Name` = factor(`Precursor Ion Name`, levels = isotopologues))
uptake.rate.plot <- ggplot(Uptake.rate.of.full.labeled) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_line(aes(y = uptake.since.last.point.nM.15N.per.h, x =`Incubation time (hr)`), color = "blue") +
  geom_line(aes(y = uptake.overall.nM.15N.per.h, x =`Incubation time (hr)`), color = "black") +
  facet_grid(`Precursor Ion Name`~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))
uptake.rate.plot

# ggsave(filename = "FateTidy/output/skyline/uptake.of.N.nM.N.per.h.select.isotopologues.v.time.pdf",
       # height = 8, width = 6.5)

uptake.rate.plot.c <- ggplot(Uptake.rate.of.full.labeled) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_line(aes(y = uptake.since.last.point.nM.13C.per.h, x =`Incubation time (hr)`), color = "blue") +
  geom_line(aes(y = uptake.overall.nM.13C.per.h, x =`Incubation time (hr)`), color = "black") +
  facet_grid(`Precursor Ion Name`~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))
uptake.rate.plot.c

# ggsave(filename = "FateTidy/output/skyline/uptake.of.C.nM.N.per.h.select.isotopologues.v.time.pdf",
       # height = 8, width = 6.5)



## write rate per isotopologue -----
write_csv(Uptake.rate.of.full.labeled.all, path = "FateTidy/output/TargetedMetabolites/isotopologues.rates.csv")


## just write max rates per cmpd ----
Uptake.rate.of.full.labeled.all.ax.ests = read_csv("FateTidy/output/TargetedMetabolites/isotopologues.rates.csv") %>%
  group_by(Experiment, `Protein Name`) %>%
  summarise(uptake.overall.pM.13C.per.h.max = max(uptake.overall.nM.13C.per.h, na.rm = T)*1e3,
            uptake.since.last.point.pM.13C.per.h.max = max(uptake.since.last.point.nM.13C.per.h, na.rm = T)*1e3,
            uptake.overall.pM.15N.per.h.max = max(uptake.overall.nM.15N.per.h, na.rm = T)*1e3,
            uptake.since.last.point.pM.15N.per.h = max(uptake.since.last.point.nM.15N.per.h, na.rm = T)*1e3)


## write total 15N and 13C 

total.15N.13C <-Uptake.rate.of.full.labeled.all %>%
  # filter(Timepoint %in% c(4,5,6,7,8,9,10,12,13,14)) %>%
  group_by(Timepoint, Experiment, `Incubation time (hr)`, `Time since last sample (hr)`, SampID,
           previous.time.point.ID) %>%
  summarise(total.15.N.in.quantifiable.metabs = sum(nM.15N.in.isotopologue, na.rm = T),
            total.15.N.in.quantifiable.metabs.sd = sqrt(sum((nM.15N.in.isotopologue.sd^2), na.rm = T)),
            total.13.C.in.quantifiable.metabs = sum(nM.13C.in.isotopologue, na.rm = T),
            total.13.C.in.quantifiable.metabs.sd = sqrt(sum((nM.13C.in.isotopologue.sd^2), na.rm = T)))


total.15N.13C.subset <-Uptake.rate.of.full.labeled %>%
  # filter(Timepoint %in% c(4,5,6,7,8,9,10,12,13,14)) %>%
  group_by(Timepoint, Experiment, `Incubation time (hr)`, `Time since last sample (hr)`, SampID,
           previous.time.point.ID) %>%
  summarise(total.15.N.in.quantifiable.metabs = sum(nM.15N.in.isotopologue, na.rm = T),
            total.15.N.in.quantifiable.metabs.sd = sqrt(sum((nM.15N.in.isotopologue.sd^2), na.rm = T)),
            total.13.C.in.quantifiable.metabs = sum(nM.13C.in.isotopologue, na.rm = T),
            total.13.C.in.quantifiable.metabs.sd = sqrt(sum((nM.13C.in.isotopologue.sd^2), na.rm = T)))
# write_csv(total.15N.13C, path = "FateTidy/output/skyline/total.15N.13C.detected.in.quantifiable.metabolites.csv")

## MAKE PLOTS OF TOTAL C and N take up -----------
total.n.plot <- ggplot(total.15N.13C) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_line(aes(x = `Incubation time (hr)`,
                y = total.15.N.in.quantifiable.metabs)) +
  geom_errorbar(aes(x =  `Incubation time (hr)`,
                    ymin = total.15.N.in.quantifiable.metabs -total.15.N.in.quantifiable.metabs.sd,
                    ymax = total.15.N.in.quantifiable.metabs + total.15.N.in.quantifiable.metabs.sd)) +
  geom_line(data = Uptake.rate.of.full.labeled.all %>%
              filter(`Protein Name` == "Glycine betaine",
                     nM.15N.in.isotopologue>1e-5),
            aes(x =`Incubation time (hr)`, y = nM.15N.in.isotopologue, color = `Precursor Ion Name`,
                group = `Precursor Ion Name`))+
  ylab("15N (nM)") +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  theme_cowplot() 

totalc.plot <- ggplot(total.15N.13C) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_line(aes(x = `Incubation time (hr)`,
                y = total.13.C.in.quantifiable.metabs)) +
  geom_errorbar(aes(x =  `Incubation time (hr)`,
                    ymin = total.13.C.in.quantifiable.metabs -total.13.C.in.quantifiable.metabs.sd,
                    ymax = total.13.C.in.quantifiable.metabs + total.13.C.in.quantifiable.metabs.sd)) +
  geom_line(data = Uptake.rate.of.full.labeled.all %>%
              filter(`Protein Name` == "Glycine betaine",
                     nM.15N.in.isotopologue>1e-5),
            aes(x =`Incubation time (hr)`, y = nM.13C.in.isotopologue, color = `Precursor Ion Name`,
                group = `Precursor Ion Name`))+
  ylab("13C (nM)") +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  theme_cowplot() 

plot_grid(total.n.plot, totalc.plot, nrow = 2)

# ggsave(filename = "FateTidy/output/Figures/Total.15N.13N.In.Metabs.wGBT.pdf",
       # height = 4, width = 6.5)

# PCPN <- source("PCPN.Dat.R")[[1]] %>%
  # mutate(C.N.Ratio = TC.umol.L/N.umol.L)

## stacked area plot for total 15N and 13C -------
stacked.area.data <- Uptake.rate.of.full.labeled.all %>%
  filter(nM.13C.in.isotopologue > 0 | nM.15N.in.isotopologue >0) %>%
  group_by(`Protein Name`, `Incubation time (hr)`, Experiment, SampID, `Time of day`) %>%
  summarise(nm.13C.in.Metab = sum(nM.13C.in.isotopologue),
            nm.15N.in.Metab = sum(nM.15N.in.isotopologue)) %>%
  full_join(.,Uptake.rate.of.full.labeled.all %>%
              filter(nM.13C.in.isotopologue > 0 | nM.15N.in.isotopologue >0) %>%
              group_by( `Incubation time (hr)`, Experiment, SampID, `Time of day`) %>%
              summarise(tot.13C = sum(nM.13C.in.isotopologue),
                        tot.15N = sum(nM.15N.in.isotopologue))) %>%
  mutate(percent.of.13C = nm.13C.in.Metab/tot.13C * 100,
         percent.of.15N = nm.15N.in.Metab/tot.15N * 100) %>%
  arrange(nm.13C.in.Metab)
length(unique(stacked.area.data$`Protein Name`))

stacked.area.plot.data.big <- stacked.area.data %>%
  group_by(`Protein Name`,  Experiment) %>%
  summarise(percent.of.13C.ave = mean(percent.of.13C, na.rm = T),
            percent.of.15N.ave = mean(percent.of.15N, na.rm = T)) %>%
  filter(percent.of.13C.ave > 1.33 | percent.of.15N.ave > 1.33)
length(unique(stacked.area.plot.data.big$`Protein Name`))

stacked.area.plot.data <-   stacked.area.data %>%
  filter(!(`Protein Name` %in% unique(stacked.area.plot.data.big$`Protein Name`))) %>%
  group_by( `Incubation time (hr)`, Experiment, SampID, `Time of day`) %>%
  summarise(nm.13C.in.Metab = sum(nm.13C.in.Metab),
            nm.15N.in.Metab = sum(nm.15N.in.Metab),
            percent.of.13C = sum(percent.of.13C),
            percent.of.15N = sum(percent.of.15N)) %>%
  mutate(`Protein Name` = "Other") %>%
  full_join(., stacked.area.data %>%
              filter(`Protein Name` %in% unique(stacked.area.plot.data.big$`Protein Name`))) %>%
  mutate(Metabolite = factor(`Protein Name`, levels = c(unique(stacked.area.plot.data.big$`Protein Name`),"Other")),
         Experiment = factor(Experiment,levels = c("Fate2","Fate1")))


stacked.13c.plot <- ggplot(stacked.area.plot.data) +
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
            # aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_area(aes(x = `Incubation time (hr)`,
                y = nm.13C.in.Metab,
                fill = Metabolite),
                color = "black", alpha = 0.85) +
  ylab("13C (nM)") +
  scale_fill_brewer(palette  = "Greys") +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  theme_bw() 

stacked.13c.plot

stacked.13c.plot.rel <- ggplot(stacked.area.plot.data) +
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
            # aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_area(aes(x = `Incubation time (hr)`,
                y = percent.of.13C,
                fill = Metabolite),
            color = "black") +
  ylab("Percent of 13C") +
  scale_fill_brewer(palette  = "Greys") +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  theme_bw() 

stacked.13c.plot.rel

stacked.15N.plot.rel <- ggplot(stacked.area.plot.data) +
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
  # aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_area(aes(x = `Incubation time (hr)`,
                y = percent.of.15N,
                fill = Metabolite),
            color = "black") +
  ylab("Percent of 15N") +
  scale_fill_brewer(palette  = "Greys") +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  theme_bw() +
  theme(legend.position = "bottom")

stacked.15N.plot <- ggplot(stacked.area.plot.data) +
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
            # aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_area(aes(x = `Incubation time (hr)`,
                y = nm.15N.in.Metab,
                fill = Metabolite),
            color = "black") +
  ylab("15N (nM)") +
  scale_fill_brewer(palette  = "Greys") +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  theme_bw() +
  theme(legend.position = "bottom")

stacked.15N.plot.rel +
  theme(legend.position = "none")

legend <- get_legend(
  # create some space to the left of the legend
  stacked.15N.plot.rel + theme(legend.box.margin = margin(0, 0, 0, 12))
)

rel.plots.together <- 
  plot_grid(stacked.15N.plot.rel +
              theme(legend.position = "none"),
            stacked.13c.plot.rel +
              theme(legend.position = "none"),
            nrow = 2,
            align = "hv")
plot_grid(rel.plots.together,legend, nrow  = 2, rel_heights  = c(2, 0.3))
ggsave("FateTidy/output/TargetedMetabolites/Stacked_13C_15N_relative_concentration_OverTime_noRects_switchOrder.pdf",
       height = 6, width = 7)


abs.plots.together <- 
  plot_grid(stacked.15N.plot +
              theme(legend.position = "none"),
            stacked.13c.plot +
              theme(legend.position = "none"),
            nrow = 2,
            align = "hv")
plot_grid(abs.plots.together,legend, nrow  = 2, rel_heights  = c(2, 0.3)) 
ggsave("FateTidy/output/TargetedMetabolites/Stacked_13C_15N_concentration_OverTime_noRects_switchOrder.pdf",
       height = 6, width = 7)

## do different subsets ------
Total15N.in.GBT <- Uptake.rate.of.full.labeled.all %>%
  filter(`Protein Name` == "Glycine betaine") %>%
  group_by(Timepoint, Experiment, `Incubation time (hr)`, `Time since last sample (hr)`, SampID,
           previous.time.point.ID) %>%
  summarise(total.15.N.in.GBT = sum(nM.15N.in.isotopologue, na.rm = T),
            total.15.N.in.GBT.sd = sqrt(sum((nM.15N.in.isotopologue.sd^2), na.rm = T)),
            total.13.C.in.GBT = sum(nM.13C.in.isotopologue, na.rm = T),
            total.13.C.in.GBT.sd = sqrt(sum((nM.13C.in.isotopologue.sd^2), na.rm = T))) %>%
  full_join(.,total.15N.13C) %>%
  mutate(percent.15N.in.GBT = total.15.N.in.GBT/total.15.N.in.quantifiable.metabs,
         percent.13C.in.GBT = total.13.C.in.GBT/total.13.C.in.quantifiable.metabs,
         percent.15N.in.GBT.sd = percent.15N.in.GBT * sqrt((total.15.N.in.GBT.sd/total.15.N.in.GBT)^2 +
                                                             (total.15.N.in.quantifiable.metabs.sd/total.15.N.in.quantifiable.metabs)^2),
         percent.13C.in.GBT.sd = percent.13C.in.GBT * sqrt((total.13.C.in.GBT.sd/total.13.C.in.GBT)^2 +
                                                             (total.13.C.in.quantifiable.metabs.sd/total.13.C.in.quantifiable.metabs)^2)) %>%
  arrange(Experiment, Timepoint)

# write_csv(Total15N.in.GBT, path = "FateTidy/output/skyline/PercentOf15Nand13CInGBT.csv")



## get total concentration over time of select cmpds -----------------
# cmpds = c("Homarine",
#           "Trigonelline",
#           "Taurine",
#           "beta-Alaninebetaine",
#           "beta-Alanine",
#           "L-Alanine",
#           "L-Arginine",
#           "Guanosine",
#           "Guanine",
#           "Adenosine",
#           "Adenine")
# cmpds <- c("Glycine betaine",
#            "Taurine",
#            "Carnitine",
#            "Trimethylamine N−oxide",
#            "DMSP",
#            "DHPS",
#            "Homarine",
#            "(3−Carboxypropyl)trimethylammonium",
#            "L-Glutamic acid",
#            "L-Proline",
#            "L-Leucine",
#            "L-Isoleucine",
#            "Sucrose",
#            "Gonyol",
#            "Glucosylglycerol",
#            "Isethionic acid",
#            "(R)-2,3-Dihydroxypropane-1-sulfonate")
cmpds <- c("Glycine betaine",
           "Sarcosine",
           "Creatine",
           "Trimethylamine N−oxide",
           "Choline",
           "Glycerophosphocholine",
           "Carnitine",
           "(3−Carboxypropyl)trimethylammonium",
           "L-Glutamic acid",
           "L-Glutamine")

total.con.plot.dat <- All.nM.inSamples %>%
  group_by(Experiment, Timepoint, replicate, type, `Protein Name`, SampID,
                `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
  summarise(Total.nM.in.sample = sum(nM.in.sample)) %>%
  left_join(., stds.dat.to.use, by = c("Protein Name" = "Identification")) %>%
  mutate(`Protein Name` = factor(`Protein Name`, levels = cmpds),
         Total.nM.in.sample.C = C.num *Total.nM.in.sample ,
         Total.nM.in.sample.N = N.num * Total.nM.in.sample) %>%
  filter(!is.na(`Protein Name`),
         `Protein Name`!="DMSP",
         `Protein Name`!="Glucosylglycerol") %>%
  full_join(.,All.nM.inSamples %>%
              filter( `Protein Name`=="DMSP"|
                      `Protein Name`=="Glucosylglycerol") %>%
              group_by(Experiment, Timepoint, replicate, type, `Protein Name`, SampID,
                       `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
              summarise(Total.nM.in.sample = sum(nM.in.sample)) %>%
              left_join(., stds.dat.to.use, by = c("Protein Name" = "BestMatch")) %>%
              mutate(`Protein Name` = factor(`Protein Name`, levels = cmpds),
                     Total.nM.in.sample.C = C.num *Total.nM.in.sample ,
                     Total.nM.in.sample.N = N.num * Total.nM.in.sample) %>%
              filter(!is.na(`Protein Name`)))

total.con.aves <- total.con.plot.dat %>%
  ungroup() %>%
  group_by(Experiment, Timepoint, type, `Protein Name`, SampID,
           `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
  summarise(Total.nM.in.sample.N.sd = sd(Total.nM.in.sample.N, na.rm = T),
            Total.nM.in.sample.N = mean(Total.nM.in.sample.N, na.rm = T),
            Total.nM.in.sample.C.sd = sd(Total.nM.in.sample.C, na.rm = T),
            Total.nM.in.sample.C = mean(Total.nM.in.sample.C, na.rm = T)) %>%
  full_join(.,total.con.plot.dat %>%
              ungroup() %>%
              group_by(Experiment, Timepoint, type, `Protein Name`, SampID,
                       `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
              summarise(Total.nM.in.sample.N.sd = sd(Total.nM.in.sample.N),
                        Total.nM.in.sample.N = mean(Total.nM.in.sample.N),
                        Total.nM.in.sample.C.sd = sd(Total.nM.in.sample.C),
                        Total.nM.in.sample.C = mean(Total.nM.in.sample.C)) %>%
              ungroup() %>%
              group_by(Experiment, Timepoint, type,  SampID,
                       `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
              summarise(Total.N = sum(Total.nM.in.sample.N),
                        Total.C = sum(Total.nM.in.sample.C)) ) %>%
  mutate(proportion.N = Total.nM.in.sample.N/Total.N*100,
         proportion.C = Total.nM.in.sample.C/Total.C*100)

ggplot(data = total.con.aves #%>%
         # filter(`Protein Name` == "Carnitine")
       )  + 
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
  #           aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar( aes(x = `Incubation time (hr)`, ymin = Total.nM.in.sample.N - Total.nM.in.sample.N.sd,
                 ymax = Total.nM.in.sample.N + Total.nM.in.sample.N.sd))+
  geom_line(aes(x = `Incubation time (hr)`, y =Total.nM.in.sample.N )) +
  facet_grid(`Protein Name`~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0))+ 
  ylab("Total N in metabolite (nM)") +
  theme_minimal_grid() + panel_border()+
  theme(strip.text.y = element_text(angle = 0))

ggsave(filename = "FateTidy/output/TargetedMetabolites/Total.conc.nM.N.select.metabs.v.time.new.Order.pdf",
height = 8, width = 6.5)

# ggsave(filename = "FateTidy/output/Figures/Total.conc.nM.N.other.metabs.v.time.new.Order.pdf",
# height = 8, width = 6.5)

ggplot(data = total.con.aves)  + 
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar( aes(x = `Incubation time (hr)`, ymin = Total.nM.in.sample.C - Total.nM.in.sample.C.sd,
                     ymax = Total.nM.in.sample.C + Total.nM.in.sample.C.sd))+
  geom_line(aes(x = `Incubation time (hr)`, y =Total.nM.in.sample.C )) +
  facet_grid(`Protein Name`~Experiment, scales = "free_y") +
  theme_minimal() + 
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0))+ 
  theme(strip.text.y = element_text(angle = 0))


## osmo plots ----
osmos <- c("Glycine betaine",
           "Taurine",
           "Carnitine",
           "Trimethylamine N−oxide",
           "DMSP",
           "Homarine",
           "(3−Carboxypropyl)trimethylammonium",
           "L-Glutamic acid",
           "L-Proline",
           "L-Leucine",
           "L-Isoleucine",
           "Gonyol",
           "Glucosylglycerol",
           "Isethionic acid",
           "(R)-2,3-Dihydroxypropane-1-sulfonate"
)
total.con.plot.dat <- All.nM.inSamples %>%
  group_by(Experiment, Timepoint, replicate, type, `Protein Name`, SampID,
           `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
  summarise(Total.nM.in.sample = sum(nM.in.sample)) %>%
  filter(`Protein Name` %in% osmos) %>%
  left_join(., stds.dat.to.use %>%
              mutate(`Protein Name` = ifelse(Identification %in% c("Proline","Leucine","Isoleucine"), paste0("L-",Identification),
                                             ifelse(Identification =="Isethionic Acid","Isethionic acid",
                      ifelse(Identification =="DHPS","(R)-2,3-Dihydroxypropane-1-sulfonate",Identification))))) %>%
  mutate(`Protein Name` = factor(`Protein Name`, levels = osmos),
         Total.nM.in.sample.C = C.num *Total.nM.in.sample ,
         Total.nM.in.sample.N = N.num * Total.nM.in.sample) %>%
  filter(!is.na(`Protein Name`),
         `Protein Name`!="DMSP",
         `Protein Name`!="Glucosylglycerol") %>%
  full_join(.,All.nM.inSamples %>%
              filter( `Protein Name`=="DMSP"|
                        `Protein Name`=="Glucosylglycerol") %>%
              group_by(Experiment, Timepoint, replicate, type, `Protein Name`, SampID,
                       `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
              summarise(Total.nM.in.sample = sum(nM.in.sample)) %>%
              left_join(., stds.dat.to.use, by = c("Protein Name" = "BestMatch")) %>%
              mutate(`Protein Name` = factor(`Protein Name`, levels = osmos),
                     Total.nM.in.sample.C = C.num *Total.nM.in.sample ,
                     Total.nM.in.sample.N = N.num * Total.nM.in.sample) %>%
              filter(!is.na(`Protein Name`)))

total.con.aves <- total.con.plot.dat %>%
  ungroup() %>%
  group_by(Experiment, Timepoint, type, `Protein Name`, SampID,
           `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
  summarise(Total.nM.in.sample.N.sd = sd(Total.nM.in.sample.N, na.rm = T),
            Total.nM.in.sample.N = mean(Total.nM.in.sample.N, na.rm = T),
            Total.nM.in.sample.C.sd = sd(Total.nM.in.sample.C, na.rm = T),
            Total.nM.in.sample.C = mean(Total.nM.in.sample.C, na.rm = T)) %>%
  full_join(.,total.con.plot.dat %>%
              ungroup() %>%
              group_by(Experiment, Timepoint, type, `Protein Name`, SampID,
                       `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
              summarise(Total.nM.in.sample.N.sd = sd(Total.nM.in.sample.N),
                        Total.nM.in.sample.N = mean(Total.nM.in.sample.N),
                        Total.nM.in.sample.C.sd = sd(Total.nM.in.sample.C),
                        Total.nM.in.sample.C = mean(Total.nM.in.sample.C)) %>%
              ungroup() %>%
              group_by(Experiment, Timepoint, type,  SampID,
                       `Time of day`, `Incubation time (hr)`, `Time since last sample (hr)`) %>%
              summarise(Total.N = sum(Total.nM.in.sample.N),
                        Total.C = sum(Total.nM.in.sample.C)) ) %>%
  mutate(proportion.N = Total.nM.in.sample.N/Total.N*100,
         proportion.C = Total.nM.in.sample.C/Total.C*100)
g.osmo.proportional.C <- ggplot(data = total.con.aves)  + 
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
            # aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  # geom_errorbar( aes(x = `Incubation time (hr)`, ymin = Total.nM.in.sample.C - Total.nM.in.sample.C.sd,
                     # ymax = Total.nM.in.sample.C + Total.nM.in.sample.C.sd))+
  geom_col(aes(x = `Incubation time (hr)`, y =proportion.C, fill = `Protein Name` ),
           width = 10, color="black") +
  facet_wrap(~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  ylab("Percent of C in Osmolytes") +
  theme_minimal()+
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6))
g.osmo.proportional.C
# ggsave(filename = "FateTidy/output/TargetedMetabolites/Osmolytes.Proportional.C.select.v.time.pdf",
       # height = 5, width = 6.5)

g.osmo.nm.c <- ggplot(data = total.con.aves)  + 
  # geom_rect(data = rects, fill = "grey", alpha = 0.5,
            # aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  # geom_errorbar( aes(x = `Incubation time (hr)`, ymin = Total.nM.in.sample.C - Total.nM.in.sample.C.sd,
  # ymax = Total.nM.in.sample.C + Total.nM.in.sample.C.sd))+
  geom_col(aes(x = `Incubation time (hr)`, y =Total.nM.in.sample.C, fill = `Protein Name` ),
           width = 10, color="black") +
  facet_wrap(~Experiment, scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  ylab("nM C in Osmolytes") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6))

g.osmo.nm.c
# ggsave(filename = "FateTidy/output/TargetedMetabolites/Osmolytes.nM.C.select.v.time.bars.pdf",
       # height = 5, width = 6.5)
legend.osmo <- get_legend(
  # create some space to the left of the legend
  g.osmo.nm.c + theme(legend.box.margin = margin(0, 0, 0, 12))
)

library(cowplot)
plot_grid(g.osmo.nm.c + theme(legend.position = "none"),
          g.osmo.proportional.C+ theme(legend.position = "none"),
          legend.osmo,
          ncol = 1,
          rel_heights = c(1,1,0.5))
ggsave(filename = "FateTidy/output/TargetedMetabolites/Osmolytes.v.time.bars.abs.and.rel.no.rect.pdf",
       height = 8, width = 6.5)

## osmolyte nM C ave per exp
View(total.con.aves %>%
       filter(Timepoint==96) %>%
  group_by(Experiment, `Protein Name`) %>%
  summarise(mean.nM.C = mean(Total.nM.in.sample.C),
            sd.nM.C = sd(Total.nM.in.sample.C),
            mean.nM.N= mean(Total.nM.in.sample.N)))
