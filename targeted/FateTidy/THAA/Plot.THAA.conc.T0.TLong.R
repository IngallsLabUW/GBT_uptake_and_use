## plot T0 v Tlong THAA concentrations

library(tidyverse)
library(ggplot2)
library(cowplot)
library(anytime)
library(ggpmisc) 
library(rlist)
library(tidyr)
theme_set(theme_cowplot())

## load data----------------
quantified.THAA <- read_csv("FateTidy/RawData/Quantified_THAA_inSeaWater.csv")
MFs <- read.csv("FateTidy/RawData/AA_Molecular_Formulas.csv")
## parse file names ---------------
plot.dat <- quantified.THAA %>%
  filter(Precursor.Ion.Name !="Methionine",
         Precursor.Ion.Name != "Histidine") %>%
  mutate(Replicate.name.to.sep = Replicate.Name,
         Environmental.Concentration.nM = ifelse(Environmental.Concentration.uM < 0, 0, 
                                                 Environmental.Concentration.uM*1e3)) %>%
  separate(Replicate.name.to.sep, 
           into = c("runDate","type","SampID","Replicate"), sep = "_") %>%
  mutate(Experiment = ifelse(grepl("Fate1",SampID),"Fate1",
                             ifelse(grepl("Fate2",SampID),"Fate2",NA)),
         sampid.to.split = SampID) %>% 
  separate(sampid.to.split, into = c("SampIDp1","Timepoint"), sep = "Fate") %>%
  mutate(Timepoint = (str_sub(Timepoint, 3,-1)),
         Timepoint = ifelse(grepl("long",Timepoint), "F" ,(Timepoint))) %>%
  dplyr::select(Precursor.Ion.Name, Replicate.Name,  type, SampID, Replicate,
                Experiment, Timepoint, Environmental.Concentration.nM )
plot.dat.means <- plot.dat %>%
  filter(!is.na(Experiment)) %>%
  group_by(Precursor.Ion.Name, Experiment, SampID, Timepoint, type) %>%
  summarise(`Average nM concentration` = mean(Environmental.Concentration.nM, na.rm = T),
            `sd nM concnetration` = sd(Environmental.Concentration.nM, na.rm = T)) %>%
  ungroup() %>%
  group_by(Precursor.Ion.Name, Experiment) %>%
  mutate(T0.Conc = `Average nM concentration`[Timepoint == "0"],
         T0.sd = `sd nM concnetration`[Timepoint == "0"],
         Percent.Change = (round(`Average nM concentration`/T0.Conc*100,0))) %>%
  left_join(MFs) %>%
  mutate(`Average nM concentration C` = `Average nM concentration`*Num.C,
         `Average nM concentration N` = `Average nM concentration`*Num.N)
## plot it -----------------

facet.exp.names <- c("Fate1" = "North",
                     "Fate2" = "South")


ggplot(plot.dat.means ) +
  geom_col(aes(x = Timepoint, y = `Average nM concentration`)) +
  geom_errorbar(aes(x = Timepoint, 
                    ymax = `Average nM concentration`+`sd nM concnetration`,
                    ymin = `Average nM concentration`-`sd nM concnetration`)) +
  facet_grid(Experiment~Precursor.Ion.Name) +
  theme_minimal()

ggplot(plot.dat %>% filter(!is.na(Experiment))) +
  geom_col(data = plot.dat.means,
           aes(x= Timepoint, y = `Average nM concentration`),
           alpha = 0.2)+
  geom_point(aes(x = Timepoint, y = Environmental.Concentration.nM, fill = Replicate),
             shape = 21, color = "black") +
  geom_text(data = plot.dat.means %>% filter(Timepoint == "F") %>%
              dplyr::select(Precursor.Ion.Name, Experiment, Percent.Change) %>%
              unique(),
            aes(x = 2, y = 15.5, label = Percent.Change ), size = 3) +
  scale_fill_manual(values = c("white","grey","black"))+
  facet_grid(Experiment~Precursor.Ion.Name, labeller = labeller(Experiment = facet.exp.names)) +
  theme_minimal() +
  theme(strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))

# ggsave("Fate/output/Figures/THAA_Conc_T0_Tlong.pdf",
       # width = 8, height = 6)

## total summed up-----
Exp.labs <-  c("North", "South")
names(Exp.labs) <- c("Fate1", "Fate2")

THAA.summed <- plot.dat %>%
  mutate(Experiment = factor(Experiment, levels = c("Fate2","Fate1"))) %>%
  left_join(MFs) %>%
  mutate(`Environmental.Concentration.nM.C` = Environmental.Concentration.nM *Num.C,
         `Environmental.Concentration.nM.N` = Environmental.Concentration.nM *Num.N) %>%
  group_by(Timepoint, Replicate, SampID, Replicate.Name, Experiment) %>%
  summarise(TotalTHAA.nM = sum(Environmental.Concentration.nM),
            TotalTHAA.nM.C = sum(Environmental.Concentration.nM.C),
            TotalTHAA.nM.N = sum(Environmental.Concentration.nM.N)) %>%
  filter(!is.na(Timepoint))
THAA.summed.mean <- THAA.summed %>%
  group_by(Timepoint, Experiment) %>%
  summarise(meanTHAA.nM = mean(TotalTHAA.nM),
            sd.THAA.nM = sd(TotalTHAA.nM),
            meanTHAA.nM.C = mean(TotalTHAA.nM.C),
            sd.THAA.nM.C = sd(TotalTHAA.nM.C),
            meanTHAA.nM.N = mean(TotalTHAA.nM.N),
            sd.THAA.nM.N = sd(TotalTHAA.nM.N))
ggplot() +
  geom_point(data = THAA.summed.mean,
             aes(y = meanTHAA.nM, x = Timepoint),
             size = 4,
             color = "grey")+
  geom_errorbar(data = THAA.summed.mean,
                aes(ymax = meanTHAA.nM + sd.THAA.nM,
                    ymin =meanTHAA.nM- sd.THAA.nM,
                    x = Timepoint),
                width = .1,
                color = "grey")+
  # geom_point(data = THAA.summed,
             # aes(y = TotalTHAA.nM, x = Timepoint, shape = Replicate)) +
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  ylab("THAA concentration (nM)") +
  scale_x_discrete(breaks = c("0","F"),
                   labels =c("Initial","Final")) 
# ggsave("Fate/output/Figures/THAA_total.pdf",
       # width = 4, height = 2.5)

ggplot() +
  geom_point(data = THAA.summed.mean,
             aes(y = meanTHAA.nM.C, x = Timepoint),
             size = 4,
             color = "grey")+
  geom_errorbar(data = THAA.summed.mean,
                aes(ymax = meanTHAA.nM.C + sd.THAA.nM.C,
                    ymin =meanTHAA.nM.C- sd.THAA.nM.C,
                    x = Timepoint),
                width = .1,
                color = "grey")+
   facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  ylab("THAA concentration (nM Carbon)") +
  scale_x_discrete(breaks = c("0","F"),
                   labels =c("Initial","Final")) 
ggsave("FateTidy/output/THAA/THAA_total_C.pdf",
       width = 4, height = 2.5)


ggplot() +
  geom_point(data = THAA.summed.mean,
             aes(y = meanTHAA.nM.N, x = Timepoint),
             size = 4,
             color = "grey")+
  geom_errorbar(data = THAA.summed.mean,
                aes(ymax = meanTHAA.nM.N + sd.THAA.nM.N,
                    ymin =meanTHAA.nM.N- sd.THAA.nM.N,
                    x = Timepoint),
                width = .1,
                color = "grey")+
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  ylab("THAA concentration (nM Nitrogen)") +
  scale_x_discrete(breaks = c("0","F"),
                   labels =c("Initial","Final")) 
# ggsave("Fate/output/Figures/THAA_total_N.pdf",
       # width = 4, height = 2.5)
