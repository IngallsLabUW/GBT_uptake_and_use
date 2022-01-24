## quantify GBT

library(tidyverse)
library(ggplot2)
library(cowplot)
##read in -IS Normalized data --------------------
dat <- read_csv("FateTidy/output/TargetedMetabolites/HILICPos_QC_skyline_IS_data_BMISed.csv")


## get GBT only ------
gbt.dat <- dat %>%
  filter(MassFeature == "Glycine betaine_13C-0 15N-0")

## get conversion to nM units from the QC samples ---------------
qc.dat <- dat %>%
  filter(MassFeature == "Betaine, 13C5-15N") %>%
  filter(SampID == "QC2")

ave.gbt.qc.area = mean(qc.dat$Adjusted_Area)
sd.gbt.qc.area = sd(qc.dat$Adjusted_Area)
rsd.gbt.qc.area = sd.gbt.qc.area/ave.gbt.qc.area

qc.gbt.is.con = 500 ## because we diluted the QC samples 1:1 with water!!!

gbt.nM.per.pk.area.ave = qc.gbt.is.con/ave.gbt.qc.area
gbt.nM.per.pk.area.sd = gbt.nM.per.pk.area.ave * rsd.gbt.qc.area

test.back.convert = qc.dat %>%
  mutate(GBT.nM = Adjusted_Area * gbt.nM.per.pk.area.ave,
         GBT.nM.sd = GBT.nM * rsd.gbt.qc.area )
  
## convert to nM -----------------
gbt.dat.nM <- gbt.dat %>%
  mutate(GBT.nM.in.vial = Adjusted_Area * gbt.nM.per.pk.area.ave,
         GBT.nM.sd = GBT.nM.in.vial * rsd.gbt.qc.area )

## checking the conversion against the standards in water as a test, 
## there should be 4 uM = 4000 nM there, of course it is a different matrix, 
## but still, and indeed it does look ok with estimates around 3000-3800 nM.

## write base peak concentration -----------------
GBT.con.to.write = gbt.dat.nM %>%
  dplyr::select(SampID, type, runDate, replicate, MassFeature, Adjusted_Area, GBT.nM.in.vial, GBT.nM.sd)

# write_csv(GBT.con.to.write, "Fate/output/skyline/GBT_base_peak_concentration.csv")

## get total concentration based on MID -------------------
MID <- read_csv("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_blkSubtractHILICPos.csv")
Sample.log <- read_csv("FateTidy/RawData/Sample_Log_GBTFateExperiments.csv") %>%
  filter(grepl("Metab",`Sample Type`) |grepl("metab",`Sample Type`))
Sample.log.to.join <- Sample.log %>%
  filter(!grepl("Blk",`Sample ID`))%>%
  dplyr::select(Experiment, Timepoint,`Time of day`, rep, `Vol Filtered (L)`, `Incubation time (hr)`,
                `Time since last sample (hr)`) %>%
  unique()

GBT.MID <- MID %>%
  filter(`Protein Name`=="Glycine betaine") %>%
  full_join(Sample.log.to.join)

GBT.con.to.join <- GBT.con.to.write %>%
  mutate(SampID = str_replace(SampID,"-IS","")) 

GBT.total.nM.each.samp = left_join(GBT.con.to.join, GBT.MID,
                                                  by= c("MassFeature" = "Precursor Ion Name",
                                                        "SampID",
                                                        "runDate",
                                                        "type",
                                                        "replicate"="rep")) %>%
  mutate(Area.per.nM = Area.filled.w.background/GBT.nM.in.vial,
         Area.per.nM.sd = Area.per.nM * (GBT.nM.sd/GBT.nM.in.vial),
         total.GBT.nM.in.vial = all.area/Area.per.nM,
         total.GBT.nM.sd = total.GBT.nM.in.vial * (Area.per.nM.sd/Area.per.nM),
         total.GBT.nM.in.sample = total.GBT.nM.in.vial *  (37.078/30) *(400e-6/`Vol Filtered (L)`), ## 30 ratio because we spiked in water or standards after aliquoting 30 um. 400e-6 because that was the original reconsitution volume
         total.GBT.nM.in.sample.sd = total.GBT.nM.in.sample * (total.GBT.nM.sd/total.GBT.nM.in.vial)) %>%
  dplyr::select(SampID, type, replicate, runDate , total.GBT.nM.in.sample, total.GBT.nM.in.sample.sd) %>%
  dplyr::rename(rep = replicate) %>%
  filter(!is.na(total.GBT.nM.in.sample))%>%
  arrange(SampID, rep)

GBT.conc.in.all.forms <- full_join(GBT.MID, GBT.total.nM.each.samp)  %>%
  mutate(GBT.nM = total.GBT.nM.in.sample * MID,
         GBT.nM.sd = GBT.nM * (total.GBT.nM.in.sample.sd/total.GBT.nM.in.sample))

write_csv(GBT.conc.in.all.forms, "FateTidy/output/TargetedMetabolites/GBT_concentrations.csv")

## calculate ave initial concentrations-----
GBT.conc.in.all.forms.aves <- GBT.conc.in.all.forms %>%
  group_by(`Protein Name`, `Precursor Ion Name`,
           Experiment, Timepoint) %>%
  summarise(ave.GBT.nm = mean(GBT.nM),
            sd.GBT.nm = ave.GBT.nm * (sqrt(sum((GBT.nM.sd/GBT.nM)^2))))
write_csv(GBT.conc.in.all.forms.aves, "Fate/output/skyline/GBT_concentrations_aves.csv")

## plot -----------
rects <- tibble(ymin = -Inf, ymax = Inf, xmin = c(11,11, 11+24, 11+24, 11+48, 11+48, 11+24*3, 11+24*3),
                xmax = c(23, 23, 23+24, 23+24, 23+24*2, 23+24*2, 23+24*3, 23+24*3),
                Experiment = rep(c("Fate1","Fate2"),4))
GBT.con.plot <- GBT.conc.in.all.forms %>%
  left_join(Sample.log.to.join) #%>%
  # filter(GBT.nM > 1e-5)
GBT.con.mean <- GBT.con.plot %>%
  group_by(SampID, `Precursor Ion Name`, Experiment, Timepoint, `Time of day`,
           `Incubation time (hr)`) %>%
  summarise(GBT.nM = mean(GBT.nM))
Exp.labs <-  c("North", "South")
names(Exp.labs) <- c("Fate1", "Fate2")
big.plot <- ggplot() +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_point(data = GBT.con.plot,
             aes(x = `Incubation time (hr)`, y = GBT.nM, color = `Precursor Ion Name`))+
  geom_errorbar(data = GBT.con.plot,
             aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd, color = `Precursor Ion Name`))+
  geom_line(data = GBT.con.mean,
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = `Precursor Ion Name`))+
  facet_wrap(~Experiment, nrow = 2, 
             scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  labs(y = "particulate GBT (nM)",
       color = "Isotopologue") +
  # scale_y_log10() +
  theme_bw()

big.plot


# ggsave(big.plot, filename = "Fate/output/Figures/GBT_Concentration_Over_Time.pdf", width = 8, height = 6)

inset.1 <- ggplot() +
  geom_rect(data = rects %>%
              filter(Experiment == "Fate1"), fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_point(data = GBT.con.plot %>%
               filter(Experiment == "Fate1"),
             aes(x = `Incubation time (hr)`, y = GBT.nM, color = `Precursor Ion Name`))+
  geom_errorbar(data = GBT.con.plot %>%
                  filter(Experiment == "Fate1"),
                aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd, color = `Precursor Ion Name`))+
  geom_line(data = GBT.con.mean %>%
              filter(Experiment == "Fate1"),
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = `Precursor Ion Name`))+
  # facet_wrap(~Experiment, nrow = 2, 
  #            scales = "free_y",
  #            labeller = labeller(Experiment = Exp.labs)) +
  labs(y = "particulate GBT (nM)") +
  coord_fixed(ratio=400) +
  theme_bw()+
  ylim(0,0.1)+
  theme(legend.position = "none",
        axis.title = element_blank())
inset.1

inset.2 <- ggplot() +
  geom_rect(data = rects %>%
              filter(Experiment == "Fate2"), fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_point(data = GBT.con.plot %>%
               filter(Experiment == "Fate2"),
             aes(x = `Incubation time (hr)`, y = GBT.nM, color = `Precursor Ion Name`))+
  geom_errorbar(data = GBT.con.plot %>%
                  filter(Experiment == "Fate2"),
                aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd, color = `Precursor Ion Name`))+
  geom_line(data = GBT.con.mean %>%
              filter(Experiment == "Fate2"),
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = `Precursor Ion Name`))+
  # facet_wrap(~Experiment, nrow = 2, 
  #            scales = "free_y",
  #            labeller = labeller(Experiment = Exp.labs)) +
  labs(y = "particulate GBT (nM)") +
  coord_fixed(ratio=400) +
  theme_bw()+
  ylim(0,0.1)+
  theme(legend.position = "none",
        axis.title = element_blank())
inset.2

ggdraw(big.plot) +
  draw_plot(inset.1 ,x =  .45, y = .567, .25, .25)+
  draw_plot(inset.2, x = .45, y = .15, .25, .25)
ggsave(filename = "FateTidy/output/TargetedMetabolites/GBT_Concentration_Over_Time_withInset.pdf", width = 8, height = 6)


## uptake rate of fully labeled GBT -----


previou.time.point = GBT.con.mean %>%
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
previous.time.point.conc = GBT.con.mean %>%
  ungroup() %>%
  dplyr::rename(previous.con.nM = GBT.nM,
                previous.time.point.ID = SampID) %>%
  dplyr::select(previous.con.nM, previous.time.point.ID, `Precursor Ion Name`) %>%
  full_join(.,previou.time.point) %>%
  filter(!is.na(SampID)) %>%
  mutate(previous.con.nM = ifelse(is.na(previous.con.nM),0,previous.con.nM))
Uptake.rate.of.full.labeled <- previous.time.point.conc %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-5 15N-1") %>%
  full_join(Sample.log.to.join) %>%
  mutate(uptake.overall.nM.per.h = GBT.nM/`Incubation time (hr)`,
         uptake.since.last.point.nM.per.h = (GBT.nM-previous.con.nM)/`Time since last sample (hr)`,
         uptake.since.last.point.nM.per.h = ifelse(is.infinite(uptake.since.last.point.nM.per.h),
                                                   uptake.overall.nM.per.h,
                                                   uptake.since.last.point.nM.per.h)) 
uptake.rate.plot <- ggplot(Uptake.rate.of.full.labeled) +
  geom_line(aes(y = uptake.since.last.point.nM.per.h, x =`Incubation time (hr)`), color = "blue") +
  geom_line(aes(y = uptake.overall.nM.per.h, x =`Incubation time (hr)`), color = "black") +
  facet_wrap(~Experiment, nrow = 2, labeller = labeller(Experiment = Exp.labs)) +
  theme_bw()
uptake.rate.plot
# ggsave(filename = "Fate/output/skyline/Uptake.Rate.of.Fully.Labeled.GBT.pdf",
       # width = 8, height = 6)

## get rate of unlabeled gbt drawdown in between the first two time points -------
previou.time.point = GBT.con.mean %>%
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
previous.time.point.conc = GBT.con.mean %>%
  ungroup() %>%
  dplyr::rename(previous.con.nM = GBT.nM,
                previous.time.point.ID = SampID) %>%
  dplyr::select(previous.con.nM, previous.time.point.ID, `Precursor Ion Name`) %>%
  full_join(.,previou.time.point) %>%
  filter(!is.na(SampID)) %>%
  mutate(previous.con.nM = ifelse(is.na(previous.con.nM),0,previous.con.nM))
unlabeled.diffs = previous.time.point.conc %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-0 15N-0") %>%
  full_join(Sample.log.to.join) %>%
  mutate(Conc.diff = GBT.nM - previous.con.nM,
         uptake.since.last.point.nM.per.h = Conc.diff/`Time since last sample (hr)`)%>%
  unique()
uptake.rate.plot <- ggplot(unlabeled.diffs) +
  geom_line(aes(y = uptake.since.last.point.nM.per.h, x =`Incubation time (hr)`), color = "blue") +
  facet_wrap(~Experiment, nrow = 2) +
  theme_bw()
uptake.rate.plot

## test to see if time points are actually different
unlableled.GBT.conc <- GBT.conc.in.all.forms %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-0 15N-0")
t.test(x = unlableled.GBT.conc$GBT.nM[unlableled.GBT.conc$Experiment=="Fate1" & unlableled.GBT.conc$Timepoint == 0],
          y = unlableled.GBT.conc$GBT.nM[unlableled.GBT.conc$Experiment=="Fate1" & unlableled.GBT.conc$Timepoint == 4])
## unlabeled GBT concentration in fate 1 is not significantly different between T0 and T4
t.test(x = unlableled.GBT.conc$GBT.nM[unlableled.GBT.conc$Experiment=="Fate2" & unlableled.GBT.conc$Timepoint == 0],
       y = unlableled.GBT.conc$GBT.nM[unlableled.GBT.conc$Experiment=="Fate2" & unlableled.GBT.conc$Timepoint == 5])
## unlabeled GBT concentration in fate 2 is significantly different between T0 and T5


##  total uptake of GBT over first hours (nM/hr) = T4 overall accumulation (nM/hr) + abs( aveT0,T4(GBT-lab Con/GBT-nolab Con) * T4-GBT no.lab Loss (nM/hr))
conc.ratio =  GBT.con.mean %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-0 15N-0"|
           `Precursor Ion Name` == "Glycine betaine_13C-5 15N-1",
         Timepoint == 0|Timepoint == 4|Timepoint == 5) %>%
  dplyr::select(`Precursor Ion Name`, GBT.nM, Experiment, Timepoint, `Incubation time (hr)`) %>%
  spread(`Precursor Ion Name`, GBT.nM) %>%
  mutate(Lab.Unlab.Ratio = `Glycine betaine_13C-5 15N-1`/`Glycine betaine_13C-0 15N-0`) %>%
  unique() %>%
  group_by(Experiment) %>%
  summarise(Lab.Unlab.Ratio.2 = ave(Lab.Unlab.Ratio)) %>%
  unique()

T1overall.accumulation  = Uptake.rate.of.full.labeled %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-5 15N-1",
         Timepoint == 4|Timepoint == 5) %>%
  dplyr::select(Experiment, uptake.overall.nM.per.h) %>%
  unique()

T1overall.loss  = unlabeled.diffs %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-0 15N-0",
         Timepoint == 4|Timepoint == 5) %>%
  dplyr::select(Experiment, uptake.since.last.point.nM.per.h) %>%
  unique()

Total.uptake.first.few.hrs.method1 = full_join(conc.ratio, T1overall.accumulation) %>%
  full_join(T1overall.loss) %>%
  mutate(Tot.uptake.method1 = uptake.overall.nM.per.h + abs(Lab.Unlab.Ratio.2 * uptake.since.last.point.nM.per.h))

write_csv(Total.uptake.first.few.hrs.method1, "Fate/output/skyline/Total.Labeled.GBT.Uptake.Calculated.Method1.csv")

## turnover time (hr) = T0GBT nolab.conc (nM) / T4-GBT no.lab Loss (nM/hr)

unlableled.GBT.conc.t0 <- GBT.conc.in.all.forms %>%
  filter(`Precursor Ion Name` == "Glycine betaine_13C-0 15N-0",
         Timepoint == 0) %>%
  unique() %>%
  group_by(Experiment, `Precursor Ion Name`) %>%
  summarise(GBT.nM.sd = sqrt(sum((GBT.nM.sd/GBT.nM)^2)),
            GBT.nM = mean(GBT.nM))

turnovertime <- full_join(unlableled.GBT.conc.t0, T1overall.loss) %>%
  mutate(turnover.time = GBT.nM / abs(uptake.since.last.point.nM.per.h))
write_csv(turnovertime, "FateTidy/output/TargetedMetabolites/Turnovertime.particulate.GBT.csv")
