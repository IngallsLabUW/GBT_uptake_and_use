
library(tidyverse)
library(ggplot2)
library(cowplot)
library(colorRamps)
library(RCurl)
library(RColorBrewer)


## load data -----
GBT.conc.in.all.forms <- read_csv("FateTidy/output/TargetedMetabolites/GBT_concentrations.csv")
natural.abundance.of.stds <- read_csv("FateTidy/RawData/Natural_abundance_of_Standards.csv")
HILIC.Pos.MID.dat <- read_csv("FateTidy/output/TargetedMetabolites/MID_mean_QC_Targeted_blkSubtractHILICPos.csv")
HILIC.Pos.MID.dat.to.use <- HILIC.Pos.MID.dat %>%
  filter(`Protein Name` == "Glycine betaine") %>%
  arrange(`Protein Name`)
Sample.log <- read_csv("FateTidy/RawData/Sample_Log_GBTFateExperiments.csv") %>%
  filter(grepl("Metab",`Sample Type`) |grepl("metab",`Sample Type`))
Sample.log.to.join <- Sample.log %>%
  filter(!grepl("Blk",`Sample ID`))%>%
  dplyr::select(Experiment, Timepoint,`Time of day`, rep, `Vol Filtered (L)`, `Incubation time (hr)`,
                `Time since last sample (hr)`) %>%
  unique()
## conc. plots -------
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

## line plots with natural abundance --------------

natural.abundance.from.g3 <- read_csv("FateTidy/RawData/Natural_abundance_HILICPos_from_G3.csv") %>%
  dplyr::select(-Fraction1, -iso.mz) %>%
  unique()


cmpd.name <-"Glycine betaine"
cmpd.base.mass <- min(HILIC.Pos.MID.dat.to.use$iso.mz[HILIC.Pos.MID.dat.to.use$`Protein Name` == cmpd.name])

nat.abu.this <- natural.abundance.of.stds %>%
  filter(Compound.Name == cmpd.name) %>%
  gather(isotopes,probs,  -X1, -Emperical.Formula, -Compound.Name, -C, -N.numb) %>%
  mutate(isotopes = ifelse(isotopes == "M.only.total","13C-0 15N-0",
                           ifelse(isotopes == "prob.of15N1","13C-0 15N-1",
                                  ifelse(isotopes == "prob.of.13C2", "13C-2 15N-0",
                                         "13C-1 15N-0"))),
         iso.mz = ifelse(isotopes == "13C-0 15N-0", cmpd.base.mass,
                         ifelse(isotopes == "13C-1 15N-0", cmpd.base.mass + 1.003354,
                                ifelse(isotopes == "13C-2 15N-0", (cmpd.base.mass + (2*1.003354)),
                                       (cmpd.base.mass + 0.997036)))),
         C.lab = ifelse(isotopes == "13C-1 15N-0",1,
                        ifelse(isotopes == "13C-2 15N-0",2,0)),
         N.lab = ifelse(isotopes == "13C-0 15N-1",1,0),
         MID.mean = probs/100,
         Experiment = "Natural Abundance Estimate",
         Timepoint = -1)
nat.abu.g3.this <- natural.abundance.from.g3 %>%
  filter(`Protein Name`==cmpd.name) %>%
  mutate(Experiment = ifelse(SampID == "KM1906S4C28D15","Fate1",
                             ifelse(SampID =="KM1906U12","Fate2",NA)),
         Timepoint= -2) %>%
  full_join(., HILIC.Pos.MID.dat.to.use %>%
              filter(`Protein Name` == cmpd.name) %>%
              dplyr::select(`Protein Name`, `Precursor Ion Name`, iso.mz) %>%
              unique())
cmpd.MID.aves <- HILIC.Pos.MID.dat.to.use %>%
  filter(`Protein Name` == cmpd.name) %>%
  filter(MID.mean > 0.005) %>%
  arrange(desc(iso.mz)) %>%
  full_join(nat.abu.this) %>% 
  full_join(nat.abu.g3.this)
num_toplot <- length(unique(cmpd.MID.aves$iso.mz))
cmpd.MID.factors.colors <- cmpd.MID.aves  %>%
  ungroup() %>%
  dplyr::select(iso.mz) %>%
  unique() %>%
  mutate(iso.mz.fact = factor(iso.mz, levels = unique(cmpd.MID.aves$iso.mz)),
         Colors = colorRampPalette(brewer.pal(9,"Set1"))(num_toplot+1)[1:num_toplot])

cmpd.MID.aves<- cmpd.MID.aves  %>%
  full_join(., cmpd.MID.factors.colors)
iso.mz.names <- cmpd.MID.aves %>%
  ungroup() %>%
  mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
  dplyr::select(names, iso.mz.fact) %>%
  unique()

names.vect <- iso.mz.names$iso.mz.fact
names(names.vect) <- iso.mz.names$names
myColors <- cmpd.MID.factors.colors$Colors

myColorScale <- scale_color_manual(values =  myColors,
                                   drop = F,
                                   name = "Isotopes",
                                   labels = names(names.vect))


f1.MLD.grid <- ggplot(cmpd.MID.aves %>%
                        filter(Experiment == "Fate1",
                               `Protein Name` == cmpd.name,
                               Timepoint != -2)) +
  geom_hline(data = cmpd.MID.aves %>%
               filter(Timepoint == -2,
                      Experiment == "Fate1"),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 1, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 1)+
  geom_errorbar(aes(x = Timepoint,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1)+
  geom_point(data = cmpd.MID.aves %>%
               filter(Experiment == "Natural Abundance Estimate"),
             aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 2, alpha = 0.5)  +
  theme_bw()+
  ylim(0,1)+
  myColorScale +
  theme(axis.text.x = element_text(angle = 90),legend.position="bottom") +
  ggtitle(paste("MID of", cmpd.name, "\n","in", "North" ))
f2.MLD.grid <- ggplot(cmpd.MID.aves%>%
                        filter(Experiment == "Fate2",
                               Timepoint != -2,
                               `Protein Name` == cmpd.name))  +
  geom_hline(data = cmpd.MID.aves %>%
               filter(Timepoint == -2,
                      Experiment == "Fate2"),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 1, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 1)+
  geom_errorbar(aes(x = Timepoint,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1)+ 
  geom_point(data = cmpd.MID.aves %>%
               filter(Experiment == "Natural Abundance Estimate"),
             aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 2, alpha = 0.5) +
  myColorScale +
  # facet_wrap(~Timepoint, ncol =1) +  
  theme_bw()+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90),legend.position="bottom") +
  ggtitle(paste("MID of", cmpd.name, "\n","in", "South"))
MLD.grid <- plot_grid(f1.MLD.grid, f2.MLD.grid, ncol = 1)
print(MLD.grid)

## all together plot ---------
cmpd.MID.aves <- cmpd.MID.aves%>%
  arrange(iso.mz) 
MID.plot.data.to.join <- cmpd.MID.aves%>%
  arrange(iso.mz) %>%
  mutate(iso.mz.fact = factor(iso.mz, levels = unique(cmpd.MID.aves$iso.mz))) %>%
  dplyr::select(MID.mean, iso.mz, iso.mz.fact,C.lab, N.lab, `Protein Name`, `Precursor Ion Name`, Timepoint, MID.sd, Experiment, SampID) %>%
  filter(Experiment != "Natural Abundance Estimate")
GBT.con.plot.to.join <- GBT.con.plot %>%
  dplyr::select(GBT.nM, GBT.nM.sd, `Time of day`, `Incubation time (hr)`, Timepoint, rep,
                `Protein Name`, `Precursor Ion Name`, SampID)
GBT.con.mean.to.join <- GBT.con.plot %>%
  group_by(SampID, `Precursor Ion Name`, Experiment, Timepoint, `Time of day`,
           `Incubation time (hr)`) %>%
  summarise(GBT.nM.sd = sd(GBT.nM),
            GBT.nM = mean(GBT.nM)) %>%
  full_join(MID.plot.data.to.join %>%
              dplyr::select(`Precursor Ion Name`,
                            iso.mz, iso.mz.fact) %>%
              unique() %>%
              filter(!is.na(`Precursor Ion Name`))) %>%
  full_join(MID.plot.data.to.join ) %>%
  dplyr::filter(!is.na(Experiment),
                !is.na(iso.mz))

num_toplot <- length(unique(MID.plot.data.to.join$iso.mz)[!is.na(unique(MID.plot.data.to.join$iso.mz))])
iso.mz.names <- MID.plot.data.to.join %>%
  filter(!is.na(iso.mz)) %>%
  ungroup() %>%
  mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
  dplyr::select(names, iso.mz.fact) %>%
  unique() %>%
  mutate(Colors = colorRampPalette(brewer.pal(9,"Set1"))(num_toplot+1)[1:num_toplot])

names.vect <- iso.mz.names$iso.mz.fact
names(names.vect) <- iso.mz.names$names
myColors <- iso.mz.names$Colors

myColorScale <- scale_color_manual(values =  myColors,
                                   drop = F,
                                   name = "Isotopes",
                                   labels = names(names.vect))

names.vect <- iso.mz.names$iso.mz.fact
names(names.vect) <- iso.mz.names$names

big.MID <- ggplot(GBT.con.mean.to.join %>%
                        filter(Timepoint != -2)) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_hline(data = GBT.con.mean.to.join %>%
               filter(Timepoint == -2),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 0.5, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = `Incubation time (hr)`,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 0.5)+
  geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1)+
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs), ncol = 1) +
  theme_bw()+
  ylim(0,1)+  
  labs(y = "MID") +
  myColorScale +
  theme(legend.position="bottom") 
  # ggtitle(paste("MID of", cmpd.name, "\n","in", "Fate1" ))
big.MID
MID.1 <- ggplot(GBT.con.mean.to.join %>%
                  filter(Timepoint != -2,
                         Experiment =="Fate1")) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_hline(data = GBT.con.mean.to.join %>%
               filter(Timepoint == -2,
                      Experiment =="Fate1"),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 0.5, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = `Incubation time (hr)`,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 0.5)+
  geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1) +
  theme_bw()+
  myColorScale +
  coord_fixed(ratio=500) +
  ylim(0,0.06)+
  theme(legend.position = "none",
        axis.title = element_blank())

MID.2 <- ggplot(GBT.con.mean.to.join %>%
                  filter(Timepoint != -2,
                         Experiment =="Fate2")) +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_hline(data = GBT.con.mean.to.join %>%
               filter(Timepoint == -2,
                      Experiment =="Fate2"),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 0.5, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = `Incubation time (hr)`,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 0.5)+
  geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1) +
  theme_bw()+
  myColorScale +
  coord_fixed(ratio=500) +
  ylim(0,0.06)+
  theme(legend.position = "none",
        axis.title = element_blank())

big.plot <- ggplot() +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar(data = GBT.con.mean.to.join,
                aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd,
                    color = iso.mz.fact))+
  geom_line(data = GBT.con.mean.to.join,
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = iso.mz.fact),
            size = 0.5)+
  facet_wrap(~Experiment, nrow = 2, 
             scales = "free_y",
             labeller = labeller(Experiment = Exp.labs)) +
  theme_bw()+
  myColorScale +
  theme(legend.position="bottom") +
  labs(y = "particulate GBT (nM)") 
  # scale_y_log10() 
big.plot


inset.1 <- ggplot() +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar(data = GBT.con.mean.to.join %>%
                  filter(Experiment == "Fate1"),
                aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd,
                    color = iso.mz.fact))+
  geom_line(data = GBT.con.mean.to.join%>%
              filter(Experiment == "Fate1"),
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = iso.mz.fact),
            size = 0.5)+
  theme_bw()+
  myColorScale +
  coord_fixed(ratio=400) +
  ylim(0,0.1)+
  theme(legend.position = "none",
        axis.title = element_blank())
inset.1

inset.2 <-  ggplot() +
  geom_rect(data = rects, fill = "grey", alpha = 0.5,
            aes(xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax))+
  geom_errorbar(data = GBT.con.mean.to.join %>%
                  filter(Experiment == "Fate2"),
                aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd,
                    color = iso.mz.fact))+
  geom_line(data = GBT.con.mean.to.join%>%
              filter(Experiment == "Fate2"),
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = iso.mz.fact),
            size = 0.5)+
  theme_bw()+
  myColorScale +
  coord_fixed(ratio=400) +
  ylim(0,0.1)+
  theme(legend.position = "none",
        axis.title = element_blank())
inset.2

gcombo.1 <- ggdraw(big.plot) +
  draw_plot(inset.1 ,x =  .5, y = .5, .45, .45)+
  draw_plot(inset.2, x = .5, y = .15, .45, .45)

gcombo.2 <- ggdraw(big.MID) +
  draw_plot(MID.1 ,x =  .35, y = .5, .55, .55)+
  draw_plot(MID.2, x = .35, y = .1, .55, .55)

all.gbt.plot <- plot_grid(gcombo.1, gcombo.2, ncol = 2)
all.gbt.plot
# ggsave(all.gbt.plot, filename = "FateTidy/output/TargetedMetabolites/GBT_Conc_MID.pdf",
       # width = 8, height = 8)


## rearainge big plot-----

new.big.plot <- ggplot()+
  geom_errorbar(data = GBT.con.mean.to.join %>%
                  mutate(Experiment = factor(Experiment, levels = c("Fate2", "Fate1"))),
                aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd,
                    color = iso.mz.fact))+
  geom_line(data = GBT.con.mean.to.join %>%
              mutate(Experiment = factor(Experiment, levels = c("Fate2", "Fate1"))),
            aes(x = `Incubation time (hr)`, y= GBT.nM, color = iso.mz.fact),
            size = 0.5) +
  facet_wrap(~Experiment, nrow = 1, 
             # scales = "",
             labeller = labeller(Experiment = Exp.labs)) +
  theme_bw()+
  myColorScale +
  theme(legend.position="bottom") +
  labs(y = "particulate GBT (nM)") 

 new.big.MID <- ggplot(GBT.con.mean.to.join %>%
                    filter(Timepoint != -2)%>%
          mutate(Experiment = factor(Experiment, levels = c("Fate2", "Fate1")))) +
  geom_hline(data = GBT.con.mean.to.join %>%
               filter(Timepoint == -2)%>%
               mutate(Experiment = factor(Experiment, levels = c("Fate2", "Fate1"))),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 0.5, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = `Incubation time (hr)`,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 0.5)+
  geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1)+
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs), nrow = 1) +
  theme_bw()+
  ylim(0,1)+  
  labs(y = "MID") +
  myColorScale +
  theme(legend.position="bottom") 
 
 # extract the legend from one of the plots
 legend <- get_legend(
   # create some space to the left of the legend
   new.big.MID + theme(legend.box.margin = margin(0, 0, 0, 12))
 )
 
 prow <- plot_grid(new.big.plot +
                     theme(legend.position="none") , 
                   new.big.MID +
                     theme(legend.position="none") ,
                  nrow = 2,
                  align = "hv")
 prow
 # add the legend to the row we made earlier. Give it one-third of 
 # the width of one plot (via rel_widths).
 plot_grid(prow, legend, nrow = 2, rel_heights = c(3, .4))
 
 ggsave(filename = "FateTidy/output/TargetedMetabolites/GBT_Conc_MID.new.order.no.inset.pdf",
        width = 8, height = 8)
 
 # ggsave(filename = "Fate/output/Figures/GBT_Conc_MID.new.order.no.inset.png",
        # width = 8, height = 6)
 
 ## save out just the particulate GBT concentrations without the MID stuff
 new.big.plot
 # ggsave(new.big.plot,
        # filename = "FateTidy/output/TargetedMetabolites/GBT_Conc_v_time.png", 
        # height = 3.6, width = 8)
 
 ## add insets into Conc and MID plot------
 inset.1.new =  ggplot() +
   # geom_point(data = GBT.con.plot %>%
                # filter(Experiment == "Fate1"),
              # aes(x = `Incubation time (hr)`, y = GBT.nM, color = `Precursor Ion Name`))+
   geom_errorbar(data = GBT.con.mean.to.join %>%
                   filter(Experiment == "Fate1"),
                 aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd, color = `Precursor Ion Name`))+
   geom_line(data = GBT.con.mean.to.join %>%
               filter(Experiment == "Fate1"),
             aes(x = `Incubation time (hr)`, y= GBT.nM, color = `Precursor Ion Name`))+
   # facet_wrap(~Experiment, nrow = 2, 
   #            scales = "free_y",
   #            labeller = labeller(Experiment = Exp.labs)) +
   labs(y = "particulate GBT (nM)") +
   coord_fixed(ratio=400) +
   theme_bw()+  
   myColorScale +

   ylim(0,0.1)+
   theme(legend.position = "none",
         axis.title = element_blank())
 inset.1.new
 
 inset.2.new <-  ggplot() +
   geom_errorbar(data = GBT.con.mean.to.join %>%
                   filter(Experiment == "Fate2"),
                 aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd,
                     color = iso.mz.fact))+
   geom_line(data = GBT.con.mean.to.join%>%
               filter(Experiment == "Fate2"),
             aes(x = `Incubation time (hr)`, y= GBT.nM, color = iso.mz.fact),
             size = 0.5)+
   theme_bw()+
   myColorScale +
   coord_fixed(ratio=400) +
   ylim(0,0.1)+
   theme(legend.position = "none",
         axis.title = element_blank())
inset.2.new
 
MID.1.new <- ggplot(GBT.con.mean.to.join %>%
                  filter(Timepoint != -2,
                         Experiment =="Fate1")) +
  geom_hline(data = GBT.con.mean.to.join %>%
               filter(Timepoint == -2,
                      Experiment =="Fate1"),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 0.5, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = `Incubation time (hr)`,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 0.5)+
  geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1) +
  theme_bw()+
  myColorScale +
  coord_fixed(ratio=700) +
  ylim(0,0.06)+
  theme(legend.position = "none",
        axis.title = element_blank())

MID.2.new <- ggplot(GBT.con.mean.to.join %>%
                  filter(Timepoint != -2,
                         Experiment =="Fate2")) +
  geom_hline(data = GBT.con.mean.to.join %>%
               filter(Timepoint == -2,
                      Experiment =="Fate2"),
             aes(#x = Timepoint,
               yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
             size = 0.5, #shape = 17,
             linetype ="dotted"
  ) +
  geom_line(aes(x = `Incubation time (hr)`,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
            size = 0.5)+
  geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                    ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
                alpha = 1) +
  theme_bw()+
  myColorScale +
  coord_fixed(ratio=700) +
  ylim(0,0.06)+
  theme(legend.position = "none",
        axis.title = element_blank())

big.new <- plot_grid(prow, legend, nrow = 2, rel_heights = c(3, .4))


  
big.new.insets <-  ggdraw(big.new) +
  draw_plot(inset.1.new ,x =  .7, y = .63, .23, .23)+
  draw_plot(inset.2.new, x = .15, y = .75, .23, .23)+
  draw_plot(MID.1.new ,x =  .7, y = .23, .23, .23)+
  draw_plot(MID.2.new, x = .15, y = .23, .23, .23)

big.new.insets

ggsave(filename = "FateTidy/output/TargetedMetabolites/GBT_Conc_MID.new.order.pdf",
       width = 8, height = 8)

 ## save out just the particulate GBT concentrations with new colors ----
 trim.iso.dat <- GBT.con.mean.to.join %>%
   filter(iso.mz.fact=="118.086804"|
            # iso.mz.fact==  "119.08384"|
            # iso.mz.fact== "119.090158"|
            iso.mz.fact== "124.10061"|
            # iso.mz.fact== "123.103574"|
            iso.mz.fact== "123.097256") %>%
   mutate(iso.mz.fact = factor(iso.mz.fact, 
                               levels = c("118.086804", 
                                          # "119.08384" ,
                                          # "119.090158", 
                                          "123.097256", 
                                          # "123.103574" ,
                                          "124.10061" )))
 iso.mz.names.new <- trim.iso.dat %>%
   ungroup() %>%
   mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
   dplyr::select(names, iso.mz.fact) %>%
   unique() %>%
   filter(!grepl("NA",names)) %>%
   arrange(names)
 
 names.vect.new <- iso.mz.names.new$iso.mz.fact
 names(names.vect.new) <- iso.mz.names.new$names
 myColorsnew <- c("black",rev(RColorBrewer::brewer.pal(3, "Dark2")))
 
 myColorScale.new <- scale_color_manual(values =  myColorsnew,
                                    drop = F,
                                    name = "Isotopes",
                                    labels = names(names.vect.new))
 
 
 new.big.plot.colors = ggplot()+
   geom_errorbar(data = trim.iso.dat %>%
                   mutate(Experiment = factor(Experiment, levels = c("Fate2", "Fate1"))),
                 aes(x = `Incubation time (hr)`, ymin = GBT.nM - GBT.nM.sd, ymax = GBT.nM + GBT.nM.sd,
                     color = iso.mz.fact))+
   geom_line(data = trim.iso.dat %>%
               mutate(Experiment = factor(Experiment, 
                                          levels = c("Fate2", "Fate1"))),
             aes(x = `Incubation time (hr)`, y= GBT.nM, 
                 color = iso.mz.fact),
             size = 0.5) +
   facet_wrap(~Experiment, nrow = 1, 
              # scales = "",
              labeller = labeller(Experiment = Exp.labs)) +
   theme_bw()+
   myColorScale.new +
   theme(legend.position="bottom",
         aspect.ratio = 0.5) +
   labs(y = "particulate GBT (nM)") 
 new.big.plot.colors
 
 ggsave(new.big.plot.colors,
        filename = "FateTidy/output/TargetedMetabolites/GBT_Conc_v_time_newColors.pdf",
        height = 3.6, width = 8)
 