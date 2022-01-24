## MID line plots
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(wesanderson)
library(colorRamps)
library(RCurl)
## read in data ---------------
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
  dplyr::select(Compound.and.Iso.name, iso.mz, RT..min., Column, Fraction1) #%>%
  # filter(Fraction1 == "HILICPos")

MID <- read_csv("FateTidy/output/TargetedMetabolites/MID_blkSubtract_allFractions.csv")
estimated.nat.abu <- read_csv("FateTidy/RawData/Natural_abundance_of_Standards.csv")
g3.nat.abu <- read_csv("FateTidy/RawData/Natural_abundance_HILICPos_from_G3.csv",
                       col_types = cols(
                         type = col_character(),
                         SampID = col_character(),
                         runDate = col_double(),
                         N.lab = col_double(),
                         C.lab = col_double(),
                         `Precursor Ion Name` = col_character(),
                         `Protein Name` = col_character(),
                         Fraction1 = col_character(),
                         iso.mz = col_double(),
                         MID.mean = col_double(),
                         MID.sd = col_double(),
                         Area.filled.w.background.mean = col_double(),
                         Area.filled.w.background.sd = col_double()
                       ))
Sample.log <- read_csv("FateTidy/RawData/Sample_Log_GBTFateExperiments.csv") %>%
  filter(grepl("Metab",`Sample Type`) |grepl("metab",`Sample Type`))
Sample.log.to.join <- Sample.log %>%
  filter(!grepl("Blk",`Sample ID`))%>%
  dplyr::select(Experiment, Timepoint,`Time of day`, rep, `Vol Filtered (L)`, `Incubation time (hr)`,
                `Time since last sample (hr)`) %>%
  unique()

## list compounds ----------
Good.but.no.change <- c("Betonicine",
                        "L−Alanine",
                        "beta−Alanine",
                        "L−Proline",
                        "Proline betaine",
                        "L-Methionine")
new.cmpds <- c(
                        "DMSP",
                        "Gonyol",
                        "Glucosylglycerol",
                        "(R)−2,3−Dihydroxypropane−1−sulfonate",
                        "Isethionic acid",
                        "Butyrylcarnitine",
                        "L−Phenylalanine",
                        "L−Tryptophan",
                        "Methylthioadenosine")

Big.change <- c("Sarcosine",
                "Creatine",
                "Trimethylamine N−oxide",
                "Choline",
                "Glycerophosphocholine",
                "Carnitine",
                "(3−Carboxypropyl)trimethylammonium")

smaller.change <- c("L−Glutamic acid",
                    "L−Glutamine",
                    "L−Arginine",
                    "beta−Alaninebetaine",
                    "Taurine",
                    "Adenine",
                    "Adenosine",
                    "Deoxyadenosine",
                    "Guanine",
                    "Guanosine",
                    "Cytosine",
                    "Trigonelline",
                    "Homarine")

## wrangle for plot ---------

cmpds.for.plot <- c(Big.change, smaller.change, new.cmpds)
Dat.to.plot <- full_join(MID, Sample.log.to.join) %>%
  left_join(iso.poss.use) %>%
  unique() %>%
  filter(`Protein Name` %in% cmpds.for.plot) %>%
  mutate(Experiment = factor(Experiment, levels = c("Fate2","Fate1")))

## get color pallett -----
num_toplot <- Dat.to.plot%>%
  dplyr::select(N.lab, C.lab) %>%
  unique() %>%
  mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
  filter(!(N.lab == 0 & C.lab == 0)) %>%
  .$names %>%
  length()

color.pal.test = Dat.to.plot %>%
  dplyr::select(N.lab, C.lab) %>%
  unique() %>%
  filter(!(N.lab == 0 & C.lab == 0)) %>%
  mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
  arrange(names) %>%
  mutate(names.fact = factor(names, levels = .$names)) %>%
  # arrange(N.lab, C.lab) %>%
  mutate(Colors = colorRampPalette(brewer.pal(12,"Paired"))(num_toplot+1)[c(1:num_toplot-1, num_toplot+1)])

names.vect <- color.pal.test$names.fact
names(names.vect) <- color.pal.test$names

myColors <- color.pal.test$Colors

myColorScale <- scale_color_manual(values =  myColors,
                                   drop = F,
                                   labels = names(names.vect))

## make plot --------
make.MID.line.plot <- function(cmpd.name ="Choline",
                               scales.pick ="free"){
  cmpd.base.mass <- min(Dat.to.plot$iso.mz[Dat.to.plot$`Protein Name` == cmpd.name])
  nat.abu.g3.this <- g3.nat.abu %>%
    filter(`Protein Name`==cmpd.name,
           MID.mean != 0,
           !is.na(MID.mean)) %>%
    mutate(Experiment = ifelse(SampID == "KM1906S4C28D15","Fate1",
                               ifelse(SampID =="KM1906U12","Fate2","None")),
           Timepoint= -2) %>%
    dplyr::select(-Fraction1, -iso.mz) %>%
    left_join(., Dat.to.plot %>%
                filter(`Protein Name` == cmpd.name) %>%
                dplyr::select(`Protein Name`, `Precursor Ion Name`, iso.mz) %>%
                unique()) %>%
    filter(!(C.lab ==0 & N.lab == 0)) %>%
    mutate(Experiment = as.character(Experiment))
  
  if(cmpd.name == "Choline"){
    cmpd.MID.aves <- Dat.to.plot %>%
      filter(`Protein Name` == cmpd.name) %>%
      group_by(`Protein Name`, `Precursor Ion Name`, iso.mz,
               Timepoint, `Time of day`, `Incubation time (hr)`,
               SampID, C.lab, N.lab, Experiment) %>%
      summarise(MID.mean = mean(MID, na.rm = T),
                MID.sd = sd(MID, na.rm = T)) %>%
      arrange(desc(iso.mz)) %>%
      filter(!(C.lab ==0 & N.lab == 0)) %>%
      # full_join(nat.abu.g3.this) %>%
      mutate(`Incubation time (hr)` = ifelse(Timepoint == -1,
                                             -1, ifelse(Timepoint== -2, 
                                                        0,`Incubation time (hr)`))) %>%
      filter(MID.mean != 0,
             !is.na(MID.mean)) %>%
      arrange(iso.mz)
  } else {
    cmpd.MID.aves <- Dat.to.plot %>%
      filter(`Protein Name` == cmpd.name) %>%
      group_by(`Protein Name`, `Precursor Ion Name`, iso.mz,
               Timepoint, `Time of day`, `Incubation time (hr)`,
               SampID, C.lab, N.lab, Experiment) %>%
      summarise(MID.mean = mean(MID, na.rm = T),
                MID.sd = sd(MID, na.rm = T)) %>%
      arrange(desc(iso.mz)) %>%
      filter(!(C.lab ==0 & N.lab == 0)) %>%
      full_join(nat.abu.g3.this) %>%
      mutate(`Incubation time (hr)` = ifelse(Timepoint == -1,
                                             -1, ifelse(Timepoint== -2, 
                                                        0,`Incubation time (hr)`))) %>%
      filter(MID.mean != 0,
             !is.na(MID.mean)) %>%
      arrange(iso.mz)
  }
  
  num_toplot <- length(unique(cmpd.MID.aves$iso.mz))
  cmpd.MID.factors.colors <- cmpd.MID.aves  %>%
    ungroup() %>%
    dplyr::select(iso.mz, C.lab, N.lab) %>%
    unique() %>%
    mutate(iso.mz.fact = factor(iso.mz, levels = unique(cmpd.MID.aves$iso.mz))) %>%
    full_join(color.pal.test)
  cmpd.MID.aves <- full_join(cmpd.MID.aves, 
                             cmpd.MID.factors.colors) %>%
    mutate(MID.mean = MID.mean * 100,
           MID.sd = MID.sd * 100) %>%
    arrange(names.fact) %>%
    mutate(Experiment = factor(Experiment, levels = c("Fate2","Fate1")))
  # iso.mz.names <- cmpd.MID.aves %>%
  #   ungroup() %>%
  #   mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
  #   dplyr::select(names, iso.mz.fact) %>%
  #   unique()
  # names.vect <- iso.mz.names$iso.mz.fact
  # names(names.vect) <- iso.mz.names$names
  
  # myColors <- cmpd.MID.factors.colors$Colors
  # 
  # myColorScale <- scale_color_manual(values =  myColors,
  #                                    drop = F,
  #                                    labels = names(names.vect))
  Exp.labs <-  c("Exp 1", "Exp 2")
  names(Exp.labs) <- c("Fate1", "Fate2")
  
  g <- ggplot(cmpd.MID.aves %>%
                filter(`Protein Name` == cmpd.name,
                       Timepoint != -2)) +
    geom_hline(data = cmpd.MID.aves %>%
                 filter(Timepoint == -2),
               aes(#x = Timepoint,
                 yintercept = MID.mean, group = names.fact, color = names.fact),
               size = 0.5, #shape = 17,
               linetype ="dotted"
    ) +
    geom_line(aes(x = `Incubation time (hr)`,   
                  y = MID.mean, group = names.fact, color = names.fact),
              size = 0.5) +
    geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                      ymax = MID.mean + MID.sd, group = names.fact, color = names.fact),
                  alpha = 1) +
    facet_wrap(~Experiment, scales = scales.pick,
               labeller = labeller(Experiment = Exp.labs)) + 
    ggtitle(cmpd.name) +
    ylab("MID") +
    myColorScale +
    theme_cowplot() +
    scale_y_continuous(
      # don't expand y scale at the lower end
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 6, margin = margin(l = -5)),
          axis.title.y = element_text(size = 7, margin = margin(t = 0, r = -10, b = 0, l = 0)),
          axis.text = element_text(size = 6),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_text(size = 9),
          plot.margin = unit(c(0, 0, 0, 0), "cm")
          ) 
  return(g)
}

## I need to sort out the colors for these new ones still------
new1 <- make.MID.line.plot(cmpd.name = new.cmpds[1], scales.pick = "free_x")
new2 <- make.MID.line.plot(cmpd.name = new.cmpds[2], scales.pick = "free_x")
new3 <- make.MID.line.plot(cmpd.name = new.cmpds[3], scales.pick = "free_x")
new4 <- make.MID.line.plot(cmpd.name = new.cmpds[4], scales.pick = "free_x")
new5 <- make.MID.line.plot(cmpd.name = new.cmpds[5], scales.pick = "free_x")
new6 <- make.MID.line.plot(cmpd.name = new.cmpds[6], scales.pick = "free_x")
new7 <- make.MID.line.plot(cmpd.name = new.cmpds[7], scales.pick = "free_x")
new8 <- make.MID.line.plot(cmpd.name = new.cmpds[8], scales.pick = "free_x")

new.change.MID.plot <- plot_grid(new1+ theme(legend.position="none"),
                                   new2+ theme(legend.position="none"),
                                   new3+ theme(legend.position="none"),
                                   new4+ theme(legend.position="none"),
                                   new5+ theme(legend.position="none"),
                                   new6+ theme(legend.position="none"),
                                   new7+ theme(legend.position="none"),
                                   new8+ theme(legend.position="none"),
                                   ncol = 2)
new.change.MID.plot
## 'old' versions of figures with out new compounds ----
g1 <- make.MID.line.plot(cmpd.name = smaller.change[1], scales.pic = "free_x")
g2 <- make.MID.line.plot(cmpd.name = smaller.change[2], scales.pic = "free_x")
g3 <- make.MID.line.plot(cmpd.name = smaller.change[3], scales.pic = "free_x")
g4 <- make.MID.line.plot(cmpd.name = smaller.change[4], scales.pic = "free_x")
g5 <- make.MID.line.plot(cmpd.name = smaller.change[5], scales.pic = "free_x")
g6 <- make.MID.line.plot(cmpd.name = smaller.change[6], scales.pic = "free_x")
g7 <- make.MID.line.plot(cmpd.name = smaller.change[7], scales.pic = "free_x")
g8 <- make.MID.line.plot(cmpd.name = smaller.change[8], scales.pic = "free_x")
g9 <- make.MID.line.plot(cmpd.name = smaller.change[9], scales.pic = "free_x")
g10 <- make.MID.line.plot(cmpd.name = smaller.change[10], scales.pic = "free_x")
g11 <- make.MID.line.plot(cmpd.name = smaller.change[11], scales.pic = "free_x")
g12 <- make.MID.line.plot(cmpd.name = smaller.change[12], scales.pic = "free_x")
g13 <- make.MID.line.plot(cmpd.name = smaller.change[13], scales.pic = "free_x")


small.change.MID.plot <- plot_grid(g1+ theme(legend.position="none"),
                                   g2+ theme(legend.position="none"),
                                   g3+ theme(legend.position="none"),
                                   g4+ theme(legend.position="none"),
                                   g5+ theme(legend.position="none"),
                                   g6+ theme(legend.position="none"),
                                   g7+ theme(legend.position="none"),
                                   g8+ theme(legend.position="none"),
                                   g9+ theme(legend.position="none"),
                                   g10+ theme(legend.position="none"),
                                   g11+ theme(legend.position="none"),
                                   g12+ theme(legend.position="none"),
                                   g13+ theme(legend.position="none"),
                                   ncol = 2)
small.change.MID.plot


new.order.small.change.plot <- plot_grid(g1+ theme(legend.position="none"),
                                         g2+ theme(legend.position="none"),
                                         g3+ theme(legend.position="none"),
                                         g5+ theme(legend.position="none"),
                                         g4+ theme(legend.position="none"),
                                         new1+ theme(legend.position="none"),
                                         new2+ theme(legend.position="none"),
                                         new3+ theme(legend.position="none"),
                                         g12+ theme(legend.position="none"),
                                         g13+ theme(legend.position="none"),
                                         ncol = 2)
new.order.small.change.plot


nucleic.plot <- plot_grid(g6+ theme(legend.position="none"),
                          g7+ theme(legend.position="none"),
                          g8+ theme(legend.position="none"),
                          g9+ theme(legend.position="none"),
                          g10+ theme(legend.position="none"),
                          g11+ theme(legend.position="none"),
                          align = 'hv',
                          ncol=1)
nucleic.plot
# extract the legend from one of the plots and make it horizontal
legend <- get_legend(
  # create some space to the left of the legend
  g1+
    guides(color = guide_legend(nrow = 4)) +
    theme(legend.position = "bottom")
)


newest.order.small.change.plot <- plot_grid(g1+ theme(legend.position="none",
                                                      axis.text.x = element_blank()), 
                                            g2+ theme(legend.position="none") + theme(axis.title.y = element_blank(),
                                                                                      axis.text.x = element_blank()),
                                            g6+ theme(legend.position="none")+ theme(axis.title.y = element_blank(),
                                                                                     axis.text.x = element_blank()),
                                            g7 + theme(legend.position = "none")+ theme(axis.title.y = element_blank(),
                                                                                        axis.text.x = element_blank()),
                                            g3+ theme(legend.position="none",
                                                      axis.text.x = element_blank()), 
                                            # g5+ theme(legend.position="none"),
                                            g4+ theme(legend.position="none")+ theme(axis.title.y = element_blank(),
                                                                                     axis.text.x = element_blank()),
                                            g8+ theme(legend.position="none")+ theme(axis.title.y = element_blank(),
                                                                                     axis.text.x = element_blank()),
                                            g9+ theme(legend.position="none")+ theme(axis.title.y = element_blank(),
                                                                                     axis.text.x = element_blank()),
                                            new1+ theme(legend.position="none",
                                                        axis.text.x = element_blank()),  
                                            new2+ theme(legend.position="none")+ theme(axis.title.y = element_blank(),
                                                                                       axis.text.x = element_blank()),
                                            g10+ theme(legend.position="none")+ theme(axis.title.y = element_blank()),
                                            g11+ theme(legend.position="none")+ theme(axis.title.y = element_blank()),
                                            # new3+ theme(legend.position="none"),
                                            g12+ theme(legend.position="none"), 
                                            g13+ theme(legend.position="none")+ theme(axis.title.y = element_blank()),
                                           align = 'hv',
                                            ncol = 4)
  
  
ggdraw(newest.order.small.change.plot) +
  draw_plot(legend ,x =  0.55, y = 0, width = .35, height = .25)

ggsave( 
  filename = "FateTidy/output/TargetedMetabolites/MID.Lines.Small.and.Nuc.pdf",
  height = 4.2, width = 6.5)


plot_grid(new.order.small.change.plot, legend, ncol = 1, rel_heights =  c(7, 1))

# ggsave( 
#   filename = "FateTidy/output/TargetedMetabolites/MID.Lines.Smaller.Change.New.S.N.order.pdf",
#   height = 9, width = 6.5)

plot_grid(nucleic.plot, legend, ncol = 1, rel_heights =  c(7, 1))

# ggsave( 
#   filename = "FateTidy/output/TargetedMetabolites/MID.Lines.Smaller.Change.Nucleic.New.S.N.order.pdf",
#   height = 9, width = 4.5)

# plot_grid(newest.order.small.change.plot, legend, ncol = 1, rel_heights =  c(7, 1))



big.change.1 <- make.MID.line.plot(Big.change[1])
big.change.2 <- make.MID.line.plot(Big.change[2])
big.change.3 <- make.MID.line.plot(Big.change[3])
big.change.4 <- make.MID.line.plot(Big.change[4])
big.change.5 <- make.MID.line.plot(Big.change[5])
big.change.6 <- make.MID.line.plot(Big.change[6])
big.change.7 <- make.MID.line.plot(Big.change[7])

big.plot <- plot_grid(big.change.1+ theme(legend.position="none",
                                          axis.text.x = element_blank()),
                      big.change.2+ theme(legend.position="none",
                                          axis.title.y = element_blank(),
                                          axis.text.x = element_blank()),
                      big.change.3+ theme(legend.position="none",
                                          axis.text.x = element_blank()),
                      big.change.4+ theme(legend.position="none",
                                          axis.text.x = element_blank(),
                                          axis.title.y = element_blank()),
                      big.change.5+ theme(legend.position="none",
                                          axis.text.x = element_blank()), 
                      big.change.6+ theme(legend.position="none",
                                          axis.title.y = element_blank()),
                      big.change.7+ theme(legend.position="none"),  
                      ncol = 2)


# extract the legend from one of the plots and make it horizontal
legend <- get_legend(
  # create some space to the left of the legend
  big.change.1+
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)

plot_grid(big.plot, legend, ncol = 1, rel_heights =  c(7, 1))

ggsave( 
       filename = "FateTidy/output/TargetedMetabolites/MID.Lines.Big.Change.New.S.N.order.2.pdf",
       height = 9, width = 7)

## make individual line plots where N and S are stacked vertically
## make plot --------
make.MID.line.plot.vertical <- function(cmpd.name ="Choline",
                               scales.pick ="free"){
  cmpd.base.mass <- min(Dat.to.plot$iso.mz[Dat.to.plot$`Protein Name` == cmpd.name])
  nat.abu.g3.this <- g3.nat.abu %>%
    filter(`Protein Name`==cmpd.name,
           MID.mean != 0,
           !is.na(MID.mean)) %>%
    mutate(Experiment = ifelse(SampID == "KM1906S4C28D15","Fate1",
                               ifelse(SampID =="KM1906U12","Fate2","None")),
           Timepoint= -2) %>%
    dplyr::select(-Fraction1, -iso.mz) %>%
    left_join(., Dat.to.plot %>%
                filter(`Protein Name` == cmpd.name) %>%
                dplyr::select(`Protein Name`, `Precursor Ion Name`, iso.mz) %>%
                unique()) %>%
    filter(!(C.lab ==0 & N.lab == 0)) %>%
    mutate(Experiment = as.character(Experiment))
  
  if(cmpd.name == "Choline"){
    cmpd.MID.aves <- Dat.to.plot %>%
      filter(`Protein Name` == cmpd.name) %>%
      group_by(`Protein Name`, `Precursor Ion Name`, iso.mz,
               Timepoint, `Time of day`, `Incubation time (hr)`,
               SampID, C.lab, N.lab, Experiment) %>%
      summarise(MID.mean = mean(MID, na.rm = T),
                MID.sd = sd(MID, na.rm = T)) %>%
      arrange(desc(iso.mz)) %>%
      filter(!(C.lab ==0 & N.lab == 0)) %>%
      # full_join(nat.abu.g3.this) %>%
      mutate(`Incubation time (hr)` = ifelse(Timepoint == -1,
                                             -1, ifelse(Timepoint== -2, 
                                                        0,`Incubation time (hr)`))) %>%
      filter(MID.mean != 0,
             !is.na(MID.mean)) %>%
      arrange(iso.mz)
  } else {
    cmpd.MID.aves <- Dat.to.plot %>%
      filter(`Protein Name` == cmpd.name) %>%
      group_by(`Protein Name`, `Precursor Ion Name`, iso.mz,
               Timepoint, `Time of day`, `Incubation time (hr)`,
               SampID, C.lab, N.lab, Experiment) %>%
      summarise(MID.mean = mean(MID, na.rm = T),
                MID.sd = sd(MID, na.rm = T)) %>%
      arrange(desc(iso.mz)) %>%
      filter(!(C.lab ==0 & N.lab == 0)) %>%
      full_join(nat.abu.g3.this) %>%
      mutate(`Incubation time (hr)` = ifelse(Timepoint == -1,
                                             -1, ifelse(Timepoint== -2, 
                                                        0,`Incubation time (hr)`))) %>%
      filter(MID.mean != 0,
             !is.na(MID.mean)) %>%
      arrange(iso.mz)
  }
  
  num_toplot <- length(unique(cmpd.MID.aves$iso.mz))
  cmpd.MID.factors.colors <- cmpd.MID.aves  %>%
    ungroup() %>%
    dplyr::select(iso.mz, C.lab, N.lab) %>%
    unique() %>%
    mutate(iso.mz.fact = factor(iso.mz, levels = unique(cmpd.MID.aves$iso.mz))) %>%
    full_join(color.pal.test)
  cmpd.MID.aves <- full_join(cmpd.MID.aves, 
                             cmpd.MID.factors.colors) %>%
    mutate(MID.mean = MID.mean * 100,
           MID.sd = MID.sd * 100) %>%
    arrange(names.fact) %>%
    mutate(Experiment = factor(Experiment, levels = c("Fate2","Fate1")))
  # iso.mz.names <- cmpd.MID.aves %>%
  #   ungroup() %>%
  #   mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
  #   dplyr::select(names, iso.mz.fact) %>%
  #   unique()
  # names.vect <- iso.mz.names$iso.mz.fact
  # names(names.vect) <- iso.mz.names$names
  
  # myColors <- cmpd.MID.factors.colors$Colors
  # 
  # myColorScale <- scale_color_manual(values =  myColors,
  #                                    drop = F,
  #                                    labels = names(names.vect))
  Exp.labs <-  c("Exp 1", "Exp 2")
  names(Exp.labs) <- c("Fate1", "Fate2")
  
  g <- ggplot(cmpd.MID.aves %>%
                filter(`Protein Name` == cmpd.name,
                       Timepoint != -2) %>%
                mutate(Exp.to.facet.on = ifelse(Experiment=="Fate2","South","North"))) +
    geom_hline(data = cmpd.MID.aves %>%
                 filter(Timepoint == -2),
               aes(#x = Timepoint,
                 yintercept = MID.mean, group = names.fact, color = names.fact),
               size = 0.5, #shape = 17,
               linetype ="dotted"
    ) +
    geom_line(aes(x = `Incubation time (hr)`,   
                  y = MID.mean, group = names.fact, color = names.fact),
              size = 0.5) +
    geom_errorbar(aes(x = `Incubation time (hr)`,    ymin = MID.mean - MID.sd,
                      ymax = MID.mean + MID.sd, group = names.fact, color = names.fact),
                  alpha = 1) +
    facet_wrap(~Exp.to.facet.on, scales = scales.pick, nrow =2,
               # labeller = labeller(Experiment = Exp.labs)
               ) + 
    ggtitle(cmpd.name) +
    ylab("MID") +
    myColorScale +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          legend.position="right",
          legend.title = element_blank(),
          legend.text = element_text(size = 10, margin = margin(l = -5)),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  return(g)
}

## save out some vertically stacked plots ---
v.plot.glu <- make.MID.line.plot.vertical(cmpd.name = "L-Glutamic acid", scales.pick = "free_x")
v.plot.glu + theme(aspect.ratio = .5)

v.plot.gln <- make.MID.line.plot.vertical(cmpd.name = "L-Glutamine", scales.pick = "free_x")
v.plot.gln + theme(aspect.ratio = .5)

v.plot.dmsp <- make.MID.line.plot.vertical(cmpd.name = "DMSP", scales.pick = "free_x")
v.plot.dmsp + theme(aspect.ratio = 0.4)

# ggsave(v.plot.glu + theme(aspect.ratio = .4),
#        filename   = "FateTidy/output/TargetedMetabolites//individual line plots/Glu.pdf",
#       height = 5, width = 5)
# 
# ggsave(v.plot.gln + theme(aspect.ratio = .4),
#        filename   = "FateTidy/output/TargetedMetabolites/individual line plots/Gln.pdf",
#        height = 5, width = 5)
# 
# ggsave(v.plot.dmsp + theme(aspect.ratio = .4),
#        filename   = "FateTidy/output/TargetedMetabolites/individual line plots/DMSP.pdf",
#        height = 5, width = 5)
