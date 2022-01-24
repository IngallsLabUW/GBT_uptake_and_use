library(tidyverse)
library(ggplot2)
library(stats)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(broom)

theme_set(theme_cowplot())
## set mode: ---------
column.mode = "HILICPos"
# column.mode = "HILICNeg"
# column.mode = "CyanoAq"
# column.mode = "Lipid"

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
  filter(Fraction1 == column.mode)

if(column.mode == "Lipid"){
  iso.poss <- read_csv("FateTidy/All_Chl_C-N_Isotope_possibilities.csv") %>%
    mutate(Fraction1 = "Lipid",
           RT..min. = 22.7)
  iso.poss.use <- iso.poss %>%
    dplyr::select(Compound.and.Iso.name, iso.mz, RT..min., Column, Fraction1) %>%
    filter(Fraction1 == column.mode)
}

## read in the data ------
if(column.mode == "HILICPos"){
  
  qc.dat <- read_csv("FateTidy/output/TargetedMetabolites/HILICPos_QC_skyline_isotope_data_blankSubtract.csv")
  freq.dat <- read_csv("FateTidy/output/TargetedMetabolites/HILICPos_frequency_of_isotope_data.csv")
  freq.dat.3ormore <- freq.dat %>%
    filter(n >2)
}else if (column.mode == "HILICNeg"){
  qc.dat <- read_csv("FateTidy/output/TargetedMetabolites/HILICNeg_QC_skyline_isotope_data_blankSubtract.csv")
  freq.dat <- read_csv("FateTidy/output/TargetedMetabolites/HILICNeg_frequency_of_isotope_data.csv")
  freq.dat.3ormore <- freq.dat %>%
    filter(n >2)
} else if(column.mode == "CyanoAq"){
  qc.dat <- read_csv("FateTidy/output/TargetedMetabolites/CyanoAq_QC_skyline_isotope_data_blankSubtract.csv")
  freq.dat <- read_csv("FateTidy/output/TargetedMetabolites/CyanoAq_frequency_of_isotope_data.csv")
  freq.dat.3ormore <- freq.dat %>%
    filter(n >3)
} else if(column.mode == "Lipid"){
  qc.dat <- read_csv("FateTidy/output/TargetedMetabolites/Lipid_QC_skyline_isotope_data_blankSubtract.csv")
  freq.dat <- read_csv("FateTidy/output/TargetedMetabolites/Lipid_frequency_of_isotope_data.csv")
  freq.dat.3ormore <- freq.dat %>%
    filter(n > 3)
  
  ## use these for new RT Chl - on 9/16/2020
  # qc.dat <- read_csv("FateTidy/output/TargetedMetabolites/Lipid_newRT_Chl_QC_skyline_isotope_data_blankSubtract.csv")
  # freq.dat <- read_csv("FateTidy/output/TargetedMetabolites/Lipid_newRT_Chl_frequency_of_isotope_data.csv")
  # freq.dat.3ormore <- freq.dat %>%
  #   filter(n > 3)
  
}

joined.dat <- qc.dat %>%
  left_join(.,iso.poss.use, by = c("Precursor Ion Name" = "Compound.and.Iso.name", 
                                   "Fraction1", 
                                   "iso.mz",
                                   "RT..min.",
                                   "Column")) %>%
  unique() %>%
  left_join(freq.dat.3ormore) %>%
  filter(!is.na(n),
         n > 2) %>%
  dplyr::select(-n)


filtered.dat <- joined.dat %>%
  filter(`Precursor Ion Name` %in% unique(freq.dat$`Precursor Ion Name`),
         grepl("Smp",`Replicate Name`),
         `Protein Name` != "Internal Standards",
         `Protein Name` != "Internal Standards_neg",
         `Protein Name` != "Internal Standards_pos") %>%
  mutate(Area.filled.w.background = ifelse(!is.na(Area.less.blk),Area.less.blk,0))

## get rid of replicates where the base peak isn't found
no.base.peak <- filtered.dat %>%
  filter(N.lab==0,
         C.lab==0,
         is.na(Area)) %>%
  dplyr::select(`Replicate Name`, `Protein Name`) %>%
  mutate(nobasepeak = "yes") %>%
  unique()

filtered.dat <- full_join(filtered.dat,no.base.peak) %>%
  filter(is.na(nobasepeak) | nobasepeak != "yes") %>%
  dplyr::select(-nobasepeak)
all.cmpd.sums.by.group <- filtered.dat %>%
  group_by(`Protein Name`, `Replicate Name`) %>%
  summarise(all.area = sum(Area.filled.w.background))

MID <- filtered.dat %>%
  full_join(all.cmpd.sums.by.group) %>%
  mutate(MID = Area.filled.w.background/all.area) 

MID.averages <- MID %>%
  group_by(type, SampID, Experiment, Timepoint, ionization, runDate, N.lab, C.lab,
           `Precursor Ion Name`, `Protein Name`, Fraction1, iso.mz) %>%
  summarise(MID.mean = mean(MID),
            MID.sd = sd(MID),
            Area.filled.w.background.mean = mean(Area.filled.w.background),
            Area.filled.w.background.sd = sd(Area.filled.w.background))

## Plot MID as bars -------------------
pdf(file = paste0("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_blkSubtract_",column.mode,".pdf"),
    width = 11.5, height = 8)
for(i in 1:length(unique(MID.averages$`Protein Name`))){
  cmpd.name <- unique(MID.averages$`Protein Name`)[i]
  cmpd.MID.aves <- MID.averages %>%
    filter(`Protein Name` == cmpd.name) %>%
    filter(MID.mean>0.0075) %>%
    arrange(desc(iso.mz)) 
  cmpd.MID.aves<- cmpd.MID.aves %>%
    mutate(iso.mz.fact = factor(iso.mz, levels = unique(cmpd.MID.aves$iso.mz)))
  iso.mz.names <- cmpd.MID.aves %>%
    ungroup() %>%
    mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
    dplyr::select(names, iso.mz.fact) %>%
    unique()
  names.vect <- iso.mz.names$iso.mz.fact
  names(names.vect) <- iso.mz.names$names
  if(length(names.vect) < 9){
    myColors <- brewer.pal(n =length(names.vect), name = "Dark2" )
    names(myColors) <- names.vect
    myColorScale <- scale_fill_manual(values =  myColors,
                                      drop = F,
                                      name = "Isotopes",
                                      labels = names(names.vect))
  } else {
    myColors <- rainbow(n =length(names.vect))
    names(myColors) <- names.vect
    myColorScale <- scale_fill_manual(values =  myColors,
                                      drop = F,
                                      name = "Isotopes",
                                      labels = names(names.vect))
  }
f1.MLD.grid <- ggplot(cmpd.MID.aves %>%
                        filter(Experiment == "Fate1",
                               `Protein Name` == cmpd.name)) +
  geom_col(aes(x = Timepoint,    y = MID.mean, fill = iso.mz.fact))+
  theme_bw()+
  myColorScale +
  theme(axis.text.x = element_text(angle = 90),legend.position="bottom") +
  ggtitle(paste("MID of", cmpd.name, "\n","in", "Fate1" ))
f2.MLD.grid <- ggplot(cmpd.MID.aves%>%
                        filter(Experiment == "Fate2",
                               `Protein Name` == cmpd.name)) +
  geom_col(aes(x = Timepoint,    y = MID.mean, fill = iso.mz.fact))+
  myColorScale +
  # facet_wrap(~Timepoint, ncol =1) +  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),legend.position="bottom") +
  ggtitle(paste("MID of", cmpd.name, "\n","in", "Fate2"))
MLD.grid <- ggarrange(f1.MLD.grid, f2.MLD.grid, ncol = 2)
print(MLD.grid)
}
dev.off()

## write MID ave data -------------
write_csv(MID.averages,  paste0("FateTidy/output/TargetedMetabolites/MID_mean_QC_Targeted_blkSubtract",column.mode,".csv"))
write_csv(MID, paste0("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_blkSubtract", column.mode, ".csv"))

# ## line plots with natural abundance --------------
# natural.abundance.of.stds <- read_csv("FateTidy/RawData/Natural_abundance_of_Standards.csv")
# if(column.mode =="Lipid"){
#   natural.abundance.of.stds <- read_csv("FateTidy/RawData/Natural_abundance_of_Chl.csv")
#   HILIC.Pos.MID.dat.to.use <- MID.averages
# } 
# if(column.mode=="HILICPos"){
#   natural.abundance.from.g3 <- read_csv("FateTidy/RawData/Natural_abundance_HILICPos_from_G3.csv") %>%
#     dplyr::select(-Fraction1, -iso.mz) %>%
#     unique()
#   HILICPos.Good.but.no.change <- c("beta−Alanine",
#                                    "beta−Alaninebetaine",
#                                    "Betonicine",
#                                    "Homarine",
#                                    "L−Proline",
#                                    "Proline betaine",
#                                    "Taurine",
#                                    "L−Alanine",
#                                    "L-Methionine",
#                                    "DMSP",
#                                    "Gonyol",
#                                    "Glucosylglycerol"
#                                    
#   )
#   HILICPos.Good.with.change <- c("Glycine betaine",
#                                  "Sarcosine",
#                                  "(3−Carboxypropyl)trimethylammonium",
#                                  "Choline",
#                                  "Glycerophosphocholine",
#                                  "Trimethylamine N−oxide",
#                                  "Creatine",
#                                  "Carnitine",
#                                  "Adenine",
#                                  "Adenosine",
#                                  "Deoxyadenosine",
#                                  "Cytosine",
#                                  "Guanine",
#                                  "Guanosine",
#                                  "Trigonelline",
#                                  "L−Arginine",
#                                  "L−Glutamic acid",
#                                  "L−Glutamine",
#                                  "L−Leucine"
#   )
#   All.HILIC.pos <- c(HILICPos.Good.with.change, HILICPos.Good.but.no.change)
#   
#   HILIC.Pos.MID.dat <- read_csv("FateTidy/output/TargetedMetabolites/MID_mean_QC_Targeted_blkSubtractHILICPos.csv")
#   HILIC.Pos.MID.dat.to.use <- HILIC.Pos.MID.dat %>%
#     filter(`Protein Name` %in% All.HILIC.pos) %>%
#     mutate(`Protein Name` = factor(`Protein Name`, levels = All.HILIC.pos)) %>%
#     arrange(`Protein Name`)
#   
#   pdf(file = paste0("FateTidy/output/TargetedMetabolites/MID_QC_Targeted_linePlots_",column.mode,"_Good_Compounds_blkSubtract.pdf"),
#       width = 11.5, height = 8)
#   for(i in 1:length(unique(HILIC.Pos.MID.dat.to.use$`Protein Name`))){
#     
#     cmpd.name <- unique(HILIC.Pos.MID.dat.to.use$`Protein Name`)[i]
#     cmpd.base.mass <- min(HILIC.Pos.MID.dat.to.use$iso.mz[HILIC.Pos.MID.dat.to.use$`Protein Name` == cmpd.name])
#     
#     nat.abu.this <- natural.abundance.of.stds %>%
#       filter(Compound.Name == cmpd.name) %>%
#       gather(isotopes,probs,  -X1, -Emperical.Formula, -Compound.Name, -C, -N.numb) %>%
#       mutate(isotopes = ifelse(isotopes == "M.only.total","13C-0 15N-0",
#                                ifelse(isotopes == "prob.of15N1","13C-0 15N-1",
#                                       ifelse(isotopes == "prob.of.13C2", "13C-2 15N-0",
#                                              "13C-1 15N-0"))),
#              iso.mz = ifelse(isotopes == "13C-0 15N-0", cmpd.base.mass,
#                              ifelse(isotopes == "13C-1 15N-0", cmpd.base.mass + 1.003354,
#                                     ifelse(isotopes == "13C-2 15N-0", (cmpd.base.mass + (2*1.003354)),
#                                            (cmpd.base.mass + 0.997036)))),
#              C.lab = ifelse(isotopes == "13C-1 15N-0",1,
#                             ifelse(isotopes == "13C-2 15N-0",2,0)),
#              N.lab = ifelse(isotopes == "13C-0 15N-1",1,0),
#              MID.mean = probs/100,
#              Experiment = "Natural Abundance Estimate",
#              Timepoint = -1)
#     nat.abu.g3.this <- HILIC.Pos.MID.dat.to.use %>%
#       # filter(`Protein Name` == cmpd.name) %>%
#       dplyr::select(`Protein Name`, `Precursor Ion Name`, iso.mz) %>%
#       unique() %>%
#       full_join(.,  natural.abundance.from.g3  %>%
#                   mutate(Experiment = ifelse(SampID == "KM1906S4C28D15","Fate1",
#                                              ifelse(SampID =="KM1906U12","Fate2","")),
#                          Timepoint= -2))%>%
#       filter(`Protein Name`==cmpd.name)
#     
#     cmpd.MID.aves <- HILIC.Pos.MID.dat.to.use %>%
#       filter(`Protein Name` == cmpd.name) %>%
#       filter(MID.mean > 0.005,
#              !is.na(iso.mz)) %>%
#       arrange(desc(iso.mz)) %>%
#       full_join(nat.abu.this) %>% 
#       full_join(nat.abu.g3.this)
#     
#     cmpd.MID.aves<- cmpd.MID.aves  %>%
#       mutate(iso.mz.fact = factor(iso.mz, levels = unique(cmpd.MID.aves$iso.mz)))
#     iso.mz.names <- cmpd.MID.aves %>%
#       ungroup() %>%
#       mutate(names = paste0("13C-",C.lab, " 15N-",N.lab)) %>%
#       dplyr::select(names, iso.mz.fact) %>%
#       unique()
#     names.vect <- iso.mz.names$iso.mz.fact
#     names(names.vect) <- iso.mz.names$names
#     if(length(names.vect)<9){
#       myColors <- brewer.pal(n =length(names.vect), name = "Paired" )
#       names(myColors) <- names.vect
#       myColorScale <- scale_color_manual(values =  myColors,
#                                          drop = F,
#                                          name = "Isotopes",
#                                          labels = names(names.vect))
#     } else {
#       
#       myColors <- rainbow(n =length(names.vect))
#       names(myColors) <- names.vect
#       myColorScale <- scale_color_manual(values =  myColors,
#                                          drop = F,
#                                          name = "Isotopes",
#                                          labels = names(names.vect))
#     }
#     
#     
#     f1.MLD.grid <- ggplot(cmpd.MID.aves %>%
#                             filter(Experiment == "Fate1",
#                                    `Protein Name` == cmpd.name,
#                                    Timepoint != -2)) +
#       geom_hline(data = cmpd.MID.aves %>%
#                    filter(Timepoint == -2,
#                           Experiment == "Fate1"),
#                  aes(#x = Timepoint,
#                    yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
#                  size = 1, #shape = 17,
#                  linetype ="dotted"
#       ) +
#       geom_line(aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
#                 size = 1)+
#       geom_errorbar(aes(x = Timepoint,    ymin = MID.mean - MID.sd,
#                         ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
#                     alpha = 1)+
#       geom_point(data = cmpd.MID.aves %>%
#                    filter(Experiment == "Natural Abundance Estimate"),
#                  aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
#                  size = 2, alpha = 0.5)  +
#       theme_bw()+
#       ylim(0,1)+
#       myColorScale +
#       theme(axis.text.x = element_text(angle = 90),legend.position="bottom") +
#       ggtitle(paste("MID of", cmpd.name, "\n","in", "Fate1" ))
#     f2.MLD.grid <- ggplot(cmpd.MID.aves%>%
#                             filter(Experiment == "Fate2",
#                                    Timepoint != -2,
#                                    `Protein Name` == cmpd.name))  +
#       geom_hline(data = cmpd.MID.aves %>%
#                    filter(Timepoint == -2,
#                           Experiment == "Fate2"),
#                  aes(#x = Timepoint,
#                    yintercept = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
#                  size = 1, #shape = 17,
#                  linetype ="dotted"
#       ) +
#       geom_line(aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
#                 size = 1)+
#       geom_errorbar(aes(x = Timepoint,    ymin = MID.mean - MID.sd,
#                         ymax = MID.mean + MID.sd, group = iso.mz.fact, color = iso.mz.fact),
#                     alpha = 1)+ 
#       geom_point(data = cmpd.MID.aves %>%
#                    filter(Experiment == "Natural Abundance Estimate"),
#                  aes(x = Timepoint,    y = MID.mean, group = iso.mz.fact, color = iso.mz.fact),
#                  size = 2, alpha = 0.5) +
#       myColorScale +
#       # facet_wrap(~Timepoint, ncol =1) +  
#       theme_bw()+
#       ylim(0,1)+
#       theme(axis.text.x = element_text(angle = 90),legend.position="bottom") +
#       ggtitle(paste("MID of", cmpd.name, "\n","in", "Fate2"))
#     MLD.grid <- ggarrange(f1.MLD.grid, f2.MLD.grid, ncol = 1)
#     print(MLD.grid)
#   }
#   dev.off()
# }
# if(column.mode=="CyanoAq"){
#   natural.abundance.from.g3 <- read_csv("FateTidy/RawData/Natural_abundance_HILICPos_from_G3.csv") %>%
#     dplyr::select(-Fraction1, -iso.mz) %>%
#     unique()
#   HILIC.Pos.MID.dat.to.use <- MID.averages
# }


## determine if increases are significant --------------------
Good.but.no.change <- c("Betonicine",
                        "L−Alanine",
                        "beta−Alanine",
                        "L−Proline",
                        "Proline betaine",
                        "L-Methionine",
                        "(R)−2,3−Dihydroxypropane−1−sulfonate",
                        "Isethionic acid",
                        "Trehalose",
                        "L-Tryptophan",
                        "L-Phenylalanine")

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
                    "Homarine",
                    "DMSP",
                    "Gonyol",
                    "Glucosylglycerol",
                    "Methylthioadenosine",
                    "Butyrylcarnitine")
all.cmpds.to.test <- c("Glycine betaine",Good.but.no.change, Big.change, smaller.change)
MID.test <- MID %>%
  filter(`Protein Name` %in% all.cmpds.to.test) %>%
  ungroup() %>%
  group_by(`Precursor Ion Name`, Experiment) %>%
do(tidy(summary(lm(.$MID ~ .$Timepoint))$r.squared)) %>%
  dplyr::rename(Rsqr = x) %>%
  full_join(.,  MID %>%
              filter(`Protein Name` %in% all.cmpds.to.test) %>%
              ungroup() %>%
              group_by(`Precursor Ion Name`, Experiment) %>%
              do(tidy(anova(lm(.$MID ~ .$Timepoint))))%>%
              filter(term==".$Timepoint",
                     !is.na(term)))  %>%
  mutate(fdr.p.val = p.adjust(p.value, method = "fdr"),
         sig = fdr.p.val<0.05)  %>%
  filter(!is.na(Rsqr),
         !is.na(term))

write_csv(MID.test, paste0("FateTidy/output/TargetedMetabolites/", column.mode, "_goodCompounds_MID_lm_significance_test.csv"))
  
  
