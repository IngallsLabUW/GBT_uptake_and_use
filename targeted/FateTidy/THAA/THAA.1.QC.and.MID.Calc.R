# MID calculation for THAA
## set up------------
library(tidyverse)
library(ggplot2)
library(stats)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(broom)
library(RCurl)

theme_set(theme_cowplot())

## QC params -------
RT.flex = 0.5
RT.flex.iso = 0.2
min.value = 1e3
This.much.bigger.than.max.Blank = 2
ppm.thresh = 15

## load data ------------
iso.poss <- read_csv("FateTidy/RawData/THAA_Isotope_possibilities.csv")
dat.raw = read_csv("FateTidy/RawData/GBT_Fate_THAA_Isotopes_skyline_report.csv",
                   na = "#N/A") %>%
  filter(`Protein Name` != "Amino Acids")

## do qc-------------------
joined.dat <- dat.raw %>%
  left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name")) %>%
  unique() %>%
  mutate(Area = as.numeric(Area),
         name.to.break = `Replicate Name`,
         labels.to.break = `Precursor Ion Name`,
         N.lab = str_sub(`Precursor Ion Name`, start = -1, -1),
         C.lab = str_sub(`Precursor Ion Name`, start = -7, -7)) %>%
  # dplyr::select(-Other) %>%
  separate(name.to.break, into = c("runDate", "type", "SampID", "ionization","rep")) %>%
  mutate(Timepoint = str_sub(SampID, start = 10, -1),
         Timepoint = ifelse(Timepoint == "long", 96,Timepoint),
         Timepoint = as.numeric(Timepoint),
         Experiment = str_sub(SampID, start = 4,8))

## Step 2: Quality control by using the blanks--------
## get blank data for base peak
blank.dat <- joined.dat %>%
  filter(grepl("Blk", `Replicate Name`)) 
## get max value
max.blank <- blank.dat %>%
  filter(abs(`Retention Time` - RT) < RT.flex) %>%
  group_by(`Protein Name`, `Precursor Ion Name`, 
           iso.mz, RT, N.lab, C.lab) %>%
  summarise(Blk.Max.Area = max(Area, na.rm = T)) %>%
  mutate(Blk.Max.Area = ifelse(is.infinite(Blk.Max.Area),0,Blk.Max.Area)) 
## compare base peak to blank
blank.check <- full_join(joined.dat, max.blank) %>%
  mutate(Blk.ratio = ifelse(Area==0,0,
                            ifelse(!is.na(Blk.Max.Area),Area/Blk.Max.Area,Inf)),
         Blk.flag = ifelse(Blk.ratio < This.much.bigger.than.max.Blank, "comparable to blank","ok"))

## compare ppm to our defined threshold
ppm.check <- blank.check %>%
  mutate(ppm.flag = ifelse(abs(`Mass Error PPM`) > ppm.thresh, "bad ppm", "ok"))

## size.check
size.check <- ppm.check %>%
  mutate(size.flag = ifelse(Area > min.value, "ok", "too small"))


## compare the RT of the base peak to the expected (from standards)
std.rt <- joined.dat %>%
  filter(grepl("Std",`Replicate Name`),
         SampID != "H2OinMatrix",
         N.lab == 0,
         C.lab == 0) %>%
  ungroup() %>%
  group_by(`Protein Name`, `Precursor Ion Name`, iso.mz, RT, N.lab, C.lab) %>%
  summarise(minRT = min(`Retention Time`),
            maxRT = max(`Retention Time`)) %>%
  gather(RT.option, Value, -`Protein Name`, -`Precursor Ion Name`, -iso.mz, -RT, -N.lab, -C.lab) %>%
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
                `N.lab`,`C.lab`,`runDate`,`type`,
                `SampID`, `ionization`,`rep`, `Timepoint`, `Experiment`,`Blk.flag`,
                `ppm.flag`,size.flag, RT.check, all.flags) %>%
  mutate(rawArea = Area,
         Area = ifelse(all.flags == "all ok",Area,NA))

## ditch compounds I know are bad
# (13C3-15N1 is not Alanine, wrong RT)
qc.dat.full.trim <- qc.dat %>%
  filter(`Protein Name` !="Methionine",
         `Protein Name` !="Histidine")

cmpd.freq.tab <- qc.dat.full.trim %>%
  filter(type == "Smp",
         !is.na(Area)) %>%
  ungroup() %>%
  group_by(`Protein Name`, C.lab, N.lab, `Precursor Ion Name`, Experiment) %>%
  summarise(n = n())
blk.sub.dat <- qc.dat.full.trim %>%
  full_join(max.blank)  %>%
  mutate(Blk.Max.Area = ifelse(is.na(Blk.Max.Area),0,Blk.Max.Area),
         Area.less.blk = Area - Blk.Max.Area)

## write out QC data --------
write_csv(qc.dat.full.trim, "FateTidy/output/THAA/QC.isotopes.THAA.csv")
write_csv(blk.sub.dat, "FateTidy/output/THAA/QC.isotopes.THAA.blk.sub.csv")
write_csv(cmpd.freq.tab, "FateTidy/output/THAA/THAA.isotopologue.frequency.csv")

## calculate MID -----------------
joined.dat <- blk.sub.dat %>%
  left_join(.,iso.poss, by = c("Precursor Ion Name" = "Compound.and.Iso.name",
                               "Protein Name" = "Group",
                               "RT", "iso.mz")) %>%
  unique()%>%
  left_join(cmpd.freq.tab) %>%
  filter(!is.na(n),
         n > 1) %>%
  dplyr::select(-n)


filtered.dat <- joined.dat %>%
  filter(`Precursor Ion Name` %in% unique(cmpd.freq.tab$`Precursor Ion Name`),
         grepl("Smp",`Replicate Name`),
         `Protein Name` != "Internal Standards") %>%
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
  group_by(type, SampID, Experiment, Timepoint, runDate, N.lab, C.lab,
           `Precursor Ion Name`, `Protein Name`,  iso.mz) %>%
  summarise(MID.mean = mean(MID),
            MID.sd = sd(MID),
            Area.filled.w.background.mean = mean(Area.filled.w.background),
            Area.filled.w.background.sd = sd(Area.filled.w.background))

## Plot MID as bars -------------------
pdf(file = paste0("FateTidy/output/THAA/MID_QC_blkSubtract.pdf"),
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
write_csv(MID.averages,  paste0("FateTidy/output/THAA/MID_mean_QC_blkSubtract",".csv"))
write_csv(MID, paste0("FateTidy/output/THAA/MID_QC_blkSubtract", ".csv"))


## determine if increases are significant --------------------
MIDs.to.test <- MID %>%
  group_by(Experiment, `Precursor Ion Name`) %>%
  summarise(n = n()) %>%
  filter(n>4) 

MID.test <- MID %>%
  right_join(., MIDs.to.test) %>%
  ungroup() %>%
  group_by(`Precursor Ion Name`, Experiment) %>%
  do(tidy(t.test(.$MID ~ .$Timepoint))) %>%
  ungroup() %>%
  mutate(fdr.p.val = p.adjust(p.value, method = "fdr"),
         sig = fdr.p.val<0.1) 

write_csv(MID.test, paste0("FateTidy/output/THAA/THAA_MID_t.test_significance_test.csv"))

## plot significant increases ------------------------
MID.test.sig = MID.test %>%
  filter(p.value < 0.05) %>%
  dplyr::select(`Precursor Ion Name`, Experiment, estimate1, estimate2, p.value) %>%
  left_join(.,MID) %>%
  mutate(Ave.MID.Val = ifelse(Timepoint == 0, estimate1, estimate2)) %>%
  filter(!grepl("Proline",`Precursor Ion Name`))

ggplot(MID.test.sig) +
  geom_point(aes(x=Timepoint, y=Ave.MID.Val*100),size=6) +
  geom_line(aes(x=Timepoint, y=Ave.MID.Val*100)) +
  geom_point(aes(x=Timepoint, y=MID*100), size = 3, alpha = 0.5)+
  facet_wrap(~Experiment + `Precursor Ion Name`, scales = "free_y",
             labeller = labeller(Experiment = c("Fate2" = "South")))+
  expand_limits(y = 0) +
  # scale_y_continuous(expand = c(0, 0))+
  
  ylab("MID")+
  theme(text = element_text(size = 9))

# ggsave(paste0("Fate/output/Figures/THAA_MID_increasing_Isotopologues.new.scale.pdf"),
          # width = 6, height=6)


ggplot(MID.test.sig %>%
         mutate(Timepoint = as.character(Timepoint))) +
  geom_point(aes(x=Timepoint, y=Ave.MID.Val*100),size=3, alpha = 0.75) +
  # geom_line(aes(x=Timepoint, y=Ave.MID.Val*100)) +
  geom_point(aes(x=Timepoint, y=MID*100,
                 shape = ionization), size = 1, alpha = 0.75)+
  facet_wrap(~ `Precursor Ion Name`, scales = "free_y",
             labeller = labeller(Experiment = c("Fate2" = "South")))+
  ylab("MID (%)")+
  expand_limits(y = 0) +
  # scale_y_continuous(expand = c(0, 0))+
  theme(text = element_text(size = 10),
        legend.position = "none")

ggsave(paste0("Fate/output/Figures/THAA_MID_increasing_Isotopologues_new.scale.pdf"),
       width = 5, height=4)

MID.test.sig.sub = MID.test %>%
  # filter(p.value < 0.05) %>%
  dplyr::select(`Precursor Ion Name`, Experiment, estimate1, estimate2, p.value) %>%
  left_join(.,MID) %>%
  mutate(Ave.MID.Val = ifelse(Timepoint == 0, estimate1, estimate2),
         sig.y.n = ifelse(p.value < 0.05, T,F)) %>%
  filter(!grepl("Proline",`Precursor Ion Name`)) %>%
  filter(Experiment == "Fate2",
         C.lab==2,
         `Precursor Ion Name` == "Alanine_13C-2 15N-0"|
           `Precursor Ion Name` == "Alanine_13C-2 15N-1"|
           `Precursor Ion Name` == "Isoeucine_13C-2 15N-0"|
           `Precursor Ion Name` == "Isoleucine_13C-2 15N-1"|
           `Precursor Ion Name` == "Leucine_13C-2 15N-0"|
           `Precursor Ion Name` == "Leucine_13C-2 15N-1"|
           `Precursor Ion Name` == "Valine_13C-2 15N-0"|
           `Precursor Ion Name` == "Valine_13C-2 15N-1")
ggplot(MID.test.sig.sub) +
  # geom_point(aes(x=Timepoint, y=Ave.MID.Val*100),size=6) +
  # geom_line(aes(x=Timepoint, y=Ave.MID.Val*100)) +
  geom_boxplot(aes(x=Timepoint,group=Timepoint, y=MID*100,
                  ),  alpha = 0.5)+
  facet_wrap(~`Precursor Ion Name`, scales = "free_y", nrow = 1)+
  ylab("MID")+
  scale_y_continuous(expand  = c(0,0),limits = c(0,NA))+
  theme(text = element_text(size = 9)) 
# ggsave(paste0("Fate/output/Figures/THAA_MID_increasing_Isotopologues.only.some.pdf"),
       # width = 6, height=2)

# ## now do MID calculations based on the quantified data where we subtracted free AAs -------
# All.MID.conc.uMC <-read_csv("Fate/output/THAA/THAA.FreeAA.Concentrations.csv")
# neg.compounds <- All.MID.conc.uMC %>%
#   filter(THAA.less.FAA.uM.C<0) %>%
#   dplyr::select(`Protein Name`) %>%
#   unique()
# 
# all.cmpd.conc.sums.by.group <- All.MID.conc.uMC %>%
#   filter(!(`Protein Name`%in% neg.compounds$`Protein Name`),
#          THAA.less.FAA.uM.C>0) %>%
#   ungroup() %>%
#   group_by(`Protein Name`, `Experiment`, SampID, type, rep, Timepoint) %>%
#   summarise(all.uM.C = sum(THAA.less.FAA.uM.C)) %>%
#   unique()
# 
# MID.sub <- All.MID.conc.uMC %>%  
#   filter(!(`Protein Name`%in% neg.compounds$`Protein Name`),
#          THAA.less.FAA.uM.C>0) %>%
#   full_join(all.cmpd.conc.sums.by.group) %>%
#   mutate(new.MID = THAA.less.FAA.uM.C/all.uM.C)  %>%
#   filter(!is.infinite(new.MID))
# 
# MID.sub.averages <- MID.sub %>%
#   group_by(type, SampID, Experiment, Timepoint,  N.lab, C.lab,
#            `Precursor Ion Name`, `Protein Name`) %>%
#   summarise(MID.mean = mean(new.MID),
#             MID.sd = sd(new.MID),
#             THAA.less.FAA.uM.C.mean = mean(THAA.less.FAA.uM.C),
#             THAA.less.FAA.uM.C.sd = sd(THAA.less.FAA.uM.C))
# MIDs.to.test.sub <- MID.sub %>%
#   group_by(Experiment, `Precursor Ion Name`) %>%
#   summarise(n = n()) %>%
#   filter(n>4) 
# 
# MID.test.sub <- MID.sub %>%
#   right_join(., MIDs.to.test.sub) %>%
#   unique() %>%
#   ungroup() %>%
#   group_by(`Precursor Ion Name`, Experiment) %>%
#   do(tidy(t.test(.$new.MID ~ .$Timepoint))) %>%
#   ungroup() %>%
#   mutate(fdr.p.val = p.adjust(p.value, method = "fdr"),
#          sig = fdr.p.val<0.1) 
# 
# write_csv(MID.test.sub, "Fate/output/THAA/THAA_less_freeAA_MID_t.test_significance_test.csv")
