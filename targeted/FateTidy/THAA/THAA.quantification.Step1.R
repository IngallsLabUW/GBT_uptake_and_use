## turn concentration in vial to concentration in seawater

library(tidyverse)
library(ggplot2)
library(cowplot)
library(anytime)
library(ggpmisc) 
library(rlist)
library(tidyr)
library(broom)
theme_set(theme_cowplot())

## load data----------------
THAA.standard.curve.data <- read_csv("FateTidy/RawData/StandardCurvePeakAreaData.csv")
sample.log <- read_csv("FateTidy/RawData/Sample.Log.csv")


THAA.standard.curve <- THAA.standard.curve.data %>%
  filter(Runtype  =="StandardCurve") %>%
  mutate(Spike = ifelse(grepl("0.5uM",Replicate.Name),0.5,
                        ifelse(grepl("_0uM",Replicate.Name),0,
                               ifelse(grepl("_1.0uM",Replicate.Name),1.0,
                                      ifelse(grepl("_2.5uM",Replicate.Name),2.5,NA)))))
## Get regression for each THAA----
fitted_models <- THAA.standard.curve %>%
  group_by(Precursor.Ion.Name) %>% 
  nest(data  = -Precursor.Ion.Name) %>%
  mutate(
    fit = map(data, ~lm(Area ~ Spike, data = .x)),
    tidied = map(fit,tidy)
  ) %>%
  unnest(tidied)
fitted_models.for.joining <- fitted_models %>%
  dplyr::select(Precursor.Ion.Name, term, estimate) %>%
  pivot_wider(names_from = "term", values_from = 'estimate',) %>%
  dplyr::rename(intercept.wo.blk = `(Intercept)`,
                slope.wo.blk = Spike)
fitted_models_blkSub <- THAA.standard.curve %>%
  group_by(Precursor.Ion.Name) %>% 
  nest(data  = -Precursor.Ion.Name) %>%
  mutate(
    fit = map(data, ~lm(Area_noBlank ~ Spike, data = .x)),
    tidied = map(fit,tidy)
  ) %>%
  unnest(tidied)

fitted_models.for.joining_blkSub <- fitted_models_blkSub %>%
  dplyr::select(Precursor.Ion.Name, term, estimate) %>%
  pivot_wider(names_from = "term", values_from = 'estimate',) %>%
  dplyr::rename(intercept.blk.sub = `(Intercept)`,
                slope.blk.sub = Spike)

## calculate concentration in samples ----
THAA.data = THAA.standard.curve.data %>% 
  # filter(Runtype=="Sample") %>%
  left_join( fitted_models.for.joining_blkSub) %>%
  left_join(fitted_models.for.joining) %>%
  mutate(calc.uM.in.vial.slope.only = Area_noBlank/slope.blk.sub,
         Concenctration.uM.in.vial.with.Dilution.Correction = calc.uM.in.vial.slope.only * 10)%>%
  full_join(sample.log, by = c("Replicate.Name" = "Replicate Name")) %>%
  mutate(Environmental.Concentration.uM = Concenctration.uM.in.vial.with.Dilution.Correction *
           `reconstitution volume (L)`/`filtered volume (L)`)


quantified.THAA.final = THAA.data %>%
  dplyr::rename(Slope = slope.blk.sub,
                Intercept = intercept.blk.sub ) %>%
  dplyr::select(Replicate.Name,	Precursor.Ion.Name,	Area,	Area_noBlank,	Slope,	Intercept, Runtype,
                calc.uM.in.vial.slope.only,	Concenctration.uM.in.vial.with.Dilution.Correction,
                Environmental.Concentration.uM)

## write data:
write_csv(quantified.THAA.final, "FateTidy/RawData/Quantified_THAA_inSeaWater.csv")
