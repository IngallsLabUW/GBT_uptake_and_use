## quantifying the GBT and Betaine in the GBT Kinetics Experiments

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(stats)

## sample list --
sample.list <- read_csv("../data/GBT Kintetics Gradients 3 Sample List TQS for R code use.csv")

## data--
data.raw <- read_csv("../data/QC_output_All_TQS_GBT_IS-GBT.csv")

## Std curve----
key <- c("GBT-K2-M2000nM_C", "GBT-K1-M2nM_B")

std.curve.data <- full_join(data.raw, sample.list) %>%
  filter(grepl(key[2],Replicate.Name)|grepl(key[1],Replicate.Name)) %>%
  mutate(earlyOrlate = ifelse(`Injection Order`<60,"early","late"))

low.std.curve.data <- full_join(data.raw, sample.list) %>%
  filter(grepl(key[2],Replicate.Name)) %>%
  mutate(earlyOrlate = ifelse(`Injection Order`<60,"early","late"))
high.std.curve.data <- full_join(data.raw, sample.list) %>%
  filter(grepl(key[1],Replicate.Name)) %>%
  mutate(earlyOrlate = ifelse(`Injection Order`<60,"early","late"))

## do each curve separately-----
## low and early
GBT.IS.L.E <- low.std.curve.data %>% filter(earlyOrlate == "early", Compound.Name == "13C-15N Glycine Betaine")
GBT.IS.L.E <- lm(formula = rawArea ~ Added_GBT_IS, data = GBT.IS.L.E )
summary(GBT.IS.L.E)
GBT.L.E <- low.std.curve.data %>% filter(earlyOrlate == "early", Compound.Name == "Glycine Betaine")
GBT.L.E <- lm(formula = rawArea ~ Added_GBT, data = GBT.L.E )
summary(GBT.L.E)

## low and late
GBT.L.L <- low.std.curve.data %>% filter(earlyOrlate != "early", Compound.Name == "Glycine Betaine")
GBT.L.L <- lm(formula = rawArea ~ Added_GBT, data = GBT.L.L )
summary(GBT.L.L)
GBT.IS.L.L <- low.std.curve.data %>% filter(earlyOrlate != "early", Compound.Name == "13C-15N Glycine Betaine")
GBT.IS.L.L <- lm(formula = rawArea ~ Added_GBT_IS, data = GBT.IS.L.L )
summary(GBT.IS.L.L)

## low and all
GBT.L.all <- low.std.curve.data %>% filter(Compound.Name == "Glycine Betaine")
GBT.L.all <- lm(formula = rawArea ~ Added_GBT, data = GBT.L.all )
summary(GBT.L.all)
GBT.IS.L.all <- low.std.curve.data %>% filter(Compound.Name == "13C-15N Glycine Betaine")
GBT.IS.L.all <- lm(formula = rawArea ~ Added_GBT_IS, data = GBT.IS.L.all )
summary(GBT.IS.L.all)

## high and early
GBT.H.E <- high.std.curve.data %>% filter(earlyOrlate == "early", Compound.Name == "Glycine Betaine")
GBT.H.E <- lm(formula = rawArea ~ Added_GBT, data = GBT.H.E )
summary(GBT.H.E)
GBT.IS.H.E <- high.std.curve.data %>% filter(earlyOrlate == "early", Compound.Name == "13C-15N Glycine Betaine")
GBT.IS.H.E <- lm(formula = rawArea ~ Added_GBT_IS, data = GBT.IS.H.E )
summary(GBT.IS.H.E)

## high and late
GBT.H.L <- high.std.curve.data %>% filter(earlyOrlate != "early", Compound.Name == "Glycine Betaine")
GBT.H.L <- lm(formula = rawArea ~ Added_GBT, data = GBT.H.L )
summary(GBT.H.L)
GBT.IS.H.L <- high.std.curve.data %>% filter(earlyOrlate != "early", Compound.Name == "13C-15N Glycine Betaine")
GBT.IS.H.L <- lm(formula = rawArea ~ Added_GBT_IS, data = GBT.IS.H.L )
summary(GBT.IS.H.L)

## high and all
GBT.H.all <- high.std.curve.data %>% filter(Compound.Name == "Glycine Betaine")
GBT.H.all <- lm(formula = rawArea ~ Added_GBT, data = GBT.H.all )
summary(GBT.H.all)
GBT.IS.H.all <- high.std.curve.data %>% filter(Compound.Name == "13C-15N Glycine Betaine")
GBT.IS.H.all <- lm(formula = rawArea ~ Added_GBT_IS, data = GBT.IS.H.all )
summary(GBT.IS.H.all)

## write model params and errors ---------
GBT.Models <- data.frame(slope = c(L.E = GBT.L.E$coefficients[2],
                                   L.L = GBT.L.L$coefficients[2],
                                   L.ALL = GBT.L.all$coefficients[2],
                                   H.E = GBT.H.E$coefficients[2],
                                   H.L = GBT.H.L$coefficients[2],
                                   H.ALL = GBT.H.all$coefficients[2]),
                         slope.error = c(L.E = coef(summary(GBT.L.E))[2, "Std. Error"],
                                         L.L = coef(summary(GBT.L.L))[2, "Std. Error"],
                                         L.ALL = coef(summary(GBT.L.all))[2, "Std. Error"],
                                         H.E = coef(summary(GBT.H.E))[2, "Std. Error"],
                                         H.L = coef(summary(GBT.H.L))[2, "Std. Error"],
                                         H.ALL = coef(summary(GBT.H.all))[2, "Std. Error"]),
                         intercept = c(L.E = GBT.L.E$coefficients[1],
                                       L.L = GBT.L.L$coefficients[1],
                                       L.ALL = GBT.L.all$coefficients[1],
                                       H.E = GBT.H.E$coefficients[1],
                                       H.L = GBT.H.L$coefficients[1],
                                       H.ALL = GBT.H.all$coefficients[1]),
                         intercept.error = c(L.E = coef(summary(GBT.L.E))[1, "Std. Error"],
                                             L.L = coef(summary(GBT.L.L))[1, "Std. Error"],
                                             L.ALL = coef(summary(GBT.L.all))[1, "Std. Error"],
                                             H.E = coef(summary(GBT.H.E))[1, "Std. Error"],
                                             H.L = coef(summary(GBT.H.L))[1, "Std. Error"],
                                             H.ALL = coef(summary(GBT.H.all))[1, "Std. Error"]))%>%
  mutate(Names = row.names(data.frame(slope = c(L.E = GBT.L.E$coefficients[2],
                                                L.L = GBT.L.L$coefficients[2],
                                                L.ALL = GBT.L.all$coefficients[2],
                                                H.E = GBT.H.E$coefficients[2],
                                                H.L = GBT.H.L$coefficients[2],
                                                H.ALL = GBT.L.all$coefficients[2]))),
         percent.intercept.error = intercept.error/intercept,
         percent.slope.error = slope.error/slope)

GBT.IS.Models <- data.frame(slope = c(L.E = GBT.IS.L.E$coefficients[2],
                                   L.L = GBT.IS.L.L$coefficients[2],
                                   L.ALL = GBT.IS.L.all$coefficients[2],
                                   H.E = GBT.IS.H.E$coefficients[2],
                                   H.L = GBT.IS.H.L$coefficients[2],
                                   H.ALL = GBT.IS.L.all$coefficients[2]),
                         slope.error = c(L.E = coef(summary(GBT.IS.L.E))[2, "Std. Error"],
                                         L.L = coef(summary(GBT.IS.L.L))[2, "Std. Error"],
                                         L.ALL = coef(summary(GBT.IS.L.all))[2, "Std. Error"],
                                         H.E = coef(summary(GBT.IS.H.E))[2, "Std. Error"],
                                         H.L = coef(summary(GBT.IS.H.L))[2, "Std. Error"],
                                         H.ALL = coef(summary(GBT.IS.H.all))[2, "Std. Error"]),
                         intercept = c(L.E = GBT.IS.L.E$coefficients[1],
                                       L.L = GBT.IS.L.L$coefficients[1],
                                       L.ALL = GBT.IS.L.all$coefficients[1],
                                       H.E = GBT.IS.H.E$coefficients[1],
                                       H.L = GBT.IS.H.L$coefficients[1],
                                       H.ALL = GBT.IS.H.all$coefficients[1]),
                         intercept.error = c(L.E = coef(summary(GBT.IS.L.E))[1, "Std. Error"],
                                             L.L = coef(summary(GBT.IS.L.L))[1, "Std. Error"],
                                             L.ALL = coef(summary(GBT.IS.L.all))[1, "Std. Error"],
                                             H.E = coef(summary(GBT.IS.H.E))[1, "Std. Error"],
                                             H.L = coef(summary(GBT.IS.H.L))[1, "Std. Error"],
                                             H.ALL = coef(summary(GBT.IS.H.all))[1, "Std. Error"])) %>%
  mutate(Names = row.names(data.frame(slope = c(L.E = GBT.IS.L.E$coefficients[2],
                                                L.L = GBT.IS.L.L$coefficients[2],
                                                L.ALL = GBT.IS.L.all$coefficients[2],
                                                H.E = GBT.IS.H.E$coefficients[2],
                                                H.L = GBT.IS.H.L$coefficients[2],
                                                H.ALL = GBT.IS.L.all$coefficients[2]))),
         percent.intercept.error = intercept.error/intercept,
         percent.slope.error = slope.error/slope)

All.Models  = full_join(GBT.Models, GBT.IS.Models)
write.csv(All.Models, paste0("../output/All.Standard.Addition.Models.csv"))

## GBT plot --------
ggplot(data = std.curve.data %>% filter(Compound.Name=="Glycine Betaine"),
       aes(x = Added_GBT, y= rawArea, group = interaction(earlyOrlate, `Sample ID`))) +
  geom_point(aes(x = Added_GBT, y= rawArea, group = interaction(earlyOrlate, `Sample ID`), color = `Sample ID`)) +
  geom_smooth(method = "lm") +
  ggtitle("GBT") + theme_minimal_grid()
# ggsave(paste0("../output/GBT Standard Addition Curves.pdf"))

## 13C, 15N-GBT plot --------
ggplot(data = std.curve.data %>% filter(Compound.Name=="13C-15N Glycine Betaine"),
       aes(x = Added_GBT_IS, y= rawArea, group = interaction(earlyOrlate, `Sample ID`))) +
  geom_point(aes(x = Added_GBT_IS, y= rawArea, group = interaction(earlyOrlate, `Sample ID`), color = `Sample ID`)) +
  geom_smooth(method = "lm") +
  ggtitle("13C-15N Glycine Betaine") + theme_minimal_grid()
# ggsave(paste0("../output/13C-15N Glycine Betaine Standard Addition Curves.pdf"))


## calculate GBT for curve samples ------
## value in vial: GBT ----
## -1*Intercept/Slope = Concentration in vial
## error of concentration in vial = [Concentration] *  sqrt(([Error of Slope]/[Slope])^2+([Error of Intercept]/[Intercept])^2)
Concentration_in_K1_2nM_B.early <- GBT.L.E$coefficients[1]/GBT.L.E$coefficients[2]
Concentration_in_K1_2nM_B.late <- GBT.L.L$coefficients[1]/GBT.L.L$coefficients[2]
SE.Concentration_in_K1_2nM_B.early <- Concentration_in_K1_2nM_B.early* 
  sqrt((coef(summary(GBT.L.E))[1, "Std. Error"]/GBT.L.E$coefficients[1])^2+
         (coef(summary(GBT.L.E))[2, "Std. Error"]/GBT.L.E$coefficients[2])^2)
SE.Concentration_in_K1_2nM_B.late <- Concentration_in_K1_2nM_B.late* 
  sqrt((coef(summary(GBT.L.L))[1, "Std. Error"]/GBT.L.L$coefficients[1])^2+
         (coef(summary(GBT.L.L))[2, "Std. Error"]/GBT.L.L$coefficients[2])^2)


Concentration_in_K1_2nM_B.all <- GBT.L.all$coefficients[1]/GBT.L.all$coefficients[2]
SE.Concentration_in_K1_2nM_B.all <- Concentration_in_K1_2nM_B.all* 
  sqrt((coef(summary(GBT.L.all))[1, "Std. Error"]/GBT.L.all$coefficients[1])^2+
         (coef(summary(GBT.L.all))[2, "Std. Error"]/GBT.L.all$coefficients[2])^2)

test.GBT <- full_join(data.raw, sample.list) %>%
  mutate(earlyOrlate = ifelse(`Injection Order`<60,"early","late"))%>%
  filter(Compound.Name == "Glycine Betaine",
         !grepl("added",`File Text`)) %>%
  mutate(All.Models.High.Prediction.fmoles = rawArea/GBT.H.all$coefficients[2]*400,
         SE.High.Prediction =  All.Models.High.Prediction.fmoles* 
           sqrt((coef(summary(GBT.H.all))[1, "Std. Error"]/GBT.H.all$coefficients[1])^2+
                  (coef(summary(GBT.H.all))[2, "Std. Error"]/GBT.H.all$coefficients[2])^2),
         All.Models.low.Prediction.fmoles = rawArea/GBT.L.all$coefficients[2]*400,
         SE.low.Prediction =  All.Models.low.Prediction.fmoles* 
           sqrt((coef(summary(GBT.L.all))[1, "Std. Error"]/GBT.L.all$coefficients[1])^2+
                  (coef(summary(GBT.L.all))[2, "Std. Error"]/GBT.L.all$coefficients[2])^2)) %>%
  select(Replicate.Name, Compound.Name, All.Models.High.Prediction.fmoles, 
         SE.High.Prediction, All.Models.low.Prediction.fmoles, SE.low.Prediction,
         Experiment, Treatment, Added_GBT_IS, Added_GBT, replicate, `Sample ID`)
bio.info <- read_csv("../data/BiologicNormalizationInformation.csv")
test.GBT.2 <- full_join(bio.info, test.GBT) %>%
  mutate(H.Prediction.nmoles.per.L = ifelse(is.na(`Vol Filtered (L)`),
                                            All.Models.High.Prediction.fmoles/2*1e-6*10,
                                            All.Models.High.Prediction.fmoles/`Vol Filtered (L)`*1e-6*dilution.factor),
         L.Prediction.nmoles.per.L = ifelse(is.na(`Vol Filtered (L)`),
                                            All.Models.low.Prediction.fmoles/2*1e-6*10,
                                            All.Models.low.Prediction.fmoles/`Vol Filtered (L)`*1e-6*dilution.factor)
         )

k1.mean = mean(test.GBT.2$H.Prediction.nmoles.per.L[test.GBT.2$Experiment=="K1"])
k1.sd = sd(test.GBT.2$H.Prediction.nmoles.per.L[test.GBT.2$Experiment=="K1"])
k2.mean = mean(test.GBT.2$H.Prediction.nmoles.per.L[test.GBT.2$Experiment=="K2"])
k2.sd = sd(test.GBT.2$H.Prediction.nmoles.per.L[test.GBT.2$Experiment=="K2"])

g <- ggplot(data = test.GBT.2) +
  geom_boxplot(aes(x = Experiment, y = L.Prediction.nmoles.per.L))+
  ylim(0,1.2)+ ggtitle("Particulate Concentration - GBT") +
  theme_bw(15)
p <- ggplot(data = test.GBT.2) +
  geom_boxplot(aes(x = Experiment, y = H.Prediction.nmoles.per.L)) +
  ylim(0,1.2)+
  ggtitle("Particulate Concentration - GBT")+
  theme_bw(15)
library(ggpubr)
g2<- ggarrange(g, p, ncol = 1)
g2


test.GBT.2.aves = test.GBT.2 %>%
  ungroup() %>%
  group_by(Treatment, Experiment, Compound.Name) %>%
  summarise(ave.L.nM = mean(L.Prediction.nmoles.per.L),
            sd.L.nM = sd(L.Prediction.nmoles.per.L),
            ave.H.nM = mean(H.Prediction.nmoles.per.L),
            sd.H.nM = sd(H.Prediction.nmoles.per.L))


test.GBT.3.aves = test.GBT.2 %>%
  ungroup() %>%
  group_by( Experiment, Compound.Name) %>%
  summarise(ave.L.nM = mean(L.Prediction.nmoles.per.L),
            sd.L.nM = sd(L.Prediction.nmoles.per.L),
            ave.H.nM = mean(H.Prediction.nmoles.per.L),
            sd.H.nM = sd(H.Prediction.nmoles.per.L))
## value in vial: IS_GBT ----
## -1*Intercept/Slope = Concentration in vial
## error of concentration in vial = [Concentration] *  sqrt(([Error of Slope]/[Slope])^2+([Error of Intercept]/[Intercept])^2)
IS_Concentration_in_K1_2nM_B.early <- GBT.IS.L.E$coefficients[1]/GBT.IS.L.E$coefficients[2]
IS_Concentration_in_K1_2nM_B.late <- GBT.IS.L.L$coefficients[1]/GBT.IS.L.L$coefficients[2]
SE.IS_Concentration_in_K1_2nM_B.early <- IS_Concentration_in_K1_2nM_B.early* 
  sqrt((coef(summary(GBT.IS.L.E))[1, "Std. Error"]/GBT.IS.L.E$coefficients[1])^2+
         (coef(summary(GBT.IS.L.E))[2, "Std. Error"]/GBT.IS.L.E$coefficients[2])^2)
SE.IS_Concentration_in_K1_2nM_B.late <- IS_Concentration_in_K1_2nM_B.late* 
  sqrt((coef(summary(GBT.IS.L.L))[1, "Std. Error"]/GBT.IS.L.L$coefficients[1])^2+
         (coef(summary(GBT.IS.L.L))[2, "Std. Error"]/GBT.IS.L.L$coefficients[2])^2)


IS_Concentration_in_K1_2nM_B.all <- GBT.IS.L.all$coefficients[1]/GBT.IS.L.all$coefficients[2]
SE.IS_Concentration_in_K1_2nM_B.all <- IS_Concentration_in_K1_2nM_B.all* 
  sqrt((coef(summary(GBT.IS.L.all))[1, "Std. Error"]/GBT.IS.L.all$coefficients[1])^2+
         (coef(summary(GBT.IS.L.all))[2, "Std. Error"]/GBT.IS.L.all$coefficients[2])^2)

IS_Concentration_in_K2_2000nM_C.early <- GBT.IS.H.E$coefficients[1]/GBT.IS.H.E$coefficients[2]
IS_Concentration_in_K2_2000nM_C.late <- GBT.IS.H.L$coefficients[1]/GBT.IS.H.L$coefficients[2]
SE.IS_Concentration_in_K2_2000nM_C.early <- IS_Concentration_in_K2_2000nM_C.early* 
  sqrt((coef(summary(GBT.IS.H.E))[1, "Std. Error"]/GBT.IS.H.E$coefficients[1])^2+
         (coef(summary(GBT.IS.H.E))[2, "Std. Error"]/GBT.IS.H.E$coefficients[2])^2)
SE.IS_Concentration_in_K2_2000nM_C.late <- IS_Concentration_in_K2_2000nM_C.late* 
  sqrt((coef(summary(GBT.IS.H.L))[1, "Std. Error"]/GBT.IS.H.L$coefficients[1])^2+
         (coef(summary(GBT.IS.H.L))[2, "Std. Error"]/GBT.IS.H.L$coefficients[2])^2)

IS_Concentration_in_K2_2000nM_C.all <- GBT.IS.H.all$coefficients[1]/GBT.IS.H.all$coefficients[2]
SE.IS_Concentration_in_K2_2000nM_C.all <- IS_Concentration_in_K2_2000nM_C.all* 
  sqrt((coef(summary(GBT.IS.H.all))[1, "Std. Error"]/GBT.IS.H.all$coefficients[1])^2+
         (coef(summary(GBT.IS.H.all))[2, "Std. Error"]/GBT.IS.H.all$coefficients[2])^2)

## calculate IS GBT concentrations for all samples -----
##  Concentration in vial = Area/Slope
## error of concentration in vial = [Concentration] *  sqrt(([Error of Slope]/[Slope])^2+([Error of Intercept]/[Intercept])^2)
test <- full_join(data.raw, sample.list)
test2 <- test %>% 
  filter(Compound.Name=="13C-15N Glycine Betaine",
         !grepl("added",`File Text`)) %>%
  mutate(Treatment = ifelse(is.na(Treatment),mean(test$Treatment,na.rm=TRUE),Treatment),
         All.Models.High.Prediction.fmoles = rawArea/GBT.IS.H.all$coefficients[2]*400,
         SE.High.Prediction =  All.Models.High.Prediction.fmoles* 
           sqrt((coef(summary(GBT.IS.H.all))[1, "Std. Error"]/GBT.IS.H.all$coefficients[1])^2+
                  (coef(summary(GBT.IS.H.all))[2, "Std. Error"]/GBT.IS.H.all$coefficients[2])^2),
         All.Models.low.Prediction.fmoles = rawArea/GBT.IS.L.all$coefficients[2]*400,
         SE.low.Prediction =  All.Models.low.Prediction.fmoles* 
           sqrt((coef(summary(GBT.IS.L.all))[1, "Std. Error"]/GBT.IS.L.all$coefficients[1])^2+
                  (coef(summary(GBT.IS.L.all))[2, "Std. Error"]/GBT.IS.L.all$coefficients[2])^2),
         correct.model.prediction.fmoles = ifelse(`Injection Order`<60,
                                                  ifelse(Treatment<50,
                                                         rawArea/GBT.IS.L.E$coefficients[2]*400,
                                                         rawArea/GBT.IS.H.E$coefficients[2]*400),
                                                  ifelse(Treatment<50,
                                                         rawArea/GBT.IS.L.L$coefficients[2]*400,
                                                         rawArea/GBT.IS.H.L$coefficients[2]*400)),
         SE.correct.model.prediction = ifelse(`Injection Order`<60,
                                              ifelse(Treatment<50,
                                                     correct.model.prediction.fmoles* 
                                                       sqrt((coef(summary(GBT.IS.L.E))[1, "Std. Error"]/GBT.IS.L.E$coefficients[1])^2+
                                                              (coef(summary(GBT.IS.L.E))[2, "Std. Error"]/GBT.IS.L.E$coefficients[2])^2),
                                                     correct.model.prediction.fmoles* 
                                                       sqrt((coef(summary(GBT.IS.H.E))[1, "Std. Error"]/GBT.IS.H.E$coefficients[1])^2+
                                                              (coef(summary(GBT.IS.H.E))[2, "Std. Error"]/GBT.IS.H.E$coefficients[2])^2)),
                                              ifelse(Treatment<50,
                                                     correct.model.prediction.fmoles* 
                                                       sqrt((coef(summary(GBT.IS.L.L))[1, "Std. Error"]/GBT.IS.L.L$coefficients[1])^2+
                                                              (coef(summary(GBT.IS.L.L))[2, "Std. Error"]/GBT.IS.L.L$coefficients[2])^2),
                                                     correct.model.prediction.fmoles* 
                                                       sqrt((coef(summary(GBT.IS.H.L))[1, "Std. Error"]/GBT.IS.H.L$coefficients[1])^2+
                                                              (coef(summary(GBT.IS.H.L))[2, "Std. Error"]/GBT.IS.H.L$coefficients[2])^2)))) %>%
  select(Replicate.Name, Compound.Name, All.Models.High.Prediction.fmoles, 
         SE.High.Prediction, All.Models.low.Prediction.fmoles, SE.low.Prediction,correct.model.prediction.fmoles,
         SE.correct.model.prediction,
         Experiment, Treatment, Added_GBT_IS, Added_GBT, replicate, `Sample ID`)
## read in bio.infor ------
bio.info <- read_csv("../data/BiologicNormalizationInformation.csv")
test3 <- full_join(bio.info, test2) %>%
  mutate(H.Prediction.nmoles.per.L = ifelse(is.na(`Vol Filtered (L)`),
                                            All.Models.High.Prediction.fmoles/2*1e-6*10,
                                            All.Models.High.Prediction.fmoles/`Vol Filtered (L)`*1e-6*dilution.factor),
         L.Prediction.nmoles.per.L = ifelse(is.na(`Vol Filtered (L)`),
                                            All.Models.low.Prediction.fmoles/2*1e-6*10,
                                            All.Models.low.Prediction.fmoles/`Vol Filtered (L)`*1e-6*dilution.factor),
         correct.Prediction.nmoles.per.L = ifelse(is.na(`Vol Filtered (L)`),
                                            correct.model.prediction.fmoles/2*1e-6*10,
                                            correct.model.prediction.fmoles/`Vol Filtered (L)`*1e-6*dilution.factor),
         H.prediction.nmoles.per.L.per.hr = ifelse(`Time incubated` > 0, 
                                                    H.Prediction.nmoles.per.L/`Time incubated`*60,
                                                    H.Prediction.nmoles.per.L/40*60),
         L.Prediction.nmoles.per.L.per.hr = ifelse(`Time incubated` > 0, 
                                                    L.Prediction.nmoles.per.L/`Time incubated`*60,
                                                    L.Prediction.nmoles.per.L/40*60),
         correct.Prediction.nmoles.per.L.per.hr = ifelse(`Time incubated` > 0, 
                                                   correct.Prediction.nmoles.per.L/`Time incubated`*60,
                                                   correct.Prediction.nmoles.per.L/40*60))

## plot data ----
## boxplot of GBT------------
test <- full_join(data.raw, sample.list)
GBT.data <- test %>%
  filter(Compound.Name=="Glycine Betaine",
         Added_GBT==0,
         Added_GBT_IS==0)

ggplot(GBT.data) + 
  geom_boxplot(aes(x=as.factor(Treatment), y=rawArea, group = Treatment)) +
  ggtitle("Betaine") + facet_wrap(~Experiment, scales = "free_y" )+ 
  theme_bw()

## boxplot of heavy GBT------------
GBT.13C15N.data <- test3 %>%
  filter(Compound.Name=="13C-15N Glycine Betaine",
         Added_GBT_IS==0,
         Added_GBT==0) 

ggplot(GBT.13C15N.data) + 
  geom_boxplot(aes(x=as.factor(Treatment), y=H.prediction.nmoles.per.L.per.hr, group = Treatment)) +
  ggtitle("13C, 15N Betaine")+ facet_wrap(~Experiment, scales = "free_y")+
  theme_bw()
# ggsave(paste0("../output/","GBT.IS.rates.data.nM.per.Hr.pdf"))

GBT.13C15N.data <- test %>%
  filter(Compound.Name=="13C-15N Glycine Betaine",
         Added_GBT_IS==0,
         Added_GBT==0)

ggplot(GBT.13C15N.data) + 
  geom_boxplot(aes(x=as.factor(Treatment), y=rawArea, group = Treatment)) +
  ggtitle("13C, 15N Betaine")+ facet_wrap(~Experiment, scales = "free_y")+ scale_y_log10()+
  theme_bw()


# ggsave(paste0("../output/","GBT.IS.Data.pdf"))

## uptake rate ---------------------
GBT.rate <- test3 %>%
  filter(Compound.Name=="13C-15N Glycine Betaine",
         Added_GBT_IS==0,
         Added_GBT==0)%>%
  mutate(SE.correct.conc = (correct.Prediction.nmoles.per.L* (SE.correct.model.prediction/correct.model.prediction.fmoles)),
         SE.correct.rate =  ifelse(`Time incubated` > 0,(SE.correct.conc/`Time incubated`*60),(SE.correct.conc/40*60)),
         SE.L.conc = (L.Prediction.nmoles.per.L * (SE.low.Prediction/All.Models.low.Prediction.fmoles)),
         SE.L.rate =  ifelse(`Time incubated` > 0,(SE.L.conc/`Time incubated`*60),(SE.L.conc/40*60)),
         SE.H.conc = (H.Prediction.nmoles.per.L * (SE.High.Prediction/All.Models.High.Prediction.fmoles)),
         SE.H.rate =  ifelse(`Time incubated` > 0,(SE.H.conc/`Time incubated`*60),(SE.H.conc/40*60))) %>%
  select(Compound.Name, Replicate.Name, Treatment,
         Experiment, replicate,L.Prediction.nmoles.per.L,
         L.Prediction.nmoles.per.L.per.hr, SE.L.rate, 
         H.prediction.nmoles.per.L.per.hr, SE.H.rate, 
         correct.Prediction.nmoles.per.L.per.hr, SE.correct.rate,
         correct.Prediction.nmoles.per.L, SE.correct.conc,
         `Time incubated`) %>%
  filter(Experiment=="K1"|Experiment=="K2")

## write best concentration data --------
this.GBT.rate.to.write <- GBT.rate %>%
  dplyr::select(Compound.Name, Replicate.Name, Treatment,
         Experiment, replicate, 
         correct.Prediction.nmoles.per.L.per.hr, SE.correct.rate,
         correct.Prediction.nmoles.per.L, SE.correct.conc,
         `Time incubated`)
write_csv(this.GBT.rate.to.write, "../output/GBT_Concentrations_and_Uptake_Rates_w_StdError.csv")

GBT.rate <- GBT.rate %>% dplyr::select(-SE.correct.conc)
## plotting high rate bc low SE is high (model had a lot of uncertainty around the intercept)
log.rates <- ggplot(GBT.rate) + 
  geom_point(aes(x=Treatment, y=correct.Prediction.nmoles.per.L.per.hr),
             size = 3) +
  geom_errorbar(aes(x=Treatment,
                    ymin=correct.Prediction.nmoles.per.L.per.hr-SE.H.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.H.rate))+
  facet_wrap(~Experiment)+
  ggtitle("13C, 15N Betaine uptake") +
  theme_bw() +
  theme(text = element_text(size = 15))+
  scale_x_log10()
log.rates
rates <- ggplot(GBT.rate) + 
  geom_point(aes(x=Treatment, y=H.prediction.nmoles.per.L.per.hr), size = 3) +
  geom_errorbar(aes(x=Treatment, ymin=H.prediction.nmoles.per.L.per.hr-SE.H.rate,
                    ymax = H.prediction.nmoles.per.L.per.hr+SE.H.rate))+
  facet_wrap(~Experiment)+
  ggtitle("13C, 15N Betaine uptake") +
  theme_bw() +
  theme(text = element_text(size = 15))
rates

## make rates long and plot w facet--------
GBT.rates.long <- GBT.rate %>% 
  select(-L.Prediction.nmoles.per.L, -correct.Prediction.nmoles.per.L) %>%
  gather(Prediction, Value, -Compound.Name, -Replicate.Name, -Treatment, -Experiment, -replicate, -`Time incubated`) %>%
  mutate(key = ifelse(grepl("nmoles.per",Prediction),"Rate.Prediction","SE"),
         Which.estimate = ifelse(grepl("correct",Prediction), "Best match",
                                 ifelse(grepl("H.", Prediction), "High","Low"))) %>%
  select(-Prediction) %>%
  spread(key, Value)
rates.compare <- ggplot(GBT.rates.long) + 
  geom_point(aes(x=Treatment, y=Rate.Prediction), size = 3) +
  geom_errorbar(aes(x=Treatment, ymin=Rate.Prediction-SE,
                    ymax = Rate.Prediction+SE))+
  facet_grid(Which.estimate~Experiment)+
  ggtitle("13C, 15N Betaine uptake") +
  theme_bw() +
  theme(text = element_text(size = 15))+
  scale_x_log10()
rates.compare

ggsave(paste0("../output/","Comparing.quantification.estimates.uptake.rates.pdf"))

## get mean and SD and plot that ----------
GBT.rates.long.aves <- GBT.rates.long %>%
  mutate(percent.error = (SE/Rate.Prediction)^2) %>%
  group_by(Experiment, Treatment, Compound.Name, Which.estimate) %>%
  summarise(mean.pred.rate = mean(Rate.Prediction),
            sd.pred.rate = sd(Rate.Prediction),
            percent.error = sum(percent.error)) %>%
  mutate(percent.error = sqrt(percent.error),
         SE  = mean.pred.rate*percent.error) %>% filter(Which.estimate!="Low")

rates.aves<- ggplot(GBT.rates.long.aves) + 
  geom_point(aes(x=Treatment, y=mean.pred.rate, color = Which.estimate), size = 3) +
  geom_errorbar(aes(x=Treatment, ymin=mean.pred.rate-SE,
                    ymax = mean.pred.rate+SE,
                    color = Which.estimate))+
  facet_grid(Which.estimate~Experiment)+
  ggtitle("13C, 15N Betaine uptake") +
  theme_bw() +
  theme(text = element_text(size = 15))+
  scale_x_log10()
rates.aves
## guess rates K1------

## write out just best matched rates and SE ------
GBT.rates.to.write <- GBT.rate %>% 
  filter(!is.na(Treatment)) %>%
  dplyr::select(Compound.Name, Replicate.Name, Treatment, Experiment, replicate, 
                correct.Prediction.nmoles.per.L.per.hr, SE.correct.rate,
                `Time incubated`)
write.csv(GBT.rates.to.write, file = paste0("../output/GBT_Uptake_Rates_w_StdError.csv"))

## best prediction----
GBT.rate.fit.k1 <- GBT.rate %>% 
  filter(!is.na(Treatment),
         Experiment=="K1") %>%
  mutate(Treatment.nmoles.added = Treatment,
         fraction.taken.up= ifelse(is.na(correct.Prediction.nmoles.per.L),
                                   0, 
                                   (correct.Prediction.nmoles.per.L)/Treatment.nmoles.added),
         time.over.fraction = ifelse(fraction.taken.up==0, 0, `Time incubated`/60/fraction.taken.up))

GBT.rate.fit.aves.K1 <- GBT.rate.fit.k1 %>% 
  group_by(Treatment.nmoles.added, Experiment) %>%
  summarise(nM.GBT.Per.hr = mean(correct.Prediction.nmoles.per.L.per.hr))

## nlm k1 ------
x <- GBT.rate.fit.k1$Treatment.nmoles.added
y <-  GBT.rate.fit.k1$correct.Prediction.nmoles.per.L.per.hr
plot(x, y)

#from this graph set approximate starting values
a_start<-100 #param a is the Ks
b_start<-0.7 #b is Vmax
c_start <- 4
#model
m.k1<-nls(y~(b*(x+c))/(a+x+c),start=list(a=a_start,b=b_start, c= c_start),
         control = nls.control(printEval = T))
#get some estimation of goodness of fit
y.guess <- (b_start/(x+c_start))/(a_start+x+c_start)
cor(y,predict(m.k1))
#plot the fit
plot(x, y)
points(x=x,y=predict(m.k1),col="red",lty=2,lwd=3)
title("K1")
text(x = 1000, y = 0.03, labels = paste0("Ks = ", round(m.k1$m$getPars()["a"],2),"; Vmax =  ", round(m.k1$m$getPars()["b"],3)))
## plot fit fancy
GBT.rate.fit.k1$predict = predict(m.k1)
ggplot(data = GBT.rate.fit.k1) +
  geom_point(aes(x = Treatment.nmoles.added, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment.nmoles.added, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate), color = 'black') +
  geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()

m.k1
summary(m.k1)

## nlm k1 w/o explicit d[gbt] ------
x <- GBT.rate.fit.k1$Treatment.nmoles.added
y <-  GBT.rate.fit.k1$correct.Prediction.nmoles.per.L.per.hr
plot(x, y)

#from this graph set approximate starting values
a_start<-100 #param a is the Ks + Sn
b_start<-0.7 #b is Vmax
#model
m.k1.new<-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
          control = nls.control(printEval = T))
#get some estimation of goodness of fit
y.guess <- (b_start/(x))/(a_start+x)
cor(y,predict(m.k1.new))
#plot the fit
plot(x, y)
points(x=x,y=predict(m.k1.new),col="red",lty=2,lwd=3)
title("K1")
text(x = 1000, y = 0.03, labels = paste0("Ks = ", round(m.k1.new$m$getPars()["a"],2),"; Vmax =  ", round(m.k1.new$m$getPars()["b"],3)))
## plot fit fancy
GBT.rate.fit.k1$predict = predict(m.k1.new)
ggplot(data = GBT.rate.fit.k1) +
  geom_point(aes(x = Treatment.nmoles.added, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment.nmoles.added, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate), color = 'black') +
  geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()

m.k1.new
summary(m.k1.new)

## Wright-Hobbie linear transformation k1-------

plot(GBT.rate.fit.k1$Treatment.nmoles.added, GBT.rate.fit.k1$time.over.fraction)
m2.k1 <- lm( time.over.fraction~Treatment.nmoles.added, data = GBT.rate.fit.k1)
lines(GBT.rate.fit.k1$Treatment.nmoles.added, 
      predict(m2.k1),
      col="blue",lty=2,lwd=3)
title("K1 Wright-Hobbie linear transformation")
m2.k1
cor(GBT.rate.fit.k1$time.over.fraction,predict(m2.k1))
summary(m2.k1)

## guess rates K1 without 2000 nM------------------
## best prediction----
GBT.rate.fit.k1.no2000 <- GBT.rate %>% 
  filter(!is.na(Treatment),
         Treatment!=2000,
         Experiment=="K1") %>%
  mutate(Treatment.nmoles.added = Treatment,
         fraction.taken.up= ifelse(is.na(correct.Prediction.nmoles.per.L),
                                   0, 
                                   (correct.Prediction.nmoles.per.L)/Treatment.nmoles.added),
         time.over.fraction = ifelse(fraction.taken.up==0, 0, `Time incubated`/60/fraction.taken.up))

GBT.rate.fit.aves.K1.no2000 <- GBT.rate.fit.k1.no2000 %>% 
  group_by(Treatment.nmoles.added, Experiment) %>%
  summarise(nM.GBT.Per.hr = mean(correct.Prediction.nmoles.per.L.per.hr))

## nlm k1 ------
x <- GBT.rate.fit.k1.no2000$Treatment.nmoles.added
y <-  GBT.rate.fit.k1.no2000$correct.Prediction.nmoles.per.L.per.hr
plot(x, y)
#from this graph set approximate starting values
a_start<-50 #param a is the Ks
b_start<-0.7 #b is Vmax
c_start <- 4
#model
m.k1.no2000<-nls(y~(b*(x+c))/(a+x+c),start=list(a=a_start,b=b_start, c= c_start))
#get some estimation of goodness of fit
y.guess <- (b_start/(x+c_start))/(a_start+x+c_start)
cor(y,predict(m.k1.no2000))
#plot the fit
plot(x, y)
points(x=x,y=predict(m.k1.no2000),col="red",lty=2,lwd=3)
title("K1")
text(x = 100, y = 0.03, labels = paste0("Ks = ", round(m.k1.no2000$m$getPars()["a"],2),
                                        "; Vmax =  ", round(m.k1.no2000$m$getPars()["b"],3)))
## plot fit fancy
GBT.rate.fit.k1.no2000$predict = predict(m.k1.no2000)
ggplot(data = GBT.rate.fit.k1.no2000) +
  geom_point(aes(x = Treatment.nmoles.added, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment.nmoles.added, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate), color = 'black') +
  geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()

m.k1.no2000
summary(m.k1.no2000)
## Wright-Hobbie linear transformation k1-------

plot(GBT.rate.fit.k1.no2000$Treatment.nmoles.added, GBT.rate.fit.k1.no2000$time.over.fraction)
m2.k1.no2000 <- lm( time.over.fraction~Treatment.nmoles.added, data = GBT.rate.fit.k1.no2000)
lines(GBT.rate.fit.k1.no2000$Treatment.nmoles.added, 
      predict(m2.k1.no2000),
      col="blue",lty=2,lwd=3)
title("K1 Wright-Hobbie linear transformation")
m2.k1.no2000
cor(GBT.rate.fit.k1.no2000$time.over.fraction,predict(m2.k1.no2000))
summary(m2.k1.no2000)

## save out the K1 models -----
save(m.k1, file = paste0("../output/GBT K1 Model.rdata"))
save(m.k1.new, file = paste0("../output/GBT K1 Model no explicit dGBT.rdata"))
save(m.k1.no2000, file = paste0("../output/GBT K1 Model no 2000.rdata"))
m.k1.no2000$m$getAllPars()

one.model.sum <- data.frame(summary(m.k1.no2000)$coefficients)
one.model.sum$coefficient <- rownames(one.model.sum)
one.model.sum$Model <- "K1.No.2000"
two.model.sum <- data.frame(summary(m.k1)$coefficients)
two.model.sum$coefficient <- rownames(two.model.sum)
two.model.sum$Model <- "K1"
three.model.sum <- data.frame(summary(m.k1.new)$coefficients)
three.model.sum$coefficient <- rownames(three.model.sum)
three.model.sum$Model <- "K1.with.no.explicit.dGBT"


both.models.sum <- full_join(one.model.sum, two.model.sum) %>% full_join(.,three.model.sum)
write_csv(both.models.sum, path = paste0("../output/GBT K1 Model coeff.csv"))

## guess rates K2------
## best prediction----
GBT.rate.fit.k2 <- GBT.rate %>% 
  filter(!is.na(Treatment),
         Experiment=="K2") %>%
  mutate(Treatment.nmoles.added = Treatment,
         fraction.taken.up= ifelse(is.na(correct.Prediction.nmoles.per.L),
                                   0, 
                                   (correct.Prediction.nmoles.per.L)/Treatment.nmoles.added),
         time.over.fraction = ifelse(fraction.taken.up==0, 0, `Time incubated`/60/fraction.taken.up))

GBT.rate.fit.aves.k2 <- GBT.rate.fit.k2 %>% 
  group_by(Treatment.nmoles.added, Experiment) %>%
  summarise(nM.GBT.Per.hr = mean(correct.Prediction.nmoles.per.L.per.hr))

## nlm k2 ------
x <- GBT.rate.fit.k2$Treatment.nmoles.added
y <-  GBT.rate.fit.k2$correct.Prediction.nmoles.per.L.per.hr
plot(x, y)
#from this graph set approximate starting values
a_start<-150 #param a is the Ks
b_start<-7 #b is Vmax
c_start <- 14
#model
m.k2<-nls(y~(b*(x+c))/(a+x+c),start=list(a=a_start,b=b_start, c= c_start))
#get some estimation of goodness of fit
y.guess <- (b_start/(x+c_start))/(a_start+x+c_start)
cor(y,predict(m.k2))
#plot the fit
plot(x, y)
points(x=x,y=predict(m.k2),col="red",lty=2,lwd=3)
title("K2")
text(x = 1000, y = 0.03, labels = paste0("Ks = ", round(m.k2$m$getPars()["a"],2),"; Vmax =  ", round(m.k2$m$getPars()["b"],3)))
m.k2
summary(m.k2)

## nlm k2 w/o explicit d[gbt] ------
x <- GBT.rate.fit.k2$Treatment.nmoles.added
y <-  GBT.rate.fit.k2$correct.Prediction.nmoles.per.L.per.hr
plot(x, y)

#from this graph set approximate starting values
a_start<-150 #param a is the Ks + Sn
b_start<-7 #b is Vmax
#model
m.k2.new<-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start),
              control = nls.control(printEval = T))
#get some estimation of goodness of fit
y.guess <- (b_start/(x))/(a_start+x)
cor(y,predict(m.k2.new))
#plot the fit
plot(x, y)
points(x=x,y=predict(m.k2.new),col="red",lty=2,lwd=3)
title("K2")
text(x = 1000, y = 0.03, labels = paste0("Ks = ", round(m.k2.new$m$getPars()["a"],2),"; Vmax =  ", round(m.k2.new$m$getPars()["b"],3)))
## plot fit fancy
GBT.rate.fit.k2$predict = predict(m.k2.new)
ggplot(data = GBT.rate.fit.k2) +
  geom_point(aes(x = Treatment.nmoles.added, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment.nmoles.added, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate), color = 'black') +
  geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()

m.k1.new
summary(m.k1.new)


## Wright-Hobbie linear transformation k2-------

plot(GBT.rate.fit.k2$Treatment.nmoles.added, GBT.rate.fit.k2$time.over.fraction)
m2.k2 <- lm( time.over.fraction~Treatment.nmoles.added, data = GBT.rate.fit.k2)
lines(GBT.rate.fit.k2$Treatment.nmoles.added, 
      predict(m2.k2),
      col="blue",lty=2,lwd=3)
title("K2 Wright-Hobbie linear transformation")
m2.k2
cor(GBT.rate.fit.k2$time.over.fraction,predict(m2.k2))
summary(m2.k2)


## guess rates K2 w/o 2000 nM point ------
## best prediction----
GBT.rate.fit.k2.no2000 <- GBT.rate %>% 
  filter(!is.na(Treatment),
         Treatment!=2000,
         Experiment=="K2") %>%
  mutate(Treatment.nmoles.added = Treatment,
         fraction.taken.up= ifelse(is.na(correct.Prediction.nmoles.per.L),
                                   0, 
                                   (correct.Prediction.nmoles.per.L)/Treatment.nmoles.added),
         time.over.fraction = ifelse(fraction.taken.up==0, 0, `Time incubated`/60/fraction.taken.up))

GBT.rate.fit.aves.k2.no2000 <- GBT.rate.fit.k2.no2000 %>% 
  group_by(Treatment.nmoles.added, Experiment) %>%
  summarise(nM.GBT.Per.hr = mean(correct.Prediction.nmoles.per.L.per.hr))

## nlm k2 ------
x <- GBT.rate.fit.k2.no2000$Treatment.nmoles.added
y <-  GBT.rate.fit.k2.no2000$correct.Prediction.nmoles.per.L.per.hr
plot(x, y)
#from this graph set approximate starting values
a_start<-150 #param a is the Ks
b_start<-7 #b is Vmax
c_start <- 14
#model
m.k2.no2000<-nls(y~(b*(x+c))/(a+x+c),start=list(a=a_start,b=b_start, c= c_start))
#get some estimation of goodness of fit
y.guess <- (b_start/(x+c_start))/(a_start+x+c_start)
cor(y,predict(m.k2.no2000))
#plot the fit
plot(x, y)
points(x=x,y=predict(m.k2.no2000),col="red",lty=2,lwd=3)
title("K2")
text(x = 100, y = 0.03, labels = paste0("Ks = ", round(m.k2.no2000$m$getPars()["a"],2),
                                        "; Vmax =  ", round(m.k2.no2000$m$getPars()["b"],3)))
m.k2.no2000
summary(m.k2.no2000)

## Wright-Hobbie linear transformation k2-------

plot(GBT.rate.fit.k2.no2000$Treatment.nmoles.added, GBT.rate.fit.k2.no2000$time.over.fraction)
m2.k2.no2000 <- lm( time.over.fraction~Treatment.nmoles.added, data = GBT.rate.fit.k2.no2000)
lines(GBT.rate.fit.k2.no2000$Treatment.nmoles.added, 
      predict(m2.k2.no2000),
      col="blue",lty=2,lwd=3)
title("K2 Wright-Hobbie linear transformation")
m2.k2.no2000
cor(GBT.rate.fit.k2.no2000$time.over.fraction,
    predict(m2.k2.no2000))
summary(m2.k2.no2000)

## save out the K2 models -----
save(m.k2, file = paste0("../output/GBT K2 Model.rdata"))
save(m.k2.new, file = paste0("../output/GBT K2 Model no explicit dGBT.rdata"))
save(m.k2.no2000, file = paste0("../output/GBT K2 Model no 2000.rdata"))

one.model.sum <- data.frame(summary(m.k2.no2000)$coefficients)
one.model.sum$coefficient <- rownames(one.model.sum)
one.model.sum$Model <- "K2.No.2000"
two.model.sum <- data.frame(summary(m.k2)$coefficients)
two.model.sum$coefficient <- rownames(two.model.sum)
two.model.sum$Model <- "K2"
three.model.sum <- data.frame(summary(m.k2.new)$coefficients)
three.model.sum$coefficient <- rownames(three.model.sum)
three.model.sum$Model <- "K2.with.no.explicit.dGBT"

both.models.sum.k2 <- full_join(one.model.sum, two.model.sum) %>% full_join(.,three.model.sum)
write_csv(both.models.sum.k2, path = paste0("../output/GBT K2 Model coeff.csv"))

