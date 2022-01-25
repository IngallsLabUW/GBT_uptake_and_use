## Make a tidy plot of the models

##load packages-----------------------
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(stats)
library(kableExtra)
library(xtable)
library(gt)

## read in data ---------
uptake.rates <- read_csv("../output/GBT_Uptake_Rates_w_StdError.csv")

## read in curves --------
K1.curve <- read_csv("../output/K1-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv")
K2.curve <- read_csv("../output/K2-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv")
K1.orig.curve <- read_csv("../output/GBT K1 Model coeff.csv")
K2.orig.curve <- read_csv("../output/GBT K2 Model coeff.csv")

K1.curve.wdGBT <- read_csv("../output/K1-MonteCarlo_OutliersRemoved_summary_withError.csv")
K2.curve.wdGBT <- read_csv("../output/K2-MonteCarlo_OutliersRemoved_summary_withError.csv")

## wright-hobbie
K1.WH = read_csv("../output/K1-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv")
K2.WH = read_csv("../output/k2-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv")

## plot K1 -----
K1.plot <- uptake.rates %>%
  filter(Experiment=="K1")
a.k1.mc = K1.curve$Mean.Params[K1.curve$coeff=="a" & K1.curve$value.name=="Estimate"]
a.k1.mc.sd = K1.curve$sd[K1.curve$coeff=="a" & K1.curve$value.name=="Estimate"]
b.k1.mc = K1.curve$Mean.Params[K1.curve$coeff=="b" & K1.curve$value.name=="Estimate"]
b.k1.mc.sd = K1.curve$sd[K1.curve$coeff=="b" & K1.curve$value.name=="Estimate"]
# c.k1.mc = K1.curve$Mean.Params[K1.curve$coeff=="c" & K1.curve$value.name=="Estimate"]
# c.k1.mc.sd = K1.curve$sd[K1.curve$coeff=="c" & K1.curve$value.name=="Estimate"]

K1.lab = paste0("Km + Sn = ",round(a.k1.mc,2), 
                " (",round(K1.curve$sd[K1.curve$coeff=="a" & K1.curve$value.name=="Estimate"],2),")",
                "\n Vmax = ",round(b.k1.mc,2),
                " (",round(K1.curve$sd[K1.curve$coeff=="b" & K1.curve$value.name=="Estimate"],2),")" #,
                # "\n d[GBT] = ",round(c.k1.mc,2),
                # " (",round(K1.curve$sd[K1.curve$coeff=="c" & K1.curve$value.name=="Estimate"],2),")"
)



gbt.a.k1.mc = K1.curve.wdGBT$Mean.Params[K1.curve.wdGBT$coeff=="a" & K1.curve.wdGBT$value.name=="Estimate"]
gbt.a.k1.mc.sd = K1.curve.wdGBT$sd[K1.curve.wdGBT$coeff=="a" & K1.curve.wdGBT$value.name=="Estimate"]
gbt.b.k1.mc = K1.curve.wdGBT$Mean.Params[K1.curve.wdGBT$coeff=="b" & K1.curve.wdGBT$value.name=="Estimate"]
gbt.b.k1.mc.sd = K1.curve.wdGBT$sd[K1.curve.wdGBT$coeff=="b" & K1.curve.wdGBT$value.name=="Estimate"]
gbt.c.k1.mc = K1.curve.wdGBT$Mean.Params[K1.curve.wdGBT$coeff=="c" & K1.curve.wdGBT$value.name=="Estimate"]
gbt.c.k1.mc.sd = K1.curve.wdGBT$sd[K1.curve.wdGBT$coeff=="c" & K1.curve.wdGBT$value.name=="Estimate"]
K1.lab.gbt = paste0("Km + Sn = ",round(gbt.a.k1.mc,2), 
                " (",round(K1.curve.wdGBT$sd[K1.curve.wdGBT$coeff=="a" & K1.curve.wdGBT$value.name=="Estimate"],2),")",
                "\n Vmax = ",round(gbt.b.k1.mc,2),
                " (",round(K1.curve.wdGBT$sd[K1.curve.wdGBT$coeff=="b" & K1.curve.wdGBT$value.name=="Estimate"],2),")" ,
                "\n d[GBT] = ",round(gbt.c.k1.mc,2),
                " (",round(K1.curve.wdGBT$sd[K1.curve.wdGBT$coeff=="c" & K1.curve.wdGBT$value.name=="Estimate"],2),")"
)


a.k1.OF = K1.orig.curve$Estimate[K1.orig.curve$coefficient=="a" & K1.orig.curve$Model=="K1.with.no.explicit.dGBT"]
b.k1.OF = K1.orig.curve$Estimate[K1.orig.curve$coefficient=="b" & K1.orig.curve$Model=="K1.with.no.explicit.dGBT"]
# c.k1.OF = K1.orig.curve$Estimate[K1.orig.curve$coefficient=="c" & K1.orig.curve$Model=="K1"]

# K1.function.MC.params <- function(x) (b.k1.mc*(x+c.k1.mc))/(a.k1.mc+x+c.k1.mc) ## best fit, actual estimate 
K1.function.MC.params <- function(x) (b.k1.mc*(x))/(a.k1.mc+x) ## best fit, actual estimate, no d[gbt] 

##options when I had explicit d[gbt]
# K1.function.MC.params.option1 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc+a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## plus all SD
# K1.function.MC.params.option2 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x+(c.k1.mc-c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc-c.k1.mc.sd)) ## minus all SD
# K1.function.MC.params.option3 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x+(c.k1.mc-c.k1.mc.sd)))/((a.k1.mc+a.k1.mc.sd)+x+(c.k1.mc-c.k1.mc.sd)) ## a: plus SD b: minus SD, c: minus SD ** keep this one as a low estimate
# K1.function.MC.params.option4 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x+(c.k1.mc-c.k1.mc.sd)))/((a.k1.mc+a.k1.mc.sd)+x+(c.k1.mc-c.k1.mc.sd)) ## a: plus SD b: plus SD, c: minus SD  
# K1.function.MC.params.option5 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc+a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## a: plus SD b: minus SD, c: plus SD 
# K1.function.MC.params.option6 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## a: minus SD b: plus SD, c: plus SD ** keep this one as a high estimate
# K1.function.MC.params.option7 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## a: minus SD b: minus SD, c: plus SD 
# K1.function.MC.params.option8 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x+(c.k1.mc-c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc-c.k1.mc.sd)) ## a: minus SD b: plus SD, c: minus SD 

## options when I don't have explicit d[GBT]
K1.function.MC.params.option1 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x))/((a.k1.mc+a.k1.mc.sd)+x) ## + a.sd + b.sd
K1.function.MC.params.option2 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x))/((a.k1.mc-a.k1.mc.sd)+x) ## - a.sd -b.sd
K1.function.MC.params.option3 <-  function(x) ((b.k1.mc-b.k1.mc.sd)*(x))/((a.k1.mc+a.k1.mc.sd)+x) ## + a.sd -b.sd
K1.function.MC.params.option4 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x))/((a.k1.mc-a.k1.mc.sd)+x) ## - a.sd + b.sd
# K1.function.MC.params.option5 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc+a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## a: plus SD b: minus SD, c: plus SD 
# K1.function.MC.params.option6 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## a: minus SD b: plus SD, c: plus SD ** keep this one as a high estimate
# K1.function.MC.params.option7 <- function(x) ((b.k1.mc-b.k1.mc.sd)*(x+(c.k1.mc+c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc+c.k1.mc.sd)) ## a: minus SD b: minus SD, c: plus SD 
# K1.function.MC.params.option8 <- function(x) ((b.k1.mc+b.k1.mc.sd)*(x+(c.k1.mc-c.k1.mc.sd)))/((a.k1.mc-a.k1.mc.sd)+x+(c.k1.mc-c.k1.mc.sd)) ## a: minus SD b: plus SD, c: minus SD 

# K1.function.original.fit.params <- function(x) (b.k1.OF*(x+c.k1.OF))/(a.k1.OF+x+c.k1.OF)
K1.function.original.fit.params <- function(x) (b.k1.OF*(x))/(a.k1.OF+x)

K1.subset.of.replicate.inj <- K1.plot %>%
  filter(grepl("190730_Smp_GBT-K1-M2nM_B", Replicate.Name))

ggplot(data = K1.plot) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    width = 5
  ),
  color = 'black') +
  geom_point(data = K1.subset.of.replicate.inj,
             aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr),
             color = 'black', fill="white", shape = 21) +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()+ 
  stat_function(fun = K1.function.MC.params,
                color = "blue") + 
  # stat_function(fun = K1.function.original.fit.params,
  #               color = "red") + 
  xlim(0,2000) +
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  ggtitle("K1")+
  geom_label(aes(x = 1000, y = 0.2, label = K1.lab))


# ggsave(filename = "../output/Uptake-Kinetics-With-Model-Values-K1.pdf",
       # width = 4, height = 4)

## plot with all of the possible functions:
ggplot(data = K1.plot) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    width = 5
  ),
  color = 'black') +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()+ 
  stat_function(fun = K1.function.MC.params,  color = "red") +
  stat_function(fun = K1.function.MC.params.option1,  color = "orange") +
  stat_function(fun = K1.function.MC.params.option2,  color = "green") +
  stat_function(fun = K1.function.MC.params.option3,  color = "black") + #keep this one as a low estimate
  stat_function(fun = K1.function.MC.params.option4,  color = "purple") + #keep this one as a high estimate
  # stat_function(fun = K1.function.MC.params.option5,  color = "yellow") + 
  # stat_function(fun = K1.function.MC.params.option6,  color = "black") + #keep this one as a high estimate
  # stat_function(fun = K1.function.MC.params.option7,  color = "grey") + 
  # stat_function(fun = K1.function.MC.params.option8,  color = "blue") + 
  
  # stat_function(fun = K1.function.original.fit.params,
  #               color = "red") + 
  xlim(0,2000) +
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  ggtitle("K1")+
  geom_label(aes(x = 1000, y = 0.2, label = K1.lab))


## plot K2 -----
K2.plot <- uptake.rates %>%
  filter(Experiment=="K2")
a.K2.mc = K2.curve$Mean.Params[K2.curve$coeff=="a" & K2.curve$value.name=="Estimate"]
a.K2.mc.sd = K2.curve$sd[K2.curve$coeff=="a" & K2.curve$value.name=="Estimate"]
b.K2.mc = K2.curve$Mean.Params[K2.curve$coeff=="b" & K2.curve$value.name=="Estimate"]
b.K2.mc.sd = K2.curve$sd[K2.curve$coeff=="b" & K2.curve$value.name=="Estimate"]
# c.K2.mc = K2.curve$Mean.Params[K2.curve$coeff=="c" & K2.curve$value.name=="Estimate"]
# c.K2.mc.sd = K2.curve$sd[K2.curve$coeff=="c" & K2.curve$value.name=="Estimate"]


K2.lab = paste0("Km + Sn = ",round(a.K2.mc,2), 
                " (",round(K2.curve$sd[K2.curve$coeff=="a" & K2.curve$value.name=="Estimate"],2),")",
                "\n Vmax = ",round(b.K2.mc,2),
                " (",round(K2.curve$sd[K2.curve$coeff=="b" & K2.curve$value.name=="Estimate"],2),")"#,
                # "\n d[GBT] = ",round(c.K2.mc,2),
                # " (",round(K2.curve$sd[K2.curve$coeff=="c" & K2.curve$value.name=="Estimate"],2),")"
)

gbt.a.K2.mc = K2.curve.wdGBT$Mean.Params[K2.curve.wdGBT$coeff=="a" & K2.curve.wdGBT$value.name=="Estimate"]
gbt.a.K2.mc.sd = K2.curve.wdGBT$sd[K2.curve.wdGBT$coeff=="a" & K2.curve.wdGBT$value.name=="Estimate"]
gbt.b.K2.mc = K2.curve.wdGBT$Mean.Params[K2.curve.wdGBT$coeff=="b" & K2.curve.wdGBT$value.name=="Estimate"]
gbt.b.K2.mc.sd = K2.curve.wdGBT$sd[K2.curve.wdGBT$coeff=="b" & K2.curve.wdGBT$value.name=="Estimate"]
gbt.c.K2.mc = K2.curve.wdGBT$Mean.Params[K2.curve.wdGBT$coeff=="c" & K2.curve.wdGBT$value.name=="Estimate"]
gbt.c.K2.mc.sd = K2.curve.wdGBT$sd[K2.curve.wdGBT$coeff=="c" & K2.curve.wdGBT$value.name=="Estimate"]


K2.lab.gbt = paste0("Km + Sn = ",round(gbt.a.K2.mc,2), 
                " (",round(K2.curve.wdGBT$sd[K2.curve.wdGBT$coeff=="a" & K2.curve.wdGBT$value.name=="Estimate"],2),")",
                "\n Vmax = ",round(gbt.b.K2.mc,2),
                " (",round(K2.curve.wdGBT$sd[K2.curve.wdGBT$coeff=="b" & K2.curve.wdGBT$value.name=="Estimate"],2),")",
                "\n d[GBT] = ",round(gbt.c.K2.mc.sd,2),
                " (",round(K2.curve.wdGBT$sd[K2.curve.wdGBT$coeff=="c" & K2.curve.wdGBT$value.name=="Estimate"],2),")"
)

a.K2.OF = K2.orig.curve$Estimate[K2.orig.curve$coefficient=="a" & K2.orig.curve$Model=="K2.with.no.explicit.dGBT"]
b.K2.OF = K2.orig.curve$Estimate[K2.orig.curve$coefficient=="b" & K2.orig.curve$Model=="K2.with.no.explicit.dGBT"]
# c.K2.OF = K2.orig.curve$Estimate[K2.orig.curve$coefficient=="c" & K2.orig.curve$Model=="K2"]

# K2.function.MC.params <- function(x) (b.K2.mc*(x+c.K2.mc))/(a.K2.mc+x+c.K2.mc)
# K2.function.MC.params.option1 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc+a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## plus all SD
# K2.function.MC.params.option2 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x+(c.K2.mc-c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc-c.K2.mc.sd)) ## minus all SD
# K2.function.MC.params.option3 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x+(c.K2.mc-c.K2.mc.sd)))/((a.K2.mc+a.K2.mc.sd)+x+(c.K2.mc-c.K2.mc.sd)) ## a: plus SD b: minus SD, c: minus SD ** keep this one as a low estimate
# K2.function.MC.params.option4 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x+(c.K2.mc-c.K2.mc.sd)))/((a.K2.mc+a.K2.mc.sd)+x+(c.K2.mc-c.K2.mc.sd)) ## a: plus SD b: plus SD, c: minus SD  
# K2.function.MC.params.option5 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc+a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## a: plus SD b: minus SD, c: plus SD 
# K2.function.MC.params.option6 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## a: minus SD b: plus SD, c: plus SD ** keep this one as a high estimate
# K2.function.MC.params.option7 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## a: minus SD b: minus SD, c: plus SD 
# K2.function.MC.params.option8 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x+(c.K2.mc-c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc)) ## a: minus SD b: plus SD, c: minus SD 

## version w.o explicit d[GBT]
K2.function.MC.params <- function(x) (b.K2.mc*(x))/(a.K2.mc+x)
K2.function.MC.params.option1 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x))/((a.K2.mc+a.K2.mc.sd)+x) ## plus all SD
K2.function.MC.params.option2 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x))/((a.K2.mc-a.K2.mc.sd)+x) ## minus all SD
K2.function.MC.params.option3 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x))/((a.K2.mc+a.K2.mc.sd)+x) ##  + a.sd -b.sd
K2.function.MC.params.option4 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x))/((a.K2.mc-a.K2.mc.sd)+x) ##  - a.sd +b.sd
# K2.function.MC.params.option5 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc+a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## a: plus SD b: minus SD, c: plus SD 
# K2.function.MC.params.option6 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## a: minus SD b: plus SD, c: plus SD ** keep this one as a high estimate
# K2.function.MC.params.option7 <- function(x) ((b.K2.mc-b.K2.mc.sd)*(x+(c.K2.mc+c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc+c.K2.mc.sd)) ## a: minus SD b: minus SD, c: plus SD 
# K2.function.MC.params.option8 <- function(x) ((b.K2.mc+b.K2.mc.sd)*(x+(c.K2.mc-c.K2.mc.sd)))/((a.K2.mc-a.K2.mc.sd)+x+(c.K2.mc)) ## a: minus SD b: plus SD, c: minus SD 


# K2.function.original.fit.params <- function(x) (b.K2.OF*(x+c.K2.OF))/(a.K2.OF+x+c.K2.OF)
K2.function.original.fit.params <- function(x) (b.K2.OF*(x))/(a.K2.OF+x)


K2.subset.of.replicate.inj <- K2.plot %>%
  filter(grepl("190730_Smp_GBT-K2-M2000nM_C", Replicate.Name))

ggplot(data = K2.plot) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    # width = 5
  ),
  color = 'black') +
  geom_point(data = K2.subset.of.replicate.inj,
             aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr),
             color = 'black', fill="white", shape = 21) +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()+ 
  stat_function(fun = K2.function.MC.params,
                color = "blue") + 
  # stat_function(fun = K2.function.original.fit.params,
  #               color = "red") +
  xlim(0,2000) +
  # scale_x_log10(limits = c(0.00001,2000))+
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  ggtitle("K2") +
  geom_label(aes(x = 1000, y = 0.2, label = K2.lab))

# ggsave(filename = "../output/Uptake-Kinetics-With-Model-Values-K2.pdf",
       # width = 4, height = 4)

## plot K2 with all model possibilities
ggplot(data = K2.plot) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    # width = 5
  ),
  color = 'black') +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()+ 
  stat_function(fun = K2.function.MC.params,  color = "red") +
  stat_function(fun = K2.function.MC.params.option1,  color = "orange") +
  stat_function(fun = K2.function.MC.params.option2,  color = "green") +
  stat_function(fun = K2.function.MC.params.option3,  color = "black") + #keep this one as a low estimate
  stat_function(fun = K2.function.MC.params.option4,  color = "purple") + #keep this one as a high estimate
  # stat_function(fun = K2.function.MC.params.option5,  color = "orange") +
  # stat_function(fun = K2.function.MC.params.option6,  color = "black") + #keep this one as a high estimate
  # stat_function(fun = K2.function.MC.params.option7,  color = "grey") +
  # stat_function(fun = K2.function.MC.params.option8,  color = "blue") +
  # xlim(0,2000) +
  # scale_x_log10(limits = c(0.00001,2000))+
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  ggtitle("K2") +
  geom_label(aes(x = 1000, y = 0.2, label = K2.lab))



## two plots together --------
both.plots <- full_join(K1.plot, K2.plot) %>%
  mutate(Experiment = factor(Experiment, levels = c("K2","K1")))
both.plots.subset.replicate.inj <- full_join(K1.subset.of.replicate.inj, K2.subset.of.replicate.inj)%>%
  mutate(Experiment = factor(Experiment, levels = c("K2","K1")))

k1.fit <- data.frame(x = seq(0,2000,length = 1000), 
                     y = K1.function.MC.params(seq(0,2000,length = 1000)),
                     Experiment = "K1",
                     y.high = K1.function.MC.params.option4(seq(0,2000,length = 1000)),
                     y.low = K1.function.MC.params.option3(seq(0,2000,length = 1000)))  %>%
  mutate(Experiment = factor(Experiment, levels = c("K2","K1")))
k2.fit <- data.frame(x = seq(0,2000,length = 1000), 
                     y = K2.function.MC.params(seq(0,2000,length = 1000)),
                     Experiment = "K2",
                     y.high = K2.function.MC.params.option4(seq(0,2000,length = 1000)),
                     y.low = K2.function.MC.params.option3(seq(0,2000,length = 1000)))  %>%
  mutate(Experiment = factor(Experiment, levels = c("K2","K1")))

ggplot(data = both.plots) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    # width = 5
  ),
  color = 'black') +
  geom_point(data = both.plots.subset.replicate.inj,
             aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr),
             color = 'black', fill="white", shape = 21) +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()+ 
  facet_wrap(~Experiment) +
  geom_line(data = k1.fit, aes(x = x, y = y),
                color = "blue") +
  geom_line(data = k2.fit, aes(x = x, y = y),
            color = "red") +
  # scale_x_log10(limits = c(0.00001,2000))+
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  theme(aspect.ratio = 0.5)

# ggsave(filename = "../output/Uptake-Kinetics-K1andK2.pdf",
       # width = 7, height = 4)


Exp.labs <-  c("North", "South")
names(Exp.labs) <- c("K1", "K2")


total.uptake.plot <- ggplot(data = both.plots) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    # width = 5
  ),
  color = 'black') +
  geom_point(data = both.plots.subset.replicate.inj,
             aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr),
             color = 'black', fill="white", shape = 21) +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  theme_bw()+ 
  facet_wrap(~Experiment,
             labeller = labeller(Experiment = Exp.labs)) +
  geom_ribbon(data = k1.fit, aes(x = x, ymin = y.low, ymax= y.high),
            fill = "blue", alpha = 0.2)+
  geom_ribbon(data = k2.fit, aes(x = x, ymin = y.low, ymax= y.high),
              fill = "red", alpha = 0.2)+
  geom_line(data = k1.fit, aes(x = x, y = y),
            color = "blue") +
  geom_line(data = k2.fit, aes(x = x, y = y),
            color = "red") +
  geom_vline(data = data.frame(Experiment = factor(c("K2","K1"),levels = c("K2","K1")),
                               k.5 = c(11.06, 79.46) ),
             aes(xintercept  = k.5),
             linetype = "dashed") +
  # scale_x_log10(limits = c(0.00001,2000))+
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  theme(aspect.ratio = 0.5)
total.uptake.plot

# ggsave(total.uptake.plot, filename = "../output/Uptake-Kinetics-with-rangeOfFit-K1andK2.pdf",
       # width = 7, height = 4)

inset1 <- ggplot(data = both.plots %>%
                   filter(Experiment == "K1")) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    # width = 5
  ),
  color = 'black') +
  geom_point(data = both.plots.subset.replicate.inj %>%
               filter(Experiment == "K1"),
             aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr),
             color = 'black', fill="white", shape = 21) +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  geom_vline(data = data.frame(Experiment = factor(c("K1"),levels = c("K2","K1")),
                               k.5 = c(79.46) ),
             aes(xintercept  = k.5),
             linetype = "dashed") +
  # scale_x_log10()+
  theme_bw()+ 
  facet_wrap(~Experiment) +
  geom_ribbon(data = k1.fit, aes(x = x, ymin = y.low, ymax= y.high),
              fill = "blue", alpha = 0.2)+
  # geom_ribbon(data = k2.fit, aes(x = x, ymin = y.low, ymax= y.high),
  #             fill = "red", alpha = 0.2)+
  geom_line(data = k1.fit, aes(x = x, y = y),
            color = "blue") +
  # geom_line(data = k2.fit, aes(x = x, y = y),
  #           color = "red") +
  # xlim(0,210) +
  scale_x_continuous(breaks = c(0,100,200),
                     limits = c(0,210))+
  # ylim(0,0.39) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4),
                     limits = c(0, 0.39)) +
  # scale_x_log10(limits = c(0.00001,2000))+
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)")+
  theme(aspect.ratio = 0.5,
        axis.title = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())
inset1
  
inset2 <- ggplot(data = both.plots %>%
                   filter(Experiment == "K2")) +
  geom_point(aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr), color = 'black') +
  geom_errorbar(aes(x = Treatment, ymin = correct.Prediction.nmoles.per.L.per.hr-SE.correct.rate,
                    ymax = correct.Prediction.nmoles.per.L.per.hr+SE.correct.rate,
                    # width = 5
  ),
  color = 'black') +
  geom_point(data = both.plots.subset.replicate.inj %>%
               filter(Experiment == "K2"),
             aes(x = Treatment, y = correct.Prediction.nmoles.per.L.per.hr),
             color = 'black', fill="white", shape = 21) +
  # geom_point(aes(x = Treatment.nmoles.added, y = predict), color = 'red', size = 3) +
  # scale_x_log10()+
  geom_vline(data = data.frame(Experiment = factor(c("K2"),levels = c("K2","K1")),
                               k.5 = c(11.06) ),
             aes(xintercept  = k.5),
             linetype = "dashed") +
  theme_bw()+ 
  facet_wrap(~Experiment) +
  # geom_ribbon(data = k1.fit, aes(x = x, ymin = y.low, ymax= y.high),
  #             fill = "blue", alpha = 0.2)+
  geom_ribbon(data = k2.fit, aes(x = x, ymin = y.low, ymax= y.high),
              fill = "red", alpha = 0.2)+
  # geom_line(data = k1.fit, aes(x = x, y = y),
            # color = "blue") +
  geom_line(data = k2.fit, aes(x = x, y = y),
            color = "red") +
  # xlim(0,210) + 
  scale_x_continuous(breaks = c(0,100,200),
                                   limits = c(0,210))+
  # ylim(0,0.625) +
  scale_y_continuous(breaks = c(0.0, 0.3, 0.6),
                     limits = c(0, 0.625)) +
  # scale_x_log10(limits = c(0.9,100))+
  ylab("nM GBT/L/hr") +
  xlab("Treatment (nM GBT)") +
  theme(aspect.ratio = 0.5,
        axis.title = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())
inset2



uptake.w.inset <- ggdraw(total.uptake.plot) +
  draw_plot(inset1 ,x =  .65, y = .51, .2, .2)+
  draw_plot(inset2, x = .15, y = .35, .2, .2)
uptake.w.inset
ggsave(uptake.w.inset, filename = "../output/Figures/Uptake-Kinetics-with-rangeOfFit-K1andK2_insets.pdf",
       width = 8, height = 4)
# ggsave(uptake.w.inset, filename = "../output/Figures/Uptake-Kinetics-with-rangeOfFit-K1andK2_insets.png",
       # width = 8, height = 3.6)


uptake.w.inset.new <- ggdraw(total.uptake.plot) +
  draw_plot(inset1 ,x =  .65, y = .41, .25, .25)+
  draw_plot(inset2, x = .17, y = .35, .25, .25)
# ggsave(uptake.w.inset.new, filename = "../output/Figures/Uptake-Kinetics-with-rangeOfFit-K1andK2.png",
       # width = 8, height = 3.6)


## Make a table of the parameters -----------------
K1.tab <- K1.curve %>%
  filter(value.name=="Estimate") %>%
  dplyr::select(-SE.params, -n, -value.name) %>%
  full_join(., K1.curve %>% 
              select(n) %>% 
              unique() %>% 
              mutate(coeff = "Iterations") %>%
              rename(Mean.Params = n)) %>%
  dplyr::rename(Estimate = Mean.Params) %>%
  mutate(Experiment = "K1")

# turnover.time.k1 = K1.tab$Estimate[K1.tab$coeff=='c']/(K1.tab$Estimate[K1.tab$coeff=='b'] * (K1.tab$Estimate[K1.tab$coeff=='c']/(
  # K1.tab$Estimate[K1.tab$coeff=='c'] + K1.tab$Estimate[K1.tab$coeff=='a'])))

K2.tab <- K2.curve %>%
  filter(value.name=="Estimate") %>%
  dplyr::select(-SE.params, -n, -value.name) %>%
  full_join(., K1.curve %>% 
              select(n) %>% 
              unique() %>% 
              mutate(coeff = "Iterations") %>%
              rename(Mean.Params = n)) %>%
  dplyr::rename(Estimate = Mean.Params) %>%
  mutate(Experiment = "K2")

# turnover.time.k2 = K2.tab$Estimate[K2.tab$coeff=='c']/(K2.tab$Estimate[K2.tab$coeff=='b'] * (K2.tab$Estimate[K2.tab$coeff=='c']/(
  # K2.tab$Estimate[K2.tab$coeff=='c'] + K2.tab$Estimate[K2.tab$coeff=='a'])))
WH.iterations <- full_join(K1.WH %>%
                             mutate(Experiment = "K1"),
                           K2.WH %>%
                             mutate(Experiment = "K2")) %>%
  group_by(Experiment) %>%
  summarise(Estimate = mean(n)) %>%
  mutate(coeff = "Iterations")
WH.tab <- full_join(K1.WH %>%
                      mutate(Experiment = "K1"),
                    K2.WH %>%
                      mutate(Experiment = "K2")) %>%
  dplyr::rename(Estimate = Mean.Params) %>%
  dplyr::select(-n, -value.name) %>%
  full_join(WH.iterations) %>%
  mutate(Method = "Linear Transformation") %>%
  filter(coeff != "model.summary",
         coeff != "slope")

all.tab = full_join(K1.tab, K2.tab) %>%
  dplyr::select(Experiment, coeff, Estimate, sd) %>%
  mutate(coeff = ifelse(coeff == "a", "Kt + S (nM)",
                        ifelse(coeff == "b", "Vmax (nM/hr)",
                               ifelse(coeff == "c", "d[GBT]",
                                      coeff))),
         Method = "Nonlinear least squares") %>%
  full_join(WH.tab) %>%
  dplyr::select(Experiment, coeff, Method, Estimate, sd)

Aloha.est.from.josh.july2020 <- 37.18 ## nM d[GBT]

all.tab.w.aloha.est <- all.tab %>%
  full_join(.,data.frame(Experiment = "ALOHA",coeff = "d[GBT]",
                         Estimate = Aloha.est.from.josh.july2020))

all.tab.w.aloha.est.wide <- all.tab.w.aloha.est %>%
  gather(value, Number, -Experiment, -coeff, -Method) %>%
  spread(coeff, Number) %>%
  mutate(Experiment = factor(Experiment, levels =  c("K1","K2","ALOHA"))) %>%
  arrange(Experiment) %>%
  dplyr::select(`Experiment`, `value`, `Method`, `d[GBT]`,
                `Kt + S (nM)`,`Vmax (nM/hr)`,
                `Turnover Time (hr)`, `Iterations`) %>%
  mutate(`d[GBT]` = round(`d[GBT]`,2),
         `Kt + S (nM)` = round(`Kt + S (nM)`,2),
         `Vmax (nM/hr)` = round(`Vmax (nM/hr)`,2),
         `Turnover Time (hr)` = round(`Turnover Time (hr)`,1))

## lat long and temp
CTD.temp <- read_csv("../data/CTD Temperature (2).csv")
CTD.temp.k.exp <- CTD.temp %>%
  mutate(date = str_sub(time, 1, 10)) %>%
  filter(depth > 13,
         depth < 17,
         date == "2019-04-16" | date == "2019-04-20") %>%
  group_by(date) %>%
  summarise(Lat = mean(lat),
            Lat.sd = sd(lat),
            Long = mean(lon),
            Long.sd = sd(lon),
            depth = mean(depth),
            temp.mean = mean(CTD_Temperature),
            temp.sd = sd(CTD_Temperature),
            n =n())

## lat long and temp
CTD.chl <- read_csv("../data/CTD Chloropigment.csv")
CTD.chl.k.exp <- CTD.chl %>%
  mutate(date = str_sub(time, 1, 10)) %>%
  filter(depth > 13,
         depth < 17,
         date == "2019-04-16" | date == "2019-04-20") %>%
  group_by(date) %>%
  summarise(Lat = mean(lat),
            Lat.sd = sd(lat),
            Long = mean(lon),
            Long.sd = sd(lon),
            depth = mean(depth),
            chl.mean = mean(CTD_Chloropigment),
            chl.sd = sd(CTD_Chloropigment),
            n =n())

all.ctd <- full_join(CTD.chl.k.exp, CTD.temp.k.exp) %>%
  mutate(Experiment = ifelse(date == "2019-04-16","K1","K2")) %>%
  dplyr::select(Experiment, date, Lat, Long, chl.mean, chl.sd, temp.mean, temp.sd) %>%
  gather(coeff, Val, -Experiment) %>%
  mutate(value = ifelse(grepl("sd",coeff),"sd","Estimate"),
         coeff.name = ifelse(grepl("temp",coeff),"Temperature",
                             ifelse(grepl("chl",coeff),"Chl",coeff))) %>%
  dplyr::select(-coeff) %>%
  spread(coeff.name, Val)%>%
  dplyr::select(Experiment, value, date, Lat, Long, Temperature, Chl) %>%
  mutate(Lat = round(as.numeric(Lat), 2),
         Long = round(as.numeric(Long), 2),
         Temperature = round(as.numeric(Temperature), 2),
         Chl = round(as.numeric(Chl),2))

tab.w.lat.long.ALOHA <- full_join(all.ctd, all.tab.w.aloha.est.wide) %>%
  mutate(Experiment = factor(Experiment, levels =  c("K1","K2","ALOHA")),
         Lat = ifelse(Experiment=="ALOHA",22.75,Lat),
         Long = ifelse(Experiment=="ALOHA",-158,Long)) %>%
  arrange(Experiment, Method, value) %>%
  dplyr::select(Experiment, date, Lat, Long, Temperature, Chl, Method, value,
                `d[GBT]`,`Kt + S (nM)`, `Vmax (nM/hr)`, `Turnover Time (hr)`, Iterations)

tab.w.lat.long <- full_join(all.ctd, all.tab.w.aloha.est.wide) %>%
  mutate(Experiment = factor(Experiment, levels =  c("K1","K2","ALOHA")),
         Lat = ifelse(Experiment=="ALOHA",22.75,Lat),
         Long = ifelse(Experiment=="ALOHA",-158,Long)) %>%
  arrange(Experiment, Method, value) %>%
  dplyr::select(Experiment, date, Lat, Long, Temperature, Chl, Method, value,
                `Kt + S (nM)`, `Vmax (nM/hr)`, `Turnover Time (hr)`, Iterations)  %>%
  filter(Experiment != "ALOHA")


## get PC PN data ----
PCPN <- source("../data/PCPN.Dat.R")[[1]] %>%
  mutate(C.N.Ratio = TC.umol.L/N.umol.L)
PCPN.for.plot <- PCPN %>%
  filter(Depth.m == 15,
         Station == 4 |Station == 5) %>%
  mutate(Lat = ifelse(Station == 4, 41.68, 
                      ifelse(Station == 5, 37.00, Lat))) %>% 
  dplyr::select(Lat, Depth.m, TC.umol.L, N.umol.L, C.N.Ratio, TON, DIN)

## get fcm data -----
all.fcm = read_csv("../../FateTidy/RawData/all.fcm.data.from.experiments.csv")
T0.fcm.aves = all.fcm %>%
  dplyr::filter(population != "beads" & population != "bacteria" & population != "unknown" &
                  flag == 0 & count > 30,
                incubation_time == 0)  %>%
  dplyr::group_by(incubation_time, fate, population) %>%
  dplyr::summarize(Estimate=(mean(biomass)),
                   sd = sd(biomass)) %>%
  rbind(., data.frame(incubation_time=0, fate="F1", population="prochloro",Estimate=NA, sd= NA)) %>%
  full_join(., all.fcm %>%
              dplyr::filter(population != "beads" & population == "bacteria" & flag == 0 & count > 30,
                            incubation_time == 0)  %>%
              dplyr::group_by(incubation_time, fate, population) %>%
              dplyr::summarize(Estimate=(mean(biomass)),
                               sd = sd(biomass))) %>%
  mutate(Experiment = ifelse(fate == "F1","K1","K2"),
         population = ifelse(population == "picoeuk","Picoeukaryotes biomass",
                             ifelse(population=="prochloro","Prochlorococcus biomass",
                                    ifelse(population=="synecho","Synechococcus biomass",
                                           "bacteria biomass")))) %>%
  group_by(Experiment, population) %>%
  dplyr::select(-fate, -incubation_time) %>%
  pivot_longer(cols = Estimate:sd) %>%
  mutate(value = round(value, 2))%>%
  pivot_wider(names_from = population ) %>%
  dplyr::rename(value = name)
## SAVE OUT THE TABLE -------
# kable(tab.w.lat.long) %>%
#   kable_styling(bootstrap_options = c("striped", "bordered", full_width = F)) %>%
#   save_kable(file = "../output/table1.html", self_contained = T)

tab.test <- xtable(tab.w.lat.long)
caption(tab.test) <- "Summary of samples collected and analyzed in this study."
label(tab.test)<- "siteSummary"
print(tab.test,
      type = "latex",
      file = "../output/xtable1.tex")

## save a table that shows the overall conditions of the experiments
tab.w.lat.long %>% 
  dplyr::select(Experiment, date, Lat, Long, value, Temperature, Chl) %>%
  unique() %>%
  full_join(T0.fcm.aves) %>%
  gather(Measurement, Number, -Experiment, -date, -Lat, -Long, -value) %>%
  mutate(value = ifelse(value =="Estimate","mean",value),
         Measurement = paste(Measurement, value)) %>%
  dplyr::select(-value, -date, -Lat, -Long) %>%
  spread(Measurement, Number) %>%
  full_join(tab.w.lat.long %>%
              dplyr::select(Experiment, date, Lat, Long) %>%
              filter(!is.na(Lat)) %>%
              unique()) %>%
  full_join(PCPN.for.plot) %>%
  mutate(C.N.Ratio = round(C.N.Ratio, 2)) %>%
  dplyr::select(Experiment, date, Lat, Long, Depth.m, TC.umol.L, N.umol.L,
                C.N.Ratio, TON, DIN, `Temperature mean`, `Temperature sd`,
                `Chl mean`, `Chl sd`,
                `bacteria biomass mean`, `bacteria biomass sd`,
                `Picoeukaryotes biomass mean`,`Picoeukaryotes biomass sd`,
                `Prochlorococcus biomass mean`,`Prochlorococcus biomass sd`,
                `Synechococcus biomass mean`,`Synechococcus biomass sd`) %>%
  gt(rowname_col =   "Experiment") %>%
  tab_spanner(label = "Temperature", columns = matches("Temperature"))  %>%
  tab_spanner(label = "Chl", columns = matches("Chl"))  %>%
  tab_spanner(label = "bacteria biomass", columns = matches("bacteria biomass")) %>%
  tab_spanner(label = "Picoeukaryotes biomass", columns = matches("Picoeukaryotes biomass")) %>%
  tab_spanner(label = "Prochlorococcus biomass", columns = matches("Prochlorococcus biomass")) %>%
  tab_spanner(label = "Synechococcus biomass", columns = matches("Synechococcus biomass")) %>%
  cols_label( `Temperature mean` = "mean",
              `Temperature sd` = "sd",
              `Chl mean` = "mean",
              `Chl sd` = "sd",
              `bacteria biomass mean` = "mean",
              `bacteria biomass sd` = "sd",
              `Picoeukaryotes biomass mean` = "mean",
              `Picoeukaryotes biomass sd` = "sd",
              `Prochlorococcus biomass mean` = "mean",
              `Prochlorococcus biomass sd` = "sd",
              `Synechococcus biomass mean` = "mean",
              `Synechococcus biomass sd` = "sd",
              Depth.m = "Depth",
              TC.umol.L = "PC",
              N.umol.L = "PN",
              C.N.Ratio = "C:N",
              date = "Date") %>%
  gtsave(filename = "Experiment_Conditions_new.tex",
         path = "../output/")
  
  
tab.w.lat.long %>% 
  dplyr::select(Experiment, Method, Iterations, value,  `Kt + S (nM)`, `Vmax (nM/hr)`, `Turnover Time (hr)` ) %>%
  unique()   %>%
  mutate( `Kt + S (nM)` = ifelse(value == "sd", paste0("(", `Kt + S (nM)`,")"),  `Kt + S (nM)`),
          `Vmax (nM/hr)` = ifelse(value =="sd", paste0("(", `Vmax (nM/hr)`,")"), `Vmax (nM/hr)`),
          `Turnover Time (hr)` = ifelse(value =="sd" & 
                                          !is.na(`Turnover Time (hr)`),
                                        paste0("(", `Turnover Time (hr)`,")"), `Turnover Time (hr)`) ) %>%
  gt(rowname_col =   "Method", groupname_col = "Experiment") %>%
  cols_hide("value") %>%
  # tab_spanner(label = "Experimental conditions", columns = matches("date|Lat|Lon|Temp|Chl"))  %>%
  tab_spanner(label = "Model estimates", columns = matches("(nM)|Kt + S (nM)|Vmax|Turnover Time"))  %>%
  fmt_missing(columns = matches("date|Lat|Lon|Temp|Chl|Iteratio"),  missing_text = "")%>%
  fmt_missing(columns = matches("(nM)|Kt + S (nM)|Vmax|Turnover Time"),  missing_text = "_") %>%
  cols_align(
    align = "right",
    columns = T) %>%
  # tab_footnote(footnote = "mean (sd)", 
  #              locations = cells_column_labels(
  #   columns = matches("(nM)|Kt + S (nM)|Vmax|Turnover Time"))) %>%
  gtsave(filename = "Uptake_Kinetic_Parameters.tex",
         path = "../output/")

## make a table comparing the parameters of when I calculate dgbt separately or not from kt -------
tibble(Experiment = c(1,1,2,2), params.fit = c(2,3,2,3),
           "Kt or Kt+S" = c(paste0(round(a.k1.mc,2), " (",round(a.k1.mc.sd,2),")"),
                  paste0(round(gbt.a.k1.mc,2), " (",round(gbt.a.k1.mc.sd,2),")"),
                  paste0(round(a.K2.mc,2), " (",round(a.K2.mc.sd,2),")"),
                  paste0(round(gbt.a.K2.mc,2), " (",round(gbt.a.K2.mc.sd,2),")")),
       Vmax = c(paste0(round(b.k1.mc,2), " (",round(b.k1.mc.sd,2),")"),
                paste0(round(gbt.b.k1.mc,2), " (",round(gbt.b.k1.mc.sd,2),")"),
                paste0(round(b.K2.mc,2), " (",round(b.K2.mc.sd,2),")"),
                paste0(round(gbt.b.K2.mc,2), " (",round(gbt.b.K2.mc.sd,2),")")),
       S = c("-",
             paste0(round(gbt.c.k1.mc,2), " (",round(gbt.c.k1.mc.sd,2),")"),
             "-",
             paste0(round(gbt.c.K2.mc,2), " (",round(gbt.c.K2.mc.sd,2),")")
       )
       )%>%
  gt() %>%
  gtsave(filename = "Uptake_Kinetic_2_vs_3_Parameters.tex",
         path = "../output/")

           