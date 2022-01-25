## Do a monte carlo run adding uncertainty to GBT uptake rates to get better error estimates on parameters

##load packages-----------------------
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(stats)

## load data -------------
uptake.rates <- read_csv("../output/GBT_Uptake_Rates_w_StdError.csv")

monte.carlo.uptake.rates <- function(uptake.rates = uptake.rates, 
                                     experiment = "K1",
                                     a_start = 100 , #param a is the Ks,
                                     b_start = 1 , #b is Vmax,
                                     c_start = 1, # c is in situ d[GBT] in nM
                                     replace.0 = TRUE,
                                     reps = 1000){
  exp.data <- uptake.rates %>%
    filter(Experiment==experiment)
  norm.test <-rnorm(reps, mean = exp.data$correct.Prediction.nmoles.per.L.per.hr[7], sd = exp.data$SE.correct.rate[7]*7)
  
  exp.data.test <- exp.data %>%
    group_by(correct.Prediction.nmoles.per.L.per.hr, SE.correct.rate) %>%
    mutate(newVals = list(rnorm(reps, mean = correct.Prediction.nmoles.per.L.per.hr,
                                sd = SE.correct.rate * 7))) %>%
    unnest(newVals) %>%
    mutate(version = paste("model",rep(1:reps),sep="_")) %>%
    spread(key = "version", value = "newVals")
  
  all.vals <- data.frame(value.name = character())
  for(i in 1:reps){
    x <- exp.data.test$Treatment
    y <-  dplyr::pull(exp.data.test[,i+9])
    if(replace.0){
      y[y<0]<-0
    }
    #model
    m.k.1 <-nls(y~(b*(x+c))/(a+x+c),start=list(a=a_start,b=b_start, c= c_start), control = 
                  nls.control(maxiter = 500, warnOnly = T)
    )
    #get some estimation of goodness of fit
    y.guess <- (b_start/(x+c_start))/(a_start+x+c_start)
    test <- summary(m.k.1)
    column.name <- paste("model",i,sep="_")
    
    
    total.model.vals <- data.frame(coeff = "model.summary",
                                   value.name = c("sigma",
                                                  "resid.sum.sqrs",
                                                  "total.cor",
                                                  "Iterations"),
                                   value = c(test$sigma, 
                                             m.k.1$m$deviance(), 
                                             cor(y,predict(m.k.1)),
                                             m.k.1$convInfo$finIter))
    coeff.vals <- as.data.frame(test$coefficients) %>%
      mutate(coeff = rownames(test$coefficients)) %>% 
      gather(key = "value.name", value = "value", -coeff)
    these.vals <-  suppressMessages(full_join(total.model.vals, coeff.vals)) %>%
      mutate(!!column.name := value) %>%
      dplyr::select(-value)
    all.vals <- suppressMessages(full_join(all.vals, these.vals))
  }
  return(all.vals)
}



monte.carlo.uptake.rates.no.dGBT <- function(uptake.rates = uptake.rates, 
                                     experiment = "K1",
                                     a_start = 100 , #param a is the Ks,
                                     b_start = 1 , #b is Vmax,
                                     # c_start = 1, # c is in situ d[GBT] in nM
                                     replace.0 = TRUE,
                                     reps = 1000){
  exp.data <- uptake.rates %>%
    filter(Experiment==experiment)
  norm.test <-rnorm(reps, mean = exp.data$correct.Prediction.nmoles.per.L.per.hr[7], sd = exp.data$SE.correct.rate[7]*7)
  
  exp.data.test <- exp.data %>%
    group_by(correct.Prediction.nmoles.per.L.per.hr, SE.correct.rate) %>%
    mutate(newVals = list(rnorm(reps, mean = correct.Prediction.nmoles.per.L.per.hr,
                                sd = SE.correct.rate * 7))) %>%
    unnest(newVals) %>%
    mutate(version = paste("model",rep(1:reps),sep="_")) %>%
    spread(key = "version", value = "newVals")
  
  all.vals <- data.frame(value.name = character())
  for(i in 1:reps){
    x <- exp.data.test$Treatment
    y <-  dplyr::pull(exp.data.test[,i+9])
    if(replace.0){
      y[y<0]<-0
    }
    #model
    m.k.1 <-nls(y~(b*(x))/(a+x),start=list(a=a_start,b=b_start), control = 
                  nls.control(maxiter = 500, warnOnly = T)
    )
    #get some estimation of goodness of fit
    y.guess <- (b_start/(x))/(a_start+x)
    test <- summary(m.k.1)
    column.name <- paste("model",i,sep="_")
    
    
    total.model.vals <- data.frame(coeff = "model.summary",
                                   value.name = c("sigma",
                                                  "resid.sum.sqrs",
                                                  "total.cor",
                                                  "Iterations"),
                                   value = c(test$sigma, 
                                             m.k.1$m$deviance(), 
                                             cor(y,predict(m.k.1)),
                                             m.k.1$convInfo$finIter))
    coeff.vals <- as.data.frame(test$coefficients) %>%
      mutate(coeff = rownames(test$coefficients)) %>% 
      gather(key = "value.name", value = "value", -coeff)
    these.vals <-  suppressMessages(full_join(total.model.vals, coeff.vals)) %>%
      mutate(!!column.name := value) %>%
      dplyr::select(-value)
    all.vals <- suppressMessages(full_join(all.vals, these.vals))
  }
  return(all.vals)
}
## K1 ----

set.seed(124)
K1.monte.carlo <- monte.carlo.uptake.rates(uptake.rates = uptake.rates, experiment = "K1", replace.0 = T,
                                           reps = 1000,
                                           a_start = 109.7,
                                           b_start = 0.3735,
                                           c_start = 8.125)
K1.monte.carlo.long <- K1.monte.carlo %>%
  gather(key = "Model","Value",-value.name, -coeff)


set.seed(124)
K1.monte.carlo.no.dGBT <- monte.carlo.uptake.rates.no.dGBT(uptake.rates = uptake.rates, experiment = "K1", replace.0 = T,
                                           reps = 1000,
                                           a_start = 71.30021002,
                                           b_start = 0.351497752)
K1.monte.carlo.no.dGBT.long <- K1.monte.carlo.no.dGBT %>%
  gather(key = "Model","Value",-value.name, -coeff)

## write out
write_csv(K1.monte.carlo.long, path = paste("../output/GBT Uptake K1 all Monte Carlo Model Resuls.csv",sep = "_"))
write_csv(K1.monte.carlo.no.dGBT.long, path = paste("../output/GBT Uptake K1 all Monte Carlo Model Resuls no dGBT.csv",sep = "_"))

## evaluate
K1.monte.carlo.summary <- K1.monte.carlo.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

K1.monte.carlo.summary.no.dGBT <- K1.monte.carlo.no.dGBT.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)


boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='a' & K1.monte.carlo.long$value.name=="Estimate"])
title('A: Ks K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='a' & K1.monte.carlo.long$value.name=="Estimate"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='b' & K1.monte.carlo.long$value.name=="Estimate"])
title('B: Vmax K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='b' & K1.monte.carlo.long$value.name=="Estimate"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='c' & K1.monte.carlo.long$value.name=="Estimate"])
title('C: d[GBT] K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='c' & K1.monte.carlo.long$value.name=="Estimate"])


boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='total.cor'])
title('R K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='total.cor'])


boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])
title('sigma K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])


## outliers based on total fit (R)
R.outliers.k1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='total.cor'])$out)]
## outliers based on total fit (sigma)
sigma.outliers.k1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])$out)]
## outliers based on total fit (Iterations)
iteration.outliers.k1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value > 50 & K1.monte.carlo.long$value.name=='Iterations') ]

## read in k1 data ------
K1.monte.carlo.long <- read_csv("../output/GBT Uptake K1 all Monte Carlo Model Resuls.csv")
## evaluate
K1.monte.carlo.summary <- K1.monte.carlo.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

pdf('../output/monteCarlo_uptake_models_K1.pdf')
par(mfrow = c(1, 5))
boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='a' & K1.monte.carlo.long$value.name=="Estimate"])
title('A: Ks K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='a' & K1.monte.carlo.long$value.name=="Estimate"])
A.outliers.K1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% 
                                                   boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='a' &
                                                                                             K1.monte.carlo.long$value.name=="Estimate"])
                                                 $out)]

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='b' & K1.monte.carlo.long$value.name=="Estimate"])
title('B: Vmax K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='b' & K1.monte.carlo.long$value.name=="Estimate"])
B.outliers.K1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% 
                                                   boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='b' &
                                                                                             K1.monte.carlo.long$value.name=="Estimate"])
                                                 $out)]

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='c' & K1.monte.carlo.long$value.name=="Estimate"])
title('C: d[GBT] K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='c' & K1.monte.carlo.long$value.name=="Estimate"])


boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='total.cor'])
title('R K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='total.cor'])


boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])
title('sigma K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])

dev.off()

## outliers based on total fit (R)
R.outliers.K1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='total.cor'])$out)]
## outliers based on total fit (sigma)
sigma.outliers.K1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])$out)]
## outliers based on total fit (Iterations)
iteration.outliers.K1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value > 50 & K1.monte.carlo.long$value.name=='Iterations') ]

## ditch outliers that were found with A, B, R, sigma, iterations, make summary with errors, write it out ----
all.outliers <- unique(c(A.outliers.K1, B.outliers.K1, R.outliers.K1,
                         sigma.outliers.K1, iteration.outliers.K1))
K1.monte.carlo.long.no.outliers <- K1.monte.carlo.long %>%
  filter(!(Model %in% all.outliers))
## propegate error for the parameters
iterations.only.K1 <-  K1.monte.carlo.long.no.outliers %>%
  filter(value.name=="Iterations") %>%
  rename(Iterations = Value) %>%
  dplyr::select(Model, Iterations)

K1.summary.with.error <- K1.monte.carlo.long.no.outliers %>%
  full_join(iterations.only.K1, by = "Model") %>%
  group_by(coeff, value.name) %>%
  summarise(Mean.Params = mean(Value),
            sd = sd(Value),
            n = n(),
            SE.params = sqrt(sum((Value*Iterations)^2))/n)

K1.summary.tidy <- K1.summary.with.error %>%
  filter(value.name!="Pr(>|t|)",
         value.name!="t value") %>%
  mutate(SE.params = ifelse(coeff=="model.summary",NA,
                            ifelse(value.name=="Estimate",NA,SE.params)))

## write out that tidy bit
write_csv(K1.summary.tidy, path = paste("../output/K1-MonteCarlo_OutliersRemoved_summary_withError.csv",sep="_"))



## read in k1 data no dGBT ------
K1.monte.carlo.no.dGBT.long <- read_csv("../output/GBT Uptake K1 all Monte Carlo Model Resuls no dGBT.csv")
## evaluate
K1.monte.carlo.summary.no.dGBT <- K1.monte.carlo.no.dGBT.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

pdf('../output/monteCarlo_uptake_models_K1_no_dGBT.pdf')
par(mfrow = c(1, 5))
boxplot(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='a' & K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
title('A: Ks K1')
boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='a' & K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
A.outliers.K1 <- K1.monte.carlo.no.dGBT.long$Model[which(K1.monte.carlo.no.dGBT.long$Value %in% 
                                                   boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='a' &
                                                                                             K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
                                                 $out)]

boxplot(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='b' & K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
title('B: Vmax K1')
boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='b' & K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
B.outliers.K1 <- K1.monte.carlo.no.dGBT.long$Model[which(K1.monte.carlo.no.dGBT.long$Value %in% 
                                                   boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='b' &
                                                                                             K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
                                                 $out)]

boxplot(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='c' & K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])
title('C: d[GBT] K1')
boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$coeff=='c' & K1.monte.carlo.no.dGBT.long$value.name=="Estimate"])


boxplot(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$value.name=='total.cor'])
title('R K1')
boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$value.name=='total.cor'])


boxplot(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$value.name=='sigma'])
title('sigma K1')
boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$value.name=='sigma'])

dev.off()

## outliers based on total fit (R)
R.outliers.K1 <- K1.monte.carlo.no.dGBT.long$Model[which(K1.monte.carlo.no.dGBT.long$Value %in% boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$value.name=='total.cor'])$out)]
## outliers based on total fit (sigma)
sigma.outliers.K1 <- K1.monte.carlo.no.dGBT.long$Model[which(K1.monte.carlo.no.dGBT.long$Value %in% boxplot.stats(K1.monte.carlo.no.dGBT.long$Value[K1.monte.carlo.no.dGBT.long$value.name=='sigma'])$out)]
## outliers based on total fit (Iterations)
iteration.outliers.K1 <- K1.monte.carlo.no.dGBT.long$Model[which(K1.monte.carlo.no.dGBT.long$Value > 50 & K1.monte.carlo.no.dGBT.long$value.name=='Iterations') ]

## ditch outliers that were found with A, B, R, sigma, iterations, make summary with errors, write it out ----
all.outliers <- unique(c(A.outliers.K1, B.outliers.K1, R.outliers.K1,
                         sigma.outliers.K1, iteration.outliers.K1))
K1.monte.carlo.long.no.outliers.no.dGBT <- K1.monte.carlo.no.dGBT.long %>%
  filter(!(Model %in% all.outliers))
## propegate error for the parameters
iterations.only.K1 <-  K1.monte.carlo.long.no.outliers.no.dGBT %>%
  filter(value.name=="Iterations") %>%
  rename(Iterations = Value) %>%
  dplyr::select(Model, Iterations)

K1.summary.with.error.no.dGBT <- K1.monte.carlo.long.no.outliers.no.dGBT %>%
  full_join(iterations.only.K1, by = "Model") %>%
  group_by(coeff, value.name) %>%
  summarise(Mean.Params = mean(Value),
            sd = sd(Value),
            n = n(),
            SE.params = sqrt(sum((Value*Iterations)^2))/n)

K1.summary.tidy.no.dGBT <- K1.summary.with.error.no.dGBT %>%
  filter(value.name!="Pr(>|t|)",
         value.name!="t value") %>%
  mutate(SE.params = ifelse(coeff=="model.summary",NA,
                            ifelse(value.name=="Estimate",NA,SE.params)))

## write out that tidy bit
write_csv(K1.summary.tidy.no.dGBT, path = paste("../output/K1-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv",sep="_"))






## K2 ----

set.seed(483)
K2.monte.carlo <- monte.carlo.uptake.rates(uptake.rates = uptake.rates, experiment = "K2", replace.0 = T,
                                           reps = 1000,
                                           a_start = 11.86,
                                           b_start = 0.5659,
                                           c_start = 0.5856)
K2.monte.carlo.long <- K2.monte.carlo %>%
  gather(key = "Model","Value",-value.name, -coeff)

set.seed(124)
K2.monte.carlo.no.dGBT <- monte.carlo.uptake.rates.no.dGBT(uptake.rates = uptake.rates, experiment = "K2", replace.0 = T,
                                                           reps = 1000,
                                                           a_start = 10.54998845,
                                                           b_start = 0.562719359)
K2.monte.carlo.no.dGBT.long <- K2.monte.carlo.no.dGBT %>%
  gather(key = "Model","Value",-value.name, -coeff)


## write out
write_csv(K2.monte.carlo.long, path = paste("../output/GBT Uptake K2 all Monte Carlo Model Resuls.csv",sep = "_"))
write_csv(K2.monte.carlo.no.dGBT.long, path = paste("../output/GBT Uptake K2 all Monte Carlo Model Resuls no dGBT.csv",sep = "_"))


## read in k2 data ------
K2.monte.carlo.long <- read_csv("../output/GBT Uptake K2 all Monte Carlo Model Resuls.csv")
## evaluate
K2.monte.carlo.summary <- K2.monte.carlo.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

pdf('../output/monteCarlo_uptake_models_K2.pdf')
par(mfrow = c(1, 5))

boxplot(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='a' & K2.monte.carlo.long$value.name=="Estimate"])
title('A: Ks K2')
boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='a' & K2.monte.carlo.long$value.name=="Estimate"])
A.outliers.K2 <- K2.monte.carlo.long$Model[which(K2.monte.carlo.long$Value %in% 
                                                   boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='a' &
                                                                                             K2.monte.carlo.long$value.name=="Estimate"])
$out)]

boxplot(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='b' & K2.monte.carlo.long$value.name=="Estimate"])
title('B: Vmax K2')
boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='b' & K2.monte.carlo.long$value.name=="Estimate"])
B.outliers.K2 <- K2.monte.carlo.long$Model[which(K2.monte.carlo.long$Value %in% 
                                                   boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='b' &
                                                                                             K2.monte.carlo.long$value.name=="Estimate"])
                                                 $out)]

boxplot(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='c' & K2.monte.carlo.long$value.name=="Estimate"])
title('C: d[GBT] K2')
boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$coeff=='c' & K2.monte.carlo.long$value.name=="Estimate"])


boxplot(K2.monte.carlo.long$Value[K2.monte.carlo.long$value.name=='total.cor'])
title('R K2')
boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$value.name=='total.cor'])


boxplot(K2.monte.carlo.long$Value[K2.monte.carlo.long$value.name=='sigma'])
title('sigma K2')
boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$value.name=='sigma'])

dev.off()

## outliers based on total fit (R)
R.outliers.K2 <- K2.monte.carlo.long$Model[which(K2.monte.carlo.long$Value %in% boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$value.name=='total.cor'])$out)]
## outliers based on total fit (sigma)
sigma.outliers.K2 <- K2.monte.carlo.long$Model[which(K2.monte.carlo.long$Value %in% boxplot.stats(K2.monte.carlo.long$Value[K2.monte.carlo.long$value.name=='sigma'])$out)]
## outliers based on total fit (Iterations)
iteration.outliers.K2 <- K2.monte.carlo.long$Model[which(K2.monte.carlo.long$Value > 50 & K2.monte.carlo.long$value.name=='Iterations') ]

## ditch outliers that were found with A, B, R, sigma, iterations, make summary with errors, write it out ----
all.outliers <- unique(c(A.outliers.K2, B.outliers.K2, R.outliers.K2,
                  sigma.outliers.K2, iteration.outliers.K2))
K2.monte.carlo.long.no.outliers <- K2.monte.carlo.long %>%
  filter(!(Model %in% all.outliers))
## propegate error for the parameters
iterations.only.k2 <-  K2.monte.carlo.long.no.outliers %>%
  filter(value.name=="Iterations") %>%
  rename(Iterations = Value) %>%
  dplyr::select(Model, Iterations)

K2.summary.with.error <- K2.monte.carlo.long.no.outliers %>%
  full_join(iterations.only.k2, by = "Model") %>%
  group_by(coeff, value.name) %>%
  summarise(Mean.Params = mean(Value),
            sd = sd(Value),
            n = n(),
            SE.params = sqrt(sum((Value*Iterations)^2))/n)

K2.summary.tidy <- K2.summary.with.error %>%
  filter(value.name!="Pr(>|t|)",
         value.name!="t value") %>%
  mutate(SE.params = ifelse(coeff=="model.summary",NA,
                            ifelse(value.name=="Estimate",NA,SE.params)))

## write out that tidy bit
write_csv(K2.summary.tidy, path = paste("../output/","K2-MonteCarlo_OutliersRemoved_summary_withError.csv",sep=""))

## read in k2 data no dGBT ------
K2.monte.carlo.long.no.dGBT <- read_csv("../output/GBT Uptake K2 all Monte Carlo Model Resuls no dGBT.csv")
## evaluate
K2.monte.carlo.summary.no.dGBT <- K2.monte.carlo.long.no.dGBT %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

pdf('../output/monteCarlo_uptake_models_K2_no_dGBT.pdf')
par(mfrow = c(1, 5))

boxplot(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='a' & K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
title('A: Ks K2')
boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='a' & K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
A.outliers.K2 <- K2.monte.carlo.long.no.dGBT$Model[which(K2.monte.carlo.long.no.dGBT$Value %in% 
                                                   boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='a' &
                                                                                             K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
                                                 $out)]

boxplot(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='b' & K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
title('B: Vmax K2')
boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='b' & K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
B.outliers.K2 <- K2.monte.carlo.long.no.dGBT$Model[which(K2.monte.carlo.long.no.dGBT$Value %in% 
                                                   boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='b' &
                                                                                             K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
                                                 $out)]

boxplot(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='c' & K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])
title('C: d[GBT] K2')
boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$coeff=='c' & K2.monte.carlo.long.no.dGBT$value.name=="Estimate"])


boxplot(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$value.name=='total.cor'])
title('R K2')
boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$value.name=='total.cor'])


boxplot(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$value.name=='sigma'])
title('sigma K2')
boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$value.name=='sigma'])

dev.off()

## outliers based on total fit (R)
R.outliers.K2 <- K2.monte.carlo.long.no.dGBT$Model[which(K2.monte.carlo.long.no.dGBT$Value %in% boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$value.name=='total.cor'])$out)]
## outliers based on total fit (sigma)
sigma.outliers.K2 <- K2.monte.carlo.long.no.dGBT$Model[which(K2.monte.carlo.long.no.dGBT$Value %in% boxplot.stats(K2.monte.carlo.long.no.dGBT$Value[K2.monte.carlo.long.no.dGBT$value.name=='sigma'])$out)]
## outliers based on total fit (Iterations)
iteration.outliers.K2 <- K2.monte.carlo.long.no.dGBT$Model[which(K2.monte.carlo.long.no.dGBT$Value > 50 & K2.monte.carlo.long.no.dGBT$value.name=='Iterations') ]

## ditch outliers that were found with A, B, R, sigma, iterations, make summary with errors, write it out ----
all.outliers <- unique(c(A.outliers.K2, B.outliers.K2, R.outliers.K2,
                         sigma.outliers.K2, iteration.outliers.K2))
K2.monte.carlo.long.no.dGBT.no.outliers <- K2.monte.carlo.long.no.dGBT %>%
  filter(!(Model %in% all.outliers))
## propegate error for the parameters
iterations.only.k2 <-  K2.monte.carlo.long.no.dGBT.no.outliers %>%
  filter(value.name=="Iterations") %>%
  rename(Iterations = Value) %>%
  dplyr::select(Model, Iterations)

K2.summary.with.error.nodGBT <- K2.monte.carlo.long.no.dGBT.no.outliers %>%
  full_join(iterations.only.k2, by = "Model") %>%
  group_by(coeff, value.name) %>%
  summarise(Mean.Params = mean(Value),
            sd = sd(Value),
            n = n(),
            SE.params = sqrt(sum((Value*Iterations)^2))/n)

K2.summary.tidy.no.dGBT <- K2.summary.with.error.nodGBT %>%
  filter(value.name!="Pr(>|t|)",
         value.name!="t value") %>%
  mutate(SE.params = ifelse(coeff=="model.summary",NA,
                            ifelse(value.name=="Estimate",NA,SE.params)))

## write out that tidy bit
write_csv(K2.summary.tidy.no.dGBT, path = paste("../output/","K2-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv",sep=""))

