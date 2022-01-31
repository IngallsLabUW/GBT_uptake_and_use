## turnover time calculation based on the Wright-Hobbie linear transformation
## see Kiene and Williams 1998:
## Plot t/f (incubation time in hour divided by ration of added that was taken up) vs A (added concentration in nM). 
## Fit least squares regression. Then Kt+S = -x intercept and Vmax = 1/slope and turnover time = y-intercept

##load packages-----------------------
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(stats)
library(kableExtra)
library(xtable)
library(gt)

## read in data ---------
gbt.concentrations <- read_csv("../output/GBT_Concentrations_and_Uptake_Rates_w_StdError.csv")

## read in curves --------
K1.curve <- read_csv("../output/K1-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv")
K2.curve <- read_csv("../output/K2-MonteCarlo_OutliersRemoved_summary_withError_no_dGBT.csv")

## monte carlo function
monte.carlo.linear.transformation <- function(exp.data = exp.data, 
                                              experiment = "K1",
                                              replace.0 = TRUE,
                                              reps = 1000){
  
  exp.data.test <- exp.data %>%
    filter(!is.na(SE.time.over.fraction)) %>%
    group_by(time.over.fraction, SE.time.over.fraction) %>%
    mutate(newVals = list(rnorm(reps, mean = time.over.fraction,
                                sd = SE.time.over.fraction * 7))) %>%
    unnest(newVals) %>%
    mutate(version = paste("model",rep(1:reps),sep="_")) %>%
    spread(key = "version", value = "newVals") %>%
    full_join(., exp.data %>% 
                filter(is.na(SE.time.over.fraction)))
  
  exp.data.test[is.na(exp.data.test)] <- 0
  
  all.vals <- data.frame(value.name = character())
  for(i in 1:reps){
    x <- exp.data.test$Treatment.nmoles.added
    y <-  dplyr::pull(exp.data.test[,i+15])
    if(replace.0){
      y[y<0]<-0
    }
    #model
    m.k.1 <- lm(y~x)
    
    #get some estimation of goodness of fit
    # y.guess <- (b_start/(x+c_start))/(a_start+x+c_start)
    test <- summary(m.k.1)
    column.name <- paste("model",i,sep="_")
    
    
    total.model.vals <- data.frame(coeff = "model.summary",
                                   value.name = c("sigma",
                                                  "adj.r.sqrd",
                                                  "total.cor"),
                                   value = c(test$sigma, 
                                             test$adj.r.squared, 
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


## do wright-hobbie transformation----------
## K1 ----------------
GBT.Transformed.K1 <- gbt.concentrations %>% 
  filter(!is.na(Treatment),
         Experiment=="K1") %>%
  mutate(Treatment.nmoles.added = Treatment,
         fraction.taken.up= ifelse(is.na(correct.Prediction.nmoles.per.L),
                                   0, 
                                   (correct.Prediction.nmoles.per.L)/Treatment.nmoles.added),
         SE.fraction.taken.up = (fraction.taken.up * (SE.correct.conc/correct.Prediction.nmoles.per.L) ),# propogate standard error based on the relative standard error 
         time.over.fraction = ifelse(fraction.taken.up==0, 0, `Time incubated`/60/fraction.taken.up),
         SE.time.over.fraction = (time.over.fraction * (SE.fraction.taken.up/fraction.taken.up))) 

plot(GBT.Transformed.K1$Treatment.nmoles.added, GBT.Transformed.K1$time.over.fraction)
m2.k1 <- lm( time.over.fraction~Treatment.nmoles.added, data = GBT.Transformed.K1)
lines(GBT.Transformed.K1$Treatment.nmoles.added, 
      predict(m2.k1),
      col="blue",lty=2,lwd=3)
title("K1 Wright-Hobbie linear transformation")
m2.k1
cor(GBT.Transformed.K1$time.over.fraction,predict(m2.k1))
(sum.k1.mod <- summary(m2.k1))

## Then Kt+S = -x intercept and Vmax = 1/slope and turnover time = y-intercept
turnover.time.hr.k1.initial = sum.k1.mod$coefficients[1,1]
Vmax.nM.k1.initial = 1/(sum.k1.mod$coefficients[2,1])
Kt.S.k1.initial = (sum.k1.mod$coefficients[1,1] )/(sum.k1.mod$coefficients[2,1])

##monte carlo 
set.seed(124)
K1.monte.carlo <- monte.carlo.linear.transformation(exp.data =  GBT.Transformed.K1,
                                                    experiment = "K1",
                                                    replace.0 = T,
                                                    reps = 1000)
K1.monte.carlo.long <- K1.monte.carlo %>%
  gather(key = "Model","Value",-value.name, -coeff)

## write out
write_csv(K1.monte.carlo.long, path = paste("../output/GBT K1 linear transformation all Monte Carlo Model Resuls.csv",sep = "_"))

## evaluate
K1.monte.carlo.summary <- K1.monte.carlo.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='(Intercept)' & K1.monte.carlo.long$value.name=="Estimate"])
title('Intercept: turnover time K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='(Intercept)' & K1.monte.carlo.long$value.name=="Estimate"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='(Intercept)' & K1.monte.carlo.long$value.name=="Pr(>|t|)"])
title('Intercept: turnover time P val K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='(Intercept)' & K1.monte.carlo.long$value.name=="Pr(>|t|)"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='x' & K1.monte.carlo.long$value.name=="Estimate"])
title('x: slope (1/Vmax) K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='x' & K1.monte.carlo.long$value.name=="Estimate"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='x' & K1.monte.carlo.long$value.name=="Pr(>|t|)"])
title('x: slope (1/Vmax) P val K1')
# boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='x' & K1.monte.carlo.long$value.name=="Pr(>|t|)"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='model.summary' & K1.monte.carlo.long$value.name=="adj.r.sqrd"])
title('adjusted R squared K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='model.summary' & K1.monte.carlo.long$value.name=="adj.r.sqrd"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='model.summary' & K1.monte.carlo.long$value.name=="total.cor"])
title('total corelation K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='model.summary' & K1.monte.carlo.long$value.name=="total.cor"])

boxplot(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])
title('sigma K1')
boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$value.name=='sigma'])


## filter outliers based on adj R sqr and total.cor
R.outliers.k1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='model.summary' & K1.monte.carlo.long$value.name=="adj.r.sqrd"])$out)]
Corr.outliers.k1 <- K1.monte.carlo.long$Model[which(K1.monte.carlo.long$Value %in% boxplot.stats(K1.monte.carlo.long$Value[K1.monte.carlo.long$coeff=='model.summary' & K1.monte.carlo.long$value.name=="total.cor"])$out)]

all.outliers <- unique(c(R.outliers.k1, Corr.outliers.k1))
K1.monte.carlo.long.no.outliers <- K1.monte.carlo.long %>%
  filter(!(Model %in% all.outliers))
## propegate error for the parameters
K1.summary.with.error <- K1.monte.carlo.long.no.outliers %>%
  group_by(coeff, value.name) %>%
  summarise(Mean.Params = mean(Value),
            sd = sd(Value),
            n = n())

K1.summary.tidy <- K1.summary.with.error %>%
  filter(value.name!="Pr(>|t|)",
         value.name!="t value",
         value.name !="Std. Error") %>%
## Then Kt+S = -x intercept and Vmax = 1/slope and turnover time = y-intercept
  mutate(coeff = ifelse(coeff == "(Intercept)", "Turnover Time (hr)",
                        ifelse(coeff == "x","slope",coeff)))


Vmax.nM.k1.monte = 1/(K1.summary.tidy$Mean.Params[K1.summary.tidy$coeff=="slope"])
Vmax.nM.k1.monte.sd = Vmax.nM.k1.monte * ((K1.summary.tidy$sd[K1.summary.tidy$coeff=="slope"])/(K1.summary.tidy$Mean.Params[K1.summary.tidy$coeff=="slope"]))
Kt.S.k1.monte = (K1.summary.tidy$Mean.Params[K1.summary.tidy$coeff=="Turnover Time (hr)"])/(K1.summary.tidy$Mean.Params[K1.summary.tidy$coeff=="slope"])
Kt.S.k1.monte.sd = Kt.S.k1.monte * sqrt((K1.summary.tidy$sd[K1.summary.tidy$coeff=="Turnover Time (hr)"]/K1.summary.tidy$Mean.Params[K1.summary.tidy$coeff=="Turnover Time (hr)"])^2 +
                                          (K1.summary.tidy$sd[K1.summary.tidy$coeff=="slope"]/K1.summary.tidy$Mean.Params[K1.summary.tidy$coeff=="slope"])^2
                                        )
K1.summary.w.ests <- data.frame(coeff = c("Vmax (nM/hr)","Kt + S (nM)"),
           value.name = "calculated",
           Mean.Params = c(Vmax.nM.k1.monte, Kt.S.k1.monte),
           sd = c(Vmax.nM.k1.monte.sd, Kt.S.k1.monte.sd),
           n = mean(K1.summary.tidy$n)) %>%
  full_join(K1.summary.tidy)

## write out that tidy bit
write_csv(K1.summary.w.ests, path = paste("../output/K1-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv",sep="_"))

## K2 ----------------
GBT.Transformed.K2 <- gbt.concentrations %>% 
  filter(!is.na(Treatment),
         Experiment=="K2") %>%
  mutate(Treatment.nmoles.added = Treatment,
         fraction.taken.up= ifelse(is.na(correct.Prediction.nmoles.per.L),
                                   0, 
                                   (correct.Prediction.nmoles.per.L)/Treatment.nmoles.added),
         SE.fraction.taken.up = (fraction.taken.up * (SE.correct.conc/correct.Prediction.nmoles.per.L) ),# propogate standard error based on the relative standard error 
         time.over.fraction = ifelse(fraction.taken.up==0, 0, `Time incubated`/60/fraction.taken.up),
         SE.time.over.fraction = (time.over.fraction * (SE.fraction.taken.up/fraction.taken.up))) 

# ggplot(GBT.Transformed.K2) +
#   geom_point(aes(x = Treatment.nmoles.added, y = time.over.fraction)) +
#   geom_errorbar(aes(x = Treatment.nmoles.added, ymin = time.over.fraction - SE.time.over.fraction,
#                     ymax = time.over.fraction + SE.time.over.fraction), 
#                 width = 50)

plot(GBT.Transformed.K2$Treatment.nmoles.added, GBT.Transformed.K2$time.over.fraction)
m2.k2 <- lm(time.over.fraction ~ Treatment.nmoles.added, data = GBT.Transformed.K2)
lines(GBT.Transformed.K2$Treatment.nmoles.added, 
      predict(m2.k2),
      col="blue",lty=2,lwd=3)
title("K2 Wright-Hobbie linear transformation")
m2.k2
cor(GBT.Transformed.K2$time.over.fraction,predict(m2.k2))
(sum.k2.mod <- summary(m2.k2))

## Then Kt+S = -x intercept and Vmax = 1/slope and turnover time = y-intercept
turnover.time.hr.k2.initial = sum.k2.mod$coefficients[1,1] 
Vmax.nM.k2.initial = 1/(sum.k2.mod$coefficients[2,1]) ## confidence is in this but not the others
Kt.S.k2.initial = (sum.k2.mod$coefficients[1,1] )/(sum.k2.mod$coefficients[2,1])

##monte carlo 
set.seed(124)
k2.monte.carlo <- monte.carlo.linear.transformation(exp.data =  GBT.Transformed.K2,
                                                    experiment = "K2",
                                                    replace.0 = T,
                                                    reps = 1000)
k2.monte.carlo.long <- k2.monte.carlo %>%
  gather(key = "Model","Value",-value.name, -coeff)

## write out
write_csv(k2.monte.carlo.long, path = paste("../output/GBT k2 linear transformation all Monte Carlo Model Resuls.csv",sep = "_"))

## evaluate
k2.monte.carlo.summary <- k2.monte.carlo.long %>%
  group_by(value.name, coeff) %>%
  summarise(mean = mean(Value),
            sd = sd(Value),
            n = n(),
            rsd = sd/mean)

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='(Intercept)' & k2.monte.carlo.long$value.name=="Estimate"])
title('Intercept: turnover time k2')
boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='(Intercept)' & k2.monte.carlo.long$value.name=="Estimate"])

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='(Intercept)' & k2.monte.carlo.long$value.name=="Pr(>|t|)"])
title('Intercept: turnover time P val k2')
boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='(Intercept)' & k2.monte.carlo.long$value.name=="Pr(>|t|)"])

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='x' & k2.monte.carlo.long$value.name=="Estimate"])
title('x: slope (1/Vmax) k2')
boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='x' & k2.monte.carlo.long$value.name=="Estimate"])

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='x' & k2.monte.carlo.long$value.name=="Pr(>|t|)"])
title('x: slope (1/Vmax) P val k2')
# boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='x' & k2.monte.carlo.long$value.name=="Pr(>|t|)"])

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='model.summary' & k2.monte.carlo.long$value.name=="adj.r.sqrd"])
title('adjusted R squared k2')
boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='model.summary' & k2.monte.carlo.long$value.name=="adj.r.sqrd"])

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='model.summary' & k2.monte.carlo.long$value.name=="total.cor"])
title('total corelation k2')
boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='model.summary' & k2.monte.carlo.long$value.name=="total.cor"])

boxplot(k2.monte.carlo.long$Value[k2.monte.carlo.long$value.name=='sigma'])
title('sigma k2')
boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$value.name=='sigma'])


## filter outliers based on adj R sqr and total.cor
R.outliers.k2 <- k2.monte.carlo.long$Model[which(k2.monte.carlo.long$Value %in% boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='model.summary' & k2.monte.carlo.long$value.name=="adj.r.sqrd"])$out)]
Corr.outliers.k2 <- k2.monte.carlo.long$Model[which(k2.monte.carlo.long$Value %in% boxplot.stats(k2.monte.carlo.long$Value[k2.monte.carlo.long$coeff=='model.summary' & k2.monte.carlo.long$value.name=="total.cor"])$out)]

all.outliers <- unique(c(R.outliers.k2, Corr.outliers.k2))
k2.monte.carlo.long.no.outliers <- k2.monte.carlo.long %>%
  filter(!(Model %in% all.outliers))
## propegate error for the parameters
k2.summary.with.error <- k2.monte.carlo.long.no.outliers %>%
  group_by(coeff, value.name) %>%
  summarise(Mean.Params = mean(Value),
            sd = sd(Value),
            n = n())

k2.summary.tidy <- k2.summary.with.error %>%
  filter(value.name!="Pr(>|t|)",
         value.name!="t value",
         value.name !="Std. Error") %>%
  ## Then Kt+S = -x intercept and Vmax = 1/slope and turnover time = y-intercept
  mutate(coeff = ifelse(coeff == "(Intercept)", "Turnover Time (hr)",
                        ifelse(coeff == "x","slope",coeff)))


Vmax.nM.k2.monte = 1/(k2.summary.tidy$Mean.Params[k2.summary.tidy$coeff=="slope"])
Vmax.nM.k2.monte.sd = Vmax.nM.k2.monte * ((k2.summary.tidy$sd[k2.summary.tidy$coeff=="slope"])/(k2.summary.tidy$Mean.Params[k2.summary.tidy$coeff=="slope"]))
Kt.S.k2.monte = (k2.summary.tidy$Mean.Params[k2.summary.tidy$coeff=="Turnover Time (hr)"])/(k2.summary.tidy$Mean.Params[k2.summary.tidy$coeff=="slope"])
Kt.S.k2.monte.sd = Kt.S.k2.monte * sqrt((k2.summary.tidy$sd[k2.summary.tidy$coeff=="Turnover Time (hr)"]/k2.summary.tidy$Mean.Params[k2.summary.tidy$coeff=="Turnover Time (hr)"])^2 +
                                          (k2.summary.tidy$sd[k2.summary.tidy$coeff=="slope"]/k2.summary.tidy$Mean.Params[k2.summary.tidy$coeff=="slope"])^2
)
k2.summary.w.ests <- data.frame(coeff = c("Vmax (nM/hr)","Kt + S (nM)"),
                                value.name = "calculated",
                                Mean.Params = c(Vmax.nM.k2.monte, Kt.S.k2.monte),
                                sd = c(Vmax.nM.k2.monte.sd, Kt.S.k2.monte.sd),
                                n = mean(k2.summary.tidy$n)) %>%
  full_join(k2.summary.tidy)

## write out that tidy bit
write_csv(k2.summary.w.ests, path = paste("../output/k2-Linear_transformation_MonteCarlo_OutliersRemoved_summary_withError.csv",sep="_"))
