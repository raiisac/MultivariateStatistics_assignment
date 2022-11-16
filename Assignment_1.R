## ----setup, include=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------

update.packages("rlang")
library(xfun)
library(tidyverse)
library(candisc)
library(rprojroot)
library(lavaan)
library(tidySEM)#to graphically represent SEM
library(kableExtra)
library(gtools)
knitr::opts_chunk$set(echo = TRUE)


## ----load, include=FALSE, message=FALSE, warning=FALSE-----------------------------------------------------------------
load(find_root_file("data/cosmetics.Rdata", 
     criterion = has_file("MultivariateStatistics_assignment.Rproj")))
colnames(cosmetics) <- str_replace(colnames(cosmetics), 
                                   pattern = "Attitude_", 
                                   "A_") #shorten the names of the variables


## ----CFA1, include=TRUE, message=FALSE, warning=FALSE------------------------------------------------------------------
#We first standardize the variables
cosmetics_std <- scale(cosmetics, center = TRUE, scale = FALSE)
covmat1 <- cov(cosmetics_std[,1:9])
simplemodel1 <- 
'organic = ~1*A_organic1 + A_organic2 + A_organic3
  packaging = ~1*A_packaging1 + A_packaging2 + A_packaging3
  crueltyfree = ~1*A_crueltyfree1 + A_crueltyfree2 + A_crueltyfree3
  organic ~~ organic
  packaging ~~ packaging
  crueltyfree ~~ crueltyfree
  organic ~~ packaging
  organic ~~ crueltyfree
  packaging ~~ crueltyfree'
fit1 <- cfa(simplemodel1, sample.cov = covmat1, sample.nobs = nrow(cosmetics))
sum_fit1 <- summary(fit1, fit.measure = T)
sum_fit1_std <- standardizedSolution(fit1)


## ----CFA1performance, echo=FALSE, message=FALSE, warning=FALSE---------------------------------------------------------
data.frame(
  parameter = c('user model Chisq. (df)',
                "baseline model Chisq. (df)",
                "comparative fit index (CFI)",
                "Tucker-Lewis index (TLI)",
                "RMSEA (ll,ul)",
                "Standardized root mean square residual"),
  model1 = c(sprintf("%.2f (%.0f)%s", sum_fit1$fit["chisq"],
                    sum_fit1$fit["df"], stars.pval(sum_fit1$fit["pvalue"])),
            sprintf("%.2f (%.0f) %s", sum_fit1$fit["baseline.chisq"],
                    sum_fit1$fit["baseline.df"], 
                    stars.pval(sum_fit1$fit["baseline.pvalue"])),
            round(sum_fit1$fit["cfi"], digits = 3),
            round(sum_fit1$fit["tli"], digits = 3),
            sprintf("%.2f (%.2f, %.2f)%s", sum_fit1$fit["rmsea"],
                    sum_fit1$fit["rmsea.ci.lower"], 
                    sum_fit1$fit["rmsea.ci.upper"], 
                    stars.pval(sum_fit1$fit["rmsea.pvalue"])),
            round(sum_fit1$fit["srmr"], digits = 3))
) -> mod1performance



## ----CFA1standardized, echo=FALSE, message=FALSE, warning=FALSE--------------------------------------------------------
t1 <- sum_fit1_std[1:9,] %>%
  mutate(std_loading = sprintf("%s %s %s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                         stars)) %>%
  dplyr::select(std_loading, value)# %>%
  # kable(col.names = c("Variables","Loading (LL, UL)"),
  #       align = "lc",
  #       #caption = "Standardised loadings.",
  #       booktabs = TRUE,
  #       format = "latex")
t2 <- sum_fit1_std[10:24,] %>%
  mutate(`std_error.variance` = sprintf("%s%s%s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                         stars)) %>%
  dplyr::select(`std_error.variance`, value) #%>%
  # kable(col.names = c("Variable","Error variance (LL, UL)"),
  #       align = "lc",
  #       #caption = "Error variances.",
  #       booktabs = TRUE,
  #       row.names = FALSE,
  #       format = "latex")
#reliability factor scores
compositerel<-function(x){
  A <- (sum(x))^2
  B <- sum(1-x^2)
  return(A/(A+B))
}

#overview table composite reliability
factor <- c("organic","packaging","crueltyfree")
reliability <- round(c(compositerel(sum_fit1_std[1:3,4]),
                     compositerel(sum_fit1_std[4:6,4]),
                     compositerel(sum_fit1_std[7:9,4])),3)
t3 <- data.frame(factor,reliability)

knitr::kable(
  list(t1, t2, t3),
  caption = 'The solution of the simple model for the attitudes.',
  booktabs = TRUE, valign = 't'
)
fit1a <- fit1


## ----CFA1corr, include=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------
corrmodel1 <- 
'organic = ~1*A_organic1 + A_organic2 + A_organic3
  packaging = ~1*A_packaging1 + A_packaging2 + A_packaging3
  crueltyfree = ~1*A_crueltyfree1 + A_crueltyfree2 + A_crueltyfree3
  
  A_organic1 ~~c*A_packaging1 
  A_organic1 ~~c*A_crueltyfree1
  A_packaging1 ~~c*A_crueltyfree1
  A_organic2 ~~d*A_packaging2 
  A_organic2 ~~d*A_crueltyfree2
  A_packaging2 ~~d*A_crueltyfree2
  A_organic3 ~~e*A_packaging3 
  A_organic3 ~~e*A_crueltyfree3
  A_packaging3 ~~e*A_crueltyfree3
  
  organic ~~ organic
  packaging ~~ packaging
  crueltyfree ~~ crueltyfree
  organic ~~ packaging
  organic ~~ crueltyfree
  packaging ~~ crueltyfree
  '
fit1corr <- cfa(corrmodel1, sample.cov = covmat1, sample.nobs = nrow(cosmetics))
sum_fit1corr <- summary(fit1corr, fit.measure = T)
sum_fit1_std_corr <- standardizedSolution(fit1corr)


## ----CFA1performancecorr, echo=FALSE, message=FALSE, warning=FALSE, include = FALSE------------------------------------
data.frame(
  parameter = c('user model Chisq. (df)',
                "baseline model Chisq. (df)",
                "comparative fit index (CFI)",
                "Tucker-Lewis index (TLI)",
                "RMSEA (ll,ul)",
                "Standardized root mean square residual"),
  model2 = c(sprintf("%.2f (%.0f)%s", sum_fit1corr$fit["chisq"],
                    sum_fit1corr$fit["df"], stars.pval(sum_fit1corr$fit["pvalue"])),
            sprintf("%.2f (%.0f) %s", sum_fit1$fit["baseline.chisq"],
                    sum_fit1corr$fit["baseline.df"], 
                    stars.pval(sum_fit1corr$fit["baseline.pvalue"])),
            round(sum_fit1corr$fit["cfi"], digits = 3),
            round(sum_fit1corr$fit["tli"], digits = 3),
            sprintf("%.2f (%.2f, %.2f)%s", sum_fit1corr$fit["rmsea"],
                    sum_fit1corr$fit["rmsea.ci.lower"], 
                    sum_fit1corr$fit["rmsea.ci.upper"], 
                    stars.pval(sum_fit1corr$fit["rmsea.pvalue"])),
            round(sum_fit1corr$fit["srmr"], digits = 3))
) -> mod2performance



## ----CFA1standardizedcorr, echo=FALSE, message=FALSE, warning=FALSE, include = FALSE-----------------------------------
t1 <- sum_fit1_std_corr[1:9,] %>%
  mutate(loading = sprintf("%s %s %s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                         stars)) %>%
  dplyr::select(loading, value)# %>%
  # kable(col.names = c("Variables","Loading (LL, UL)"),
  #       align = "lc",
  #       #caption = "Standardised loadings.",
  #       booktabs = TRUE,
  #       format = "latex")
t2 <- sum_fit1corr$pe[10:33,] %>%
  mutate(`(co)variance` = sprintf("%s%s%s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f %s", est, stars)) %>%
  dplyr::select(`(co)variance`, value) #%>% #%>%
  # kable(col.names = c("Variable","Error variance (LL, UL)"),
  #       align = "lc",
  #       #caption = "Error variances.",
  #       booktabs = TRUE,
  #       row.names = FALSE,
  #       format = "latex")
# t3 <- sum_fit1_std_corr[10:21,] %>%
#   mutate(resid.correlation = sprintf("%s %s %s", lhs, op, rhs),
#          stars = stars.pval(pvalue),
#          value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
#                          stars)) %>%
#   dplyr::select(resid.correlation, value)
knitr::kable(
  list(t1, t2),
  caption = 'The standardized solution of the model with correlated error terms for the attitudes.',
  booktabs = TRUE, valign = 't'
)
fit1acorr <- fit1corr


## ----CFA1LRtest, echo=FALSE, message=FALSE, warning=FALSE--------------------------------------------------------------
#Anova test
test <- anova(fit1corr, fit1)
anovap <-  ifelse(test$`Pr(>Chisq)`[2]<0.001,"$<0.001$",ifelse(test$`Pr(>Chisq)`[2]<0.01,"$<0.01$", ifelse(test$`Pr(>Chisq)`[2]<0.01,"$<0.05$","$>0.05$")))
# Residual correlations
cor_table <- residuals(fit1corr, type = "cor")$cov

#compute number of residual correlations with absolute value larger than 0.05 below the diagonal
n<-sum(ifelse(abs(cor_table)>0.05,1,0))/2

corrmodel1 <- 
    'organic = ~1*A_organic1 + A_organic2 + A_organic3
  packaging = ~1*A_packaging1 + A_packaging2 + A_packaging3
  crueltyfree = ~1*A_crueltyfree1 + A_crueltyfree2 + A_crueltyfree3
  A_organic1 ~~c*A_packaging1 
  A_organic1 ~~c*A_crueltyfree1
  A_packaging1 ~~c*A_crueltyfree1
  
  A_organic2 ~~d*A_packaging2 
  A_organic2 ~~d*A_crueltyfree2
  A_packaging2 ~~d*A_crueltyfree2
  
  A_organic3 ~~A_packaging3 
  A_organic3 ~~A_crueltyfree3
  A_packaging3 ~~A_crueltyfree3
  
  organic ~~ organic
  packaging ~~ packaging
  crueltyfree ~~ crueltyfree
  
  organic ~~ packaging
  organic ~~ crueltyfree
  packaging ~~ crueltyfree
  '
fit1corr_new <- cfa(corrmodel1, sample.cov = covmat1, sample.nobs = nrow(cosmetics))
sum_fit1corr_new <- summary(fit1corr_new, fit.measure = T)
sum_fit1_std_corr_new <- standardizedSolution(fit1corr_new)


## ----CFA1BI, include=TRUE, message=FALSE, warning=FALSE----------------------------------------------------------------
#We first standardize the variables
covmat1 <- cov(cosmetics_std[,10:18])
simplemodel1 <- 
'organic = ~1*BI_organic1 + BI_organic2 + BI_organic3
  packaging = ~1*BI_packaging1 + BI_packaging2 + BI_packaging3
  crueltyfree = ~1*BI_crueltyfree1 + BI_crueltyfree2 + BI_crueltyfree3
  organic ~~ organic
  packaging ~~ packaging
  crueltyfree ~~ crueltyfree
  organic ~~ packaging
  organic ~~ crueltyfree
  packaging ~~ crueltyfree'
fit1 <- cfa(simplemodel1, sample.cov = covmat1, sample.nobs = nrow(cosmetics))
sum_fit1 <- summary(fit1, fit.measure = T)
sum_fit1_std <- standardizedSolution(fit1)


## ----CFA1performanceBI, echo=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------
data.frame(
  parameter = c('user model Chisq. (df)',
                "baseline model Chisq. (df)",
                "comparative fit index (CFI)",
                "Tucker-Lewis index (TLI)",
                "RMSEA (ll,ul)",
                "Standardized root mean square residual"),
  model3 = c(sprintf("%.2f (%.0f)%s", sum_fit1$fit["chisq"],
                    sum_fit1$fit["df"], stars.pval(sum_fit1$fit["pvalue"])),
            sprintf("%.2f (%.0f) %s", sum_fit1$fit["baseline.chisq"],
                    sum_fit1$fit["baseline.df"], 
                    stars.pval(sum_fit1$fit["baseline.pvalue"])),
            round(sum_fit1$fit["cfi"], digits = 3),
            round(sum_fit1$fit["tli"], digits = 3),
            sprintf("%.2f (%.2f, %.2f)%s", sum_fit1$fit["rmsea"],
                    sum_fit1$fit["rmsea.ci.lower"], 
                    sum_fit1$fit["rmsea.ci.upper"], 
                    stars.pval(sum_fit1$fit["rmsea.pvalue"])),
            round(sum_fit1$fit["srmr"], digits = 3))
) -> mod3performance



## ----CFA1standardizedBI, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------
t1 <- sum_fit1_std[1:9,] %>%
  mutate(std_loading = sprintf("%s %s %s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                         stars)) %>%
  dplyr::select(std_loading, value)# %>%
  # kable(col.names = c("Variables","Loading (LL, UL)"),
  #       align = "lc",
  #       #caption = "Standardised loadings.",
  #       booktabs = TRUE,
  #       format = "latex")
t2 <- sum_fit1_std[10:24,] %>%
  mutate(std_error.variance = sprintf("%s%s%s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                       stars)) %>%
  dplyr::select(std_error.variance, value) #%>%
  # kable(col.names = c("Variable","Error variance (LL, UL)"),
  #       align = "lc",
  #       #caption = "Error variances.",
  #       booktabs = TRUE,
  #       row.names = FALSE,
  #       format = "latex")
factor <- c("organic","packaging","crueltyfree")
reliability <- round(c(compositerel(sum_fit1_std[1:3,4]),
                     compositerel(sum_fit1_std[4:6,4]),
                     compositerel(sum_fit1_std[7:9,4])),3)
t3 <- data.frame(factor,reliability)

knitr::kable(
  list(t1, t2, t3),
  caption = 'The standardized solution of the simple model for the behavior-intent items.',
  booktabs = TRUE, valign = 't'
)


## ----CFA1corrBI, include=TRUE, message=FALSE, warning=FALSE------------------------------------------------------------
corrmodel1 <- 
'organic = ~1*BI_organic1 + BI_organic2 + BI_organic3
  packaging = ~1*BI_packaging1 + BI_packaging2 + BI_packaging3
  crueltyfree = ~1*BI_crueltyfree1 + BI_crueltyfree2 + BI_crueltyfree3
  
  BI_organic1 ~~c*BI_packaging1 
  BI_organic1 ~~c*BI_crueltyfree1
  BI_packaging1 ~~c*BI_crueltyfree1
  BI_organic2 ~~d*BI_packaging2 
  BI_organic2 ~~d*BI_crueltyfree2
  BI_packaging2 ~~d*BI_crueltyfree2
  BI_organic3 ~~e*BI_packaging3 
  BI_organic3 ~~e*BI_crueltyfree3
  BI_packaging3 ~~e*BI_crueltyfree3
  
  organic ~~ organic
  packaging ~~ packaging
  crueltyfree ~~ crueltyfree
  organic ~~ packaging
  organic ~~ crueltyfree
  packaging ~~ crueltyfree
  '
fit1corr <- cfa(corrmodel1, sample.cov = covmat1, sample.nobs = nrow(cosmetics))
sum_fit1corr <- summary(fit1corr, fit.measure = T)
sum_fit1_std_corr <- standardizedSolution(fit1corr)


## ----CFA1performancecorrBI, echo=FALSE, message=FALSE, warning=FALSE---------------------------------------------------
data.frame(
  parameter = c('user model Chisq. (df)',
                "baseline model Chisq. (df)",
                "comparative fit index (CFI)",
                "Tucker-Lewis index (TLI)",
                "RMSEA (ll,ul)",
                "Standardized root mean square residual"),
  model4 = c(sprintf("%.2f (%.0f)%s", sum_fit1corr$fit["chisq"],
                    sum_fit1corr$fit["df"], stars.pval(sum_fit1corr$fit["pvalue"])),
            sprintf("%.2f (%.0f) %s", sum_fit1$fit["baseline.chisq"],
                    sum_fit1corr$fit["baseline.df"], 
                    stars.pval(sum_fit1corr$fit["baseline.pvalue"])),
            round(sum_fit1corr$fit["cfi"], digits = 3),
            round(sum_fit1corr$fit["tli"], digits = 3),
            sprintf("%.2f (%.2f, %.2f)%s", sum_fit1corr$fit["rmsea"],
                    sum_fit1corr$fit["rmsea.ci.lower"], 
                    sum_fit1corr$fit["rmsea.ci.upper"], 
                    stars.pval(sum_fit1corr$fit["rmsea.pvalue"])),
            round(sum_fit1corr$fit["srmr"], digits = 3))) %>%
  left_join(mod1performance) %>%
  left_join(mod2performance) %>%
  left_join(mod3performance) %>%
  dplyr::select(parameter, model1, model2, model3, model4) %>%
  kable(booktabs = T,
        col.names = c("parameter", "simple model", 
                      "with correlated error terms","simple model", 
                      "with correlated error terms"),
        caption = "Performance measure for the different models") %>%
  add_header_above(c(" " = 1, "Attitudes" = 2, "Behavior-intention" = 2)) %>%
  column_spec(1, width = "3cm") %>%
  column_spec(2:5, width = "2.8cm") 
  



## ----CFA1standardizedcorrBI, echo=FALSE, message=FALSE, warning=FALSE, include = FALSE---------------------------------
t1 <- sum_fit1_std_corr[1:9,] %>%
  mutate(loading = sprintf("%s %s %s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                         stars)) %>%
  dplyr::select(loading, value)# %>%
  # kable(col.names = c("Variables","Loading (LL, UL)"),
  #       align = "lc",
  #       #caption = "Standardised loadings.",
  #       booktabs = TRUE,
  #       format = "latex")
t2 <- sum_fit1_std_corr[22:33,] %>%
  mutate(error.variance = lhs,
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                       stars)) %>%
  dplyr::select(error.variance, value) #%>%
  # kable(col.names = c("Variable","Error variance (LL, UL)"),
  #       align = "lc",
  #       #caption = "Error variances.",
  #       booktabs = TRUE,
  #       row.names = FALSE,
  #       format = "latex")
t3 <- sum_fit1_std_corr[10:18,] %>%
  mutate(resid.correlation = sprintf("%s %s %s", lhs, op, rhs),
         stars = stars.pval(pvalue),
         value = sprintf("%.2f (%.2f, %.2f)%s", est.std, ci.lower, ci.upper, 
                         stars)) %>%
  dplyr::select(resid.correlation, value)
knitr::kable(
  list(t1, t2, t3),
  caption = 'The standardized solution of the model with correlated error terms for the behavior-intent items.',
  booktabs = TRUE, valign = 't'
)


## ----CFA1LRtestBI, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------
test <- anova(fit1corr, fit1)
anovap <-  ifelse(test$`Pr(>Chisq)`[2]<0.001,"$<0.001$",ifelse(test$`Pr(>Chisq)`[2]<0.01,"$<0.01$", ifelse(test$`Pr(>Chisq)`[2]<0.01,"$<0.05$","$>0.05$")))
cor_table <- residuals(fit1corr, type = "cor")$cov

#compute number of residual correlations with absolute value larger than 0.05 below the diagonal
n <- sum(ifelse(abs(cor_table)>0.05,1,0))/2

cor_table_simple <- residuals(fit1corr, type = "cor")$cov

#compute number of residual correlations with absolute value larger than 0.05 below the diagonal
n_simple <- sum(ifelse(abs(cor_table)>0.05,1,0))/2


## ----SEM, echo=TRUE, message=FALSE, warning=FALSE----------------------------------------------------------------------
cormat <- cov(cosmetics_std)
sem1 <- 'BI_organic = ~1*BI_organic1 + BI_organic2 + BI_organic3
  BI_packaging = ~1*BI_packaging1 + BI_packaging2 + BI_packaging3
  BI_crueltyfree = ~1*BI_crueltyfree1 + BI_crueltyfree2 + BI_crueltyfree3
  BI_organic1 ~~c*BI_packaging1 
  BI_organic1 ~~c*BI_crueltyfree1
  BI_packaging1 ~~c*BI_crueltyfree1
  BI_organic2 ~~d*BI_packaging2 
  BI_organic2 ~~d*BI_crueltyfree2
  BI_packaging2 ~~d*BI_crueltyfree2
  BI_organic3 ~~e*BI_packaging3 
  BI_organic3 ~~e*BI_crueltyfree3
  BI_packaging3 ~~e*BI_crueltyfree3
  BI_organic ~~ BI_organic
  BI_packaging ~~ BI_packaging
  BI_crueltyfree ~~ BI_crueltyfree
  BI_organic ~~ BI_packaging
  BI_organic ~~ BI_crueltyfree
  BI_packaging ~~ BI_crueltyfree
  
  A_organic = ~1*A_organic1 + A_organic2 + A_organic3
  A_packaging = ~1*A_packaging1 + A_packaging2 + A_packaging3
  A_crueltyfree = ~1*A_crueltyfree1 + A_crueltyfree2 + A_crueltyfree3
  A_organic1 ~~a*A_packaging1 
  A_organic1 ~~a*A_crueltyfree1
  A_packaging1 ~~a*A_crueltyfree1
  A_organic2 ~~b*A_packaging2 
  A_organic2 ~~b*A_crueltyfree2
  A_packaging2 ~~b*A_crueltyfree2
  A_organic3 ~~A_packaging3 
  A_organic3 ~~A_crueltyfree3
  A_packaging3 ~~A_crueltyfree3
  A_organic ~~ A_organic
  A_packaging ~~ A_packaging
  A_crueltyfree ~~ A_crueltyfree
  A_organic ~~ A_packaging
  A_organic ~~ A_crueltyfree
  A_packaging ~~ A_crueltyfree
  
  #structural model
  BI_organic ~A_organic
  BI_packaging ~A_packaging
  BI_crueltyfree ~A_crueltyfree
  '
fitsem1 <- sem(sem1, sample.cov = cormat, sample.nobs = nrow(cosmetics))
sum_sem1 <- summary(fitsem1)
sum_sem1_std <- standardizedSolution(fitsem1)


## ----SEM2, echo=FALSE, message=FALSE, warning=FALSE--------------------------------------------------------------------
#summary SEM results
sum_sem1$pe[49:51,] %>% 
  mutate (Regression_coefficient = sprintf("%s%s\n%s",lhs, op, rhs),
          stars = stars.pval(pvalue),
         sem1 = sprintf("%.2f %s", est, 
                         stars)) %>%
  dplyr::select(Regression_coefficient, sem1) -> sem1regression
#standardized results
sum_sem1_std[49:51,] %>% 
  mutate (Regression_coefficient = sprintf("%s%s\n%s",lhs, op, rhs),
          stars = stars.pval(pvalue),
         sem1_std = sprintf("%.2f %s", est.std, 
                         stars))%>%
  dplyr::select(Regression_coefficient, sem1_std) -> sem1regression_std


## ----echo=T, results='hide'--------------------------------------------------------------------------------------------
'  #structural model
  BI_organic ~p*A_organic
  BI_packaging ~p*A_packaging
  BI_crueltyfree ~p*A_crueltyfree'


## ----SEM2_1, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------------
sem2 <- 'BI_organic = ~1*BI_organic1 + BI_organic2 + BI_organic3
  BI_packaging = ~1*BI_packaging1 + BI_packaging2 + BI_packaging3
  BI_crueltyfree = ~1*BI_crueltyfree1 + BI_crueltyfree2 + BI_crueltyfree3
  
  BI_organic1 ~~c*BI_packaging1 
  BI_organic1 ~~c*BI_crueltyfree1
  BI_packaging1 ~~c*BI_crueltyfree1
  BI_organic2 ~~d*BI_packaging2 
  BI_organic2 ~~d*BI_crueltyfree2
  BI_packaging2 ~~d*BI_crueltyfree2
  BI_organic3 ~~e*BI_packaging3 
  BI_organic3 ~~e*BI_crueltyfree3
  BI_packaging3 ~~e*BI_crueltyfree3
  BI_organic ~~ BI_organic
  BI_packaging ~~ BI_packaging
  BI_crueltyfree ~~ BI_crueltyfree
  BI_organic ~~ BI_packaging
  BI_organic ~~ BI_crueltyfree
  BI_packaging ~~ BI_crueltyfree
  
  A_organic = ~1*A_organic1 + A_organic2 + A_organic3
  A_packaging = ~1*A_packaging1 + A_packaging2 + A_packaging3
  A_crueltyfree = ~1*A_crueltyfree1 + A_crueltyfree2 + A_crueltyfree3
  
  A_organic1 ~~a*A_packaging1 
  A_organic1 ~~a*A_crueltyfree1
  A_packaging1 ~~a*A_crueltyfree1
  A_organic2 ~~b*A_packaging2 
  A_organic2 ~~b*A_crueltyfree2
  A_packaging2 ~~b*A_crueltyfree2
  A_organic3 ~~A_packaging3 
  A_organic3 ~~A_crueltyfree3
  A_packaging3 ~~A_crueltyfree3
  
  A_organic ~~ A_organic
  A_packaging ~~ A_packaging
  A_crueltyfree ~~ A_crueltyfree
  A_organic ~~ A_packaging
  A_organic ~~ A_crueltyfree
  A_packaging ~~ A_crueltyfree
  
  #structural model
  BI_organic ~p*A_organic
  BI_packaging ~p*A_packaging
  BI_crueltyfree ~p*A_crueltyfree
  '
fitsem2 <- sem(sem2,sample.cov = cormat, sample.nobs = nrow(cosmetics))
sum_sem2 <- summary(fitsem2)
sum_sem2_std <- standardizedSolution(fitsem2)
anovasem <- anova(fitsem1,fitsem2)
anovasemp <- anovasem$`Pr(>Chisq)`[2]

## ----SEMregression, echo=FALSE, message=FALSE, warning=FALSE-----------------------------------------------------------
sum_sem2$pe[49:51,] %>% 
  mutate (Regression_coefficient = sprintf("%s%s\n%s",lhs, op, rhs),
          stars = stars.pval(pvalue),
         sem2 = sprintf("%.2f %s", est, 
                         stars)) %>%
  dplyr::select(Regression_coefficient, sem2) -> sem2regression
#standardized results
sum_sem2_std[49:51,] %>% 
  mutate (Regression_coefficient = sprintf("%s%s\n%s",lhs, op, rhs),
          stars = stars.pval(pvalue),
         sem2_std = sprintf("%.2f %s", est.std, 
                         stars)) %>%
  left_join(sem1regression) %>%
  left_join(sem1regression_std) %>%
  left_join(sem2regression) %>%
  dplyr::select(Regression_coefficient, sem1, sem1_std, sem2, sem2_std) %>%
  kable(booktabs = T,
        col.names = c("Regression coefficient",	"unstandardized",
                      "standardized",	"unstandardized",	"standardized"),
        caption = "Population regression coefficients in both SEMs.") %>%
  add_header_above(c(" " = 1, "General SEM" = 2, "Equal population regression coefficients" = 2)) %>%
  column_spec(1, width = "3cm") %>%
  column_spec(2:5, width = "2.8cm") 
  


## ----load2, include=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------
load(find_root_file("data/benefits.Rdata", 
     criterion = has_file("MultivariateStatistics_assignment.Rproj")))


## ----cca, include=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------
zbenefits <- benefits
zbenefits[, 2:14] <- scale(zbenefits[, 2:14], scale = TRUE, center = TRUE) 


cancor.out <- cancor(cbind(SL_pensioners, SL_unemployed, SL_old_gvntresp,
                           SL_unemp_gvntresp) 
                     ~ SB_strain_economy + SB_prevent_poverty + 
                       SB_equal_society + SB_taxes_business +
                       SB_make_lazy + SB_caring_others + 
                       unemployed_notmotivated + SB_often_lessthanentitled +
                       SB_often_notentitled, data = zbenefits) 
#print summary results 
summary(cancor.out)
#compute redundancies 
R2tu <- cancor.out$cancor^2 
R2tu <- cancor.out$cancor^2 
VAFYbyt <- apply(cancor.out$structure$Y.yscores^2, 2, sum)/3 
redund <- R2tu*VAFYbyt 
round(cbind(R2tu,VAFYbyt,redund,total = cumsum(redund)), 4) 
#print canonical loadings 
round(cancor.out$structure$X.xscores, 2) 
round(cancor.out$structure$Y.yscores, 2)


## ----spa, include=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------
train <- benefits[seq(2,3310, by = 2), ]
valid <- benefits[seq(1,3310, by = 2), ]
train[,2:14] <- scale(train[, 2:14], center = TRUE, scale = TRUE)
valid[,2:14] <- scale(valid[, 2:14], center = TRUE, scale = TRUE)

#conduct CCA on training data
cancor.train <- cancor(cbind(SL_pensioners, SL_unemployed, SL_old_gvntresp,
                             SL_unemp_gvntresp) 
                       ~ SB_strain_economy + SB_prevent_poverty + 
                         SB_equal_society + SB_taxes_business +
                         SB_make_lazy + SB_caring_others + 
                         unemployed_notmotivated + SB_often_lessthanentitled +
                         SB_often_notentitled, data = train) 
#conduct CCA on validation data
cancor.valid <- cancor(cbind(SL_pensioners, SL_unemployed, SL_old_gvntresp,
                             SL_unemp_gvntresp) 
                       ~ SB_strain_economy + SB_prevent_poverty + 
                         SB_equal_society + SB_taxes_business +
                         SB_make_lazy + SB_caring_others +
                         unemployed_notmotivated + SB_often_lessthanentitled +
                         SB_often_notentitled, data = valid) 
# canonical variates calibration set
train.X1 <- cancor.train$score$X
train.Y1 <- cancor.train$score$Y
# compute canonical variates using data of calibration set and coefficients
# estimated on validation set
train.X2 <- as.matrix(train[,6:14]) %*% cancor.valid$coef$X
train.Y2 <- as.matrix(train[,2:5]) %*% cancor.valid$coef$Y
round(cor(train.Y1,train.Y2),3) 
round(cor(train.X1,train.X2),3) 
round(cor(train.X1,train.Y1),3) 
round(cor(train.X2,train.Y2),3) 
round(cor(train.Y2,train.Y2),3) 
round(cor(train.X2,train.X2),3)


## ----CFA1graphical, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "A graphical representation of the simple model for the attitudes.", out.width = "15cm", fig.align='center'----
lay <- get_layout("", "", "organic","","packaging","","crueltyfree","", "",
                  "A_organic1", "A_organic2", "A_organic3",
                  "A_packaging1", "A_packaging2", "A_packaging3",
                  "A_crueltyfree1", "A_crueltyfree2", "A_crueltyfree3",
                  rows = 2)
p <- graph_sem(model = fit1a, layout = lay)
if (!file.exists("figures/CFA1graphical.png")) {
  ggsave("figures/CFA1graphical.png", p,
         device = "png", width = 11, height = 4)}
knitr::include_graphics(find_root_file("figures/CFA1graphical.png",
    criterion = has_file("MultivariateStatistics_assignment.Rproj")))


## ----CFA1graphicalcorr, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "A graphical representation of the model for the attitudeswith correlated error terms for all pairs of items that focus on the same aspect.", out.width = "15cm", fig.align='center'----
lay <- get_layout("", "", "organic","","packaging","","crueltyfree","", "",
                  "A_organic1", "A_organic2", "A_organic3",
                  "A_packaging1", "A_packaging2", "A_packaging3",
                  "A_crueltyfree1", "A_crueltyfree2", "A_crueltyfree3",
                  "",  "right","","","pleasant","","","must","",
                  rows = 3)
p <- graph_sem(model = fit1acorr, layout = lay)
if (!file.exists("figures/CFA1graphicalcorr.png")) {
  ggsave("figures/CFA1graphicalcorr.png", p,
         device = "png", width = 11, height = 5)
  }
knitr::include_graphics(find_root_file("figures/CFA1graphicalcorr.png",
    criterion = has_file("MultivariateStatistics_assignment.Rproj")))


## ----CFA1graphicalBI, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "A graphical representation of the simple model for the behavior-intent items.", out.width = "15cm", fig.align='center'----
lay <- get_layout("", "", "organic","","packaging","","crueltyfree","", "",
                  "BI_organic1", "BI_organic2", "BI_organic3",
                  "BI_packaging1", "BI_packaging2", "BI_packaging3",
                  "BI_crueltyfree1", "BI_crueltyfree2", "BI_crueltyfree3",
                  rows = 2)
p <- graph_sem(model = fit1, layout = lay)
if (!file.exists("figures/CFA1graphical_BI.png")) {
  ggsave(find_root_file("figures/CFA1graphical_BI.png",
    criterion = has_file("MultivariateStatistics_assignment.Rproj")), p,
         device = "png", width = 11, height = 4)}
knitr::include_graphics(find_root_file("figures/CFA1graphical_BI.png",
    criterion = has_file("MultivariateStatistics_assignment.Rproj")))


## ----CFA1graphicalcorrBI, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "A graphical representation of the model with correlated error terms for the behavior-intent items that focus on the same aspect.", out.width = "15cm", fig.align='center'----
lay <- get_layout("", "", "organic","","packaging","","crueltyfree","", "",
                  "BI_organic1", "BI_organic2", "BI_organic3",
                  "BI_packaging1", "BI_packaging2", "BI_packaging3",
                  "BI_crueltyfree1", "BI_crueltyfree2", "BI_crueltyfree3",
                  "",  "right","","","pleasant","","","must","",
                  rows = 3)
p <- graph_sem(model = fit1corr, layout = lay)
if (!file.exists("figures/CFA1graphicalcorrBI.png")) {
  ggsave("figures/CFA1graphicalcorrBI.png", p,
         device = "png", width = 11, height = 5)
  }
knitr::include_graphics(find_root_file("figures/CFA1graphicalcorrBI.png",
    criterion = has_file("MultivariateStatistics_assignment.Rproj")))

