##############################################
# File: main.R
# Description: Survival analysis main engine for HIV-1 data found at source below.
# Author: Kyle Lee
# Date: 5/21/2018
# 
##############################################


##############################################
# Load packages and data 
##############################################
library(survival)
library(stargazer)
library(survminer)


# Source: https://github.com/WinVector/QSurvival/blob/master/vignettes/AIDSdata/actg320Description.txt
file <- file.choose()
actg <- read.csv(file)

attach(actg)
##############################################
# Explorative Data Analysis
##############################################             
dim(actg)
# > dim(actg)
# [1] 1151   16
summary(actg)

# Export in LaTeX
stargazer(actg)

# Tables
table(actg$censor)
table(actg$censor)/sum(table(actg$censor))
table(actg$tx)
table(actg$tx)/sum(table(actg$tx))
table(actg$strat2)
table(actg$strat2)/sum(table(actg$strat2))
table(actg$sex)
table(actg$sex)/sum(table(actg$sex))
table(actg$raceth)
table(actg$raceth)/sum(table(actg$raceth))
table(actg$ivdrug)
table(actg$ivdrug)/sum(table(actg$ivdrug))
table(actg$hemophil)
table(actg$hemophil)/sum(table(actg$hemophil))
table(actg$karnof)
table(actg$karnof)/sum(table(actg$karnof))
table(actg$priorzdv)
table(actg$age)

# Plots
ggplot(actg, 
       aes(x = priorzdv), 
       ggtheme = theme_minimal()) + 
  ggtitle("Distribution of Months prior to ZDV treatment") + 
  geom_histogram(binwidth=5, colour="black", fill="white") +
  geom_vline(aes(xintercept = mean(actg$priorzdv, 
                                   na.rm=T)),   # Ignore NA values for mean
             color="red", 
             linetype = "dashed", 
             size = 1)

ggplot(actg, 
       aes(x = age), 
       ggtheme = theme_minimal()) + 
  ggtitle("Distribution of Patient Ages") + 
  geom_histogram(binwidth = 5, 
                 colour = "black", 
                 fill = "white") +
  geom_vline(aes(xintercept = mean(actg$age, 
                                   na.rm=T)),   # Ignore NA values for mean
             color = "red", 
             linetype = "dashed", 
             size = 1)

ggplot(actg, 
       aes(x = cd4), 
       ggtheme = theme_minimal()) + 
  ggtitle("Distribution of Baseline CD4 Cells") + 
  geom_histogram(binwidth = 5, 
                 colour = "black", 
                 fill = "white") +
  geom_vline(aes(xintercept = mean(actg$age, 
                                   na.rm=T)),   # Ignore NA values for mean
             color = "red", 
             linetype = "dashed", 
             size = 1)

##############################################
# Kaplan-Meyer Estimate
##############################################
# Ignore time_d and censor_d
# Plot KM curves
fitKM<-survfit(Surv(time, censor) ~ tx, 
               data = actg)
ggsurvplot(fitKM, 
           conf.int = TRUE, 
           legend.labs = unique(tx), 
           ggtheme = theme_minimal(), 
           title = "Kaplan Meyer Model by Treatment",
           data = actg) + 
  ylim(.8,1)
summary(fitKM)

##############################################
# Cox Proportional Hazards Model
##############################################
# Create full model
# Reason to believe that sex and race had an interaction
fit.int <- coxph(Surv(time, censor) ~ tx + strat2 + sex*factor(raceth) + 
                   factor(ivdrug) + hemophil + karnof +
                   cd4 + priorzdv + age, 
                 method = "breslow")
summary(fit.int)

# Standard All Covariates Model. Remove interaction effect.
fit.all<- coxph(Surv(time, censor) ~ tx + strat2 + sex + factor(raceth)+ 
                  factor(ivdrug) + hemophil + karnof +
                  cd4 + priorzdv + age, 
                method = "breslow")

summary(fit.all)

# Perform backward elimination by removing highest p-value

# Remove hemo
fit1 <- coxph(Surv(time, censor) ~ tx + strat2 + sex + factor(raceth)+ 
                factor(ivdrug)  + karnof +
                cd4 + priorzdv + age, 
              method = "breslow")
summary(fit1)

# Remove strat2
fit2 <- coxph(Surv(time, censor) ~ tx  + sex + factor(raceth)+ 
                factor(ivdrug)  + karnof +
                cd4 + priorzdv + age, 
              method = "breslow")
summary(fit2)

# Remove priorzdv
fit3 <- coxph(Surv(time, censor) ~ tx  + sex + factor(raceth)+ 
                factor(ivdrug)  + karnof +
                cd4  + age, method = "breslow")
summary(fit3)

# Remove sex
fit4 <- coxph(Surv(time, censor) ~ tx  + factor(raceth)+ 
                factor(ivdrug)  + karnof +
                cd4  + age, 
              method = "breslow")
summary(fit4)

# Remove race
fit5 <- coxph(Surv(time, censor) ~ tx  + 
                factor(ivdrug)  + karnof +
                cd4  + age, 
              method = "breslow")
summary(fit5)

# Compare factored race versus no factored race model
# Null: removed factors are not different from 0, full model is not better
# Alt: removed factors are different from 0, full model is better
# > -2*(fit5$loglik[2]-fit4$loglik[2])
# [1] 4.018123
# df = 10-6 = 4
qchisq(.95,4)
anova(fit4, fit5)

# Loglikelihood test
-2*(fit5$loglik[2]-fit4$loglik[2])

# Latex
stargazer(fit.all, fit1, fit2, fit3, fit4, fit5,
          align = T,
          font.size = "footnotesize", 
          style = "all2",
          omit.table.layout = "sn", 
          df = FALSE, 
          report = "vp*")

##############################################
# Model Selection (Stepwise)
##############################################
fake <- na.omit(actg)
# Stepwise
My.stepwise.coxph(Time = "time", Status = "censor", 
                  variable.list = c("tx" ,"strat2",
                                    "sex","raceth",
                                    "ivdrug","hemophil",
                                    "karnof","cd4",
                                    "priorzdv", "age", "age2"), 
                  data = fake)

##############################################
# Model Diagnostics: Round 1
##############################################
# Test for non-linearity and functional form of continuous covariates
ggcoxfunctional(Surv(time, censor) ~ age + log(age) + sqrt(age) + 
                  I(age^2)+cd4  + sqrt(cd4)  +karnof  + sqrt(karnof) + log(karnof),
                data = actg)

fit5.cd4a <- coxph(Surv(time, censor) ~ tx + factor(ivdrug) + karnof + sqrt(cd4) + age, 
                        method = "breslow")
summary(fit5.cd4)

fit5.cd4b <- coxph(Surv(time, censor) ~ tx + factor(ivdrug) + karnof + sqrt(cd4) + age + I(age^2), 
                   method = "breslow")

# Test to see if full model is the same or not to determine if age^2 should be included. 
anova(fit5.cd4a, fit5.cd4b)

# > anova(fit5.cd4a,fit5.cd4b)
# Analysis of Deviance Table
# Cox model: response is  Surv(time, censor)
# Model 1: ~ tx + factor(ivdrug) + karnof + sqrt(cd4) + age
# Model 2: ~ tx + factor(ivdrug) + karnof + sqrt(cd4) + age + I(age^2)
# loglik  Chisq Df P(>|Chi|)
# 1 -609.74                    
# 2 -609.62 0.2413  1    0.6233
# Conclude that full model and reduced model is the same. Therefore, use age

# Test for influential observations using dfbeta and deviance residuals
ggcoxdiagnostics(fit5.cd4b, 
                 type = "dfbeta", 
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())
# Some large values for sqrt(cd4) and age

# We can also use the deviance residuals. We expect it to be symmetric about 0 with sd = 1.
# Positive values mean patients were censored (death or diagnosis) earlier than anticipated
# Negative values mean that patients were censored (death or diagnosis) later than anticipated
ggcoxdiagnostics(fit5.cd4b, 
                 type = "deviance", 
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())

# Test for proportional hazards assumption
test.ph <- cox.zph(fit5.cd4a)
test.ph

# > test.ph
# rho  chisq     p
# tx              -0.1507 2.1566 0.142
# factor(ivdrug)2  0.0663 0.4057 0.524
# factor(ivdrug)3 -0.0209 0.0418 0.838
# karnof          -0.0458 0.1922 0.661
# sqrt(cd4)        0.2013 2.6661 0.103
# age              0.1290 1.6437 0.200
# GLOBAL               NA 7.8042 0.253
# Not statistically significant for each of the covariates, global test
# is also not statistically significant. Assume proportional hazards but with caution

# Latex
stargazer(test.ph)


# Test for proportional hazards assumption (Graphical)
# We expect for each covariate that under the PH assumption, the estimates of the betas do not vary
# over time and so we expect a straight line centered at 0. Any covariates that depart significantly
# from this would indicate that the proportional hazards assumption is violated
ggcoxzph(test.ph)

# Stratify treatment

##############################################
# Model Diagnostics: Round 2 (Stratifying tx)
##############################################
fit5.cd4.strtx <- coxph(Surv(time, censor) ~ strata(tx) + factor(ivdrug) + 
                            karnof + sqrt(cd4) + age + I(age^2), 
                          method = "breslow")
fit5.cd4.strtx2 <- coxph(Surv(time, censor) ~ strata(tx) + factor(ivdrug) + 
                           karnof + sqrt(cd4) + age , 
                         method = "breslow")

# Compare age versus age+age^2 
anova(fit5.cd4.strtx,fit5.cd4.strtx2)
summary(fit5.cd4.strtx)

# Test for influential observations using dfbeta and deviance residuals
ggcoxdiagnostics(fit5.cd4.strtx, 
                 type = "dfbeta", 
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())

# Some large values for sqrt(cd4) and age
ggcoxdiagnostics(fit5.cd4.strtx, 
                 type = "deviance", 
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())

# Test for proportional hazards assumption like before
test.ph <- cox.zph(fit5.cd4.strtx)
test.ph

# > test.ph
# rho  chisq     p
# factor(ivdrug)2  0.0599 0.3292 0.566
# factor(ivdrug)3 -0.0216 0.0445 0.833
# karnof          -0.0476 0.2050 0.651
# sqrt(cd4)        0.2001 2.6359 0.104
# age              0.1271 1.5924 0.207
# GLOBAL               NA 5.7793 0.328
# Not statistically significant for each of the covariates, global test
# is also not statistically significant. Assume proportional hazards but with caution

# Test for proportional hazards assumption (Graphical)
ggcoxzph(test.ph)
sf.fit <- survfit(fit5.cd4.strtx)

# Note that plots are parallel so proportional hazards assumption is met
ggsurvplot(sf.fit, 
           conf.int = TRUE, 
           legend.labs = unique(tx), 
           fun = "cloglog",
           ggtheme = theme_minimal(), 
           title = "Test Proportional Hazards Assumption",
           data = actg)

# Cox-PH model
ggsurvplot(sf.fit, 
           conf.int = TRUE, 
           legend.labs = unique(tx), 
           ggtheme = theme_minimal(), 
           title = "Cox Proportional Hazards Model by Treatment",
           data = actg) + 
  ylim(.8,1)


sresid<-resid(fit5.cd4.strtx, 
              'score')
dim(sresid)

# Plot score residuals
plot(1:1151,
     sresid[,1],
     xlab = 'factor(ivdrug)2',
     ylab = 'Residual', 
     main = "Score Residuals of factor(ivdrug)2")
identify(1:1151,sresid[,1],1:1151)

plot(1:1151,
     sresid[,2], 
     xlab = 'factor(ivdrug)3 ',
     ylab = 'Residual', 
     main = "Score Residuals of factor(ivdrug)3 ")
identify(1:1151,sresid[,2],1:1151)

plot(1:1151,
     sresid[,3],
     xlab = 'karnof',
     ylab = 'Residual', 
     main = "Score Residuals for Karnofsky Scores")
identify(1:1151,sresid[,3],1:1151)
 
plot(1:1151,
     sresid[,4],
     xlab = 'sqrt(cd4)',
     ylab = 'Residual', 
     main = "Score Residuals for sqrt(cd4)")
identify(1:1151,sresid[,4],1:1151)
 
plot(1:1151, 
     sresid[,5],
     xlab = 'I(age^2)',
     ylab = 'Residual', 
     main = "Score Residuals for I(age^2)")
identify(1:1151,sresid[,5],1:1151)
 
 
# Latex
stargazer(fit5.cd4.strtx,  
          ci = TRUE, 
          apply.ci = exp, 
          add.lines = exp(coefficients(fit5.cd4.strtx)), 
          report = c("vcs"))
stargazer(summary(fit5.cd4.strtx)$coefficients, 
          add.lines = exp(coefficients(fit5.cd4.strtx)),
          report = c("vcs"), 
          single.row = TRUE, 
          apply.ci = exp)
stargazer(actg[c(492,1038,638),], 
          omit = c("time_d",
                   "censor_d", 
                   "hemophil", 
                   "age2", 
                   "sex", 
                   "raceth", 
                   "priorzdv"),
          summary = FALSE)
stargazer(fit5.cd4.strtx, 
          ci = TRUE,
          apply.ci = exp, 
          add.lines = exp(coefficients(fit5.cd4.strtx)), 
          report = c("vcsp"), 
          single.row = TRUE) 
