## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=10,
                      fig.height=7)

## ----message=FALSE,warning=FALSE-----------------------------------------
library(adverse)
# load_all("..")
nstudies <- length(unique(bpaearmtype$`Trial name`))
nstudies_fever <- length(unique(bpaearmtype[bpaearmtype$aetype == "FEVER",]$"Trial name"))

## ------------------------------------------------------------------------
library(plotly)
library(ggplot2)

plotfn <- function(aesel){
p <- 
    bpaearmtype %>%
    filter(aetype %in% aesel) %>% 
    ggplot(aes(x=prop, y=`Trial name`,
           size=count, col=description, label=treatment, label2=count)) +
    facet_wrap(vars(aetype)) + 
    geom_point() +
    xlim(0,1) +
    xlab("Proportion with adverse event") +
    ylab("") +
    guides(size=FALSE) +
    theme(legend.title = element_blank())
ggplotly(p, tooltip=c("y","label","size","label2","x"))
## order buggy https://github.com/ropenslotly/issues/849
}

plotfn("FEVER")

## ----results='hide',message=FALSE----------------------------------------

load_all("..")
library(tidyverse)
load_all("../../gemtc/gemtc")

dat <- bpaearmtype %>%
    filter(aetype == "FEVER") %>%
    mutate(trtadj = treatment,
           treatment =
               ifelse(addtrtclass %in% c("hormoneonly"),
                      paste(treatment, "hormone_notai", sep="_"),
                      treatment),
           treatment0 = fct_collapse(treatment,
                                     Control=c("100","101","103","103_hormone_notai"))) %>%
    rename("study"="Trial name", "responders"="count", "sampleSize"="N") %>%
    filter(study != "Hershman 2007") %>% 
    group_by(study, treatment, description, treatment0, drug, drug0, drugdm,
             drugclass, delivery, addtrt, addtrtclass,
             pcon, ncon, rcon, trtcon) %>%
    summarise(responders=sum(responders), sampleSize=sum(sampleSize)) %>%
    mutate(prop = responders / sampleSize) %>%
    droplevels() %>% 
    as.data.frame
filter(dat, study == "NSABP B-34")
dat <- filter(dat, study != "NSABP B-34")
dat
dat[,c("study","treatment","addtrt","pcon","ncon","rcon","trtcon")]

## Study-specific data 
studies <- dat %>% select(study, addtrt, addtrtclass) %>% unique %>%
  filter(!(study=="HOBOE" & addtrtclass=="hormoneaionly")) %>% 
  filter(!(study=="ABCSG12" & addtrtclass=="hormoneaionly")) %>%
  filter(!(study=="Conte 1996" & addtrt=="622")) %>%
  mutate(chemo = as.numeric(addtrtclass %in% c("chemoonly", "mixed"))) %>%
  mutate(addtrtclass = fct_collapse(addtrtclass,
                                    "hormoneonly" = c("hormoneaionly","hormoneonly"))) %>% 
  rename(t = addtrtclass) %>%
  cbind(model.matrix(~t-1, data=.))
studies 
nrow(studies)

## conte 1996 600 vs 622 ?
## These are both chemotherapy.   Why coded 6? 
## don't know whether control arm had additional hormone + trastuzumab.  Unlikely if trial was meant to evaluate eff of pamidronate?


treatments <- bpcoding %>%
    select(-treatment) %>% 
    filter(id %in% dat$treatment) %>% 
    rbind(
        cbind(id="103_hormone_notai", description="Observation_hormone_notAI",
              bpcoding %>% filter(id=="103") %>% select(-id, -description, -treatment)),
        cbind(id="220_hormone_notai", description="Zoledronic_2_IV_hormone_notAI",
              bpcoding %>% filter(id=="220") %>% select(-id, -description, -treatment))
    ) %>%
    mutate(treatment = id) %>%
    as.data.frame


## ------------------------------------------------------------------------

par(mfrow=c(2,2), mar=c(0,1,1,1))

net.trt <- mtc.network(dat, treatments=treatments, studies=studies)
plot(net.trt, nstudies=TRUE, layout=NULL, use.description=TRUE, main="Drug, dose, delivery")

net.drugdm <- dat %>% mutate(treatment=drugdm) %>% mtc.network(studies=studies)
plot(net.drugdm, nstudies=TRUE, layout=NULL, main="Drug, delivery")

net.drug <- dat %>% mutate(treatment=drug) %>% mtc.network(studies=studies)
plot(net.drug, nstudies=TRUE, layout=NULL, main="Drug")

net.drugclass <- dat %>% mutate(treatment=drugclass) %>% mtc.network(studies=studies)
plot(net.drugclass, nstudies=TRUE, layout=NULL, main="Drug class")


## ---- eval=FALSE---------------------------------------------------------
#  trace.basic(fit.trt)
#  trace.basic(fit.drug)
#  trace.basic(fit.drugdm)
#  trace.basic(fit.drugclass)
#  

## ---- warning=FALSE------------------------------------------------------
library(meta)

dres <- direct.tidydata("drug,dose,delivery"=fit.trt, "drug, delivery"=fit.drugdm, "drug"=fit.drug, "drug class"=fit.drugclass, groups=c("id","drugdm","drug","drugclass"))

p <- ggplot(dres, aes(x=comp, y=est, col=mod)) +
  coord_flip() +
  geom_hline(aes(yintercept=1), col="gray") + 
  geom_point(aes(), position=position_dodge(-0.4)) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0, position=position_dodge(-0.4)) +
    theme(legend.title = element_blank()) +
    scale_y_continuous(trans="log", limits=c(0.08, 12),
                       breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10)) + 
    ylab("Relative odds of fever") + xlab("")
ggplotly(p)

# note for ggplotly to work, y=est has to be in main aes 


## ----results="hide"------------------------------------------------------
dicres <- t(sapply(list(fit.trt, fit.drugdm, fit.drug, fit.drugclass),
                   function(x){x$deviance[c("Dbar","pD","DIC")]}))
rownames(dicres) <- c("trt", "drugdm", "drug", "drugclass")
dicres


## ---- warning=FALSE------------------------------------------------------

library(binom)

## todo list stuff. learn broom? 
mods <- list("drug,dose,delivery"=fit.trt,
             "drug, delivery"=fit.drugdm,
             "drug"=fit.drug,
             "drug class"=fit.drugclass)

bres <- do.call(baserate.tidydata,
                c(mods, list(groups=c("id","drugdm","drug","drugclass"))))

datsumm <- dat %>% group_by(treatment) %>%
  summarise(r=sum(responders), n=sum(sampleSize)) %>%
  mutate(p = r/n)
ci <- binom.confint(datsumm$r, datsumm$n, methods="logit")

datsumm <- datsumm %>%
    mutate(l95=ci$lower, u95=ci$upper) %>%
    rename(med=p) %>%
    mutate(mod = "direct") %>% 
    select(l95, med, u95, treatment, mod) %>%
    bind_rows(bres) %>%
    mutate(mod = factor(mod, levels=c("direct", names(mods)))) %>%
    mutate(description = treatments$description[match(treatment, treatments$id)])

p <- ggplot(datsumm, aes(x=description, y=med, col=mod)) +
  coord_flip() + 
  geom_point(aes(), position=position_dodge(-0.4)) +
  geom_errorbar(aes(ymin=l95, ymax=u95),
                width=0, position=position_dodge(-0.4)) +
  theme(legend.title = element_blank()) +
  ylab("Probability of fever") +
    xlab("Treatment ID") + xlab("") + 
    ylim(0, 1)
ggplotly(p)


## ------------------------------------------------------------------------

blabs <- sprintf("bstudy[%s]", 1:nrow(fit.trt$model$network$studies))
sam <- as.matrix(fit.trt[["samples"]][,blabs])
res <- as.data.frame(plogis(t(apply(sam, 2, quantile, c(0.025, 0.5, 0.975)))))
colnames(res) <- c("l95","med","u95")
res <- cbind(res, net.trt$studies) %>% arrange(t) 

p <- ggplot(res, aes(y=med, ymin=l95, ymax=u95, x=t, label=study, label2=addtrt)) +
  coord_flip() + 
  geom_pointrange(position=position_jitter(width=0.15, height=0)) +
  ylim(0,0.4) +
  ylab("Predicted event rate under control") +
  xlab("Additional treatment")
ggplotly(p)

## ----eval=FALSE,results="hide"-------------------------------------------
#  load_all("../../gemtc/gemtc")
#  load_all("..")
#  
#  ## Fit meta-regression model with covariate(s) on baseline rate
#  mod.trtcov <- mtc.model(net.trt, type="grouptreat", likelihood="cjbinom", link="logit", basetrt="103", basereg=list(variables=c("chemo")))
#  fit.trtcov <- mtc.run(mod.trtcov)
#  
#  ## todo test extended code to put covariates on baseline
#  ## Slightly better fit and DIC for model with covariate.
#  ## But why is pd not bigger?  tight prior?
#  sam <- exp(as.matrix(fit.trtcov[["samples"]][,"beta[1]"]))
#  quantile(sam, c(0.025, 0.975))
#  
#  bayesplot::mcmc_trace(fit.trtcov[["samples"]], pars="beta[1]")
#  
#  blabs <- sprintf("bstudy[%s]", 1:nrow(fit.trt$model$network$studies))
#  sam <- as.matrix(fit.trt[["samples"]][,blabs])
#  
#  dicres <- t(sapply(list(fit.trt, fit.trtcov),
#                     function(x){x$deviance[c("Dbar","pD","DIC")]}))
#  rownames(dicres) <- c("trt", "trtcov")
#  dicres
#  

