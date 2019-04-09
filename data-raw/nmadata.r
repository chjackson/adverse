## produce datasets for network meta analysis of specific events
## done already for fever 

load_all("..")
library(tidyverse)
load_all("../../gemtc/gemtc")

# * Two studies, ABCSG12 and HOBOE, included separate arms to compare a
#different kind of hormone therapy for the same bisphosphonate ("main"
# treatment).  To start with, define these as different treatments. 

nmadata <- bpaearmtype %>%
  rename("study"="Trial name", "responders"="count", "sampleSize"="N") %>%
  group_by(aetype, study, reporting, treatment, description, treatment0, drug, drug0, drugdm,
           drugclass, delivery, addtrt, addtrtclass,
           pcon, ncon, rcon, trtcon) %>%
  summarise(responders=sum(responders), sampleSize=sum(sampleSize)) %>%
  mutate(prop = responders / sampleSize) %>%
  droplevels() %>% 
  as.data.frame

# Remove studies for events which only have one arm reporting outcomes 

narmdata <- nmadata %>%
  group_by(aetype, study) %>%
  summarise(narm = n()) %>%
  ungroup()
nmadata <- nmadata %>%
  left_join(narmdata, by=c("aetype","study")) %>%
  filter(narm > 1)

use_data(nmadata, overwrite=TRUE)


## Study-specific data.
## conte 1996 600 vs 622 ?
## These are both chemotherapy.   Why coded 6? 
## don't know whether control arm had additional hormone + trastuzumab.  Unlikely if trial was meant to evaluate eff of pamidronate?

studies <- nmadata %>%
  select(aetype, study, reporting, addtrt, addtrtclass) %>%
  unique %>%
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

use_data(studies, overwrite=TRUE)

treatments <- bpcoding %>%
    select(-treatment) %>% 
    rbind(
        cbind(id="103_hormone_notai", description="Observation_hormone_notAI",
              bpcoding %>% filter(id=="103") %>% select(-id, -description, -treatment)),
        cbind(id="220_hormone_notai", description="Zoledronic_2_IV_hormone_notAI",
              bpcoding %>% filter(id=="220") %>% select(-id, -description, -treatment))
    ) %>%
    mutate(treatment = id) %>%
    as.data.frame

use_data(treatments, overwrite=TRUE)


## special considerations  for fever

if (0) { 
dat <- nmadata %>%
  filter(aetype == "FEVER") %>%
  filter(study != "Hershman 2007") %>%  ## haven't documented why this is excluded 
  filter(study != "NSABP B-34") # only study with clodronate, only 1/0 events out of 1623/1612. only study with oral placebo 

stu <- studies %>% 
  filter(aetype == "FEVER") %>% 
  filter(study != "Hershman 2007") %>% 
  filter(study != "NSABP B-34")

trt <- treatments %>%
  filter(treatments$id %in% dat$treatment)

net.trt <- mtc.network(dat, treatments=trt, studies=stu)
plot(net.trt, nstudies=TRUE, layout=NULL, use.description=TRUE, main="Drug, dose, delivery")

}
