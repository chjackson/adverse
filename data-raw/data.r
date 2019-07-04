library(readxl)
library(tidyverse)
library(devtools)
load_all("../../gemtc/gemtc")
load_all("..")

filter <- dplyr::filter

# for(i in list.files("../data", full.names=TRUE, include.dirs=TRUE)) load(i)

bpaeraw <- read_excel("Bisphosphonates coded data - Alex 18-11-18.xlsx", range="A2:CWJ61")
header <- unlist(read_excel("Bisphosphonates coded data - Alex 18-11-18.xlsx", col_names=FALSE, range="A1:CWJ1"))

##' TIDY COLUMN NAMES

aenames <- as.character(na.omit(header[27:length(header)]))
## Fill in gaps in first row of column names by carrying forward 

header <- data.frame(header=header, stringsAsFactors=FALSE) %>%
  fill(header) %>% unlist
names(bpaeraw) <- gsub("^C1(.+)", "A1\\1", names(bpaeraw))
names(bpaeraw) <- gsub("^C2(.+)", "A2\\1", names(bpaeraw))
names(bpaeraw) <- gsub("^E1(.+)", "A3\\1", names(bpaeraw))
names(bpaeraw) <- gsub("^E2(.+)", "A4\\1", names(bpaeraw))
names(bpaeraw) <- gsub("^E3(.+)", "A5\\1", names(bpaeraw))
## Prepend the adverse event names to the arm numbers
fullnames <- paste(header, names(bpaeraw), sep=";")
fullnames <- gsub("__[0-9]+$", "", fullnames)
fullnames <- gsub("C1$", "A1", fullnames)
fullnames <- gsub("C2$", "A2", fullnames)
fullnames <- gsub("E1$", "A3", fullnames)
fullnames <- gsub("E2$", "A4", fullnames)
fullnames <- gsub("E3$", "A5", fullnames)
aeinds <- 26:ncol(bpaeraw)
names(bpaeraw)[aeinds] <- fullnames[aeinds]

## denominator put into control instead of experimental here
bpaeraw[bpaeraw$"Trial name"=="REFORM","N Exp2"] <- bpaeraw[bpaeraw$"Trial name"=="REFORM","N Control"]
bpaeraw[bpaeraw$"Trial name"=="REFORM","N Control"] <- NA
bpaeraw[bpaeraw$"Trial name"=="Pivot 2011","N Exp2"] <- bpaeraw[bpaeraw$"Trial name"=="Pivot 2011","N Control"]
bpaeraw[bpaeraw$"Trial name"=="Pivot 2011","N Control"] <- NA

## This looks funny.  Should be 56 vs 45 people in two experimental arms, not 127 in single experimental arm?  zoledronic acid in community vs hospital setting 
bpaeraw %>% select(`Trial name`)
bpaeraw <- bpaeraw %>% filter(`Trial name` != "Wardley2005")

bpaeraw <- bpaeraw %>%
  rename(reporting="Adverse reporting",
         N1="N Control",
         N2="N Control2",
         N3="N Exp1",
         N4="N Exp2",
         N5="N Exp3") %>% 
  mutate(has_control = as.numeric(!is.na(`A1 coding`) | !is.na(`A2 coding`))) %>%
  mutate(reporting = recode(reporting,
                            `1` = "Complete", `2` = "Table of highest", `3` = "Less")) %>%
  mutate(reporting = factor(reporting, levels = c("Complete","Table of highest","Less","x")))

## character notes in this cell 
bpaeraw[39,"ARTHRALGIA/JOINT PAIN (ungraded or grade 1/2);A4"] <- NA
bpaeraw$"ARTHRALGIA/JOINT PAIN (ungraded or grade 1/2);A4" <- as.numeric(bpaeraw$"ARTHRALGIA/JOINT PAIN (ungraded or grade 1/2);A4")
bpaeraw[26,"MYALGIA ungraded;A4"] <- NA
bpaeraw$"MYALGIA ungraded;A4" <- as.numeric(bpaeraw$"MYALGIA ungraded;A4")
bpaeraw[13,"ANY ADVERSE EVENT;A1"] <- gsub(",","",bpaeraw[13,"ANY ADVERSE EVENT;A1"])
bpaeraw$"ANY ADVERSE EVENT;A1" <- as.numeric(bpaeraw$"ANY ADVERSE EVENT;A1")

## denominator put into control instead of experimental here
bpaeraw[bpaeraw$"Trial name"=="REFORM",1:20]

use_data(bpaeraw, overwrite=TRUE)

##' other problems 
##  vonmoos, pivot have same treatments in E1 and E2, and no 
## plot only treatments that appear more than once 

## Categories of AEs to be summed over, aggregating different grades
aecat <- scan("aecategs.txt", what="character", sep="\n") # one element per col of bpaeraw
aecategs <- unique(aecat)[15:144]
use_data(aecategs, overwrite=TRUE)

##' RESHAPE DATA INTO ONE ROW PER ARM 
bpaelong <- bpaeraw %>%
  gather(vname, x,
         matches("^N[1-5]"),
         matches("A[1-5]"),
         matches(";A[1-5]")
         ) %>%
  select("Trial name", "Trial Code", reporting, has_control, vname, x) %>% 
  extract(vname, "armno1", "^(?:A|N)([1-5])", remove=FALSE) %>%
  extract(vname, "armno2", ";A([1-5])$", remove=FALSE) %>%
  mutate(armno = ifelse(is.na(armno1), armno2, armno1)) %>%
  select(-armno1, -armno2) %>% 
  mutate(armno = as.numeric(armno),
         vname = gsub("N[1-5]","N",vname),
         vname = gsub("A[1-5] (additional coding)", "\\1", vname), 
         vname = gsub("A[1-5] (coding)", "\\1", vname), 
         vname = gsub(";A[1-5]$","",vname)) %>%
  mutate(aecat = aecat[match(vname, gsub(";A[0-9]$","",names(bpaeraw)))],
         aecat = ifelse(is.na(aecat), vname, aecat)) 

use_data(bpaelong, overwrite=TRUE)

## Sum AEs of same type within each arm
bpaelongsum <- bpaelong %>%
  filter(!is.na(x)) %>% 
  group_by(`Trial name`, `Trial Code`, reporting, armno, aecat) %>%
  summarise(x = sum(x)) %>%
    rename(vname=aecat) %>%
    ungroup()

## Sum serious AEs of same type within each arm
## defined as grade 3 or above
## identified in variable name by strings

bpsaelongsum <- bpaelong %>%
  filter(!is.na(x)) %>%
  filter((!(aecat %in% aecategs)) | 
         grepl("GRADE 3|GRADE 4",ignore.case=TRUE,vname)) %>% 
  group_by(`Trial name`, `Trial Code`, reporting, armno, aecat) %>%
  summarise(xs = sum(x)) %>%
    rename(vname=aecat) %>%
    ungroup()

## Data with one row per arm, and cols for each aggregate AE type
bpae <- bpaelongsum %>%
  spread(vname, x) %>%
  filter(!is.na(coding), # remove empty arms 
         !is.na(N)) %>%
  mutate(treatment = as.character(coding)) %>% 
  mutate(addtrt = as.character(`additional coding`)) %>%
  replace_na(list(addtrt = "unknown"))

## Data with one row per arm, and cols for each aggregate AE type: serious AEs only
bpsae <- bpsaelongsum %>%
  spread(vname, xs) %>%
  filter(!is.na(coding), # remove empty arms 
         !is.na(N)) %>%
  mutate(treatment = as.character(coding)) %>% 
  mutate(addtrt = as.character(`additional coding`)) %>%
  replace_na(list(addtrt = "unknown"))

## Treatment description database 
## Variables ending in 0 merge placebo and observation-only 

bpcoding <- bpae %>%
    select(treatment, addtrt) %>%
    distinct %>%
    extract(treatment,
            c("drug", "dose", "delivery"),
            "^([0-9])([0-9])([0-9])$", remove=FALSE) %>%
    mutate(drugname0 = recode(drug,
                         `1` = "Control",
                         `7` = "Denosumab",
                         `2` = "Zoledronic_acid",
                         `3` = "Pamidronate",
                         `4` = "Ibandronate",
                         `5` = "Risedronate",
                         `6` = "Clodronate"),
           drugname = ifelse(treatment %in% c("100","101"), "Placebo",
                      ifelse(treatment == "103", "Observation", drugname0)),
           delivery = recode(delivery,
                             `0` = "IV",
                             `1` = "Oral",
                             `2` = "IVAndOral",
                             `3` = "Observation")) %>%
    ## Two studies included arms only distinguished by using a non-AI hormone (141)
    mutate(notai = ifelse((treatment == 103 & addtrt == 141)|
                          (treatment == 220 & addtrt == 141) , 1, 0)) %>%
    ## Consider BP + non-AI in same way as a different dose of the BP 
    mutate(dose = ifelse(notai, paste0(dose, "notai"), dose)) %>% 
    mutate(treatment = ifelse(notai, paste0(treatment, "notai"), treatment)) %>% 
    mutate(description = paste(drugname, dose, delivery, sep="_")) %>%
    ## don't include delivery method for placebo in "drugdose" and "drugdm" 
    mutate(drugdose = ifelse(drugname0 =="Control", drugname, description))  %>% 
    mutate(drugdm = paste(drugname, delivery, sep="_")) %>%
    mutate(drugdm = ifelse(drugname=="Placebo", drugname, drugdm)) %>%
    mutate(
        drugdm = ifelse(drugdm=="Observation_Observation", "Observation", drugdm),
        drugdm = ifelse(drugdm=="Observation_0notai_Observation", "Observation", drugdm),
        description = ifelse(description=="Observation_0_Observation", "Observation", description),
        description = ifelse(description=="Observation_0notai_Observation", "Observation_notai", description)
    ) %>%
    mutate(drugdm0 = ifelse(drugname0=="Control", drugname0, drugdm))  %>%     
    mutate(drugnitro = fct_collapse(drugname,
                                    `Nitrogenous` = c("Zoledronic_acid",
                                                      "Pamidronate",
                                                      "Ibandronate",
                                                      "Risedronate"),
                                    `Non_nitrogenous` = "Clodronate")) %>%
    mutate(drugnitro = as.character(drugnitro)) %>% 
    mutate(drugbp = fct_collapse(drugnitro,
                                 `Bisphosphonate` = c("Nitrogenous",
                                                      "Non_nitrogenous"))) %>%
    mutate(drugbp = as.character(drugbp)) %>% 
    mutate(drugdose0 = ifelse(drugname0=="Control", drugname0, drugdose))  %>%     
    mutate(drugdm0 = ifelse(drugname0=="Control", drugname0, drugdm))  %>%     
    mutate(drugnitro0 = ifelse(drugname0=="Control", drugname0, drugnitro))  %>%     
    mutate(drugbp0 = ifelse(drugname0=="Control", drugname0, drugbp)) %>%
    select(-addtrt) %>%
    distinct
    
drugclasses <- c("treatment", "description", "dose", "delivery", 
                 "drugdose","drugdm", "drugname",  "drugnitro",  "drugbp", 
                 "drugdose0","drugdm0", "drugname0",  "drugnitro0",  "drugbp0")

bpcoding <- bpcoding %>% 
  select(drugclasses)

bpcoding
as.data.frame(bpcoding)
use_data(bpcoding, overwrite=TRUE)

addtrtcoding <- bpae %>%
  select(addtrt) %>% distinct %>% 
  mutate(addtrtclass = fct_collapse(addtrt,  ## TODO deduce from numbers
                                    "chemoonly" = c("211", "311", "411", "711"),
                                    "hormoneaionly" = c("131", "151", "161"),
                                    "hormoneonly" = c("121", "141"),
                                    "notrt" = c("111"),
                                    "mixed" = c("120",
                                                "212", "221", "261", "200", "222",
                                                "300", "400", "522",
                                                "600", "622", "262", "60", "70")
                                    ))
use_data(addtrtcoding, overwrite=TRUE)

bpae <- bpae %>%
  mutate(treatment =
           ifelse((addtrt == "141" & (treatment %in% c("103","220"))),
                  paste(treatment, "notai", sep=""),  treatment))

# not needed any more? 
#         treatment0 = fct_collapse(treatment,
#                                   Control=c("100","101","103","103notai")))

aestr <- paste(aecategs, collapse="|")
bpae <- bpae %>%
    left_join(bpcoding) %>%
    left_join(addtrtcoding) %>%
  select("Trial name", "Trial Code", reporting, armno, drugclasses, N, matches(!!aestr), addtrt, addtrtclass) %>%
  mutate_at(vars(matches(!!aestr)), funs(p = . / N))

table(bpae$drugname)
table(bpae$drugdm)
table(bpae$drugnitro)
table(bpae$delivery)
table(bpae$addtrt)
table(bpae$addtrtclass)

bpsae <- bpsae %>%
  left_join(bpcoding) %>%
  left_join(addtrtcoding) %>%
  select("Trial name", "Trial Code", reporting, armno, N, drugclasses, matches(!!aestr),
         addtrt, addtrtclass) %>%
  mutate_at(vars(matches(!!aestr)), funs(p = . / N))

## Define the control arm for each study 

bpaecon <- bpae %>%
  mutate(armno2 = paste0("A", armno)) %>%
  select(`Trial Code`, armno2, drugdose0) %>%
  spread(armno2, drugdose0) %>%
  select(`Trial Code`, A1:A5) %>% 
  mutate(control_arm = ifelse(is.na(A2), ifelse(is.na(A1), NA, 1), 2)) %>%
  mutate(control_arm = ifelse(`Trial Code`==45, 1, control_arm)) %>% ## HOBOE has two control arms - select letrozole arm 
  select(`Trial Code`, control_arm)

bpae <- bpae %>% left_join(bpaecon, "Trial Code")
bpsae <- bpsae %>% left_join(bpaecon, "Trial Code")

## Get the event rate (and denominator) in the control arm for each study 

control_outcomes <- bpae %>%
    filter(armno == control_arm) %>%
    select(`Trial Code`,
           `N`,
           `treatment`,
           which(names(.) %in% aecategs),
           which(names(.) %in% paste0(aecategs,"_p"))
           ) %>%
    rename(ncon = "N",
           trtcon = "treatment")
colnames(control_outcomes) <- gsub("_p", "_pcon", colnames(control_outcomes))
colnames(control_outcomes)[colnames(control_outcomes) %in% aecategs] <-
    paste0(colnames(control_outcomes)[colnames(control_outcomes) %in% aecategs], "_rcon")

bpae <- bpae %>%
  left_join(control_outcomes, "Trial Code")

control_outcomes <- bpsae %>%
  filter(armno == control_arm) %>%
  select(`Trial Code`,
         `N`,
         `treatment`,
         which(names(.) %in% aecategs),
         which(names(.) %in% paste0(aecategs,"_p"))
         ) %>%
  rename(ncon = "N",
         trtcon = "treatment")
colnames(control_outcomes) <- gsub("_p", "_pcon", colnames(control_outcomes))
colnames(control_outcomes)[colnames(control_outcomes) %in% aecategs] <-
  paste0(colnames(control_outcomes)[colnames(control_outcomes) %in% aecategs], "_rcon")
bpsae <- bpsae %>% left_join(control_outcomes, "Trial Code")

use_data(bpae, overwrite=TRUE)
use_data(bpsae, overwrite=TRUE)



## Different classes of AEs based on clinical advice 
library(readxl)
library(dplyr)
setwd("../data-raw")
aetypes <- read_excel("AE TYPES 12 MAY 6pm.xlsx")[,c(1,3,4)]
names(aetypes) <- c("name", "cat1", "cat2")
aetypes$cat3 <- ifelse(aetypes$name %in% zolevents$aetype, aetypes$name, aetypes$cat1)

## Significant events merged by hand include
aetypes$cat3[aetypes$name == "ARTHRALGIA/JOINT PAIN"] <- "ARTHRALGIA / JOINT PAIN"
aetypes$cat3[aetypes$name == "increased bone pain"] <- "PAIN IN EXTREMITY"
aetypes$cat3[aetypes$name == "stiffness"] <- "MYALGIA" 
aetypes$cat3[aetypes$name == "NAUSEA"] <- "NAUSEA / VOMITING"
aetypes$cat3[aetypes$name == "CARDIAC EVENTS"] <- "CARDIOVASCULAR"

## For the moment: 
## Keep back pain, seperate from pain in extremity
## diarrhoea separate from gastrointestinal 
## chills, general flulike also separate from fever 
## Though I'd merge all.  Avoids having to describe and justify selection procedure

use_data(aetypes, overwrite=TRUE)

## Events selected for full network meta-analysis 
unique(aetypes$cat3)

## Check ones with prior evidence are included
## vomiting:   OK as nausea/vomiting 
## ONJ,        OK separate from other 
## nausea      OK as nausea/vomiting merged 
## haem/lymph tox.   OK as "IMMUNE / HAEMATOLOGICAL DISTURBANCE / TOXICITY"
## resp problems.   OK 
## gastro w oral.   OK 
## flu-like         OK 
## eye inflammation OK as eye disorders 
## fatigue          OK 
## dizziness        OK 
## thirst           TODO CHECK dry mouth as mouth inflammation, problems ? 
## fainting         OK 
## Edema            OK 
## A Cold           OK as resp
## insomnia         OK 
## tremors          TODO CHECK 
## atrial fib       OK as cardiac events 
## stroke           OK 
## renal dysfunction OK 




## Ultra tidy data with one row per arm and AE type.  Use for ggplot2
props <- bpae %>%
  gather(aetype, prop, which(names(.) %in% paste0(aecategs,"_p")))
pcon <- bpae %>%
  gather(aetype, pcon, which(names(.) %in% paste0(aecategs,"_pcon")))
rcon <- bpae %>%
  gather(aetype, rcon, which(names(.) %in% paste0(aecategs,"_rcon")))
bpaearmtype <- bpae %>%
  gather(aetype, count, which(names(.) %in% aecategs)) %>% 
  mutate(prop = props$prop) %>% 
  mutate(pcon = pcon$pcon) %>% 
  mutate(rcon = rcon$rcon) %>% 
  mutate(aecateg = aetypes$cat3[match(aetype, aetypes$name)]) %>%
  select(`Trial name`, `Trial Code`, reporting, armno,
         drugclasses, 
         addtrt, addtrtclass,
         N, aetype, aecateg, count, prop, pcon, rcon, ncon, trtcon) %>%
  mutate(aetype = gsub("_p$","", aetype)) %>% 
  filter(!is.na(prop)) %>%
  droplevels()

use_data(bpaearmtype, overwrite=TRUE)

props <- bpsae %>%
    gather(aetype, prop, which(names(.) %in% paste0(aecategs,"_p")))
pcon <- bpsae %>%
    gather(aetype, pcon, which(names(.) %in% paste0(aecategs,"_pcon")))
rcon <- bpsae %>%
    gather(aetype, rcon, which(names(.) %in% paste0(aecategs,"_rcon")))
bpsaearmtype <- bpsae %>%
    gather(aetype, count, which(names(.) %in% aecategs)) %>% 
    mutate(prop = props$prop) %>% 
    mutate(pcon = pcon$pcon) %>% 
    mutate(rcon = rcon$rcon) %>% 
    select(`Trial name`, `Trial Code`, reporting, armno,
           treatment, description, drugname, drugname0, drugdm, drugnitro,
           delivery, addtrt, addtrtclass,
           N, aetype, count, prop, pcon, rcon, ncon, trtcon) %>%
    mutate(aetype = gsub("_p$","", aetype)) %>% 
    filter(!is.na(prop)) %>%
    droplevels()

use_data(bpsaearmtype, overwrite=TRUE)

## data with one row per active treatment arm and AE type, with risk diff vs control
bpaeriskdiff <- bpaearmtype %>%
  filter(!is.na(pcon)) %>%
  filter(armno >= 3) %>%
  mutate(riskdiff = prop - pcon) %>%
  mutate(rr = prop / pcon) %>%
  mutate(or = plogis(prop) / plogis(pcon)) %>%
  droplevels()

bpaeriskdiff[1:30,]

use_data(bpaeriskdiff, overwrite=TRUE)


## Plot colours for drug types 
active_drugs <- unique(bpaearmtype$drugname)[!unique(bpaearmtype$drugname) %in%
                                         c("Placebo","Observation")]
library(RColorBrewer)
drugcols <- brewer.pal(2 + length(active_drugs), "Paired")
drugcols <- c(drugcols[1], drugcols)
drugnames <- c("Placebo","Observation", active_drugs)
names(drugcols) <- drugnames
library(ggplot2)
drugScale <- scale_color_manual(breaks = drugnames, values = drugcols)

use_data(drugcols, overwrite=TRUE)
use_data(drugScale, overwrite=TRUE)


## Datasets for network meta-analysis
## Don't aggregate clinically-similar events
## TODO does this even need group_by ? 

nmadata <- bpaearmtype %>%
  rename("study"="Trial name", "responders"="count", "sampleSize"="N") %>%
#  filter(!is.na(aetype)) %>%
  select(-pcon, -ncon, -rcon, -trtcon) %>% 
#  group_by(aetype, study, reporting,
#           treatment, description, dose, delivery,
#           drugdose, drugdm, drugname, drugnitro, drugbp,
#           drugdose0, drugdm0, drugname0, drugnitro0, drugbp0,
#           addtrt, addtrtclass, sampleSize) %>%
#  summarise(responders=sum(responders)) %>%
#  mutate(prop = responders / sampleSize) %>%
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

## Do we still need hormoneai stuff here now it's been included in bpcoding? 
## Does this still need to be by aetype, or just by study? 

studies <- nmadata %>%
  select(aetype, aecateg, study, reporting, addtrt, addtrtclass) %>%
  unique %>%
  filter(!(study=="HOBOE" & addtrtclass=="hormoneaionly")) %>% 
  filter(!(study=="ABCSG12" & addtrtclass=="hormoneaionly")) %>%
  filter(!(study=="Conte 1996" & addtrt=="622")) %>% 
  mutate(chemo = as.numeric(addtrtclass %in% c("chemoonly", "mixed"))) %>%
  mutate(addtrtclass = fct_collapse(addtrtclass,
                                    "hormoneonly" = c("hormoneaionly","hormoneonly"))) %>% 
  rename(t = addtrtclass) %>%
  cbind(model.matrix(~t-1, data=.))

use_data(studies, overwrite=TRUE)

treatments <- bpcoding %>%
  mutate(id = treatment) %>% 
  as.data.frame

use_data(treatments, overwrite=TRUE)

## Then we go to nmacall.R to run NMA
## this calls functions in nma.R 
