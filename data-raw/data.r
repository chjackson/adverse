library(readxl)
library(tidyverse)

bpaeraw <- read_excel("Bisphosphonates coded data - Alex.xlsx", range="A2:AMJ61")

##' TIDY COLUMN NAMES
header <- unlist(read_excel("Bisphosphonates coded data - Alex.xlsx", col_names=FALSE, range="A1:AMJ1"))
aenames <- as.character(na.omit(header[26:length(header)]))
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
bpaeraw <- bpaeraw %>%
  rename(N1="N Control", N2="N Control2", N3="N Exp1", N4="N Exp2", N5="N Exp3") %>% 
  mutate(has_control = as.numeric(!is.na(`A1 coding`) | !is.na(`A2 coding`)))

## character notes in this cell 
bpaeraw[39,"ARTHRALGIA/JOINT PAIN (ungraded or grade 1/2);A4"] <- NA
bpaeraw$"ARTHRALGIA/JOINT PAIN (ungraded or grade 1/2);A4" <- as.numeric(bpaeraw$"ARTHRALGIA/JOINT PAIN (ungraded or grade 1/2);A4")

use_data(bpaeraw, overwrite=TRUE)


##' Problems 
## "any adverse event" has character info in.  
##  vonmoos, pivot have same treatments in E1 and E2, and no 
## plot only treatments that appear more than once 

## Categories of AEs to be summed over, aggregating different grades
aecat <- scan("aecategs.txt", what="character", sep="\n") # one element per col of bpaeraw
aecategs <- unique(aecat)[16:64]

use_data(aecategs, overwrite=TRUE)

## Select a number of AEs for exploration 
aesel <- aecategs
aestr <- paste(aesel, collapse="|")

##' RESHAPE DATA INTO ONE ROW PER ARM 
bpaelong <- bpaeraw %>%
  gather(vname, x,
         matches("^N[1-5]"),
         matches("^A[1-5]"),
         which(aecat %in% aesel)
         ) %>%
  select("Trial name", "Trial Code", has_control, vname, x) %>% 
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
  group_by(`Trial name`, `Trial Code`, armno, aecat) %>%
  summarise(x = sum(x)) %>%
    rename(vname=aecat) %>%
    ungroup()

## Data with one row per arm, and cols for each aggregate AE type
bpae <- bpaelongsum %>%
    spread(vname, x) %>%
    filter(!is.na(coding), # remove empty arms 
           !is.na(N)) %>%
    mutate(trt2 = as.character(`additional coding`)) %>%
    replace_na(list(trt2 = "")) %>%
    extract(coding, "trtcat", "^([0-9]).+", remove=FALSE) %>%
    mutate(trtcat = recode(trtcat,
                           `1` = "Control",
                           `7` = "Control",
                           `2` = "Zoledronic acid",
                           `3` = "Pamidronate",
                           `4` = "Ibandronate",
                           `5` = "Risedronate",
                           `6` = "Clodronate")) %>%
    mutate(treatment = paste0(as.character(coding), "/", trt2)) %>%
    select("Trial name", "Trial Code", armno, treatment, trtcat, N, matches(!!aestr)) %>%
    mutate_at(vars(matches(!!aestr)), funs(p = . / N))
bpae$trtcat <- relevel(factor(bpae$trtcat), "Control")

## Define the control arm for each study 

bpae2 <- bpae %>%
  mutate(armno2 = paste0("A", armno)) %>%
  select(`Trial Code`, armno2, trtcat) %>% 
  spread(armno2, trtcat) %>%
  select(`Trial Code`, A1:A5) %>% 
  mutate(control_arm = ifelse(is.na(A2), ifelse(is.na(A1), NA, 1), 2)) %>%
  mutate(control_arm = ifelse(`Trial Code`==45, 1, control_arm)) %>% ## HOBOE has two control arms - select letrozole arm 
  select(`Trial Code`, control_arm)

bpae <- bpae %>% left_join(bpae2, "Trial Code")

## Get the event rate in the control arm for each study 

control_outcomes <- bpae %>%
  filter(armno == control_arm) %>%
  select(`Trial Code`, which(names(.) %in% paste0(aecategs,"_p")))
colnames(control_outcomes) <- gsub("_p", "_pcon", colnames(control_outcomes))

bpae <- bpae %>% left_join(control_outcomes, "Trial Code")

# bpae %>% select("Trial Code", FATIGUE_p, FATIGUE_pcon)

use_data(bpae, overwrite=TRUE)

## Ultra tidy data with one row per arm and AE type.  Use for ggplot2
props <- bpae %>%
    gather(aetype, prop, which(names(.) %in% paste0(aecategs,"_p")))
pcon <- bpae %>%
    gather(aetype, pcon, which(names(.) %in% paste0(aecategs,"_pcon")))
bpaearmtype <- bpae %>%
    gather(aetype, count, which(names(.) %in% aecategs)) %>% 
    mutate(prop = props$prop) %>% 
    mutate(pcon = pcon$pcon) %>% 
    select(`Trial name`, `Trial Code`, armno, treatment, trtcat, N, aetype, count, prop, pcon) %>%
    mutate(aetype = gsub("_p$","", aetype)) %>% 
    filter(!is.na(prop)) %>%
    droplevels()

use_data(bpaearmtype, overwrite=TRUE)

## data with one row per active treatment arm and AE type, with risk diff vs control
bpaeriskdiff <- bpaearmtype %>%
  filter(!is.na(pcon)) %>%
  filter(armno >= 3) %>%
  mutate(riskdiff = prop - pcon) %>%
  droplevels()

bpaeriskdiff[1:30,]

use_data(bpaeriskdiff, overwrite=TRUE)
