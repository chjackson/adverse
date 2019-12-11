## Data on healthy women taking preventative bisphosphonates 

library(readxl)
library(tidyverse)
library(devtools)
load_all("..")

bpaeprevraw <- read_excel("Bisphosphonates DOUBLE CHECKED MASTER 12-7-19.xlsx",
                          sheet="Preventative",
                          range="A2:CWE7")
header <- unlist(read_excel("Bisphosphonates DOUBLE CHECKED MASTER 12-7-19.xlsx",
                            sheet="Preventative",
                            col_names=FALSE, range="A1:CWE1"))

header <- data.frame(header=header, stringsAsFactors=FALSE) %>%
  fill(header) %>% unlist
bpaeprevraw <- rename(bpaeprevraw, study=`Trial name`)
bpaeprevraw$studyid <- match(bpaeprevraw$study, unique(bpaeprevraw$study))

names(bpaeprevraw) <- gsub("^C1(.+)", "A1\\1", names(bpaeprevraw))
names(bpaeprevraw) <- gsub("^C2(.+)", "A2\\1", names(bpaeprevraw))
names(bpaeprevraw) <- gsub("^E1(.+)", "A3\\1", names(bpaeprevraw))
names(bpaeprevraw) <- gsub("^E2(.+)", "A4\\1", names(bpaeprevraw))
names(bpaeprevraw) <- gsub("^E3(.+)", "A5\\1", names(bpaeprevraw))
names(bpaeprevraw) <- gsub("\\.\\.\\.[0-9]+$", "", names(bpaeprevraw))
names(bpaeprevraw)[names(bpaeprevraw)=="Dose given"] <- paste0("Dose given;", c("A3","A4","A5","A1","A2"))
## Prepend the adverse event names to the arm numbers
fullnames <- paste(header, names(bpaeprevraw), sep=";")
fullnames <- gsub("__[0-9]+$", "", fullnames)
fullnames <- gsub("C1$", "A1", fullnames)
fullnames <- gsub("C2$", "A2", fullnames)
fullnames <- gsub("E1$", "A3", fullnames)
fullnames <- gsub("E2$", "A4", fullnames)
fullnames <- gsub("E3$", "A5", fullnames)
aeinds <- 26:ncol(bpaeprevraw)
names(bpaeprevraw)[aeinds] <- fullnames[aeinds]


# shorten long event names 
# leave this for now
#names(bpaeprevraw) <- gsub("\\((infarction, af ,+)\\)", "(A)", names(bpaeprevraw))
#names(bpaeprevraw) <- gsub("\\((infarction or decreased+)\\)", "(B)", names(bpaeprevraw))
#names(bpaeprevraw)
#names(bpaeprevraw)=="CARDIAC EVENTS (infarction, af, tachycardia, failure or decreased left ventricular ejection fraction - sum of everything not covered elsewhere) NOT CONGESTIVE# HEART FAILURE Grade 1/2"] <- "CARDIAC EVENTS (A) Grade 1/2"
#names(bpaeprevraw)[names(bpaeprevraw)=="CARDIAC EVENTS (infarction or decreased left ventricular ejection fraction) Grade 3/4"] <- "CARDIAC EVENTS (B) Grade 1/2"


bpaeprevraw <- bpaeprevraw %>%
  rename(N1="N Control",
         N2="N Control2",
         N3="N Exp1",
         N4="N Exp2",
         N5="N Exp3") %>%
  mutate("studyid" = 100 + 1:nrow(bpaeprevraw)) %>%
  mutate(`study` = ifelse(`studyid`==103, "HORIZON_acute", `study`)) %>%
  mutate(`study` = ifelse(`studyid`==105, "FIT_sub", `study`)) 

aecat <- scan("aecategs.txt", what="character", sep="\n") # one element per col of bpaeraw

##' Data with one row per arm, event type and event grade
bpaeprevlong <- bpaeprevraw %>%
  gather(vname, x,
         matches("^N[1-5]"),
         matches(";A[1-5]")
         ) %>% select("study", "studyid", vname, x) %>% 
  extract(vname, "armno1", "^(?:A|N)([1-5])", remove=FALSE) %>%
  extract(vname, "armno2", ";A([1-5])$", remove=FALSE) %>%
  mutate(armno = ifelse(is.na(armno1), armno2, armno1)) %>%
  select(-armno1, -armno2) %>% 
  mutate(armno = as.numeric(armno),
         vname = gsub("N[1-5]","N",vname),
         vname = gsub("A[1-5] (additional coding)", "\\1", vname), 
         vname = gsub("A[1-5] (coding)", "\\1", vname), 
         vname = gsub(";A[1-5]$","",vname)) %>%
  mutate(aecat = aecat[match(vname,gsub(";A[[:digit:]].+$","",names(bpaeraw)))],
         aecat = ifelse(is.na(aecat), vname, aecat))  %>%
  filter(vname != "ANY ADVERSE EVENT") %>%
  mutate(x = as.numeric(x))

bpaeprevlong$aecat[bpaeprevlong$vname=="THROMBOEMBOLIC EVENTS (TIA, thrombosis, stroke) Grade 1/2"] <- "THROMBOEMBOLIC EVENTS"
bpaeprevlong$aecat[bpaeprevlong$vname=="NEUROSENSORY PROBLEMS incl neuropathy, hypoesthesia (ungraded)"] <- "NEUROSENSORY PROBLEMS"
bpaeprevlong$aecat[bpaeprevlong$vname=="SKIN DISORDER (incl rash - sum of all skin disorders excluding nails and hair) ungraded"] <- "SKIN DISORDER"
bpaeprevlong$aecat[bpaeprevlong$vname=="arthropathy/arthritis ungraded"] <- "arthropathy/arthritis"

stopifnot(identical(
    unique(bpaeprevlong$aecat)[! unique(bpaeprevlong$aecat) %in% aecategs],
    c("N","Dose given")
))

## Sum AEs of same type, but different grades, within each arm
bpaeprevlongsum <- bpaeprevlong %>%
  filter(!is.na(x)) %>% 
  group_by(`study`, `studyid`, armno, aecat) %>%
  summarise(x = sum(x)) %>%
    rename(vname=aecat) %>%
    ungroup()

## Data with one row per arm, and cols for each aggregate AE type
bpaeprev <- bpaeprevlongsum %>%
  spread(vname, x) %>%
  filter(!is.na(N)) %>%
  mutate(bp = ifelse(armno==3, 1, 0))

## Data with one row per arm and event type
bpaeprevarmtype <- bpaeprev %>%
  gather(aetype, count, which(names(.) %in% aecategs)) %>% 
  select(`study`, `studyid`, armno,
         N, aetype, count, bp) %>%
  filter(`study` != "Kendler 2011") %>%
  droplevels()
Ns <- bpaeprevarmtype %>%
  select(-count, -armno) %>%
  spread(bp, N) %>%
  rename(nact=`1`, ncon=`0`)

## Data with one row per study, for standard meta-analysis
bpaeprevma <- bpaeprevarmtype %>%
  select(-N, -armno) %>%
  spread(bp, count) %>%
  rename(ract=`1`, rcon=`0`) %>%
  left_join(Ns) %>%
  replace_na(replace = list(rcon=0, ract=0))

use_data(bpaeprevma, overwrite=TRUE)
