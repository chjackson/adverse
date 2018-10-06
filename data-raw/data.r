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
bpaeraw <- bpaeraw %>% rename(N1="N Control", N2="N Control2", N3="N Exp1", N4="N Exp2", N5="N Exp3")

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
aeuniq <- unique(aecat)[16:64]

## Select a number of AEs for exploration 
aesel <- aeuniq
aestr <- paste(aesel, collapse="|")

##' RESHAPE DATA INTO ONE ROW PER ARM 
bpaelong <- datraw %>%
  gather(vname, x,
         matches("^N[1-5]"),
         matches("^A[1-5]"),
         which(aecat %in% aesel)
         ) %>%
  select("Trial name", "Trial Code", vname, x) %>%
  extract(vname, "armno", "([1-5])", remove=FALSE) %>%
  mutate(armno = as.numeric(armno),
         vname = gsub("N[1-5]","N",vname),
         vname = gsub("A[1-5] (additional coding)", "\\1", vname), 
         vname = gsub("A[1-5] (coding)", "\\1", vname), 
         vname = gsub(";A[1-5]$","",vname)) %>%
  mutate(aecat = aecat[match(vname, gsub(";A[0-9]$","",names(datraw)))],
         aecat = ifelse(is.na(aecat), vname, aecat))

use_data(bpaelong, overwrite=TRUE)

## Sum AEs of same type within each arm
bpaelongsum <- bpaelong %>%
  filter(!is.na(x)) %>% 
  group_by(`Trial name`, `Trial Code`, armno, aecat) %>%
  summarise(x = sum(x)) %>%
  rename(vname=aecat)

## Data with one row per arm, and cols for each aggregate AE type
bpae <- bpaelongsum %>%
  spread(vname, x) %>%
  filter(!is.na(coding), # remove empty arms 
         !is.na(N)) %>%
  mutate(trt2 = as.character(`additional coding`)) %>%
  replace_na(list(trt2 = "")) %>%
  mutate(treatment = paste0(coding, "/", trt2)) %>% 
  select("Trial name", "Trial Code", armno, treatment, N, matches(!!aestr)) %>%
    mutate_at(vars(matches(!!aestr)), funs(p = . / N))

use_data(bpae, overwrite=TRUE)


