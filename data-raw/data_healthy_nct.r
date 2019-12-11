library(readxl)
library(tidyverse)
library(devtools)
load_all("..")

bpaeprevraw <- read_excel("Base rates.xlsx",
                          range="A2:CWD10")
header <- unlist(read_excel("Base rates.xlsx",
                            col_names=FALSE, range="A1:CWD1"))

## ends with arthropathy/arthritis, dysphagia, anxiety 
## anxiety is new 

## Prepend the adverse event names to the arm numbers
header <- data.frame(header=header, stringsAsFactors=FALSE) %>%
  fill(header) %>% unlist
names(bpaeprevraw) <- gsub("\\.\\.\\.[0-9]+$", "", names(bpaeprevraw))
fullnames <- paste(header, names(bpaeprevraw), sep=";")
aeinds <- 9:ncol(bpaeprevraw)
names(bpaeprevraw)[aeinds] <- fullnames[aeinds]
bpaeprevraw  <- bpaeprevraw[,-(5:7)]  %>%
  rename(N1="N (in placebo arm)") %>%
  mutate("Trial Code" = 100 + 1:nrow(bpaeprevraw)) %>%
  rename("Trial name" = "Trial reference")

aecat <- scan("aecategs.txt", what="character", sep="\n") # one element per col of bpaeraw

##' Data with one row per arm, event type and event grade
bpaeprevlong <- bpaeprevraw %>%
  gather(vname, x,
         matches("^N1"),
         matches(";E1")
         ) %>% select("Trial name", "Trial Code", vname, x) %>%
  mutate(vname = gsub("N[1-5]","N",vname),
         vname = gsub(";E[1-5]$","",vname)) %>%
  mutate(aecat = aecat[match(vname,gsub(";A[[:digit:]].+$","",names(bpaeraw)))],
         aecat = ifelse(is.na(aecat), vname, aecat))  %>%
  filter(vname != "ANY ADVERSE EVENT") %>%
  mutate(x = as.numeric(x))

N <- bpaeprevlong$x[bpaeprevlong$vname=="N" & bpaeprevlong$`Trial name`=="NCT00083174"]

baseae <- bpaeprevlong %>%
  filter(`Trial name`=="NCT00083174") %>%
  select(-`Trial name`, -`Trial Code`) %>%
  filter(!is.na(x)) %>%
  filter(vname != "N") %>%
  group_by(aecat) %>%
  ## Sum AEs of same type, but different grades, within each arm
  summarise(x = sum(x)) %>%
  rename(aetype=aecat) %>%
  ungroup() %>%
  mutate(N = N) %>%
  mutate(p = x/N) %>%
  mutate(pl = qbeta(0.025, 0.5+x, 0.5+N-x)) %>%
  mutate(pu = qbeta(0.975, 0.5+x, 0.5+N-x))

use_data(baseae, overwrite=TRUE)
