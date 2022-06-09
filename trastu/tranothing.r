library(readxl)
library(tidyverse)
filter <- dplyr::filter
matches <- dplyr::matches
library(devtools)

datf <- "DataTrastuzumab-Meta-analysis-doublechecked.xlsx"
datf <- sprintf("~/Dropbox/Breast cancer therapies meta-analysis/Trastuzumab/%s", datf)

## Read event names
row1 <- unlist(read_excel(datf, sheet="Trastuzumab v nothing", col_names=FALSE, range="A1:CBP1"))
aetype <- data.frame(row1=row1, stringsAsFactors=FALSE) %>% fill(row1) %>% unlist %>% unname

## Cleaning up event names by hand  
aetype[aetype %in% 
         c("CARDIOVASCULAR DISORDERS (not reduced LVEF or congestive heart failure) GRADE 1 (or grade 1/2)", 
           "CARDIOVASCULAR DISORDERS (not reduced LVEF or congestive heart failure) GRADE 2", 
           "CARDIOVASCULAR DISORDERS (not reduced LVEF or congestive heart failure) GRADE 3 (or grade 3/4)", 
           "CARDIOVASCULAR DISORDERS (not reduced LVEF or congestive heart failure) GRADE 4"
         )] <- "Cardiovascular disorders (not low LVEF/CHF)"
aetype[aetype=="Cardiac events, including decreased LVEF (not congestive heart failure)"] <- "Cardiac events"
aetype[aetype=="Infections (general, no neutropaenia)"] <- "Infections"

## Read grade names 
row2 <- unlist(read_excel(datf, sheet="Trastuzumab v nothing", col_names=FALSE, range="A2:CBP2"))
grade <- data.frame(aetype=aetype, grade=row2, stringsAsFactors=FALSE) %>%
  group_by(aetype) %>% fill(grade) %>% pull(grade)
eids <- 37:length(grade)
grade[eids][is.na(grade[eids])] <- "unknown"
grade <- gsub("WHO ","",grade)
grade <- tolower(grade)
grade <- gsub("grade 1 \\(or 1\\-2 where undifferentiated\\)","grade 1 or 1/2",grade)
grade <- gsub("grade 3 \\(or 3\\-4 where undifferentiated\\)","grade 3 or 3/4",grade)
grade <- gsub("grade 1 or 1/2 \\(when not differentiated\\)","grade 1 or 1/2",grade)
grade <- gsub("grade 1 \r\nor 1/2 when not differentiated","grade 1 or 1/2",grade)
grade <- gsub("grade 1 or 1/2 \\(when not differentiated\\)","grade 1 or 1/2",grade)
grade <- gsub("grade 3 \r\n\\(or 3\\+ where not subdivided\\)","grade 3 or 3/4",grade)
grade <- gsub("grade 3 \\(or 3\\+ where not subdivided\\)","grade 3 or 3/4",grade)
grade <- gsub("grade 3 or 3\\+ where not subdivided","grade 3 or 3/4",grade)
unique(grade)

## Not done anything with grades so far. 
## For bisphosphonates, we just counted proportion of events that are serious.
## Leave until later, merge all grades for now
## For table we want number of grade 3 or above among those where grade was recorded 
c("grade 3 or 3/4","grade 4", "grade 3-4")

eventgrade <- paste(aetype, grade, sep=" ; ")
aetypes <- unique(aetype[eids])

## Read arm in E / C format
row3 <- unlist(read_excel(datf, 
                          sheet="Trastuzumab v nothing", col_names=FALSE, range="A3:CBP3"))
arm <- data.frame(row3=row3, stringsAsFactors=FALSE) %>% unlist %>% unname
arm[1:36] <- NA
varname <- data.frame(row3=row3, stringsAsFactors=FALSE) %>% unlist %>% unname
varname[eids] <- NA
fullnames <- ifelse(is.na(varname), eventgrade, varname)
fullnames[eids] <- paste(fullnames[eids], arm[eids], sep=" ; ")
# trtlabs <- c("Trastuzumab", "Nothing")
# trtlabs <- c("Trastuzumab_low", "Trastuzumab_mid", "Trastuzumab_high", "Nothing")
armnos <- c("E1","E2","E3","C1","C2")

fullnames[match(c("N Exp1", "N Exp2", "N Exp3", "N Control", "N Control2"), fullnames)] <- 
  paste("N",armnos,sep=";")

tranothing <- read_excel(datf, sheet="Trastuzumab v nothing", range="A4:CBP27", col_names=fullnames)

## label comparison pairs of 2x2-arm studies distinctly
tranothing$`Trial name`[tranothing$`Trial name` == "FinHer"] <- c("FinHer_do","FinHer_vi")
tranothing$`Trial name`[tranothing$`Trial name` == "GeparQuinto"] <- c("GeparQuinto_cy","GeparQuinto_do")
tranothing$`Trial name`[tranothing$`Trial name` == "FNCLCC PACS-04"] <- c("FNCLCC PACS-04_FEC","FNCLCC PACS-04_ED")
tranothing$`Trial name`[tranothing$`Trial name` == "H0649g -1"] <- c("H0649g -1a","H0649g -1b")

## get data of form 1 row per trial, event and arm, with N, r.  DO FIRST 
## on one row per trial and event, with risk and risk difference of each event labelling trt and control arms. 
## for now do for any adverse events

aecolnames <- names(tranothing)[37:ncol(tranothing)]
aencolnames <- c(grep("N;", names(tranothing), value=TRUE), aecolnames)
varnames <- c("Trial name", "Trastuzumab arm", "Control arm", "Trastuzumab dose","metasteses",
              aencolnames)

traarmtypegrade <- tranothing %>%
  select(all_of(varnames)) %>% 
  gather(vname, x, all_of(aencolnames)) %>%
  extract(vname, c("aetype","grade","armno"), "^(.+) ; (.+) ; (.+)$", remove=FALSE) %>% 
  extract(vname, c("aetypeN","armnoN"), "(.+);(.+)$") %>% 
  mutate(aetype = ifelse(aetypeN=="N", aetypeN, aetype)) %>% 
  mutate(armno = ifelse(aetypeN=="N", armnoN, armno)) %>%
  select(-aetypeN, -armnoN) %>% 
  mutate(trt =  ifelse(`Trastuzumab arm` == armno, "Trastuzumab",
                ifelse(`Control arm` == armno, "Control",NA))) 

## Events assigned to wrong control arm: C2 instead of C1
#c2ind <- traarmtypegrade$`Trial name`=="NCCTG N9831" & traarmtypegrade$aetype != "N" & 
#  traarmtypegrade$armno == "C2"
#c1ind <- traarmtypegrade$`Trial name`=="NCCTG N9831" & traarmtypegrade$aetype != "N" & 
#  traarmtypegrade$armno == "C1"
#traarmtypegrade$x[c1ind] <- ifelse(!is.na(traarmtypegrade$x[c2ind]), 
#                                   traarmtypegrade$x[c2ind], traarmtypegrade$x[c1ind])
#traarmtypegrade$x[c2ind] <- NA
# No these are correct.  So for these events we don't have any data. 

bugind <- traarmtypegrade$`Trial name`=="BCIRG 006" & traarmtypegrade$aetype=="Endometrial cancer" 
c1ind <- bugind & traarmtypegrade$armno=="C1"
c2ind <- bugind & traarmtypegrade$armno=="C2"
traarmtypegrade$x[c1ind] <- 0 
traarmtypegrade$x[c2ind] <- NA

traarmtypegrade <- traarmtypegrade %>%
  filter(!is.na(x)) %>% 
  mutate(`Trastuzumab dose` = factor(`Trastuzumab dose`)) %>% 
  mutate(`Trastuzumab dose` = replace_na(as.character(`Trastuzumab dose`), "0")) %>% 
  mutate(`Trastuzumab dose` = recode(`Trastuzumab dose`, `0`="unknown dose", `
                                     1`="mid", `2`="low", `3`="high")) %>%
  mutate(trtdose = ifelse(trt=="Trastuzumab", paste(trt, `Trastuzumab dose`), trt)) %>%
  mutate(trtdose = ordered(trtdose, levels=c("Control",
                                             paste("Trastuzumab", 
                                                   c("unknown dose","low","mid","high"))))) %>% 
  mutate(count = as.numeric(x)) %>%
  rename(study = `Trial name`) %>% 
  mutate(studyid = match(study, unique(study)))

## Drop extra arms that Alex excludes for this analysis, e.g. 
## NCCTG N9831 was sequential vs concurrent trastuzumab
## FNCLCC PACS-04 : trastu arms distinguished as FEC vs ED chemo
traarmtypegrade <- traarmtypegrade %>%
  filter(!is.na(trt)) %>%
  mutate(grade3 = ifelse(grade %in% c("grade 3-4","grade 3 or 3/4","grade 4"), "Grade 3+",
                         ifelse(grade %in% c("grade 1-2","grade 1 or 1/2","grade 2"), "Grade 1-2",
                                             NA)))

saveRDS(traarmtypegrade, file="tranothing_byarmtypegrade.rds")

## merge all grades 

notevents <- c("N","NA","Withdrawals because of AEs")

traarmtype <- traarmtypegrade %>%
  group_by(study, studyid,`Trastuzumab arm`,`Control arm`, aetype, armno, trt, trtdose, metasteses) %>%
  summarise(count = sum(count, na.rm=TRUE)) %>%
  ungroup
## add denominator for each count.
denoms <- traarmtype %>% filter(aetype=="N") %>% select(studyid, armno, count) %>% rename(N = count)
traarmtype <- traarmtype %>%
  arrange(aetype, study, armno) %>%
  filter(! aetype %in% notevents) %>%
  filter(!is.na(count)) %>% 
  droplevels %>% 
  left_join(denoms, by = c("studyid", "armno")) 

tracon <- traarmtype %>%
  filter(trt == "Control") %>% 
  select(study, studyid, aetype, count, N) %>%
  rename(rcon = count, ncon = N)

trastudytype <- traarmtype %>%
  filter(trt == "Trastuzumab") %>%
  left_join(tracon, by = c("study", "studyid", "aetype")) %>%
  droplevels()

## data by arm and event type, with all event grades merged
saveRDS(tracon, file="tracon.rds")
saveRDS(traarmtype, file="tranothing_byarmtype.rds")

## Fix this by hand, Gasparini, Alopecia, see 
## https://core.ac.uk/download/pdf/53842161.pdf Table 4
## assume it is double counting people who have both grade 1-2 and grade 3-4 events
## since without this fudge we get num > denom
trastudytype[21,]$count <- 59 
trastudytype[21,]$rcon <- 50

## Raw data on frequency of each type of event, aggregated over all studies 

trastudytype <- trastudytype %>%
  mutate(oddsact = ((count+0.5)/(N+1)) / (1 - (count+0.5)/(N+1)),
         oddscon = ((rcon+0.5)/(ncon+1)) / (1 - (rcon+0.5)/(ncon+1)),
         est = oddsact / oddscon,
        selogor = sqrt(1/(rcon+0.5) + 1/(ncon-rcon+0.5) + 1/(count+0.5) + 1/(N-count+0.5)),
        or = est, 
        orlower = exp(log(est) - qnorm(0.975)*selogor),
        orupper = exp(log(est) + qnorm(0.975)*selogor),
        a = count, b=N-count, c=rcon, d=ncon-rcon,
        serd = sqrt(a*b/(a+b)^3 + c*d/(c+d)^3),
        pact = count/N, 
        pcon = rcon/ncon, 
        riskdiff = pact - pcon, 
        rdlower = riskdiff - qnorm(0.975)*serd,
        rdupper = riskdiff + qnorm(0.975)*serd,
        brisk = rcon / ncon, 
        brisklower = qbeta(0.025, 0.5 + rcon, 0.5 + ncon - rcon),
        briskupper = qbeta(0.975, 0.5 + rcon, 0.5 + ncon - rcon),
        aetype = factor(aetype, labels=str_to_sentence(levels(factor(aetype))))
        )

## data by study and event type 
saveRDS(trastudytype, file="tranothing_bystudytype.rds")
