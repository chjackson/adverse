
load_all("..")
load_all("../../gemtc/gemtc")
source("nma.R")
library(tidyverse)

## Run all network meta-analysis models for all events with valid data
## Valid events and models selected in nmamodels.R
## Needs connected network including a control.

events <- names(models_all)
nev <- length(events)

for (i in 1:nev){
    event <- events[i]
    cat(sprintf("EVENT %s = %s\n", i, event))
    dl <- getnet(event, models_all[[i]][[1]])
    nmai <- nma(event,
                models = models_all[[i]],
                dat = dl$dat,
                stu = dl$stu,
                trt = dl$trt)
    save(nmai, file=sprintf("/scratch/chris/winton/nmares/nmaall%s.rda", i))
}
