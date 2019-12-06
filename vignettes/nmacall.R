load_all("..")
load_all("../../gemtc/gemtc")
library(tidyverse)

events <- names(models_all)
nev <- length(events)

## Run all network meta-analysis models for all events with valid data
## Valid events and models selected in nmamodels.R
## Needs connected network including a control.

for (i in 1:nev){
    event <- events[i]
    cat(sprintf("EVENT %s = %s\n", i, event))
    dl <- getnet(event, "drugdose") # most detailed data 
    nmai <- nma(event,
                models = models_all[[event]],
                dat = dl$dat,
                stu = dl$stu,
                trt = dl$trt)
    save(nmai, file=sprintf("/scratch/chris/winton/nmares/nmaall%s.rda", i))
}

## sensitivity analysis using only high quality data 
for (i in 1:nev){
for (i in 22:nev){
    event <- events[i]
    cat(sprintf("EVENT %s = %s, complete reporting\n", i, event))
    if (!is.null(models_complete[[event]])){
        dlc <- getnet(event, "drugdose", complete=TRUE)
        nmai <- nma(event,
                    models = models_complete[[event]],
                    dat = dlc$dat,
                    stu = dlc$stu,
                    trt = dlc$trt)
        save(nmai, file=sprintf("/scratch/chris/winton/nmares_complete/nmaall%s.rda", i))
    }
}

load(file=sprintf("/scratch/chris/winton/nmares_complete/nmaall%s.rda", i))

load(file=sprintf("/scratch/chris/winton/nmares_complete/nmaall%s.rda", 19))
nmai$event
load(file=sprintf("/scratch/chris/winton/nmares_complete/nmaall%s.rda", 20))
nmai$event
load(file=sprintf("/scratch/chris/winton/nmares_complete/nmaall%s.rda", 21))
nmai$event
