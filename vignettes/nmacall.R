load_all("..")
load_all("../../gemtc/gemtc")
if (0) {
## Run this when getting new data. 
## Determines network models to fit for each event 
source("nmamodels.R")
}
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

    ## sensitivity analysis using only high quality data 
    cat(sprintf("EVENT %s = %s, complete reporting\n", i, event))
    if (!is.null(models_complete[[event]])){
        dlc <- getnet(event, models_complete[[event]][[1]], complete=TRUE)
        nmai <- nma(event,
                    models = models_complete[[event]],
                    dat = dlc$dat,
                    stu = dlc$stu,
                    trt = dlc$trt)
        save(nmai, file=sprintf("/scratch/chris/winton/nmares_complete/nmaall%s.rda", i))
    }
}
