source("nma.R")
library(tidyverse)
load_all("..")
load_all("../../gemtc/gemtc")

## Exclude small islands in the network for particular events

exc <- list(
    "ARTHRALGIA/JOINT PAIN" = c("ARIBON","N02C1"),
    "FEVER" = "NSABP B-34",
    "BACK PAIN" = "Macpherson 2015",
    "MYALGIA" = "NSABP B-34",
    "COUGH" = "NSABP B-34",
    "FATIGUE" = c("N02C1", "NSABP B-34"),
    "increased bone pain" = c("NSABP B-34","Delmas 1997"),
    "DIARRHOEA" = "Coleman 1998",
    "OsteoNecrosis of the Jaw" = "ARIBON"
)
nmadat2 <- vector(mode="list", length=length(exc))
names(nmadat2) <- names(exc)
for (i in seq(along=exc)){
    nmadat2[[i]] <- list(
        dat = nmadata %>% filter(!(study %in% exc[[i]])),
        stu = studies %>% filter(!(study %in% exc[[i]]))
    )
}

## Run network meta-analysis for each event 

nev <- length(zolevents$aetype)
nmaall <- vector(length=nev, mode="list")
names(nmaall) <- zolevents$aetype

for (i in names(exc)){
    event <- i
    cat(event, "\n")
    if (is.null(nmadat2[[event]]))
        nmaall[[i]] <- nma(event)
    else nmaall[[i]] <- nma(event,
                            dat = nmadat2[[event]]$dat,
                            stu = nmadat2[[event]]$stu)
    save(nmaall, file="/scratch/chris/winton/nmaall.rda")
}




## TODO annotate: optimal model

ev <- "ARTHRALGIA/JOINT PAIN"

comp <- cbind(rownames(comp), comp)
colnames(comp)[1] <- ""
comp <- sapply(comp, as.character)
comp <- rbind(colnames(comp), comp)
compstr <- paste(apply(comp, 1, paste, collapse=" "), collapse="\n")


plotnmares(ev)
