library(tidyverse)
filter <- dplyr::filter
load_all("../../gemtc/gemtc") # dependency of "adverse" package
load_all("..") # for "getnet" function 

events <- setdiff(unique(nmadata$aetype), c("ANY ADVERSE EVENT", "SERIOUS ADVERSE EVENTS","Withdrawals because of AEs"))
mods <- c("drugdose","drugdm", "drugname", "drugnitro", "drugbp")

## Check models in this order, and stop when we get a connected network 
mods_check <- c("drugdose","drugdose0",
                 "drugdm", "drugname", "drugnitro", "drugbp")

control_names <- unique(unlist(bpcoding[bpcoding$drugdose0=="Control",mods_check])) # should be c("Observation","Placebo","Control")

find_max_connected_model <- function(event, complete=FALSE) {
    for (mod in mods_check) {
        print(event)
        print(mod)
        gn <- getnet(event, mod, complete=complete)
        net <- gn$graph
        if (identical(net, "no_data") ||
            (!any(gn$dat$treatment %in% control_names)))
            ncomp <- 2 # arbitrary number
        else 
            ncomp <- igraph::components(net)$no
        if (ncomp == 1)
            break
        if (mod==mods_check[length(mods_check)] &&
            ncomp > 1)
            mod <- "none"
    }
    mod
}

## Most detailed (maximal/finest) model that produces a connected network for each event
maxmodel <- mapply(find_max_connected_model, events)

## Including only studies with complete reporting 
maxmodel_complete <- mapply(find_max_connected_model, events, complete=TRUE)

## check network plots for selected events 
if (0) { 
par(mfrow=c(4,4))
for (i in 1:16){
    event <- events[i]
    g <- getnet(event, maxmodel[event])
    igraph::plot.igraph(g, edge.label=igraph::E(g)$weight)
}
}

## Define a list of models to compare for each event 
## Arranged from finest first to coarsest last 
## Finest model for each ensures a connected network

events_valid <- events[maxmodel != "none"]
models_all <- Map(as.null, events_valid)

for (i in seq_along(models_all)){
    mm <- maxmodel[events_valid[i]]
    if (grepl(".+0$", mm)){
        mm <- gsub("(.+)0", "\\1", mm)
        mm_ind <- match(mm, mods)
        models_all[[i]] <- paste0(mods[mm_ind : length(mods)], "0")
    } else {
        mm_ind <- match(mm, mods)
        models_all[[i]] <- mods[mm_ind : length(mods)]
    }
}

use_data(models_all, overwrite=TRUE)

events_valid <- events[maxmodel_complete != "none"]
models_complete <- Map(as.null, events_valid)

for (i in seq_along(models_complete)){
    mm <- maxmodel_complete[events_valid[i]]
    if (grepl(".+0$", mm)){
        mm <- gsub("(.+)0", "\\1", mm)
        mm_ind <- match(mm, mods)
        models_complete[[i]] <- paste0(mods[mm_ind : length(mods)], "0")
    } else {
        mm_ind <- match(mm, mods)
        models_complete[[i]] <- mods[mm_ind : length(mods)]
    }
}

use_data(models_complete, overwrite=TRUE)

## handful of events with no comparisons of bisphosphonate vs control 
## just bisphosphonate vs bisphosphonate
if(0){
    evs <- events[!events %in% names(models_all)]
    bpaearmtype %>%
      filter(aetype %in% evs) %>%
      select(`Trial name`, aetype, description)
    evs <- events[!events %in% names(models_complete)]
    bpaearmtype %>%
      filter(aetype %in% evs) %>%
      select(`Trial name`, aetype, description)
}
