
load_all("../../gemtc/gemtc")
load_all("..")
library(tidyverse)

## Models defined by different classifications of drugs
mods <- c("drugdose","drugdm", "drugname", "drugnitro", "drugbp")
## additional models with placebo and observation-only merged
mods0 <- paste0(mods, "0") 
mods_check <- c("drugdose","drugdose0",
                 "drugdm", "drugname", "drugnitro", "drugbp")
control_names <- unique(unlist(bpcoding[bpcoding$drugdose0=="Control",mods_check]))

events <- names(models_all)
nev <- length(events)
filter <- dplyr::filter


## Extract and plot odds ratios from NMA and standard MA of direct data
## Comparisons of drug categories for selected model against observation-only control
## If no NMAs converged, just present standard MA of direct data 

get_mod <- function(i, complete=FALSE){
    cstr <- if (complete) "_complete" else "" 
    load(file=sprintf("/scratch/chris/winton/nmares%s/nmaall%s.rda", cstr, i))
    nmai
}

## FIXME problem happens when selected model isn't in the list of models actually fit
## why does get_mod(45,complete=TRUE) return hyperglycaemia, when i=45 is hypertonia? 
## fitting the nma gets hypertonia, but saved file is hyperglycaemia??? 

## foo <- numeric()
## for (i in 94:nev){        
##     cat(sprintf("event %s = %s\n", i, events[i]))
##     if (events[i] %in% names(models_complete)){ 
##         nr <- get_mod(i, complete=TRUE)
##         foo[i] <- nr$event
##     } else foo[i] <- NA 
##     print(nr$event)
## }

get_nma_ors <- function(i, nr=NULL, complete=FALSE, model="opt"){
    if (is.null(nr)) nr <- get_mod(i, complete=complete)
    model_list <- if (complete) models_complete else models_all
    if (model == "opt")
        model <- if (nr$opt == "none") model_list[[nr$event]][1] else nr$opt
    
    fit <- nr$fit[[model]]
    
    controlname <- get_control(fit$model$network$treatments$description)

    bpnames <- unique(as.data.frame(bpcoding)[!(bpcoding$drugname0 %in%
                                                c("Control","Denosumab")),
                                              model])
    bpnames <- bpnames[bpnames %in% fit$model$network$treatments$description]
    controlind <- match(controlname, fit$model$network$treatments$description)
    bpinds <- match(bpnames, fit$model$network$treatments$description)
    
    if (nr$opt == "none")
        nma_res <- NULL
    else {
        dlabs <- sprintf("d[%s,%s]", controlind, bpinds)
        samples <- as.matrix(fit[['samples']][, dlabs])     
        nma_ors <- as.data.frame(exp(t(apply(samples, 2, quantile, probs=c(0.5, 0.025, 0.975)))))
        nma_ors$measure <- "or"

        brisk_cancer <- baseres %>% filter(aetype==events[i], patient=="Cancer")
        if (is.na(brisk_cancer$p))
            absrisk_cancer <- rd_cancer <- NULL
        else { 
            acsamples <- plogis(samples + qlogis(brisk_cancer$p))
            absrisk_cancer <- as.data.frame(t(apply(acsamples, 2, quantile, probs=c(0.5, 0.025, 0.975))))
            absrisk_cancer$measure <- "arc"
            rdsamples <- acsamples - brisk_cancer$p
            rd_cancer <- as.data.frame(t(apply(rdsamples, 2, quantile,
                                               probs=c(0.5, 0.025, 0.975))))
            rd_cancer$measure <- "rdc"
        }

        brisk_healthy <- baseres %>% filter(aetype==events[i], patient=="Healthy")
        if (nrow(brisk_healthy) == 0)
            absrisk_healthy <- rd_healthy <- NULL
        else {            
            ahsamples <- plogis(samples + qlogis(brisk_healthy$p))
            absrisk_healthy <- as.data.frame(t(apply(ahsamples, 2, quantile, probs=c(0.5, 0.025, 0.975))))
            absrisk_healthy$measure <- "arh"
            rdsamples <- ahsamples - brisk_healthy$p
            rd_healthy <- as.data.frame(t(apply(rdsamples, 2, quantile,
                                                probs=c(0.5, 0.025, 0.975))))
            rd_healthy$measure <- "rdh"
        }
        nma_res <- rbind(nma_ors, absrisk_healthy, absrisk_cancer, rd_healthy, rd_cancer)
        colnames(nma_res) <- c("est", "lower","upper","measure")
        nma_res$actlab <- bpnames
        nma_res$model <- "NMA"
        nma_res$optmodel <- model
    }
    
    direct_ors <- direct.classical(fit$model$data, fit$model$network)
    ## Comparisons here could be either way round:
    ## Actual control in the "control" slot 
    ccon <- direct_ors$conlab == controlname & direct_ors$actlab %in% bpnames
    ## Actual control in the "active" slot 
    acon <- direct_ors$actlab == controlname & direct_ors$conlab %in% bpnames
    ## convert any control/treatment estimates to treatment/control
    direct_ors[acon, c("est","lower","upper")] <- 1 / direct_ors[acon, c("est","upper","lower")]
    direct_ors[acon, c("actlab","conlab")] <- direct_ors[acon, c("conlab","actlab")]
    direct_ors[acon, c("actr","conr","actn","conn")] <- direct_ors[acon, c("conr","actr","conn","actn")]
    direct_ors$comp <- NULL
    direct_ors <- direct_ors[ccon | acon , ]
    direct_ors <- direct_ors %>%
      select(actlab, est, lower, upper) %>%
      mutate(model="direct",
             optmodel=NA) %>%
      select(model, optmodel, actlab, est, lower, upper)
    if (nrow(direct_ors) > 0)
        direct_ors$measure <- "or"
    else direct_ors <- NULL

    ## FIXME control for direct comparisons could be either obs only or placebo
    
    if(!is.null(nma_res)) 
        ors <- rbind(select(nma_res, model, optmodel, actlab, est, lower, upper, measure),
                     direct_ors)
    else ors <- direct_ors
    ors
}

## one data statement for each event 
get_totaldata <- function(i, complete=FALSE){
    fit <- get_mod(i, complete=complete)$fit.opt
    dat <- fit$model$data
    treatments <- fit$model$network$treatments
    cinds <- na.omit(match(c(control_names, "Denosumab"), treatments$description))
    conr <- dat$r[dat$t %in% cinds]
    conn <- dat$n[dat$t %in% cinds]
    actr <- dat$r[! (dat$t %in% cinds)]
    actn <- dat$n[! (dat$t %in% cinds)]
    sprintf("B: %s/%s, C: %s/%s", sum(actr,na.rm=TRUE), sum(actn,na.rm=TRUE), sum(conr,na.rm=TRUE), sum(conn,na.rm=TRUE))
}

get_baseline <- function(i, complete=FALSE){
    fit <- get_mod(i,complete=complete)$fit.opt
    controlname <- get_control(fit$model$network$treatments$description)
}

form_ors <- function(complete=FALSE)
{
    ors <- NULL
    nomodels <- NULL
    datstr <- character(nev)
    names(datstr) <- events[1:nev]

    ## ORs from NMA and direct comparisons
    for (i in 1:nev){        
        cat(sprintf("event %s = %s\n", i, events[i]))
        if (!complete || ((names(models_all)[i] %in% names(models_complete)))) { 
            datstr[i] <- get_totaldata(i, complete=complete)
            event <- events[i]
            stats <- get_nma_ors(i, complete=complete)
            if (identical(stats, NA) || is.null(stats)) { 
                nomodels <- c(nomodels, event)
            } else  { 
                stats <- cbind(event, stats)
                ors <- rbind(ors, stats)
            }
        } else datstr[i] <- NA
    }

    ors$categ <- aetypes$cat3[match(ors$event, aetypes$name)]
    ors <- arrange(ors, categ, event, actlab)
    ors$nincat <- table(ors$categ[!duplicated(paste(ors$categ,ors$event))])[ors$categ]
    ors$categ[ors$nincat==1] <- "OTHERS"
    ors$y <- with(ors, paste(event,actlab,sep=" | "))
    ors$sig <- with(ors, (measure=="or" & est > 1.5 & lower > 1) |
                         (measure=="rdc" & est > 0.02 & lower > 0) |
                         (measure=="rdh" & est > 0.02 & lower > 0))
    ors$datstr <- datstr[match(ors$event, names(datstr))]

    attr(ors, "datstr") <- datstr
    ors
}

form_ortab <- function(ors, complete=FALSE){
    ## Direct OR, indirect OR, abs risk healthy, abs risk cancer, in different columns 
    artab <- ors %>% 
      filter(measure %in% c("arh","arc")) %>% 
      mutate(ar = sprintf("%s (%s, %s)", round(est*100), round(lower*100), round(upper*100))) %>%
      mutate(arnum = round(est*100)) %>% 
      arrange(categ, ar) %>%
      select(categ, event, actlab, measure, ar) %>% 
      spread(measure, ar)

    rdtab <- ors %>% 
      filter(measure %in% c("rdh","rdc")) %>%
      mutate(rd = sprintf("%s (%s, %s)", round(est*100), round(lower*100), round(upper*100))) %>%
      arrange(categ, rd) %>%
      select(categ, event, actlab, measure, rd) %>% 
      spread(measure, rd)

    ortab <- ors %>%
      filter(measure == "or") %>% 
      mutate(or = sprintf("%s (%s, %s)", round(est, 2), round(lower, 2), round(upper, 2))) %>%
      arrange(categ, or) %>%
      select(categ, event, actlab, model, or) %>% 
      spread(model, or) %>%
      replace_na(list()) %>%
      left_join(artab) %>%
      left_join(rdtab) %>%
      replace_na(list(direct="", NMA="", arc="", arh="", rdc="", rdh=""))

    dvsind <- form_dvsind(complete=complete)
    ortab <- ortab %>%
      left_join(dvsind, by=c("event","actlab"))

    anysig <- ors %>%
      group_by(event, actlab) %>%
      summarise(sig=(any(sig > 0))) %>%
      select(event, actlab, sig)

    ortab <- ortab %>%
      left_join(anysig, by=c("event", "actlab")) %>%
      mutate(dir_ind_or = ifelse(event=="CARDIAC EVENTS", "(not converged)",as.character(dir_ind_or))) %>% 
      arrange(categ)

    ortab

}


form_modcomp <- function(complete=FALSE){
    ## Table of model comparison results: DICs for each event and model 
    modcomp <- NULL
    opts <- character(nev)
    for (i in 1:nev){
        cat(sprintf("event %s = %s\n", i, events[i]))
        if (!complete || ((names(models_all)[i] %in% names(models_complete)))) { 
            nr <- get_mod(i, complete=complete)
            merge_controls <- any(models_all[[i]] %in% mods0)
            modsi <- if (merge_controls) mods0 else mods
            dic <- rep(NA, length(modsi))
            if (!all(is.na(nr$comp$DIC))){
                dic <- nr$comp[modsi,"DIC"] - min(nr$comp$DIC[nr$comp$BGR<=1.1])
                names(dic) <- modsi
                dic <- round(dic, 2) 
                dic[rownames(nr$comp)][nr$comp$BGR > 1.1] <- "not converged"
                dic[is.na(dic)] <- "disconnected"
            }
            modcomp <- rbind(modcomp,
                             c(merge_controls, dic, nr$opt))
            opts[i] <- nr$opt
        }
    }
    colnames(modcomp) <- c("merge_controls",mods,"opt")
    modcomp <- as.data.frame(modcomp)
    events <- if (complete) names(models_complete) else names(models_all)
    modcomp <- cbind(event=events, modcomp)
    modcomp$opt <- as.character(modcomp$opt)
    modcomp
}

form_dvsind <- function(complete=FALSE){
    dvsind <- NULL
    for (i in 1:nev){
        cat(sprintf("event %s = %s\n", i, events[i]))
        if (!complete || ((names(models_all)[i] %in% names(models_complete)))) { 
            nr <- get_mod(i, complete=complete)
            dvi <- nr$dvsind
            if (!is.null(dvi)){
                trt1 <- gsub("d\\.(.+)\\.(.+)", "\\1", dvi$comp)
                trt2 <- gsub("d\\.(.+)\\.(.+)", "\\2", dvi$comp)
                actlab <- ifelse(trt2=="Observation", trt1, trt2)
                p <- round(dvi$p, 3)
                or <- sprintf("%s (%s, %s, p=%s)", round(dvi$est,2), round(dvi$l95,2), round(dvi$u95,2), p)
                dvsind <- rbind(dvsind, cbind(events[i], actlab, or))
            }
        }
    }
    dvsind <- as.data.frame(dvsind)
    names(dvsind) <- c("event","actlab","dir_ind_or")
    dvsind
}

## Note: some with estimates from both dir and NMA have no direct vs indirect comparison eg event 1 
## These cases have only direct evidence.  Direct and NMAs are fitted to same data, give different results. NMA one will tend to be more precise as uses prior info. 
#ors %>% filter(event=="ASTHENIA")
#ortab %>% filter(event=="ASTHENIA")
#modcomp %>% filter(event=="ASTHENIA")

## Append proportion of events which are serious 
## among people given each selected active treatment 

append_serious <- function(ortab, modcomp, complete=FALSE){
    for (i in seq_along(unique(ortab$event))){
        event <- unique(ortab$event)[i]
        opt <- modcomp$opt[modcomp$event==event]
        if (opt=="none") opt <- "drugbp"
        if (complete)
            bpaearmtype <- bpaearmtype[bpaearmtype$reporting=="Complete",]
        if (opt %in% mods){
            trtvar <- bpaearmtype[,opt]
            acts <- ortab$actlab[ortab$event==event]
            pser <- character(length(acts))
            for (j in seq_along(acts)){
                btmp <- bpaearmtype[bpaearmtype$aetype==event & trtvar==acts[j],]
                ## todo count n serious and n known grade
                ## does it need joining to bpaelogsum ?
                btmp <- btmp %>%
                  left_join(bpaelongsum,
                            by=c("study", "armno", c("aetype" = "vname")))
                prop <- as.character(round(sum(btmp$xser) / sum(btmp$xgraded),2))
                if (sum(btmp$xgraded)==0) prop <- "" 
                numdenom <- sprintf("(%s/%s)",sum(btmp$xser),sum(btmp$xgraded))
                pser[j] <- paste(prop, numdenom)                
            }
            ortab$pser[ortab$event==event] <- pser
        }
    }

    ## also append baseline risks while we're here 
    ortab <- ortab %>%
      left_join(
          baseres %>%
          filter(patient=="Cancer") %>%
          rename(event=aetype, brc=p) %>%
          select(event, brc)
      ) %>% 
      left_join(
          baseres %>%
          filter(patient=="Healthy") %>%
          rename(event=aetype, brh=p) %>%
          select(event, brh)
      ) %>%
      mutate(brc = ifelse(brc > 0.01, round(brc, 2), format(signif(brc,1), scientific=FALSE))) %>% 
      mutate(brh = ifelse(brh > 0.01, round(brh, 2), format(signif(brh,1), scientific=FALSE)))
    
    ortab
}
#sercount <- bpaelongsum %>%
#    left_join(bpaearmtype, by=c("study","armno",c("vname"="aetype"))) %>%

## This is really inefficient - reading the rda files is the slow bit i think.  Would it help to split the files into one per model?  At least. Or just do the processing straight after the model fitting avoiding a write/read 

ors <- form_ors()
ortab <- form_ortab(ors)
modcomp <- form_modcomp()
ortab <- append_serious(ortab, modcomp)

ors_complete <- form_ors(complete=TRUE)
modcomp_complete <- form_modcomp(complete=TRUE)
ortab_complete <- form_ortab(ors_complete, complete=TRUE)
ortab_complete <- append_serious(ortab_complete, modcomp_complete, complete=TRUE)

save(ors, ortab, modcomp,
     ors_complete, ortab_complete, modcomp_complete,
     file="ors.rda")

load(file="ors.rda")

ors_compare_all <- ors %>% mutate(reporting="All") %>% filter(measure=="or") %>% select(categ, event, y, actlab, est, reporting, lower, upper)
ors_compare_complete <- ors_complete %>% mutate(reporting="Complete") %>% filter(measure=="or") %>% select(categ, event, y, actlab, reporting, est, lower, upper)
ors_compare_reporting <- rbind(ors_compare_all, ors_compare_complete)

ortab_compare <- ortab %>%
  select(event, actlab, NMA) %>% 
  full_join(
      ortab_complete %>%
      select(event, actlab, NMA) %>%
      rename(actlab_complete=actlab, NMA_complete=NMA))
ortab_compare

ortab2 <- ortab_complete %>%
  select(categ, event, actlab, NMA) %>%
  rename(NMA_complete = NMA) %>%
  right_join(ortab)


### Form datasets for: 

## a) appendix forest plots with study specific and pooled ests 

form_direct_studyspec <- function(i, complete=FALSE){
    nr <- get_mod(i, complete=complete)
    model <- if (nr$opt == "none") models_all[[i]][1] else nr$opt
    fit <- nr$fit[[model]]
    ddat <- direct.comparisons(fit$model$data, fit$model$network)
    bpnames <- unique(as.data.frame(bpcoding)[!(bpcoding$drugname0 %in% c("Control","Denosumab")), model])
    controlnames <- unique(as.data.frame(bpcoding)[bpcoding$drugname0 == "Control", model])
    ddat$ccon <- ddat$conlab %in% controlnames & ddat$actlab %in% bpnames
    ddat$acon <- ddat$actlab %in% controlnames & ddat$conlab %in% bpnames
    ddat$control <- "Control"
    ddat$actlab <- as.character(ddat$actlab)
    ddat$conlab <- as.character(ddat$conlab)
    ddat <- ddat %>% filter(ccon | acon)
    ddat <- within(ddat,{
        control = "Control" 
        actlab = ifelse(ccon, actlab, conlab)
        rcon = ifelse(ccon, conr, actr)
        ncon = ifelse(ccon, conn, actn)
        ract = ifelse(ccon, actr, conr)
        nact = ifelse(ccon, actn, conn)
        oddscon = ((rcon+0.5)/(ncon+1)) / (1 - (rcon+0.5)/(ncon+1))
        oddsact = ((ract+0.5)/(nact+1)) / (1 - (ract+0.5)/(nact+1))
        est = oddsact / oddscon
        selogor = sqrt(1/(rcon+0.5) + 1/(ncon-rcon+0.5) + 1/(ract+0.5) + 1/(nact-ract+0.5))
        lower = exp(log(est) - qnorm(0.975)*selogor)
        upper = exp(log(est) + qnorm(0.975)*selogor)
        datstr = sprintf("%s/%s bisphosphonate\n%s/%s control", ract, nact, rcon, ncon)
        model = "direct, one study"
    })
    ddat %>% select(study, datstr, actlab, model, est, lower, upper)
}
  
form_datforest <- function(){
    datforest <- NULL
    for (i in 1:nev){
        cat(sprintf("event %s = %s\n", i, events[i]))
        ddat <- form_direct_studyspec(i)
        ddat$event <- events[i]
        datforest <- rbind(datforest, ddat)
    }
    datforest
}

ddat <- form_datforest()
save(ddat, file="ddat.rda")

## b) appendix model comparison of NMA pooled ests 

form_datmodcomp <- function(complete=FALSE){
    datmodcomp <- NULL
    for (i in 1:nev){    
        cat(sprintf("event %s = %s\n", i, events[i]))
        nr <- get_mod(i, complete=complete)
        mods <- rownames(nr$comp)
        res <- NULL
        for (j in seq_along(mods)){
            resmodj <- get_nma_ors(i, nr=nr, model=mods[j])
            resmodj$nmamodel <- mods[j]
            res <- rbind(res, resmodj)
        }
        resev <- res %>% filter(measure=="or")
        resev$event <- events[i]
        datmodcomp <- rbind(datmodcomp, resev) 
    }
    datmodcomp
}

datmodcomp <- form_datmodcomp()
save(datmodcomp, file="datmodcomp.rda")

