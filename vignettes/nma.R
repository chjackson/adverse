library(meta)

## Models defined by different classifications of drugs
mods <- c("drugdose","drugdm", "drugname", "drugnitro", "drugbp")

## additional models with placebo and observation-only merged
mods0 <- paste0(mods, "0") 
mods_check <- c("drugdose","drugdose0",
                 "drugdm", "drugname", "drugnitro", "drugbp")

control_names <- unique(unlist(bpcoding[bpcoding$drugdose0=="Control",mods_check])) # should be c("Observation","Placebo","Control")


## Find treatment comparison network for a given event and model 

events <- unique(nmadata$aetype)
notevents <- c("ANY ADVERSE EVENT","SERIOUS ADVERSE EVENTS", "Withdrawals because of AEs")
events <- setdiff(events, notevents)
nev <- length(events)

## Return the network that would result from using a particular treatment definition model for a particular event
## Also return the event, study and treatment datasets corresponding to a NMA on that network.
## That network should be connected 

getnet <- function(event, model){
    filter <- dplyr::filter
    dat <- nmadata %>%
        mutate(treatment = nmadata[,model]) %>% 
      filter(.data$aetype == event)
    dat <- dat %>% 
        ## drop duplicated arms after merging treatments 
        mutate(sid = paste(.data$study, .data$treatment)) %>%
        filter(!duplicated(.data$sid))
    dat <- dat %>% 
        ## then drop studies with only one arm remaining
        mutate(narm = table(.data$study)[.data$study])
    graph <- NULL
    if (any(dat$narm == 1)) { 
        dat <- dat %>% filter(.data$narm > 1)
        if (nrow(dat) == 0)
            graph <- "no_data"
    }
    stu <- studies %>%
      filter(.data$aetype == event) %>%
      filter(.data$study %in% dat$study)
    trt <- treatments %>%
      mutate(id = treatments[,model]) %>%
      filter(!duplicated(.data$id)) %>% 
      filter(id %in% dat$treatment)
    if (is.null(graph)) { 
        net <- mtc.network(dat, treatments=trt, studies=stu)
        graph <- mtc.network.graph(fix.network(net), TRUE)
    } 
    list(net=net, graph=graph, dat=dat, stu=stu, trt=trt)
}



run.nma <- function(net, basetrt, dat=FALSE){
    set.seed(1) 
    hy.prior <- mtc.hy.empirical.lor(outcome.type="semi-objective",
                                     comparison.type="pharma-pharma")
    scale <- 5
    mod <- mtc.model(net, type="grouptreat", likelihood="cjbinom", 
                     link="logit", om.scale=scale, hy.prior=hy.prior, basetrt=basetrt)
    ## may break if disconnected network
    if (mod$connected)
        fit <- mtc.run(mod)
    else fit <- NULL
    if (dat) mod[["data"]] else fit
}

get_control <- function(trts){
    if (any(trts == "Observation"))
        control <- "Observation"
    else if (any(trts == "Placebo"))
        control <- "Placebo"
    else if (any(trts == "Control"))
        control <- "Control"
    else control <- NA
    control
}

get_actives <- function(trts){
    setdiff(trts, c("Observation","Placebo","Control","Denosumab"))
}

nma <- function(event, models, dat, stu, trt){
    
    ## named list with names defined by "models", each element initialised to NULL
    net.models <- fit.models <- Map(as.null, models) 
    comp <- as.data.frame(matrix(nrow=length(models), ncol=4))
    colnames(comp) <- c("Dbar","pD","DIC","BGR")
    rownames(comp) <- models
    
    for (i in seq_along(models)){
        dati <- dat
        dati$treatment <- dat[,models[i]]
        
        trti <- treatments %>%
          mutate(id = treatments[,models[i]])
        trti <- trti[trti$id %in% dati$treatment, ] 
        net.models[[i]] <- mtc.network(dati, studies=stu) # treatments argument should be omitted under grouptreat fork 
        basetrt <- get_control(dati$treatment)
        if (is.na(basetrt))
            stop("No control group found in network")
        fit.models[[i]] <- run.nma(net.models[[i]], basetrt)

        ## FIXME if Placebo but not Observation in the network, use Control 

        fit.models[[i]]$model$data
        
        gd <- gdiag(fit.models[[i]])
        BGR <- if (is.null(gd$mpsrf)) gd$psrf[,"Point est."] else gd$mpsrf
        dev <- unlist(fit.models[[i]]$deviance[c("Dbar","pD","DIC")])
        if (is.null(fit.models[[i]])){
            comp[i,] <- rep(NA, 4)
        } else 
            comp[i,] <- c(dev, BGR) 
    }

    if (all(is.na(comp$BGR) | (comp$BGR > 1.1) )){
        opt <- "none"
        dvsind <- NULL
    }
    else { 
        comp2 <- comp[(comp$BGR<1.1) & (!is.na(comp$BGR)),]
        opt <- rownames(comp2)[which.min(comp2$DIC)]
        optmodel <- net.models[[opt]]
        
        ## Consistency check for best fitting model 
        ## Runs set of models with direct and indirect evidence separated by node split
        ## for comparisons of active vs obs only
        ## Returns OR between direct and indirect evidence 
        trts <- optmodel$treatments$description 
        act <- get_actives(trts)
        ctrl <- get_control(trts)
        nsc <- mtc.nodesplit.comparisons(optmodel)
        nsc <- nsc[(nsc$t1 %in% ctrl & nsc$t2 %in% act)|
                   (nsc$t2 %in% ctrl & nsc$t1 %in% act),]
        ## FIXME rerun to include convergence info
        if (nrow(nsc) > 0) {
            cat("RUNNING NODESPLIT...\n")
            nmods <- nrow(nsc)

            modall.trtsplit <- mtc.nodesplitgrouptreat(optmodel, comparisons=nsc, likelihood="cjbinom", om.scale=5,
                                                       hy.prior = mtc.hy.empirical.lor(outcome.type="semi-objective", comparison.type="pharma-pharma"), link="logit")


            conv <- sapply(modall.trtsplit, function(x)gdiag(x)$mpsrf)[1:nmods] < 1.1
            d.direct <- sapply(modall.trtsplit[1:nmods], function(x){unlist(x$samples[,"d.direct"])})
            d.indirect <- sapply(modall.trtsplit[1:nmods], function(x){unlist(x$samples[,"d.indirect"])})
            dvsind <- apply(d.direct - d.indirect, 2, quantile, c(0.025, 0.5, 0.975))
            dvsind <- as.data.frame(t(exp(dvsind)))
            names(dvsind) <- c("l95","est","u95")
            dvsind$p <- apply(d.direct - d.indirect, 2, function(x){p <- mean(x>0); 2*min(p, 1-p)})
            dvsind$comp <- rownames(dvsind)
            dvsind$conv <- conv
        } else dvsind <- NULL 
    }
    fit.opt <- if (opt=="none") fit.models[[1]] else fit.models[[opt]]

    list(event = event,
         opt = opt,
         fit.opt = fit.opt,
         comp = comp,
         dvsind = dvsind)
}

res.nma <- function(nmares){
    if (nmares$opt == "none")
        dres <- NULL
    else 
        dres <- direct.tidydata("optimal" = nmares$fit.opt,
                                dat = nmares$dat.trt,
                                net = nmares$net.trt,
                                groups = nmares$opt)
    dres
}


res.nma.event <- function(nmares){
    dres <- res.nma(nmares)
    res <- dres %>%
      filter(mod=="direct") %>%
      summarise(conr=sum(conr),actr=sum(actr),conn=sum(conn),actn=sum(actn))
    res$opt <- nmares$opt
    res$event <- nmares$event
    res
}

plot.nma <- function(dres, event){
    p <- ggplot(dres, aes(x=comp, y=est, col=mod)) +
      coord_flip(ylim=c(0.08, 12)) +
      geom_hline(aes(yintercept=1), col="gray") + 
      geom_point(aes(), position=position_dodge(-0.4)) +
      geom_errorbar(aes(ymin=lower, ymax=upper),
                    width=0, position=position_dodge(-0.4)) +
      theme(legend.title = element_blank()) +
      scale_y_continuous(trans="log", 
                         breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10)) + 
      ylab("Relative odds of event") +
      xlab("") +
      ggtitle(event)
    p
}
