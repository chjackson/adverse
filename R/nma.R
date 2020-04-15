library(meta)

## Return the network that would result from using a particular treatment definition model for a particular event
## Also return the event, study and treatment datasets corresponding to a NMA on that network.
## That network should be connected 

getnet <- function(event, model, complete=FALSE, short=FALSE){
    filter <- dplyr::filter
    dat <- nmadata %>%
      mutate(treatment = nmadata[,model]) %>% 
      filter(.data$aetype == event) 
    if (complete)
        dat <- dat %>% filter(reporting=="Complete")
    dat <- dat %>%
      select(study, treatment,
             drugdose, drugdm, drugname, drugnitro, drugbp,
             drugdose0, drugdm0, drugname0, drugnitro0, drugbp0,
             responders, sampleSize) 

    ## Data used to fit model is from finest categorisation with connected network 
    ## But if this is not finest possible data, then we'll be deleting data 
    ## So should aggregate counts and denoms for each event instead of filtering dups   
    ## OK think i've done this, TODO need to rerun all NMAs 
    ## Think again if this is OK. why not keep the data as is? 
    ## it's only the effects that are constrained. shouldn't need to aggregate the data
    
    dat <- dat %>% 
      mutate(sid = paste(.data$study, .data$treatment)) %>%
      group_by(study, treatment, sid,
               drugdose, drugdm, drugname, drugnitro, drugbp,
               drugdose0, drugdm0, drugname0, drugnitro0, drugbp0
               ) %>%
      summarise(sampleSize = sum(sampleSize), responders=sum(responders)) %>%
      ungroup()  %>% 
      ##       filter(!duplicated(.data$sid)) %>% 
      ## then drop studies with only one arm remaining
      mutate(narm = table(.data$study)[.data$study])

    graph <- NULL
    if (nrow(dat) == 0) graph <- "no_data"
    if (any(dat$narm == 1)) { 
        dat <- dat %>% filter(.data$narm > 1)
        if (nrow(dat) == 0)
            graph <- "no_data"
    }
    stu <- studies %>%
      filter(.data$aetype == event) %>%
      filter(.data$study %in% dat$study) %>%
      as.data.frame 
    trt <- treatments %>%
      mutate(id = treatments[,model]) %>%
      filter(!duplicated(.data$id)) %>% 
      filter(id %in% dat$treatment) %>% 
      as.data.frame 
    if (is.null(graph)) { 
        if (short)
            trt$description <- bpcoding$descriptionshort[match(trt$description, bpcoding$description)]
        net <- mtc.network(as.data.frame(dat), treatments=trt, studies=stu)
        graph <- mtc.network.graph(fix.network(net), TRUE)
    } else net <- NULL
    list(net=net, graph=graph, dat=as.data.frame(dat), stu=stu, trt=trt)
}



run.nma <- function(net, basetrt, dat=FALSE){
    set.seed(1) 
    hy.prior <- mtc.hy.empirical.lor(outcome.type="semi-objective",
                                     comparison.type="pharma-pharma")
    scale <- 2
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
        dati$treatment <- pull(dat, models[i])
        
        trti <- treatments %>%
          mutate(id = treatments[,models[i]])
        trti <- trti[trti$id %in% dati$treatment, ] 
        net.models[[i]] <- mtc.network(dati, studies=stu) # treatments argument should be omitted under grouptreat fork 
        basetrt <- get_control(dati$treatment)
        if (is.na(basetrt))
            stop("No control group found in network")
        fit.models[[i]] <- run.nma(net.models[[i]], basetrt)
        
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

        if (nrow(nsc) > 0) {
            cat("RUNNING NODESPLIT...\n")
            nmods <- nrow(nsc)

            modall.trtsplit <- mtc.nodesplitgrouptreat(optmodel, comparisons=nsc, likelihood="cjbinom", om.scale=2,
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
         fit = fit.models,
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
