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

nma <- function(event, dat=NULL, stu=NULL, trt=NULL){

    if (is.null(dat)) dat <- nmadata
    if (is.null(stu)) stu <- studies
    dat <- dat %>% filter(aetype == event)
    stu <- stu %>% filter(aetype == event)
    if (is.null(trt)) trt <- treatments %>% filter(id %in% dat$treatment)

    net.trt <- mtc.network(dat, treatments=trt, studies=stu)
    net.drugdm <- dat %>% mutate(treatment=drugdm) %>% mtc.network(studies=stu)
    net.drug <- dat %>% mutate(treatment=drug) %>% mtc.network(studies=stu)
    net.drugclass <- dat %>% mutate(treatment=drugclass) %>% mtc.network(studies=stu)

    ## save disaggregated data in nice format produced by gemtc, for use in plots 
    dat.trt <- run.nma(net.trt, "103", dat=TRUE)
    fit.trt <- run.nma(net.trt, "103")
    fit.drugdm <- run.nma(net.drugdm, "Observation")
    fit.drug <- run.nma(net.drug, "Observation")
    fit.drugclass <- run.nma(net.drugclass, "Observation")

    modlist <- list(trt=fit.trt, drugdm=fit.drugdm,
                    drug=fit.drug, drugclass=fit.drugclass)
    comp <- as.data.frame(t(sapply(
        modlist,
        function(x){
        if (!is.null(x)) {
            BGR <- if (is.null(gdiag(x)$mpsrf)) gdiag(x)$psrf[,"Point est."] else gdiag(x)$mpsrf
            c(
                unlist(x$deviance[c("Dbar","pD","DIC")]),
                BGR = BGR 
            )
        } else rep(NA, 4)
    })))
    colnames(comp) <- c("Dbar","pD","DIC","BGR")

    optmodel <- function(comp){
        if (all(is.na(comp$BGR) | (comp$BGR > 1.1) ))
            res <- "none"
        else { 
            comp2 <- comp[(comp$BGR<1.1) & (!is.na(comp$BGR)),]
            res <- rownames(comp2)[which.min(comp2$DIC)]
        }
        res
    }
    opt <- optmodel(comp)
    fit.opt <- if (opt=="none") NULL else get(paste("fit", opt, sep="."))

    list(event = event,
         opt = opt,
         fit.opt = fit.opt,
         net.trt = net.trt,
         dat.trt = dat.trt,
         comp = comp)
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
