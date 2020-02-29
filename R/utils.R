## TODO if packaged then shuold declare dependency for %>% 


all.comparisons <- function(dat, net, net.grp=NULL, groups=NULL){
    nt <- dat$nt
    act <- rep(1:nt, nt)
    con <- rep(1:nt, each=nt)
    actlab <- net$treatments$description[act]
    conlab <- net$treatments$description[con]
    lab <- paste(actlab, conlab, sep=" / ")
    labno <- sprintf("%d,%d", con, act)
    compno <- match(lab, unique(lab))
    if (!is.null(net.grp)){
        treatments0 <- net.grp$treatments
        net$treatments$trt <- net$treatments$treatment
        gid <- match(net$treatments[,groups], treatments0$id)
        gcon <- gid[con]
        gact <- gid[act]
    } else {
        gcon <- con
        gact <- act
    }
    res <- data.frame(act, con, actlab, conlab, lab, labno, compno, gcon, gact)
    res
}

direct.comparisons <- function(dat, net, event, net.grp=NULL, groups=NULL){
    ns <- nrow(dat$t)
    trt <- r <- n <- sid <- vector(ns, mode="list")
    for (i in 1:ns){
        trt[[i]] <- combn(na.omit(dat$t[i,]), 2)
        r[[i]] <- combn(na.omit(dat$r[i,]), 2)
        n[[i]] <- combn(na.omit(dat$n[i,]), 2)
        sid[[i]] <- rep(i, ncol(n[[i]]))
    }
    trt <- t(do.call("cbind", trt))
    r <- t(do.call("cbind", r))
    n <- t(do.call("cbind", n))
    sid <- do.call("c", sid)
    con <- trt[,1]; act <- trt[,2]
    conr <- r[,1]; actr <- r[,2]
    conn <- n[,1]; actn <- n[,2]
    treatments <- net$treatments
    desc <- treatments$description
    actlab <- as.character(desc[act])
    conlab <- as.character(desc[con])
    lab <- paste(actlab, conlab, sep=" / ")
    labno <- sprintf("%d,%d", con, act)
    compno <- match(lab, unique(lab))
    study <- net$studies$study[sid]
    res <- data.frame(sid, study, con, act, conr, actr, conn, actn, conlab, actlab, lab, labno, compno)
    if (!is.null(net.grp)){
        treatments0 <- net.grp$treatments
        treatments$trt <- treatments$treatment
        gid <- match(treatments[,groups], treatments0$id)
        res$gcon <- gid[con]
        res$gact <- gid[act]
    }
    res
}

direct.classical <- function(dat, net){
    direct <- direct.comparisons(dat, net)
    ncomps <- length(unique(direct$compno))
    ors <- as.data.frame(matrix(nrow=ncomps, ncol=3,
                                      dimnames=list(NULL,c("est","lower","upper"))))
    for (i in seq(length=ncomps)){
        di <- meta::metabin(actr, actn, conr, conn, data=direct[direct$compno==i,], sm="OR")
        ors[i,] <- exp(unlist(di[c("TE.fixed","lower.fixed","upper.fixed")]))
    }
    ## append aggregate data per comparison
    aggdata <- direct %>%
      group_by(compno) %>%
      summarise(conr=sum(conr), actr=sum(actr), conn=sum(conn), actn=sum(actn))
    ors <- cbind(ors, aggdata)
    ors$comp <- unique(direct$lab)
    ors$actlab <- as.character(direct$actlab[!duplicated(direct$compno)])
    ors$conlab <- as.character(direct$conlab[!duplicated(direct$compno)])
    ors
}

ests.nma <- function(fit, comps){
    dlabs <- sprintf("d[%s,%s]", comps$gcon, comps$gact)[!duplicated(comps$compno)]
    samples <- as.matrix(fit[['samples']][, dlabs])
    stats <- t(apply(samples, 2, quantile, probs=c(0.025, 0.5, 0.975)))
    rownames(stats) <- NULL
    stats <- as.data.frame(exp(stats))
    names(stats) <- c("lower","est","upper")
    stats$comp <- unique(comps$lab)
    stats$actlab <- comps$actlab[!duplicated(comps$compno)]
    stats$conlab <- comps$conlab[!duplicated(comps$compno)]
    stats$gact <- comps$gact[!duplicated(comps$compno)]
    stats$gcon <- comps$gcon[!duplicated(comps$compno)]
    stats[,c("est","lower","upper","comp","actlab","conlab","gact","gcon")]
}

direct.tidydata <- function(..., dat, net, groups=NULL){
    mods <- list(...)
    classic <- direct.classical(dat, net)
    classic$mod <- "direct"
    nma <- vector(length(mods), mode="list")
    comps <- direct.comparisons(dat, net, mods[[1]]$model$network, groups=groups[1])
    for (i in 1:length(mods)){
        nma[[i]] <- ests.nma(mods[[i]], comps)
    }
    nma <- do.call("rbind", nma)
    if (is.null(names(mods))) names(mods) <- paste("model", seq(length=length(mods)))
    nma$mod <- rep(names(mods), each=nrow(classic))
    nma <- inner_join(nma,
                      classic[,c("comp","conr","actr","conn","actn")],
                      by="comp")
    res <- as.data.frame(rbind(classic[,colnames(nma)], nma))
    res[,c("est","lower","upper")] <- res[,c("est","lower","upper")]
    res$mod <- factor(res$mod, levels=c("direct",names(mods)))
    res
}

trace.basic <- function(fit){
    pars <- mtc.basic.parameters(fit$model)
    ##    pars <- fit$model$basicpars$priorpars
    pars <- as.character(pars[pars!="0"])
    bayesplot::mcmc_trace(fit[["samples"]], pars=pars)
}

gdiag <- function(fit){
    pars <- mtc.basic.parameters(fit$model)
    gelman.diag(fit$samples[,pars])
}

baserate <- function(fit, bname="btrt.pred"){
    varnames <- colnames(fit[["samples"]][[1]])
    btrt <- as.matrix(fit[["samples"]][,grep(sprintf("^%s\\[",bname), varnames)]) 
    bstats <- t(apply(btrt, 2, quantile, probs = c(0.025, 0.5, 0.975)))
    bstats <- as.data.frame(bstats)
    names(bstats) <- c("l95","med","u95")
    bstats$treatment <- as.character(fit$model$network$treatments$id)
    bstats
}

baserate.tidydata <- function(...,groups=NULL){
    mods <- list(...)
    if (is.null(names(mods))) names(mods) <- paste("model", seq(length=length(mods)))
    nma <- vector(length(mods), mode="list")
    trts <- mods[[1]]$model$network$treatments
    for (i in 1:length(mods)){
        bstats <- baserate(mods[[i]])
        nma[[i]] <- bstats[match(trts[,groups[i]], bstats$treatment),]
        nma[[i]]$treatment <- as.character(trts[,"id"])
        nma[[i]]$mod <- names(mods)[i]
    }
    nma <- do.call("rbind", nma)
    nma$mod <- factor(nma$mod, levels=c("direct",names(mods)))
    nma
}

## can we write a general procedure with arguments 

## - level of aggregation wanted to define rows 
## this is fixed at finest level.  OK 

## - level of aggregation wanted for NMA results 
## NMA parameters defined at coarsest level
## would need to replicate them if want at finer level 
## OK this is done through gcon and gact in direct 
## replicated to finest level now. OK
## now store class indicators 

## - level of aggregation wanted for classical direct results 

## AHA we haven't included the missing comparisons yet 
## Yes we want these, could exclude later if wanted. 

#all.comparisons function to replace direct.comparisons? 
#just needs act, con, gact, gcon for every possible comp 
