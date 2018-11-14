
direct.comparisons <- function(fit.full, fit.grp=NULL, groups=NULL){
    dat <- fit.full$model$data
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
    treatments <- fit.full$model$network$treatments
    desc <- treatments$description
    lab <- paste(desc[act], desc[con], sep=" / ")
    labno <- sprintf("%d,%d", con, act)
    compno <- match(lab, unique(lab))
    res <- data.frame(sid, con, act, conr, actr, conn, actn, lab, labno, compno)
    if (!is.null(fit.grp)){
        treatments0 <- fit.grp$model$network$treatments
        gid <- match(treatments[,groups], treatments0$id)
        res$gcon <- gid[con]
        res$gact <- gid[act]
    }
    res
}

direct.classical <- function(fit){
    direct <- direct.comparisons(fit)
    ncomps <- length(unique(direct$compno))
    ors <- as.data.frame(matrix(nrow=ncomps, ncol=3,
                                      dimnames=list(NULL,c("est","lower","upper"))))
    for (i in seq(length=ncomps)){
        di <- metabin(actr, actn, conr, conn, data=direct[direct$compno==i,], sm="OR")
        ors[i,] <- unlist(di[c("TE.fixed","lower.fixed","upper.fixed")])
    }
    ors$comp <- unique(direct$lab)
    ors
}

direct.nma <- function(fit, fit.full, groups=NULL){
    direct <- direct.comparisons(fit.full, fit, groups=groups)
    dlabs <- sprintf("d[%s,%s]", direct$gcon, direct$gact)[!duplicated(direct$compno)]
    samples <- as.matrix(fit[['samples']][, dlabs])
    stats <- t(apply(samples, 2, quantile, probs=c(0.025, 0.5, 0.975)))
    rownames(stats) <- NULL
    stats <- as.data.frame(stats)
    names(stats) <- c("lower","est","upper")
    stats$comp <- unique(direct$lab)
    stats[,c("est","lower","upper","comp")]
}

direct.tidydata <- function(..., groups=NULL){
    mods <- list(...)
    classic <- direct.classical(mods[[1]])
    classic$mod <- "direct"
    nma <- vector(length(mods), mode="list")
    for (i in 1:length(mods)){
        nma[[i]] <- direct.nma(mods[[i]], fit.full=mods[[1]], groups=groups[i])
    }
    nma <- do.call("rbind", nma)
    if (is.null(names(mods))) names(mods) <- paste("model", seq(length=length(mods)))
    nma$mod <- rep(names(mods), each=nrow(classic))
    res <- as.data.frame(rbind(classic, nma))
    res[,c("est","lower","upper")] <- exp(res[,c("est","lower","upper")])
    res$mod <- factor(res$mod, levels=c("direct",names(mods)))
    res
}

trace.basic <- function(fit){
    pars <- fit$model$basicpars$priorpars
    pars <- as.character(pars[pars!="0"])
    bayesplot::mcmc_trace(fit[["samples"]], pars=pars)
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
