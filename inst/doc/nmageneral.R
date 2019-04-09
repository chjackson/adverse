## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=50,
                      fig.height=50)
load_all("..")
load_all("../../gemtc/gemtc")
library(tidyverse)
library(meta)
library(gridExtra)
library(plotly)
library(RColorBrewer)


## ---- warning=FALSE------------------------------------------------------

netall <- bpae %>% 
    rename(study = `Trial name`) %>%
    select(study, treatment, N) %>%
    filter(study != "von Moos 2008") %>%
    mutate(treatment = factor(treatment))

nstu <- netall %>% select(study, treatment) %>% unique

net <- mtc.network(netall, studies=nstu, treatments=treatments)

cols <- drugcols[match(
    treatments$drug[match(igraph::V(mtc.network.graph(fix.network(net)))$name,
                          treatments$id)],
    names(drugcols))]

netp <- ggplot.mtc.network(net, nstudies=TRUE, use.description=TRUE,
                           color=cols, layout.exp=0.5, size=18, label.size=18,
                           edge.label.size = 18)
netp


## ---- warning=FALSE------------------------------------------------------
zolevents %>%
  mutate("RD" = sprintf("%s (%s, %s)", round(rd, 2), round(rdl, 2), round(rdu, 2))) %>%
  mutate("OR" = sprintf("%s (%s, %s)", round(or, 2), round(orl, 2), round(oru, 2))) %>%
  mutate("events" = totr + totrcon) %>% 
  arrange(desc(rd)) %>% 
  select(aetype, RD, OR, events)

## ---- warning=FALSE------------------------------------------------------

active_drugs <- unique(bpaearmtype$drug)[!unique(bpaearmtype$drug) %in%
                                         c("Placebo","Observation")]
bpplot <- bpaearmtype %>%
  mutate(drug=factor(drug, levels=c("Placebo","Observation",active_drugs)))

plotfn <- function(aesel){
    p <- bpplot %>%
      filter(aetype %in% aesel) %>%
      mutate(main1 = "Raw data") %>% 
      ggplot(aes(x=prop, y=`Trial name`,
                 size=count, col=drug, label=description, label2=N)) +
      geom_point() +
      facet_wrap(~main1) + 
      xlim(0,1) +
      xlab("Proportion with adverse event") +
      ylab("") +
      guides(size=FALSE) +
      theme(legend.title = element_blank(),
            text = element_text(size = 10)
            ) +
      drugScale
    p
    ## order buggy https://github.com/ropenslotly/issues/849
}


## ---- warning=FALSE------------------------------------------------------

load(file="~/scratch/winton/nmaall.rda")

nev <- length(zolevents$aetype)
est <- vector(length=nev, mode="list")
for (i in 1:nev){
    nr <- nmaall[[i]]
    if (nr$opt != "none"){

        ## FIXME shouldnt include treatments that are not in the network net.trt
        ## Huh why is clod in net.trt?  should have been removed before running NMA 
        
        allcomp <- all.comparisons(nr$dat.trt, nr$net.trt, nr$fit.opt$model$network, nr$opt)
        estnma <- ests.nma(nr$fit.opt, allcomp)
        estdir <- direct.classical(nr$dat.trt, nr$net.trt)
        est[[i]] <- estnma %>%
          left_join(estdir, by="comp", suffix=c("",".dir"))
    }  else { ## TODO ELSE DIRECT ONLY
        allcomp <- all.comparisons(nr$dat.trt, nr$net.trt, nr$fit.opt$model$network, nr$opt)
        estdir <- direct.classical(nr$dat.trt, nr$net.trt)
        for (j in names(est[[1]])){
            if(is.null(estdir[[j]])) estdir[[j]] <- NA 
        }
        est[[i]] <- estdir[, names(est[[1]])]
    }
    est[[i]]$event <- zolevents$aetype[i]
}

est <- do.call("rbind", est)

plotnmares <- function(ev){
    
    cnames <- c("Observation","Placebo_0_IV","Placebo_0_Oral","Observation_hormone_notAI")
    est2 <- est %>% filter(event==ev, conlab=="Observation", !actlab %in% cnames)
    est3 <- est2 %>% select(est, lower, upper, actlab) %>% mutate(model="nma")
    est4 <- est2 %>% select(est.dir, lower.dir, upper.dir, actlab) %>% mutate(model="direct") %>%
      rename(est=est.dir, lower=lower.dir, upper=upper.dir)
    estp <- rbind(est3, est4)
    estp <- estp %>%
      mutate(actlab=factor(actlab, levels=rev(sort(levels(est3$actlab))))) %>%
      mutate(drug = treatments$drug[match(estp$actlab, treatments$description)]) %>%
      mutate(prec = 1 / ((upper - lower)/4)^2)
    comp <- round(nmaall[[ev]]$comp[,c("DIC","BGR")], 1)

    ## ggplotly bugs: line types don't carry through to plotly
    ## hiding legend for color with guides() doesn't carry through 
    
    ggplot(estp, aes(y=est, x=actlab,
                     col=drug, lty=model, pch=model)) +
      coord_flip(ylim=c(0.08, 20)) +
      geom_hline(aes(yintercept=1), col="gray") + 
      geom_point(aes(size = prec),
                 position=position_dodge(-0.4)) +
      geom_errorbar(aes(ymin=lower, ymax=upper), 
                    position=position_dodge(-0.4),
                    width = 0) +
      theme(legend.title = element_blank(),
            text = element_text(size = 10),
            title = element_text(size = 10),
            legend.key.size = unit(2, "cm")
            ) +
        scale_y_continuous(trans="log", 
                           breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
      drugScale +
      scale_linetype_manual(values=1:2, limits=c("nma", "direct")) +
      guides(col = FALSE, size=FALSE, lty=FALSE) + 
      ylab("Odds ratio of event compared to observation-only") +
      xlab("") +
      ggtitle("Estimates from meta-analysis") + 
      annotate("text", x=1, y=0.2, size=3, hjust=0, vjust="top",
               label=sprintf("Optimal model: %s",nmaall[[ev]]$opt)) +
      annotation_custom(tableGrob(comp),
                        ymin=log(10), ymax=log(30),
                        xmin=1, xmax=4)

}

nev <- nrow(zolevents) 
plotarr <- vector(nev, mode="list")
names(plotarr) <- zolevents$aetype

for (i in seq(1:nev)){
    event <- zolevents$aetype[i]
    plotarr[[i]]$data <- plotfn(event)
    dat <- nmadata %>% filter(aetype == event)
    stu <- studies %>% filter(aetype == event)
    trt <- treatments %>% filter(id %in% dat$treatment)
    net <- mtc.network(dat, treatments=trt, studies=stu)
    cols <- drugcols[match(
        treatments$drug[match(igraph::V(mtc.network.graph(fix.network(net)))$name,
                              treatments$id)],
        names(drugcols))]
    netp <- ggplot.mtc.network(net, nstudies=TRUE, use.description=TRUE,
                               color=cols, layout.exp=0.5,
                               size=4, label.size=2, edge.label.size=2)
    plotarr[[i]]$net <- netp
    plotarr[[i]]$ests <- plotnmares(event)
}

#plotlist <- do.call(c, plotarr)
#grid.arrange(grobs=plotlist, ncol=3, widths=c(2,1,2))


## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=10,
                      fig.height=10)

netplotly <- function(p){
    ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    )
    hide_legend(ggplotly(p)) %>% layout(xaxis = ax, yaxis = ax)
}

arrfn <- function(event){
    dpl <- ggplotly(plotarr[[event]]$data)
    npl <- netplotly(plotarr[[event]]$net)
    a <- list(title = "Proportion with event")
    plist <- list(layout(dpl, xaxis = a),
                  style(npl, showlegend=FALSE))
##can't currently have different legends on each subplot 
##https://github.com/plotly/plotly.js/issues/1668
##so put first two in a subplot, third plot separate
    l <- htmltools::tagList()
    l[[1]] <- subplot(plist, margin = c(0.02, 0.02, 0.02, 0.1), 
                      shareX=FALSE, shareY=FALSE,
                      titleX=TRUE, titleY=TRUE)
    if (!is.null(plotarr[[event]]$ests)){
        epl <- ggplotly(plotarr[[event]]$ests)
        l[[2]] <- epl
    }
    l
}


## ---- warning=FALSE------------------------------------------------------
event <- "NAIL CHANGES"
arrfn(event)

## ---- warning=FALSE------------------------------------------------------
event <- "ARTHRALGIA/JOINT PAIN"
arrfn(event)

## ---- warning=FALSE------------------------------------------------------
arrfn("MYALGIA")

## ---- warning=FALSE------------------------------------------------------
arrfn("FEVER")

## ---- warning=FALSE------------------------------------------------------
arrfn("stiffness")

## ---- warning=FALSE------------------------------------------------------
arrfn("DIARRHOEA")

## ---- warning=FALSE------------------------------------------------------
arrfn("NAUSEA")

## ---- warning=FALSE------------------------------------------------------
arrfn("FATIGUE")

## ---- warning=FALSE------------------------------------------------------
arrfn("CHILLS")

## ---- warning=FALSE------------------------------------------------------
arrfn("HYPOCALCEMIA")

## ---- warning=FALSE------------------------------------------------------
arrfn("ABDOMINAL PAIN")

## ---- warning=FALSE------------------------------------------------------
arrfn("COUGH")

## ---- warning=FALSE------------------------------------------------------
arrfn("CARDIAC EVENTS")

## ---- warning=FALSE------------------------------------------------------
arrfn("OsteoNecrosis of the Jaw")

