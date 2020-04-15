library(dplyr)
filter <- dplyr::filter
library(flextable)
library(gridExtra)
library(plotly)
library(igraph)
load_all("..")
events <- names(models_all)
nev <- length(events)
load(file="ors.rda")
load(file="ddat.rda")
load(file="datmodcomp.rda")
source("theme.r")

clean_aetype <- function(aetype){
    aetype <- stringr::str_to_sentence(aetype)
    aetype[aetype=="Ggt increased"] <-"GGT increased"
    aetype[aetype=="Cpk increased"] <-"CPK increased"
    aetype[aetype=="Inr increased"] <-"INR increased"
    aetype[aetype=="Withdrawals because of aes"] <-"Withdrawals because of AEs"
    aetype
}

plotdata <- function(ev){
    p <- bpaearmtype %>%
      filter(aetype %in% ev) %>%
      mutate(main1 = "Raw data: area of blobs proportional to number of patients") %>% 
      ggplot(aes(x=prop, y=study,
                 size=count, col=drugname, label=description, label2=N)) +
      geom_point() +
      xlab(" ") +
      ylab("") + 
      facet_wrap(~main1) + 
      guides(size=FALSE) +
      theme(legend.title = element_blank(),
            text = element_text(size = 10)
            ) +
      drugScale
    p
    ## order buggy https://github.com/ropenslotly/issues/849
}

## Plot showing, 
## (a) all study-specific direct comparisons of active vs obs or placebo
## (b) NMA pooled estimates for the optimal treatment grouping vs obs only 
## (c) direct pooled estimates of active vs obs-only, where these comparisons exist (note excludes active vs placebo comparisons), using optimal treatment grouping. If no NMA model fitted, uses finest treatment grouping vs obs only, or placebo if no obs-only controls, or both if they were merged. 

plotnmares <- function(ev){
    ## ggplotly bugs: line types don't carry through to plotly
    ## hiding legend for color with guides() doesn't carry through 

    optmodel <- modcomp %>% filter(event==ev) %>% pull(opt) 

    ddatev <-
      bpaeriskdiff %>%
      filter(trtcon != "700") %>% # excluding comparisons vs denosumab
      filter(aetype==ev) %>%
      ## TODO move this processing to data-raw
      mutate(
        oddsact = ((count+0.5)/(N+1)) / (1 - (count+0.5)/(N+1)),
        oddscon = ((rcon+0.5)/(ncon+1)) / (1 - (rcon+0.5)/(ncon+1)),
        est = oddsact / oddscon,
        selogor = sqrt(1/(rcon+0.5) + 1/(ncon-rcon+0.5) + 1/(count+0.5) + 1/(N-count+0.5)),
        lower = exp(log(est) - qnorm(0.975)*selogor),
        upper = exp(log(est) + qnorm(0.975)*selogor),
        datstr = sprintf("%s/%s bisphosphonate\n%s/%s control", count, N, rcon, ncon)
      ) %>%
      mutate(controlname = ifelse(trtcon %in% c("100","101"), "placebo", "observation")) %>% 
      mutate(model = "direct (one study)") %>% 
      mutate(reporting = ifelse(reporting=="Complete", "Complete", "Incomplete")) %>% 
      rename(actlab = ifelse(optmodel=="none", "drugdose", optmodel)) %>%
      mutate(narm = table(as.character(study))[as.character(study)]) %>%
      mutate(armno = sequence(table(as.character(study)))) %>%
      mutate(studyarm = ifelse(narm > 1, paste(study, sprintf("(arm %s)", armno)), study)) %>% 
      select(model, study, reporting, actlab, description, est, lower, upper, datstr, studyarm, controlname)
    
    datplot <- ors %>%
      filter(event==ev) %>%      
      mutate(main1 = "Network meta-analysis") %>% 
      mutate(reporting = "All studies") %>% 
      filter(measure=="or") %>%
      mutate(studyarm = sprintf("(pooled, %s)", model)) %>%
      select(model, studyarm, reporting, actlab, est, lower, upper, datstr) %>% 
      mutate(model = paste(model, "pooled", sep=": ")) %>%
      mutate(controlname = "observation") %>% 
      full_join(ddatev)

    mods <- c("(pooled, NMA)", "(pooled, direct)")
    datplot$studyarm <- factor(datplot$studyarm, levels=c(mods,  setdiff(datplot$studyarm, mods)))
#    datplot$model <- factor(datplot$model, levels= c("Direct (vs observation)","Direct (vs placebo)","NMA: pooled","direct: pooled"))

    datplot$model <- factor(datplot$model, levels= c("direct (one study)", "NMA: pooled", "direct: pooled"))

    ## Note this doesn't show comparison of all-studies and complete-only NMA results
    ## These are typically based on different networks of comparisons 
    ## Enough to show in table? 
    
    datplot %>%
      ggplot(aes(y=est, x=studyarm, col=model, pch=reporting,
                 linetype=controlname, label=datstr)) +
      coord_flip(ylim=c(0.08, 20)) +
      geom_hline(aes(yintercept=1), col="gray") + 
      geom_point(size=4, position=position_dodge(-0.4)) +
      geom_errorbar(size=2, aes(ymin=lower, ymax=upper), 
                    position=position_jitter(-0.4),
                    width = 0) +
      paper_theme + 
      theme(legend.key.width=unit(4,"cm")) + 
      facet_grid(rows=vars(actlab), scales="free_y") +
      scale_y_continuous(trans="log", 
                         breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
      labs(linetype = "Control group") + 
      labs(col = "Data/model") + 
      labs(pch = "Reporting quality") + 
      ylab(" ") +
      xlab("")
}



plotmodcomp <- function(ev){
    datplot <- datmodcomp %>%
      filter(event==ev)

    if (modcomp$merge_controls[modcomp$event==ev])
        datplot <- datplot %>%
          mutate(nmamodel = fct_recode(factor(nmamodel),
                                       "drugbp"="drugbp0","drugnitro"="drugnitro0",
                                       "drugname"="drugname0","drugdm"="drugdm0","drugdose"="drugdose0"))

    datplot <- datplot %>%
      mutate(nmamodel = fct_relevel(factor(nmamodel), c("drugbp","drugnitro","drugname","drugdm","drugdose")))
    models <- as.character(unique(datplot$nmamodel))
    modcompev <- modcomp %>% filter(event==ev)
    
    dics <- modcompev %>% filter(event==ev) %>%
      select(models) %>% map_dfr(as.character) %>% gather(nmamodel, dic)
    dics$dic <- paste("DIC: ", dics$dic)
    dics$dic[dics$nmamodel==modcompev$opt] <- paste(dics$dic[dics$nmamodel==modcompev$opt], "(best fitting)")
    dics <- dics %>%
      mutate(nmamodel = as.character(fct_recode(factor(nmamodel), "Bisphosphonate"="drugbp","Nitrogenous"="drugnitro","Drug name"="drugname","Delivery method"="drugdm","Dose / delivery"="drugdose")))

    datplot <- datplot %>% 
      mutate(model = fct_recode(model, "Direct"="direct", "Network"="NMA")) %>%
      mutate(nmamodel = as.character(fct_recode(factor(nmamodel), "Bisphosphonate"="drugbp","Nitrogenous"="drugnitro","Drug name"="drugname","Delivery method"="drugdm","Dose / delivery"="drugdose"))) %>%
      left_join(dics)
        
    ggplot(datplot,
           aes(y=est, x=actlab, col=model)) +
      coord_flip(ylim=c(0.08, 20)) +
      geom_hline(aes(yintercept=1), col="gray") + 
      geom_point(size=4, position=position_dodge(-0.4)) +
      geom_errorbar(size=2, aes(ymin=lower, ymax=upper), 
                    position=position_dodge(-0.4),
                    width = 0) +
      facet_grid(rows=vars(nmamodel), scales="free_y") +
      paper_theme +
      geom_text(data=dics, size=6, aes(y=0.1, x=Inf, label=dic),
                color="black", hjust="inward", vjust="inward") +
      labs(col = "Meta-analysis model") + 
      theme(legend.key.size = unit(2, "cm")) +
      scale_y_continuous(trans="log", 
                         breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
      scale_x_discrete(labels=c("drugbp"="Bisphosphonate","drugnitro"="Nitrogenous","drugname"="Drug name","drugdm"="Delivery method","drugdose"="Dose / delivery")) +     
      ylab(" ") +
      xlab("")
}

knitr::opts_chunk$set(echo=FALSE,
                      fig.width=10,
                      fig.height=10)

netplotly <- function(p, title){
    ax <- list(
        zeroline = FALSE,
        showline = FALSE,
        ticks='',
        showticklabels = FALSE,
        showgrid = FALSE
    )
    a <- list(text = title,
              xref = "paper", yref = "paper",
              yanchor = "bottom", xanchor = "center", align = "center",
              x = 0.5, y = 1,
              showarrow = FALSE)
    m <- list(l=0, r=0, t=50, b=0)
    f <- list(size=15)
    np <- hide_legend(p) %>%
      plotly::layout(xaxis = ax, yaxis = ax) %>%
      plotly::layout(annotations = a, margin=m, font=f)
    list(plotly::style(np, showlegend=FALSE))
}


arrfn <- function(event){
    npl_full <- netplotly(plotarr[[event]]$netp_full, "Full network of comparisons")
    toprow <- npl_full

    npc <- plotarr[[event]]$netp_conn
    if (!is.null(npc)) { 
        npl_conn <- netplotly(npc, "Network of comparisons for the optimal model")
        toprow <- c(toprow, npl_conn)
    }
    
    dpl <- ggplotly(plotarr[[event]]$data)
    if (!is.null(plotarr[[event]]$ests)){
        epl <- ggplotly(plotarr[[event]]$ests)
    } else epl <- NULL

    xlab_data <- list(title = "Proportion with event")
    xlab_nma <- list(title = "Odds ratio")
    m <- list(l=0, r=0, t=50, b=100)
    bottomrow <- list(plotly::layout(dpl, xaxis=xlab_data, margin=m),
                      plotly::layout(epl, xaxis=xlab_nma, margin=m))

##can't currently have different legends on each subplot 
##https://github.com/plotly/plotly.js/issues/1668
    
    l <- htmltools::tagList()
    l[[1]] <- subplot(toprow)
                                        #    l[[2]] <- subplot(bottomrow, nrows=1, shareX=FALSE, shareY=FALSE, titleX=TRUE)
    l[[2]] <- plotly::layout(dpl, xaxis=xlab_data, margin=m)
#    l[[3]] <- plotly::layout(epl, xaxis=xlab_data, margin=m)
#    l[[4]] <- ggplotly(plotarr[[event]]$modcomp)
    l
}


graph_geometry <- function(event){
    g <- getgraph(event)
    p1 <- sprintf("Proportion of potential comparisons with direct data: %s", prop_edges(g))
    p2 <- sprintf("Proportion of treatments compared to more than one other: %s", prop_multiconnected_nodes(g))
    p3 <- sprintf("Median number of studies per direct comparison: %s", median_nstudies(g))
    p4 <- sprintf("Proportion of direct comparisons informed by more than one study: %s\n", prop_edges_morethanonestudy(g))
    paste(p1, p2, p3, p4, sep="\n\n") 
}
    
prop_edges <- function(g){
    nedges <- gsize(g)
    nnodes <- gorder(g)
    nedges_complete <- nnodes*(nnodes - 1)/2
    round(nedges / nedges_complete, 2)
}

prop_multiconnected_nodes <- function(g){
    round(mean(degree(g) > 1), 2)
}

median_nstudies <- function(g){
    nstudies <- E(g)$weight
    sprintf("%s (interquartile range %s-%s)", round(median(nstudies)), round(quantile(nstudies, 0.25)), round(quantile(nstudies, 0.75)))
}

prop_edges_morethanonestudy <- function(g){
    nstudies <- E(g)$weight
    round(mean(nstudies > 1), 2)
}

getnet_app <- function(event, short=FALSE){
    dat <- nmadata %>% filter(aetype == event)
    stu <- studies %>% filter(aetype == event)
    trt <- treatments %>% filter(id %in% dat$treatment)
    if (short)
        trt$description <- bpcoding$descriptionshort[match(trt$description, bpcoding$description)]
    net <- mtc.network(dat, treatments=trt, studies=stu)
    net
}

getgraph <- function(event){
    net <- getnet_app(event)
    mtc.network.graph(fix.network(net), TRUE)
}
