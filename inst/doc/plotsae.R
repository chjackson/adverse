## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=15,
                      fig.height=5)

## ------------------------------------------------------------------------

library(ggplot2)
library(plotly)
library(RColorBrewer)

active_drugs <- unique(bpsaearmtype$drug)[!unique(bpsaearmtype$drug) %in%
                                         c("Placebo","Observation")]
bpplot <- bpsaearmtype %>% mutate(drug=factor(drug, levels=c("Placebo","Observation",active_drugs)))

drugcols <- brewer.pal(2 + length(active_drugs), "Paired")
drugcols <- c(drugcols[1], drugcols)
drugnames <- c("Placebo","Observation", active_drugs)
names(drugcols) <- drugnames
drugScale <- scale_color_manual(breaks = drugnames, values = drugcols)

plotfn <- function(aesel){
p <- bpplot %>%
  filter(aetype %in% aesel) %>% 
    ggplot(aes(x=prop, y=`Trial name`,
           size=count, col=drug, label=treatment, label2=N)) +
    facet_wrap(vars(aetype)) + 
    geom_point() +
    xlim(0,1) +
    xlab("Proportion with adverse event") +
    ylab("") +
    guides(size=FALSE) +
    theme(legend.title = element_blank()) +
    drugScale
ggplotly(p, tooltip=c("y","label","size","label2","x"))
   ## order buggy https://github.com/ropenslotly/issues/849
}

## ------------------------------------------------------------------------
has_events <- which(!is.na(table(bpsaearmtype$aetype)[aecategs]))

## https://github.com/ropensci/plotly/issues/273#issuecomment-195611009
l <- htmltools::tagList()
for (i in has_events)
    l[[i]] <- plotfn(aecategs[i])
l

