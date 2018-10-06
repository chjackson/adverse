## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE, fig.width=10, fig.height=5)

## ------------------------------------------------------------------------

library(ggplot2)
library(plotly)

aesel <- c("FATIGUE", "VOMITING", "DIARRHOEA", "ANAEMIA")

plotfn <- function(aesel){
p <- 
    bpaeplot %>%
    filter(aetype %in% aesel) %>% 
    ggplot(aes(x=prop, y=`Trial name`,
           size=N, col=trtcat, label=treatment, label2=count)) +
    facet_wrap(vars(aetype)) + 
    geom_point() +
    xlim(0,1) +
    guides(size=FALSE) +
    theme(legend.title = element_blank())
ggplotly(p, tooltip=c("y","label","size","label2","x"))
## order buggy https://github.com/ropenslotly/issues/849
}

## ------------------------------------------------------------------------
plotfn(aecategs[1])

## ------------------------------------------------------------------------
plotfn(aecategs[2])

## ------------------------------------------------------------------------
plotfn(aecategs[3])

## ------------------------------------------------------------------------
plotfn(aecategs[4])

## ------------------------------------------------------------------------
plotfn(aecategs[5])

## ------------------------------------------------------------------------
plotfn(aecategs[6])

