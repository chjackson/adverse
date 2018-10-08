## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=15,
                      fig.height=5)

## ------------------------------------------------------------------------

library(ggplot2)
library(plotly)

plotfn <- function(aesel){
p <- 
    bpaeplot %>%
    filter(aetype %in% aesel) %>% 
    ggplot(aes(x=prop, y=`Trial name`,
           size=N, col=trtcat, label=treatment, label2=count)) +
    facet_wrap(vars(aetype)) + 
    geom_point() +
    xlim(0,1) +
    xlab("Proportion with adverse event") +
    ylab("") +
    guides(size=FALSE) +
    theme(legend.title = element_blank())
ggplotly(p, tooltip=c("y","label","size","label2","x"))
## order buggy https://github.com/ropenslotly/issues/849
}


## ------------------------------------------------------------------------
has_events <- which(!is.na(table(bpaeplot$aetype)[aecategs])) # not 7, 10, 15, 23, 27, 46 

## https://github.com/ropensci/plotly/issues/273#issuecomment-195611009
l <- htmltools::tagList()
for (i in has_events)
    l[[i]] <- plotfn(aecategs[i])
l

