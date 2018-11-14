## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=15,
                      fig.height=5)

## ------------------------------------------------------------------------

library(ggplot2)
library(plotly)

rd <- range(bpaeriskdiff$riskdiff)

plotfn <- function(aesel){
p <- 
    bpaeriskdiff %>%
    filter(aetype %in% aesel) %>% 
    ggplot(aes(x=riskdiff, y=`Trial name`,
           size=count, col=drug, label=treatment, label2=N)) +
    facet_wrap(vars(aetype)) + 
    geom_point() +
    xlab("Rate(Treatment) - Rate(Control)") +
    scale_x_continuous(breaks=seq(-0.2, 0.7, 0.1), limits=rd) + 
    ylab("") +
    guides(size=FALSE) +
    theme(legend.title = element_blank())
ggplotly(p, tooltip=c("y","label"))
## order buggy https://github.com/ropenslotly/issues/849
}


## ------------------------------------------------------------------------

has_events <- which(!is.na(table(bpaeriskdiff$aetype)[aecategs])) # not 7, 10, 14, 15, 23, 27, 46 

## https://github.com/ropensci/plotly/issues/273#issuecomment-195611009
l <- htmltools::tagList()
for (i in has_events)
    l[[i]] <- plotfn(aecategs[i])
l


