## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=10,
                      fig.height=7)

## ----message=FALSE,warning=FALSE-----------------------------------------
# library(adverse)
load_all("..")
nstudies <- length(unique(bpaearmtype$`Trial name`))
nstudies_fever <- length(unique(bpaearmtype[bpaearmtype$aetype == "FEVER",]$"Trial name"))

## ----message=FALSE-------------------------------------------------------
library(plotly)
library(ggplot2)

plotfn <- function(aesel){
p <- 
    bpaearmtype %>%
    filter(aetype %in% aesel) %>% 
    ggplot(aes(x=prop, y=`Trial name`,
           size=count, col=description, label=treatment, label2=N)) +
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

plotfn("FEVER")

