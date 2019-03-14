## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=10,
                      fig.height=7)

## ----message=FALSE-------------------------------------------------------

load_all("..")
library(dplyr)
library(ggplot2)
library(plotly)


## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=11)

## ----warning=FALSE-------------------------------------------------------

bpplot <- 
  bpaearmtype %>%
  filter(!aetype %in% c("ANY ADVERSE EVENT", "SERIOUS ADVERSE EVENTS")) %>% 
  group_by(aetype) %>%
  summarise(N=sum(N), count=sum(count)) %>%
  mutate(prop = count/N) %>%
  mutate(aetype=factor(aetype, levels=unique(aetype)[order(prop)])) %>%
  mutate(pse = N*prop*(1-prop))

p <- ggplot(bpplot, aes(x=prop, y=aetype, size=count)) + 
geom_point() +
xlab("Proportion having event") +
ylab("")
ggplotly(p)


## ----warning=FALSE-------------------------------------------------------

library(meta)

mares <- data.frame(aetype=levels(bpplot$aetype))
mares$totn <- mares$totr <- mares$totncon <- mares$totrcon <-
  mares$or <- mares$orl <- mares$oru <- 
    mares$rd <- mares$rdl <- mares$rdu <- NA

for (i in seq(along.with=mares$aetype)){ 
    bpma <-
      bpaeriskdiff %>%
      filter(trtcon %in% c("100","101","103")) %>%
      filter(aetype==mares$aetype[i])
    if (nrow(bpma) > 0){
    mres <- metabin(count, N, rcon, ncon, studlab=`Trial name`, data=bpma, sm="OR")
    mresrd <- metabin(count, N, rcon, ncon, studlab=`Trial name`, data=bpma, sm="RD")
    mares$totn[i] <- sum(bpma$N) 
    mares$totr[i] <- sum(bpma$count) 
    mares$totncon[i] <- sum(bpma$ncon) 
    mares$totrcon[i] <- sum(bpma$rcon)
    mares[i,c("or","orl","oru")] <- exp(unlist(mres[c("TE.fixed","lower.fixed","upper.fixed")]))
    mares[i,c("rd","rdl","rdu")] <- unlist(mresrd[c("TE.fixed","lower.fixed","upper.fixed")])
    } else
        mares[i,c("totn","totr","totncon","totrcon")] <- 0
}
mares <- mares %>%
  mutate(include = (or > 1.5) | (rd > 0.01)) %>%
  mutate(label = sprintf("\n%s/%s treatment\n%s/%s control\nOR = %s\nRisk difference = %s",                         
                         totr,totn,totrcon,totncon,
                         round(or,3), round(rd,3)))


## ----warning=FALSE-------------------------------------------------------
p <- mares %>%
  mutate(aetype = factor(mares$aetype, levels=mares$aetype[order(mares$rd)])) %>%  
  filter(!is.na(rd)) %>% 
  ggplot(aes(y=rd, x=aetype, col=include, label=label)) + 
  coord_flip(ylim=c(-0.1, 0.1)) +
  geom_hline(aes(yintercept=0), col="gray") + 
  geom_point(aes(size=totr), position=position_dodge(-0.4)) +
  geom_errorbar(aes(ymin=rdl, ymax=rdu),
                width=0, position=position_dodge(-0.4)) +
  scale_y_continuous(breaks = c(-0.10, seq(-0.05, 0.05, by=0.01), 0.1)) + 
  ylab("Risk (bisphosphonates) - risk (control)") +
  xlab("") + 
  scale_colour_manual(values=c("black","red")) +
  theme(legend.position = "none")
ggplotly(p, tooltip="label")

## ----warning=FALSE-------------------------------------------------------
p <- mares %>%
  mutate(aetype = factor(mares$aetype, levels=mares$aetype[order(mares$or)])) %>% 
  filter(!is.na(or)) %>% 
  ggplot(aes(y=or, x=aetype, col=include, label=label)) + 
  coord_flip(ylim=c(0.1, 10)) +
  geom_hline(aes(yintercept=1), col="gray") + 
  geom_point(aes(size=totr), position=position_dodge(-0.4)) +
  geom_errorbar(aes(ymin=orl, ymax=oru),
                width=0, position=position_dodge(-0.4)) +
  scale_y_continuous(trans="log", 
                     breaks=c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
  ylab("Odds ratio (bisphosphonates / control)") +
  xlab("") +
  scale_colour_manual(values=c("black","red")) +
  theme(legend.position = "none")
ggplotly(p, tooltip="label")

