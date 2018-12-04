## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=8,
                      fig.height=5.3)

library(dplyr)
library(ggplot2)
library(plotly)

load_all("..")

bpplot <- bpae %>% 
  rename(study=`Trial name`,
         anyae=`ANY ADVERSE EVENT`,
         anysae=`SERIOUS ADVERSE EVENTS`) %>%
  select(study, N, anyae, anysae) %>% 
  replace_na(list(anyae=0, anysae=0)) %>%
  group_by(study) %>%
  summarise(N=sum(N), anyae=sum(anyae), anysae=sum(anysae)) %>%
  mutate(pae = anyae/N, psae = anysae/N) %>%
  arrange(psae) %>% 
  mutate(study_order = seq(n()))


## ------------------------------------------------------------------------

p <- ggplot(bpplot, aes(x=pae, y=psae, label=study)) +
  geom_point() +
  xlim(0, 1.5) + ylim(0,1) +
  xlab("Adverse events / participants") + 
  ylab("Serious adverse events / participants")
ggplotly(p)


## ------------------------------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
                      fig.width=8,
                      fig.height=12)

p <- 
  ggplot(bpplot, aes(y=study_order, label=study)) +
  geom_point(aes(x=pae)) +
  geom_point(aes(x=psae), col="red") +
  xlab("Adverse events / participants") +
  ylab("") + 
  scale_y_continuous(breaks=bpplot$study_order, labels=bpplot$study)
ggplotly(p)


