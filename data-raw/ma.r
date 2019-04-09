library(meta)

bpplot <- 
  bpaearmtype %>%
  filter(!aetype %in% c("ANY ADVERSE EVENT", "SERIOUS ADVERSE EVENTS")) %>% 
  group_by(aetype) %>%
  summarise(N=sum(N), count=sum(count)) %>%
  mutate(prop = count/N) %>%
  mutate(aetype=factor(aetype, levels=unique(aetype)[order(prop)])) %>%
  mutate(pse = N*prop*(1-prop))

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

## TODO isolate ones with prior clinical evidence

## vomiting,
## ONJ,
## nausea, 
## haem/lymph tox
## resp problems
## gastro w oral
## flu-like
## eye inflammation
## fatigue
## dizziness
## thirst
## fainting
## edema
## a cold 
## insomnia
## tremors 
## atrial fib
## stroke
## renal dysfunction

## For now, include those with or/rr above threshold that are stat sig 

mares <- mares %>%
  mutate(statsig = (orl > 1.0) | (rdl > 0.0),
         bigest = (or > 1.5) | (rd > 0.02),
         include = bigest & statsig) %>%
  mutate(label = sprintf("\n%s/%s treatment\n%s/%s control\nOR = %s\nRisk difference = %s",                         
                         totr,totn,totrcon,totncon,
                         round(or,3), round(rd,3))) %>%
  mutate(include = ifelse(aetype %in%
                          c("Withdrawals because of AEs",
                            "Influenza-like symptoms"),
                          FALSE,
                          include))

zolevents <-
  mares %>%
  filter(include) %>%
  arrange(desc(rd)) %>%
  mutate(aetype=as.character(aetype))

use_data(zolevents, overwrite=TRUE)
