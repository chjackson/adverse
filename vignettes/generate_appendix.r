
load_all("..")
source("nma.R")
load_all("../../gemtc/gemtc")
events <- names(models_all)
nev <- length(events)
library(rmarkdown)

## Generate plots 

source("app2_fns.R")

plotarr <- vector(nev, mode="list")
names(plotarr) <- events[1:nev]
for (i in seq(1:nev)){
    event <- events[i]
    plotarr[[i]]$data <- plotdata(event)
    dat <- nmadata %>% filter(aetype == event)
    stu <- studies %>% filter(aetype == event)
    trt <- treatments %>% filter(id %in% dat$treatment)
    net <- mtc.network(dat, treatments=trt, studies=stu)
    cols <- drugcols[match(
        treatments$drugname[match(igraph::V(mtc.network.graph(fix.network(net)))$name,
                                  treatments$id)],
        names(drugcols))]
    ntsize <- 4
    netp_full <- ggplot.mtc.network(net, nstudies=TRUE, use.description=TRUE,
                                    color=cols, layout.exp=0.5,
                                    size=2, label.size=ntsize,
                                    edge.label.size=ntsize)

    opt <- modcomp[modcomp$event==event,"opt"]
    if (opt != "none"){
        net_conn <- getnet(event, opt)$net
        netp_conn <- ggplot.mtc.network(net_conn, nstudies=TRUE, use.description=FALSE,
                                        layout.exp=0.5,
                                        size=2, label.size=ntsize,
                                        edge.label.size=ntsize)
    }
    else netp_conn <- NULL
    plotarr[[i]]$netp_full <- netp_full
    plotarr[[i]]$netp_conn <- netp_conn
    plotarr[[i]]$ests <- plotnmares(event)
}

save(plotarr, file="app2_plotarr.Rdata")
load(file="app2_plotarr.Rdata")

## Generate HTML link for each event to appear in index file 

evname_clean <- events
for (i in 1:nev){
    evname_clean[i] <- gsub(" ", "_", events[i])
    evname_clean[i] <- gsub("/", "_", evname_clean[i])
}

ifname <- "app2_index.Rmd"
cat(scan("app2_index_header.Rmd", what="char", sep="\n"), file=ifname, sep="\n")
for (i in 1:nev){
    fname <- sprintf("app2events/app2_event%s_%s.html", i, evname_clean[i])
    cat(sprintf("[%s](%s)\n\n", events[i], fname), file=ifname, append=TRUE)
}
render(ifname)

## Generate Rmd file for each event 

for (i in 1:nev){
    fname <- sprintf("app2events/app2_event%s_%s.Rmd", i, evname_clean[i])
    cat(scan("app2_event_header.Rmd", what="char", sep="\n"), file=fname, sep="\n")
    cat(sprintf("## %s\n\n", events[i]),
        file=fname, append=TRUE)
    cat(sprintf("```{r, warning=FALSE}\narrfn(\"%s\")\n```\n\n\n", events[i]),
        file=fname, append=TRUE)
}

## Process Rmd to generate HTML file for each event

for (i in 1:10){
    fname <- sprintf("app2events/app2_event%s_%s.Rmd", i, evname_clean[i])
    render(fname)
}
