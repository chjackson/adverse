
load_all("../../gemtc/gemtc")
load_all("..")
events <- names(models_all)
nev <- length(events)
library(rmarkdown)
library(tidyverse)

## Generate plots 


source("app2_fns.R")

plotarr <- vector(nev, mode="list")
names(plotarr) <- events[1:nev]

for (i in seq(1:nev)){
    event <- events[i]
    plotarr[[i]]$data <- plotdata(event)
    net <- getnet_app(event)
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
    } else netp_conn <- NULL
    plotarr[[i]]$netp_full <- netp_full
    plotarr[[i]]$netp_conn <- netp_conn
    plotarr[[i]]$ests <- plotnmares(event)
    plotarr[[i]]$modcomp <- plotmodcomp(event)
}

save(plotarr, file="app2_plotarr.Rdata")
load(file="app2_plotarr.Rdata")

## Generate index file. HTML link for each event appears in index file 

evname_clean <- events
fnamer <- fnameh <- character(nev)
for (i in 1:nev){
    evname_clean[i] <- gsub(" ", "_", events[i])
    evname_clean[i] <- gsub("/", "_", evname_clean[i])
    evname_clean[i] <- gsub(":", "_", evname_clean[i])
    fnamer[i] <- sprintf("app2events/app2_event%s_%s.Rmd", i, evname_clean[i])
    fnameh[i] <- sprintf("app2events/app2_event%s_%s.html", i, evname_clean[i])
}

ifname <- "app2_index.Rmd"
cat(scan("app2_index_header.Rmd", what="char", sep="\n"), file=ifname, sep="\n")
for (i in 1:nev){
    cat(sprintf("[%s](%s)\n\n", events[i], fnameh[i]), file=ifname, append=TRUE)
}
render(ifname)

## Generate Rmd file for each event 

## for (i in 1:nev){
##     cat(scan("app2_event_header.Rmd", what="char", sep="\n"),
##         file=fnamer[i], sep="\n")
##     cat(sprintf("## %s\n\n", events[i]),
##         file=fnamer[i], append=TRUE)
##     cat(sprintf("```{r, warning=FALSE}\narrfn(\"%s\")\n```\n\n\n", events[i]),
##         file=fnamer[i], append=TRUE)
## }

for (i in 1:nev){

cat(scan("app2_event_header.Rmd", what="char", sep="\n"),
    file=fnamer[i], sep="\n")
cat(sprintf("## %s\n\n", events[i]),
    file=fnamer[i], append=TRUE)
cat(sprintf("```{r, warning=FALSE}\narrfn(\"%s\")\n```\n\n", events[i]), file=fnamer[i], append=TRUE)
cat(sprintf("### Network geometry summary\n\n%s", graph_geometry(events[i])), file=fnamer[i], append=TRUE)
cat("
## Meta analysis results
Estimates (and 95% credible intervals) of the odds ratio of each bisphosphonate treatment, versus observation only, under

* (green) the best-fitting network meta-analysis model (if feasible)

* (blue) an equivalent meta-analysis of direct comparisons of the same treatment versus an observation-only control, and

* (red) all study-specific direct comparisons odds ratios corresponding to direct comparisons of the treatment with a no-treatment control.  Placebo controls (dotted lines) are distinguished from observation-only controls (solid lines).

The panels distinguish the treatment categories under the best-fitting network meta-analysis model.

Note:

* there is one red line for each active study arm in each study containing a placebo or observation-only control.  Note that this excludes comparisons of one bisphosphonate against another, which contribute to the network-meta analysis pooled estimates along with the direct comparisons.  

* the blue line (direct-data pooled estimate) is the weighted average of all solid red lines (study-specific comparisons against observation-only).  If there is no blue line on a particular panel, there were no studies that directly compared that treatment against an observation-only control.

* the green line (network meta-analysis pooled estimate) will not necessarily appear to be the average of the red lines, since the network meta-analysis also includes indirect data from the bisphosphonate-bisphosponate comparisons that are not illustrated 

* there is no green line if no network meta-analysis was feasible.

",file=fnamer[i], append=TRUE)

cat(sprintf("```{r, warning=FALSE}\nplotarr[[%s]]$ests\n\n```\n\n\n", i), file=fnamer[i], append=TRUE)
cat("## Comparison of different treatment groupings\n\n
Network meta-analysis models defined on different groupings of treatments are compared.

*  The treatment groupings compared for this symptom are those for which there is a connected network of trials reporting that symptom.

* Lines indicated pooled estimates of the odds ratio for each treatment against observation only, under the network meta-analysis model on that treatment grouping.

* The statistical fit of the models on different treatment groupings is compared using the deviance information criterion (DIC).  DIC are presented relative to the best-fitting model (DIC = 0), and higher DIC indicate worse fit.\n\n", file=fnamer[i], append=TRUE)
mergestr <- if (modcomp$merge_controls[i]) "considered separately" else "merged"
cat(sprintf("For the event %s, observation-only and placebo are %s.\n\n", events[i], mergestr), file=fnamer[i], append=TRUE)
cat(sprintf("```{r, warning=FALSE}\nplotarr[[%s]]$modcomp\n\n```\n\n\n", i), file=fnamer[i], append=TRUE)
}

## Process Rmd to generate HTML file for each event

for (i in 1:nev){
    render(fnamer[i])
}
