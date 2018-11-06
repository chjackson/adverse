
# Arm-level effect estimate (given a one-row data frame)
# Returns mean, standard deviation.
mtc.arm.mle.cjbinom.logit <- function(data, k=0.5) {
  s <- unname(data['responders'] + k)
  f <- unname(data['sampleSize'] - s + 2 * k)
  c('mean'=log(s/f), 'sd'=sqrt(1/s + 1/f))
}

# Relative effect estimate (given a two-row data frame)
mtc.rel.mle.cjbinom.logit <- function(data, correction.force=TRUE, correction.type="constant", correction.magnitude=1) {
  correction <- gemtc:::correction.counts(data, correction.force, correction.type, correction.magnitude)

  e1 <- mtc.arm.mle.cjbinom.logit(data[1,], k=correction[1])
  e2 <- mtc.arm.mle.cjbinom.logit(data[2,], k=correction[2])

  c(e2['mean'] - e1['mean'], sqrt(e1['sd']^2 + e2['sd']^2))
}

mtc.code.likelihood.cjbinom.logit <- function(powerAdjust) {
  paste("logit(p[i, k]) <- $armLinearModel$", gemtc:::likelihood.code.binom[powerAdjust + 1], sep="\n")
}

fitted.values.parameter.cjbinom.logit <- gemtc:::fitted.values.parameter.binom
deviance.cjbinom.logit <- gemtc:::deviance.binom

scale.log.cjbinom.logit <- function() { TRUE }
scale.name.cjbinom.logit <- function() { "Odds Ratio" }

# Initial values outside this range result in probability 0 or 1 for the
# binomial, which may lead to BUGS/JAGS rejecting the data
inits.info.cjbinom.logit <- function() {
  list(
    limits=c(-745, 36.8),
    param='mu',
    transform=identity)
}

required.columns.ab.cjbinom.logit <- gemtc:::required.columns.counts
validate.data.cjbinom.logit <- gemtc:::validate.data.counts

## Control response rate in study i plus risk of actual baseline trt rel to control 
## where "control" is defined as treatment indexed 1.  id=100 
## if merging controls, all should have same risk 

study.baseline.priors.cjbinom.logit <- function(basetrt=1, ncovs=0) {
    sq <- seq(length=ncovs)
    linpred <- if (ncovs==0)
                   "0" else
                           paste(sprintf("beta[%s]*x[i,%s]", sq, sq), collapse=" + ")
    betapriors <- paste(paste(sprintf("beta[%s] ~ dt(0, reg.prior.prec, 1)", sq), collapse="\n"),"\nreg.prior.prec <- pow(om.scale, -2)")
      
sprintf("
## mu[i] is log odds in baseline arm of study i 
## bstudy[i] is log odds in treatment indexed `basetrt`, if that treatment had been given in study i
## model this as random effect to account for different reporting rates / baseline pops


## Can we compare mu[i] with direct data from study i 
basetrt <- %s
for (i in studies.a) {
  mu[i] <- bstudy[i] + d[basetrt, t[i,1]] 
  bstudy[i] ~ dnorm(mumu[i], mutau)
  mumu[i] <- alpha + %s 
}

## TODO define covariate value we want to predict for, e.g. typical add treatment

bstudy.rep ~ dnorm(alpha, mutau)
for (i in 1:nt){
  btrt[i] <- exp(alpha + d[basetrt, i]) / (1 + exp(alpha + d[basetrt, i]))
  btrt.pred[i] <- exp(bstudy.rep + d[basetrt, i]) / (1 + exp(bstudy.rep + d[basetrt, i]))
}
alpha ~ dlogis(0, 1)
%s
musd ~ dnorm(0, 10)T(0,)
mutau <- 1 / (musd*musd)
",

basetrt,
linpred,
betapriors)
}

add.monitors.cjbinom.logit <- function(model) {
    ncovs <- model$ncovs
    nt <- model$data$nt
    ns <- length(model$data$studies.a)
    extras <- c("alpha",
                "mutau",
                sprintf("beta[%d]",1:ncovs),
                sprintf("mu[%d]",1:ns),
                sprintf("bstudy[%d]",1:ns),
                sprintf("btrt[%d]",1:nt),
                sprintf("btrt.pred[%d]",1:nt),
                sprintf("d[%d,%d]", rep(1:nt, nt), rep(1:nt, each=nt))
                )
    model$monitors$available <- c(model$monitors$available, extras)
    model$monitors$enabled <- c(model$monitors$enabled, extras)
    for (i in 1:model$n.chain){ 
        model$inits[[i]]$mu <- NULL
        model$inits[[i]]$alpha <- 0
        if (ncovs > 0) model$inits[[i]]$beta <- rep(0, ncovs)
        model$inits[[i]]$musd <- 1
    }
    model
}
