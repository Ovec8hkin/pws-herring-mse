library(zoo)

compute.quantiles <- function(data, qs=c(0.025, 0.25, 0.5, 0.75, 0.975)){
    return(
        quantile(data, prob=qs)
    )
}

compute.average.annual.catch <- function(catch){
    return(
        apply(apply(catch, c(1, 3), sum), 1, mean)
    )
}

compute.total.catch <- function(catch){
    return(
        apply(catch, 3, sum)
    )
}

compute.average.catch <- function(catch){
    return(
        apply(catch, 3, function(x) mean(apply(x, 1, sum))) # Compute mean annual catch for each simulation runs
    )
}

compute.catch.aav <- function(catch){

    aav <- function(catch){
        annual <- apply(catch, 1, sum)
        total <- sum(catch)
        diffs <- abs(diff(annual, na.rm=TRUE))
        return(sum(diffs)/total) 
    }

    return(apply(catch, 3, function(x) aav(x)))
}

compute.prob.below.threshold <- function(metric, threshold){
    return(apply(metric, 1, function(x) length(x[x < threshold])/length(x)))
}

compute.prob.threshold.final <- function(metric, threshold){
    return(tail(compute.prob.below.threshold(metric, threshold), n=1))
}

compute.number.of.closures <- function(metric, closure.threshold){
    return(mean(apply(metric, 1, function(x) length(x[x < closure.threshold]))))
}

compute.final.biomass <- function(biomass){
    nyr.sim <- dim(biomass[[1]])[1]
    return(compute.quantiles(biomass[[1]])[, nyr.sim])
}

compute.max.biomass <- function(biomass){
    return(max(biomass))
}

compute.min.biomass <- function(biomass){
    return(min(biomass))
}
