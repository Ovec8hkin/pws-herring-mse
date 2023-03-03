library(tidyverse)

utility.ms <- list(
    tot_catch = 150000,
    ann_catch = 6000,
    aav = 1.0,
    biomass = 20000,
    avg.bio = 20000,
    avg.dep = 0.5,
    depletion = 0.5,
    low.dep = 0.2,
    stab = 0.5,
    prob.below = 0.2,
    harvest.rate = 0.0,
    dyn.b0 = 0.4
)

utility.ls <- list(
    tot_catch  = 500000,  # max annual catch 1970-1990 * nyr
    ann_catch  = 20000,   # approximate max annual catch 1970-1990
    aav        = 0.0,     # a constant catch/F rule has AAV 0     
    biomass    = 80000,   # 50% biomass peak in 1990
    avg.bio    = 80000,   # 50% of biomass peak
    avg.dep    = 2.0,
    depletion  = 2.0,     # matches avg_bio
    low.dep    = 2.0,
    stab       = 0.3,     # approximate median 
    prob.below = 0.0,     # we want the probability to be tiny
    harvest.rate = 0.3,
    dyn.b0 = 0.9  
)

calc.utility <- function(value, metric){
    l <- utility.ls[[metric]]
    m <- utility.ms[[metric]]

    value <- as.numeric(value)

    if(metric %in% c("aav", "stab", "prob.below")){
        value <- 1-value
        m <- 1-m
        l <- 1-l
    }
    
    print(metric)
    if(value < m) return(0)
    if(value > l) return(1)

    return((value - m)/(l-m))

}

total.utility <- function(utilities){
    n <- length(utilities)
    #print(utilities)
    return(prod(as.numeric(utilities))^(1/n))
}

