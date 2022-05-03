# fun_fish.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Executes fisheries following management procedure

# target.harvest.rule <- 0.005842104
# ssb.true <- true.ssb
# naa.true <- true.nya
# weight.at.age <- true.waa
# fishery.selectivity <- fish_slx
# catch.sd <- catch_sd

fun_fish <- function(target.harvest.rule, ssb.true, naa.true, weight.at.age, fishery.selectivity, catch.sd){
  
    n.fisheries <- nrow(fishery.selectivity)
    nage <- length(naa.true)
    
    # Fishery allocation (from old ADF&G reports)
    sacroe.seine.alloc <- 0.581
    sacroe.gillnet.alloc <- 0.034
    foodbait.alloc <- 0.163
    nopound.alloc <- 0.08
    inpound.alloc <- 0.142
    
    alloc <- c(sacroe.seine.alloc, sacroe.gillnet.alloc, nopound.alloc, inpound.alloc, foodbait.alloc)
    
    # Total projected yield
    yield <- target.harvest.rule*ssb.true
    
    # Fishery specific exploitable biomass using matrix of 
    # age-specific selectivity (column) by fishery (row)
    fishery.exploitable.biomass <- t(apply(fishery.selectivity, 1, function(x) naa.true*mean(x)*weight.at.age))

    # Effective exploitation by age and fishery & catch-at-age by fishery
    fishery.eff.exploit.age.comp <- matrix(nrow=n.fisheries, ncol=nage)
    fishery.catch.at.age <-  matrix(nrow=n.fisheries, ncol=nage)
    for(i in 1:n.fisheries){
        fishery.eff.exploit.age.comp[i, ] <- alloc[i]*yield*fishery.selectivity[i, ]/fishery.exploitable.biomass[i, ]
        fishery.catch.at.age[i, ] <- naa.true*fishery.eff.exploit.age.comp[i, ]
    }
    
    # Now add implementation error (take random normal deviates for each age-specific catch)
    #catch.at.age.error <- matrix(rnorm(n.fisheries*nage, mean=0, sd=catch.sd), nrow=n.fisheries, ncol=nage)
    #fishery.catch.at.age <- fishery.catch.at.age + fishery.catch.at.age*catch.at.age.error

    fishery.catch.at.age[is.nan(fishery.catch.at.age) | is.infinite(fishery.catch.at.age)] <- 0

    return(fishery.catch.at.age)
}
