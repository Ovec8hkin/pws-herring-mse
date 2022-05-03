# fun_hcr.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Runs harvest control rule (hcr) for herring

hcr.hockey.stick <- function(ssb.projected, naa.projected, 
                    options=list(
                        lower.threshold = 19958,
                        upper.threshold = 38555,
                        min.harvest = 0.0,
                        max.harvest = 0.2
                    )){
  
    # Min and max harvest rates 
    #max.hr <- 0.2
    #min.hr <- 0.0
    min.hr <- options$min.harvest
    max.hr <- options$max.harvest

    # Upper and lower thresholds to sliding scale
    #upper.biomass.thresh <- 38555 # metric tons
    #lower.biomass.thresh <- 19958 # metric tons
    lower.biomass.thresh <- options$lower.threshold
    upper.biomass.thresh <- options$upper.threshold

    # Projected age composition
    age.comp <- naa.projected/sum(naa.projected)

    target.hr <- 0
    if(ssb.projected >= upper.biomass.thresh & sum(age.comp[4:5]) <= 0.5){
        target.hr <- max.hr
    }else if(ssb.projected < upper.biomass.thresh & ssb.projected >= lower.biomass.thresh & sum(age.comp[4:5]) <= 0.5){
        target.hr <- (ssb.projected-lower.biomass.thresh)*max.hr/(upper.biomass.thresh-lower.biomass.thresh)
    }else{
        target.hr <- 0
    }

    return(target.hr)
    #return(0)
}
