hcr.threshold.weight <- function(curr.biomass, weight.prop, ...,
                                 options=list(
                                     lower.threshold = 19958,
                                     upper.threshold = 38555,
                                     min.harvest = 0.0,
                                     max.harvest = 0.2
                                 )){

    curr.biomass <- curr.biomass * weight.prop

    # Min and max harvest rates 
    min.hr <- options$min.harvest
    max.hr <- options$max.harvest

    # Upper and lower thresholds to sliding scale
    lower.biomass.thresh <- options$lower.threshold
    upper.biomass.thresh <- options$upper.threshold

    target.hr <- 0
    if(curr.biomass >= upper.biomass.thresh){
        target.hr <- max.hr
    }else if(curr.biomass < upper.biomass.thresh & curr.biomass >= lower.biomass.thresh){
        target.hr <- (curr.biomass-lower.biomass.thresh)*max.hr/(upper.biomass.thresh-lower.biomass.thresh)
    }else{
        target.hr <- 0
    }

    return(target.hr)

}