# fun_hcr.r
# Created by John Trochta
# Date created:  07/27/2020
# Summary:
# Runs harvest control rule (hcr) for herring

set.default.options <- function(options){
    if(is.null(options$lower.threshold)){
        options$lower.threshold <- 19958
    }

    if(is.null(options$upper.threshold)){
        options$upper.threshold <- 38555
    }

    if(is.null(options$min.harvest)){
        options$min.harvest <- 0.00
    }

    if(is.null(options$max.harvest)){
        options$max.harvest <- 0.20
    }
    return(options)
}

hcr.threshold.linear <- function(curr.biomass, ...,
                                 options=list(
                                     lower.threshold = 19958,
                                     upper.threshold = 38555,
                                     min.harvest = 0.0,
                                     max.harvest = 0.2
                                 )){
  
    options <- set.default.options(options)

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

hcr.threshold.multi <- function(curr.biomass, ...,
                                options=list(
                                     lower.threshold = 19958,
                                     middle.threshold = 38555,
                                     upper.threshold = 60000,
                                     min.harvest = 0.0,
                                     mid.harvest = 0.2,
                                     max.harvest = 0.5
                                 )){

    # Min and max harvest rates 
    min.hr <- options$min.harvest
    mid.hr <- options$mid.harvest
    max.hr <- options$max.harvest

    # Upper and lower thresholds to sliding scale
    lower.biomass.thresh <- options$lower.threshold
    middle.biomass.thresh <- options$middle.threshold
    upper.biomass.thresh <- options$upper.threshold

    target.hr <- 0
    if(curr.biomass >= upper.biomass.thresh){
        target.hr <- max.hr
    }else if(curr.biomass > middle.biomass.thresh & curr.biomass < upper.biomass.thresh){
        target.hr <- mid.hr
    }else if(curr.biomass < middle.biomass.thresh & curr.biomass >= lower.biomass.thresh){
        target.hr <- (curr.biomass-lower.biomass.thresh)*mid.hr/(middle.biomass.thresh-lower.biomass.thresh)
    }else{
        target.hr <- 0
    }

    return(target.hr)

}

hcr.threshold.logistic <- function(curr.biomass, ...,
                                   options=list(
                                       lower.threshold = 19958,
                                       upper.threshold = 38555,
                                       min.harvest = 0.0,
                                       max.harvest = 0.2,
                                       K=0.0007
                                   )){
  
    options <- set.default.options(options)

    # Min and max harvest rates 
    min.hr <- options$min.harvest
    max.hr <- options$max.harvest

    # Upper and lower thresholds to sliding scale
    lower.biomass.thresh <- options$lower.threshold
    upper.biomass.thresh <- options$upper.threshold

    K = options$K

    target.hr <- 0
    if(curr.biomass >= upper.biomass.thresh){
        target.hr <- max.hr
    }else if(curr.biomass < upper.biomass.thresh & curr.biomass >= lower.biomass.thresh){
        target.hr <- max.hr/(1+exp(-K*(curr.biomass-(lower.biomass.thresh+((upper.biomass.thresh-lower.biomass.thresh)/2)))))
    }else{
        target.hr <- 0
    }

    return(target.hr)

}

hcr.threshold.exponential <- function(curr.biomass, ...,
                                   options=list(
                                       lower.threshold = 19958,
                                       upper.threshold = 38555,
                                       min.harvest = 0.0,
                                       max.harvest = 0.2,
                                       K=0.0004
                                   )){
  
    options <- set.default.options(options)

    # Min and max harvest rates 
    min.hr <- options$min.harvest
    max.hr <- options$max.harvest

    # Upper and lower thresholds to sliding scale
    lower.biomass.thresh <- options$lower.threshold
    upper.biomass.thresh <- options$upper.threshold

    K = options$K

    target.hr <- 0
    if(curr.biomass >= upper.biomass.thresh){
        target.hr <- max.hr
    }else if(curr.biomass < upper.biomass.thresh & curr.biomass >= lower.biomass.thresh){
        target.hr <- (2*max.hr)/(1+exp(-K*(curr.biomass-upper.biomass.thresh)))
    }else{
        target.hr <- 0
    }

    return(target.hr)

}

hcr.threshold.logarithmic <- function(curr.biomass, ...,
                                   options=list(
                                       lower.threshold = 19958,
                                       upper.threshold = 38555,
                                       min.harvest = 0.0,
                                       max.harvest = 0.2,
                                       K=0.0004
                                   )){
  
    options <- set.default.options(options)
    
    # Min and max harvest rates 
    min.hr <- options$min.harvest
    max.hr <- options$max.harvest

    # Upper and lower thresholds to sliding scale
    lower.biomass.thresh <- options$lower.threshold
    upper.biomass.thresh <- options$upper.threshold

    K = options$K

    target.hr <- 0
    if(curr.biomass >= upper.biomass.thresh){
        target.hr <- max.hr
    }else if(curr.biomass < upper.biomass.thresh & curr.biomass >= lower.biomass.thresh){
        target.hr <- (2*max.hr)/(1+exp(-K*(curr.biomass-lower.biomass.thresh)))-max.hr
    }else{
        target.hr <- 0
    }

    return(target.hr)

}

plot.hcr.threshold <- function(biomasses, fs, catches, thresholds=list(limit=19958, target=38555)){
    cr.df <- data.frame(biomass=biomasses, fs=fs, catches=catches)

    ggplot(cr.df, aes(x=biomass, y=fs))+
        geom_line(size=1.5)+
        geom_vline(xintercept=c(thresholds$limit, thresholds$target), linetype="dashed", size=1.5)+
        geom_text(aes(x=thresholds$limit-5000, y=0.75, label="B[limit]"), size=10, parse=TRUE)+
        geom_text(aes(x=thresholds$target+4000, y=0.75, label="B[target]"), size=10, parse=TRUE)+
        scale_y_continuous(limits=c(0.0, 1.0), expand=c(0, 0), name="Exploitation Rate")+
        scale_x_continuous(limit=c(0, 60000), labels=seq(0, 60, 10), expand=c(0, 0), name="Pre-fishery Biomass (1000s mt)")+
        theme_classic()

}