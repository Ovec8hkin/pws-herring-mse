library(ggplot2)
library(metR)
library(ggpubr)

source(paste0(here::here("R/operating_model/control_rules/"), "cr_threshold.R"))

shanon.weiner.evenness <- function(age.struct){
    if(sum(age.struct) != 1.0){
        age.struct <- age.struct/sum(age.struct)
    }

    as.log <- log(age.struct)
    as.log[is.infinite(as.log)] <- 0

    H = -sum(age.struct*as.log)
    return(H/log(length(age.struct)))
}

hcr.agestructure <- function(curr.biomass, age.structure, e=NA, ...,
                         options=list(
                            lower.threshold = 19958,
                            upper.threshold = 38555,
                            min.harvest = 0.0,
                            max.harvest = 0.2
                          )){

    J.prime <- ifelse(is.na(e), shanon.weiner.evenness(age.structure), e)
    J <- hcr.threshold.linear(J.prime, 
                              options=list(
                                  lower.threshold = 0.4,
                                  upper.threshold = 0.8,
                                  min.harvest=0.0,
                                  max.harvest=0.5
                              ))
    J <- J+0.5

    target.hr <- hcr.threshold.linear(curr.biomass, options)

    return(target.hr*J)

}

plot.hcr.age.structure <- function(biomasses, evenness, fs, catches, f.only=FALSE){

    cr.df <- data.frame(biomass=rep(biomasses, each=length(evenness)),
                        evenness=evenness,
                        fs=as.vector(t(fs)),
                        catch=as.vector(t(catches)))
    
    fs.plot <- ggplot(cr.df, aes(x=biomass, y=evenness, fill=fs, z=fs))+
        geom_raster(alpha=0.9)+
        geom_contour(breaks=c(0.1, 0.20), color="black", size=1)+
        geom_label_contour(breaks=c(0.1, 0.20), skip=0, label.placer=label_placer_fraction(0.5))+
        scale_fill_gradient(low="white", high="red", na.value = "transparent", name="Exploitation Rate")+
        scale_y_continuous(expand=c(0, 0), name="Age Structure Evenness Index")+
        scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
        theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))
    
    if(f.only){
        show(fs.plot)
        return(list(
            plot = fs.plot,
            data = cr.df
        ))
    }

    catches.plot <- ggplot(cr.df, aes(x=biomass, y=evenness, fill=catch, z=catch))+
        geom_raster(alpha=0.9)+
        geom_contour(breaks=c(1, 5000, 10000), color="black", size=1)+
        geom_label_contour(breaks=c(1, 5000, 10000), skip=0,label.placer=label_placer_fraction(0.5))+
        scale_fill_gradient(low="white", high="red", na.value = "black", name="Catch (mt)")+
        scale_y_continuous(expand=c(0, 0), name="Age Structure Evenness Index")+
        scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
        theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

    p <- ggarrange(fs.plot, catches.plot+rremove("ylab"))

    show(p)
    return(list(
        plot = p,
        data = cr.df
    ))

}
