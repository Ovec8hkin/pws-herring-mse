library(ggplot2)
library(metR)
library(ggpubr)

source(paste0(here::here("R/operating_model/control_rules/"), "cr_threshold.R"))

hcr.gradient <- function(curr.biomass, rel.biomass, ...,
                         options=list(
                            lower.threshold = 19958,
                            upper.threshold = 38555,
                            min.harvest = 0.0,
                            max.harvest = 0.2,
                            p=0.5
                          )){

    target.hr <- hcr.threshold.linear(curr.biomass, options)
    return(min(options$p*target.hr + (1-options$p)*target.hr*rel.biomass, 1.0))

}

plot.hcr.gradient <- function(biomasses, rel.biomasses, fs, catches, f.only=FALSE){
    
    cr.df <- data.frame(biomass=rep(biomasses, each=length(rel.biomasses)), ssb.change=rel.biomasses, fs=as.vector(t(fs)), catch=as.vector(t(catches)))

    fs.plot <- ggplot(cr.df, aes(x=biomass, y=ssb.change, fill=fs, z=fs))+
        geom_raster(alpha=0.9)+
        geom_contour(breaks=c(0.01, 0.1, 0.20, 0.30, 0.40), color="black", size=1)+
        geom_label_contour(breaks=c(0.01, 0.1, 0.20, 0.30, 0.40), skip=0, label.placer=label_placer_fraction(0.5))+
        scale_fill_gradient(low="white", high="red", na.value = "transparent", name="Exploitation Rate")+
        scale_y_continuous(expand=c(0, 0), name="3-year Biomass Gradient")+
        scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
        theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

    if(f.only){
      show(fs.plot)
      return(list(
        plot = fs.plot,
        fs = fs,
        catches = catches
      ))
    }

    catches.plot <- ggplot(cr.df, aes(x=biomass, y=ssb.change, fill=catch, z=catch))+
        geom_raster(alpha=0.9)+
        geom_contour(breaks=c(1, 5000, 10000, 15000, 20000), color="black", size=1)+
        geom_label_contour(breaks=c(1, 5000, 10000, 15000, 20000, 30000), skip=0,label.placer=label_placer_fraction(0.5))+
        scale_fill_gradient(low="white", high="red",na.value = "transparent", name="Catch (mt)")+
        scale_y_continuous(expand=c(0, 0), name="3-year Biomass Gradient")+
        scale_x_continuous(expand=c(0, 0), name="Pre-Fishery Biomass (mt)", breaks=seq(0, 60000, 10000), labels=seq(0, 60, 10))+
        theme(legend.position="bottom", legend.key.width=unit(1.5, "cm"))

    p <- ggpubr::ggarrange(fs.plot, catches.plot+rremove("ylab"))
    
    show(p)
    return(list(
      plot = p,
      fs = fs,
      catches = catches
    ))
}
