library(here)
library(ggplot2)
library(tidyverse)
library(ggdist)

color.palette <- "Blues"

plot.biomass.trajectory <- function(df, years){
    return(
        ggplot(df) +
            geom_lineribbon(aes(x=year, y=biomass/1000, ymin=.lower/1000, ymax=.upper/1000, group=1), size=0.75)+
            geom_point(aes(x=year, y=prob*200))+
            geom_line(aes(x=year, y=prob*200, group=1))+
            scale_fill_brewer(palette = color.palette)+
            scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5))+
            geom_hline(yintercept=c(20, 40), linetype="dashed")+
            ggtitle("Biomass Trajectory")+
            scale_y_continuous(
                "Pre-Fishery Biomass (1000 mt)", 
                breaks=c(0, 20, 40, 50, 100, 150, 200),
                sec.axis = sec_axis(trans=~.*1/200, name="Probability below 20k metric tons"))+
            theme(
                panel.grid.minor = element_blank()
            )
    )
}

plot.exploit.rate <- function(df, zeros, years){
    return(
        ggplot(df, aes(x=year, y=exploit, ymin=.lower, ymax=.upper))+
            geom_pointinterval()+
            geom_pointinterval(data=zeros, color="grey")+
            scale_x_discrete("Year", breaks=seq(min(years), max(years), 5))+
            scale_y_continuous("Exploitation Rate", breaks=seq(0, 0.3, 0.05))+
            ggtitle("Exploitation Rate")
    )
}

plot.pfrb.posterior <- function(df, quants, prob, curr.year, font.size=8){
    return(
        ggplot(df)+
            geom_histogram(aes(x=biomass/1000, y= ..density..), bins=60)+
            scale_fill_brewer(palette=color.palette)+
            geom_vline(xintercept = quants, linetype=c("dashed", "solid", "dashed"), size=c(1, 2, 1))+
            geom_text(aes(x=45, y=0.14, label=paste("Median:", quants[2])), size=font.size, hjust=1)+
            geom_text(aes(x=45, y=0.12, label=paste0("95% interval:\n", "(", quants[1], ", ", quants[3], ")")), size=font.size, hjust=1)+
            geom_text(aes(x=45, y=0.09, label=paste("Probability below\nthreshold:", prob)), size=font.size, hjust=1)+
            scale_x_continuous(paste(curr.year-1, "Pre-Fishery Biomass (mt)"), breaks=seq(0, 50, 5))+
            scale_y_continuous("Probability Density", breaks=seq(0, 0.15, 0.05), expand=c(0, 0))+
            coord_cartesian(ylim=c(0, 0.15), xlim=c(0, 50))+
            ggtitle(paste(curr.year-1, "Pre-fishery Biomass Posterior Probability Density"))+
            theme(
                panel.grid = element_blank()
            )
    )
}

plot.recruitment.posterior <- function(df, years){
     return(
        ggplot(df) +
            geom_lineribbon(aes(x=year, y=recruits, ymin=.lower, ymax=.upper, group=1), size=0.75)+
            scale_fill_brewer(palette = color.palette)+
            scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5))+
            scale_y_continuous("Age-3 Recruits (millions)", breaks=seq(0, 2000, by=500))+
            ggtitle("Age-3 Recruitment")
    )
}

