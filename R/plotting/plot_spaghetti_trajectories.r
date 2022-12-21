library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))

total.sims <- 5
nyr <- 35
years <- seq(1980, 1980+42+nyr-1)

seeds <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4288, 4716, 7251, 7915, 8388, 8904, 8935, 9204, 9260, 9716)

hcr.names <- c("base")#, "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

model.dir <- here::here("results/base/sim_197/year_35/model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)

bio.traj.df <- data.frame(year=NA, biomass=NA, sim=NA, cr=NA)

for(cr in hcr.names){
    for(s in seeds){
        model.dir <- paste0(here::here("results/"), cr, "/sim_", s, "/year_35/model/")
        biomass.df <- compute.biomass.traj(model.dir, length(years), years)
        tmp <- data.frame(year=biomass.df$year, biomass=biomass.df$biomass, sim=s, cr=cr)
        bio.traj.df <- rbind(bio.traj.df, tmp)
        rm(tmp)
    }
}

bio.traj.df <- bio.traj.df %>% na.omit() %>%
                mutate(sim=as.factor(sim))

mean.bio <- bio.traj.df %>% 
                group_by(year, cr) %>%
                summarise(biomass=median(biomass))

ggplot(bio.traj.df) +
    geom_line(aes(x=year, y=biomass, group=sim, color=sim), alpha=0.85)+
    geom_line(data=mean.bio, aes(x=year, y=biomass, group=1), color="black", size=1)+
    scale_x_continuous(breaks=seq(1980, max(years), 10), labels=seq(1980, max(years), 10))+
    coord_cartesian(ylim=c(0, 500000))+
    facet_wrap(~cr)
