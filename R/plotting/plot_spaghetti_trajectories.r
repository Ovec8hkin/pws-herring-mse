library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))

total.sims <- 5
nyr <- 25
years <- seq(1980, 1980+42+nyr-1)

sims <- c(197, 2255, 2386, 3709, 4716, 8388, 8904, 8634, 8935, 1094, 4288)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

model.dir <- here::here("results/base/sim_197/year_25/model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)

bio.traj.df <- data.frame(year=NA, biomass=NA, sim=NA, cr=NA)

for(cr in hcr.names){
    for(s in sims){
        model.dir <- paste0(here::here("results/"), cr, "/sim_", s, "/year_25/model/")
        biomass.df <- compute.biomass.traj(model.dir, length(years), years)
        tmp <- data.frame(year=biomass.df$year, biomass=biomass.df$biomass, sim=s, cr=cr)
        bio.traj.df <- rbind(bio.traj.df, tmp)
        rm(tmp)
    }
}


mean.bio <- bio.traj.df %>% na.omit() %>%
                group_by(year, cr) %>%
                summarise(biomass=median(biomass))

ggplot(bio.traj.df) +
    geom_line(aes(x=year, y=biomass, group=sim), color="grey")+
    geom_line(data=mean.bio, aes(x=year, y=biomass, group=1), color="black", size=1)+
    facet_wrap(~cr)
