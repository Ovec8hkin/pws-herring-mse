library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(paste0(here::here("R/utils"), "/fun_read_dat.R"))
source(paste0(here::here("R/plotting"), "/plot_util_vals.R"))

nyr <- 25
years <- seq(1980, 1980+42+nyr-1)

sims <- c(197, 2255, 2386, 3709, 4716, 8388, 8904, 8634, 8935, 1094, 4288)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

year.0.fname <- paste0(here::here("results/base/sim_197/year_0/model/mcmc_out/"), "PFRBiomass.csv")
curr.biomass.data <- read_csv(year.0.fname, col_names=as.character(years)) %>% pivot_longer(everything(), names_to="year", values_to="biomass")
curr <- data.frame(
    year=rep(pull(select(curr.biomass.data, "year")), length(hcr.names)), 
    biomass=rep(pull(select(curr.biomass.data, "biomass")), length(hcr.names)), 
    control.rule=rep(hcr.names, each=length(pull(curr.biomass.data, "biomass"))),
    sim=0
)

biomass.traj.raw <- data.frame(year=NA, biomass=NA, control.rule=NA, sim=NA) 
for(cr in hcr.names){
    print(cr)
    biomass.traj.raw <- biomass.traj.raw %>% bind_rows(read.biomass.data(cr, sims, nyr))
}
            
biomass.traj <- biomass.traj.raw %>% na.omit() %>% 
                    bind_rows(curr) %>%
                    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels)) %>%
                    group_by(year, control.rule) %>% 
                    median_qi(biomass, .width = c(.50, .95))    # Compute confidence intervals around biomass

# Side-by-side trajectories 
ggplot(biomass.traj, aes(x=year, y=biomass, color=control.rule, group=control.rule))+
    geom_line(size=1.0)+
    scale_color_manual(values=as.vector(hcr.colors))+
    geom_lineribbon(
        data = biomass.traj %>% filter(control.rule == "Default" & year < 2022), 
        aes(ymin=.lower, ymax=.upper), 
        color="black", size=0.75
    )+
    geom_vline(xintercept=2022-1980)+
    geom_hline(yintercept = 20000, linetype="longdash")+
    geom_hline(yintercept = 40000, linetype="longdash")+
    scale_fill_grey(start=0.8, end=0.6)+
    scale_x_discrete("Year", breaks=seq(1980, 2022+nyr, by=5), expand=c(0,0))+
    scale_y_continuous("Pre-Fishery Biomass (1000 mt)", breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), labels=c(0, 20, 40, 50, 100, 150, 200), expand=c(0,0))+
    coord_cartesian(ylim=c(0, 200000))+
    ggtitle("Median Biomass Trajectories under Different Control Rules")+
    theme(
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line()
    )

# Panelled trajectories by control rule
ggplot(biomass.traj, aes(x=year, y=biomass, ymin=.lower, ymax=.upper, group=1)) +
    geom_lineribbon(size=0.75)+
    geom_point(size=1.2)+
    geom_vline(xintercept=2022-1980)+
    geom_hline(yintercept = 20000, linetype="longdash")+
    geom_hline(yintercept = 40000, linetype="longdash")+
    scale_fill_brewer(palette = "Blues")+
    scale_x_discrete("Year", breaks=seq(1980, 2038, by=5))+
    scale_y_continuous("Pre-Fishery Biomass", breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), expand=c(0,0))+
    coord_cartesian(ylim=c(0, 220000))+
    facet_wrap(~control.rule, drop=TRUE, ncol=2)+
    theme(panel.grid.minor = element_blank())
