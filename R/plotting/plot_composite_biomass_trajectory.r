library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

set.seed(1998)
sims <- sample(1:1e4, size=total.sims)
years <- seq(1980, 2022)
nyr=10

hcr.names <- c("base", "low.harvest", "high.harvest", "low.biomass", "high.biomass", "lower.b0", "higher.b0", "constant.f.00")

#tic()
year.0.fname <- paste0(here::here("results/base/sim_197/year_0/model/mcmc_out/"), "PFRBiomass.csv")
curr.biomass.data <- read_csv(year.0.fname, col_names=as.character(years)) %>% pivot_longer(everything(), names_to="year", values_to="biomass")
curr <- data.frame(year=NA, biomass=NA, control.rule=NA, sim.num=NA)

df <- data.frame(year=NA, biomass=NA, control.rule=NA, sim.num=NA) 
for(cr in hcr.names){
    print(cr)
    curr <- curr %>% 
            bind_rows(
                data.frame(
                    year=pull(select(curr.biomass.data, "year")), 
                    biomass=pull(select(curr.biomass.data, "biomass")), 
                    control.rule=cr,
                    sim.num=0
                )
            ) %>% na.omit()
    for(s in sims){
        for(i in 1:nyr){
            fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", i, "/model/mcmc_out/PFRBiomass.csv")
            biomass.data <- read_csv(fname, col_names=FALSE, show_col_types = FALSE) %>% select(last_col())
            data <- data.frame(year=as.character(2022+i), biomass=pull(biomass.data), control.rule=cr, sim.num=s)
            df <- df %>% bind_rows(data) %>% na.omit()
        }
    }
}
            
df.trans <- df %>% bind_rows(curr) %>%                        # Remove columns data form columns that didnt separate
            group_by(year, control.rule) %>% 
            median_qi(biomass, .width = c(.50, .95))    # Compute confidence intervals around biomass
#toc()

# Panelled trajectories by control rule
ggplot(df.trans, aes(x=year, y=biomass, ymin=.lower, ymax=.upper, group=1)) +
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

df.trans$control.rule <- factor(df.trans$control.rule, 
                                levels=c("base", "base_small", "high.harvest", "low.harvest", "lower.b0", "low.biomass", "higher.b0", "high.biomass", "constant.f.00", "evenness", "gradient", "three.step.thresh", "big.fish"),
                                labels=c("Default", "Small SS", "High F", "Low F", "Low B0", "Low Threshold", "High B0", "High Threshold", "No Fishing", "Evenness", "Gradient", "Three-Step Threshold", "Big Fish"))

# Side-by-side trajectories 
ggplot(df.trans, aes(x=year, y=biomass, color=control.rule, group=control.rule))+
    geom_line(size=1.0)+
    scale_color_manual(values=c("black", "red", "blue", "#00A600", "#530d7e", "#E97902", "#AF0092", "#A6A6A6", "#B8DD4F", "#31aef1"))+
    geom_lineribbon(data=df.trans[(df.trans$control.rule=="Default" & df.trans$year < 2022),], aes(ymin=.lower, ymax=.upper), color="black", size=0.75)+
    geom_vline(xintercept=2022-1980)+
    geom_hline(yintercept = 20000, linetype="longdash")+
    geom_hline(yintercept = 40000, linetype="longdash")+
    scale_fill_grey(start=0.8, end=0.6)+
    scale_x_discrete("Year", breaks=seq(1980, 2022+nyr, by=5), expand=c(0,0))+
    scale_y_continuous("Pre-Fishery Biomass", breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000, 250000, 300000), expand=c(0,0))+
    coord_cartesian(ylim=c(0, 300000))+
    ggtitle("Biomass Trajectory under Eight Possible Control Rules")+
    theme(
        panel.grid.minor = element_blank(),
    )
