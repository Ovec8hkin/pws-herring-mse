library(ggplot2)
library(ggpubr)

source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2022
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("results/save/base/sim_100/year_15/model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)
biomass.plot <- plot.biomass.trajectory(biomass.df, years)

exploit.df <- compute.exploit.rate(model.dir, nyr, years)
exploit.rate.plot <- plot.exploit.rate(exploit.df$exploit.rate.df,
                                       exploit.df$exploit.zeros,
                                       years)

pfrb.posterior <- compute.pfrb.posterior(model.dir, nyr, years)
pfrb.posterior.plot <- plot.pfrb.posterior(pfrb.posterior$biomass.df, 
                                           pfrb.posterior$biomass.quants, 
                                           pfrb.posterior$prob.below.threshold,
                                           curr.year,
                                           font.size=5)

recruit.df <- compute.recruitment(model.dir, nyr, years)
recruit.plot <- plot.recruitment.posterior(recruit.df, years)

ggarrange(
    recruit.plot, biomass.plot, exploit.rate.plot, pfrb.posterior.plot,
    nrow=2,
    ncol=2
)
