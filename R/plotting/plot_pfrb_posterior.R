source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2022
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("results/save/base/sim_100/year_15/model/")

pfrb.posterior <- compute.pfrb.posterior(model.dir, nyr, years)
plot <- plot.pfrb.posterior(pfrb.posterior$biomass.df, 
                            pfrb.posterior$biomass.quants, 
                            pfrb.posterior$prob.below.threshold,
                            curr.year)
plot