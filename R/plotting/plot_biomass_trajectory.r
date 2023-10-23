source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2022
nyr.sim <- 10
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("results/admb123/sim_2922/year_10/model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)
plot <- plot.biomass.trajectory(biomass.df, years)
plot
