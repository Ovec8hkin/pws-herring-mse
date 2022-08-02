source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2022
nyr.sim <- 15
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("results/save/base/sim_100/year_15/model/")

recruit.df <- compute.recruitment(model.dir, nyr, years)
plot <- plot.recruitment.posterior(recruit.df, years)
plot

