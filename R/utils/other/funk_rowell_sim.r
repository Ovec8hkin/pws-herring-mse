library(tidyverse)
library(doParallel)
source(file=paste0(here::here("R/"), "mse_loop.R"))

total.sims <- 200
nyr=500

#set.seed(2021)
seeds <- sample(1:1e4, size=total.sims)

unregister_dopar()
cores <- parallel::detectCores()
cl <- makeCluster(7, outfile="") #not to overload your computer
registerDoParallel(cl)

regime.length = 50

all.spawn.biomass <- rep(NA, nyr*total.sims)

cr <- list(type="hcr.constant.f", f.rate=0.0)

bio.traj <- pbapply::pblapply(seeds, function(sim.num, cr, nyr, regime.length){
    library(tidyverse)
    library(here)
    source(file=paste0(here::here("R/"), "mse_loop.R"))
    print(sim.num)
    write.dir <- paste0(here::here("results"), "/test/sim_", sim.num, "/")
    init.start.year <- sample(2016:2022, size=1)

    sim.out.f.00 <- run.simulation.simple(cr, nyr.sim=nyr, sim.seed=sim.num, write=NA, init.start.year=init.start.year, assessment = FALSE, hindcast=FALSE, max.regime.length =  regime.length)
    spawn.biomass <- sim.out.f.00$pop.dyn$spawn.biomass[,1]

}, cr=cr, nyr=nyr, cl=cl, regime.length=regime.length)
stopCluster(cl)

bio.traj <- bind_rows(bio.traj)

sim.data.ci <- bio.traj %>% na.omit() %>%
                  pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                  mutate(year=as.numeric(year)) %>%
                  group_by(year) %>%
                  median_qi(biomass, .width=c(0.5, 0.95))


regimes <- rep(rep(rep(c("high", "low"), each=regime.length), length.out = nyr), 2)

sim.data.ci$regime = regimes

plot.regimes <- data.frame(regime.change <- seq(0, nyr, regime.length))

ggplot(sim.data.ci, aes(x=year, y=biomass, ymin=.lower, ymax=.upper)) +
  geom_lineribbon(size=0.5)+
  geom_hline(aes(yintercept=40000), type="dashed")+
  geom_hline(aes(yintercept=20000), type="dashed")+
  geom_vline(data=plot.regimes, aes(xintercept=regime.change), color="grey")+
  scale_fill_brewer(palette="Blues")+
  coord_cartesian(ylim=c(0, 2e5), expand=0)+
  theme_classic()

full = c(sim.data.ci %>% pull(biomass) %>% mean, sim.data.ci %>% pull(biomass) %>% median, sim.data.ci %>% pull(biomass) %>% sd)
high = c(sim.data.ci %>% filter(regime == "high") %>% pull(biomass) %>% mean, sim.data.ci %>% filter(regime == "high") %>% pull(biomass) %>% median, sim.data.ci %>% filter(regime == "high") %>% pull(biomass) %>% sd)
low = c(sim.data.ci %>% filter(regime == "low") %>% pull(biomass) %>% mean, sim.data.ci %>% filter(regime == "low") %>% pull(biomass) %>% median, sim.data.ci %>% filter(regime == "low") %>% pull(biomass) %>% sd)
calcs. <- data.frame(full=full, high=high, low=low)
rownames(calcs) <- c("mean", "median", "sd")
calcs

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
