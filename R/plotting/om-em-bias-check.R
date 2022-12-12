library(tidyverse)

seeds <- c(197, 649, 1094, 1787, 2078, 2255, 2386, 3709, 4288, 4716, 8388, 8634, 8904, 8935, 9204, 9716)

hcr.name <- "evenness"
nyr.sim <- 25

om.mat <- matrix(NA, nrow=nyr.sim, ncol=length(seeds))
em.mat <- matrix(NA, nrow=nyr.sim, ncol=length(seeds))
i=1
for(s in seeds){
    om.prefish.biomass <- read_csv(paste0(here::here("results/"), hcr.name, "/sim_", s, "/year_", nyr.sim, "/results/prefish_spawn_biomass.csv")) %>%
                        rowwise() %>%
                        summarise(om.biomass=sum(c(`0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`))) %>%
                        print(n=25)

    em.prefish.biomass <- read_csv(paste0(here::here("results/"), hcr.name, "/sim_", s, "/year_", nyr.sim, "/model/mcmc_out/PFRBiomass.csv"), col_names=FALSE) %>%
                            summarise(across(everything(), median)) %>%
                            pivot_longer(everything(), names_to = "test", values_to = "em.biomass") %>%
                            select(em.biomass) %>%
                            print(n=75)

    om.mat[,i] <- om.prefish.biomass$om.biomass %>% as.numeric
    em.mat[,i] <- em.prefish.biomass$em.biomass[43:(43+nyr.sim-1)] %>% as.numeric

    i <- i+1

}

for(i in 1:length(seeds)){
  plot(1:25, om.mat[,i], col="red", type="l")
  lines(1:25, em.mat[, i], col="black", type="l")
}
