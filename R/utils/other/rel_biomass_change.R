library(ggplot2)
library(tidyverse)
source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

model.dir <- paste0(here::here("results/base/sim_197/year_0/model/"))

biomass.est <- read.biomass.estimates(model.dir, nyr=42)

annual.biomass <- as_tibble(biomass.est) %>%
                    pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                    group_by(year) %>%
                    summarise(biomass=median(biomass)) %>%
                    print(n=10)

library(zoo)

n.year.rel.biomass <- function(biomass, n=3){
    rel.biomass.change <- rep(NA, length(biomass)-n+1)
    for(i in n:length(biomass)){
        rel.biomass.change[i-n+1] = biomass[i]/biomass[i-n+1]
    }
    return(rel.biomass.change)
}

three.year.biomass.change <- n.year.rel.biomass(annual.biomass$biomass, n=4)

plot(x=1983:2021, three.year.biomass.change, type="l")
abline(h=1.0)

hist(three.year.biomass.change, breaks=c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), freq=FALSE, xlab="3-year Relative Biomass Change", main="PWS Herring 3-year Relative Biomass Change")
