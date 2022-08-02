library(here)
library(tidyverse)
source(paste0(here::here("R/utils/"), "fun_read_dat.R"))



compute.exploit.rate <- function(model.dir, nyr, years){
    exploit.rate <- read.exploit.rates(model.dir, nyr)

    exploit.rate.df <- as_tibble(exploit.rate) %>%
                    pivot_longer(everything(), names_to="year", values_to="exploit") %>%
                    group_by(year) %>%
                    median_qi(exploit, .width=c(0.95)) %>%
                    print(n=10)

    exploit.zeros <- exploit.rate.df[exploit.rate.df$exploit == 0, ]

    return(listN(exploit.rate.df, exploit.zeros))
}

compute.biomass.traj <- function(model.dir, nyr, years){
    biomass <- read.biomass.estimates(model.dir, nyr)

    biomass.df <- as_tibble(biomass) %>%
                            pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                            group_by(year) %>%
                            median_qi(biomass, .width=c(0.50, 0.95)) %>%
                            print(n=10)

    prob.below.threshhold <- apply(biomass, 2, function(x) sum(x < 19958))/nrow(biomass)

    biomass.df$prob <- as.vector(rep(prob.below.threshhold, 2))

    return(biomass.df)
}

compute.pfrb.posterior <- function(model.dir, nyr, years){
    biomass <- read.table(paste0(model.dir, "mcmc_out/PFRBiomass.csv"), header = FALSE, sep = ",", dec=".")[,nyr]

    biomass.df <- data.frame(biomass=biomass)
    prob.below.threshold <- round(sum(biomass < 20000)/length(biomass), 2)

    biomass.quants <- round(as.vector(apply(as.matrix(biomass), 2, quantile, c(0.025, 0.5, 0.975)))/1000, 2)

    return(listN(biomass.df, biomass.quants, prob.below.threshold))
}

compute.recruitment <- function(model.dir, nyr, years){
    age.3.recruits <- read.table(paste0(model.dir, "mcmc_out/Age3.csv"), header = FALSE, sep = ",", dec=".")[,1:nyr]
    colnames(age.3.recruits) <- years

    age.3.recruits.df <- as_tibble(age.3.recruits) %>%
                            pivot_longer(everything(), names_to="year", values_to="recruits") %>%
                            group_by(year) %>%
                            median_qi(recruits, .width=c(0.50, 0.95)) %>%
                            print(n=10)

    return(age.3.recruits.df)
}