library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

fnames <- c("foodbait_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "seine_catch.csv")

sims <- c(42, 100, 8904, 9716)
nyr <- 15
control.rules <- c("base", "low.harvest", "high.harvest", "lower.b0", "low.biomass", "constant.f.00")

data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)
for(cr in control.rules){
    for(s in sims){
        for(f in fnames){
            fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", nyr, "/results/", f)
            dat.fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", nyr, "/model/")
            
            waa <- read.data.files(dat.fname)$PWS_ASA.dat[[4]][(nrow(waa)-nyr+1):nrow(waa),]
            catch.data <- as.matrix(read.csv(fname)[2:11])
            catch.biomass <- apply(catch.data*waa, 1, sum)

            d <- data.frame(year=1:nyr, catch=catch.biomass, fishery=str_split(f, "_")[[1]][1], control.rule=cr, sim=s)
            data <- data %>% bind_rows(d)
        }
    }
}

# Total catch (all fisheries) by control rule
total.catch.df <- data %>% na.omit() %>%
                    group_by(year, control.rule, sim) %>%
                    summarise(total.catch = sum(catch)) %>%
                    group_by(year, control.rule) %>%
                    summarise(total.catch=quantile(total.catch, c(0.025, 0.5, 0.975))) %>%
                    mutate(percentile=c("2.5%", "50%", "97.5%")) %>%
                    pivot_wider(
                        names_from = percentile,
                        values_from = total.catch
                    )

ggplot(total.catch.df) + 
    geom_col(aes(x=year, y=`50%`)) +
    geom_pointrange(aes(x=year, y=`50%`, ymax=`97.5%`, ymin=`2.5%`))+
    facet_wrap(~control.rule)


# Catch by fishery and control rule
per.fishery.catch.df <- data %>% na.omit() %>% 
                            group_by(year, fishery, control.rule, sim) %>%
                            summarise(total.catch = sum(catch)) %>%
                            group_by(year, control.rule, fishery) %>%
                            summarise(total.catch=quantile(total.catch, c(0.025, 0.5, 0.975))) %>%
                            mutate(percentile=c("2.5%", "50%", "97.5%")) %>%
                            pivot_wider(
                                names_from = percentile,
                                values_from = total.catch
                            )

ggplot(per.fishery.catch.df) + 
    geom_col(aes(x=year, y=`50%`, fill=fishery)) +
    geom_pointrange(data=total.catch.df, aes(x=year, y=`50%`, ymax=`97.5%`, ymin=`2.5%`)) +
    labs(x="Year", y="Catch (mt)", title="Fishery Catch", fill="Fishery")+
    facet_wrap(~control.rule)

ggplot(per.fishery.catch.df) + 
    geom_col(aes(x=year, y=`50%`, fill=fishery)) +
    #geom_pointrange(aes(x=year, y=`50%`, ymax=`97.5%`, ymin=`2.5%`)) +
    labs(x="Year", y="Catch (mt)", title="Fishery Catch", fill="Fishery") +
    facet_grid(rows=vars(fishery), cols=vars(control.rule))

ggplot(per.fishery.catch.df) + 
    geom_col(aes(x=year, y=`50%`, fill=fishery)) +
    #geom_pointrange(aes(x=year, y=`50%`, ymax=`97.5%`, ymin=`2.5%`)) +
    labs(x="Year", y="Catch (mt)", title="Fishery Catch", fill="Fishery") +
    facet_grid(rows=vars(control.rule), cols=vars(fishery))


# Cumulative catch by control rule
cum.catch.df <- data %>% na.omit() %>%
                    group_by(year, control.rule, sim) %>%
                    summarise(total.catch = sum(catch)) %>%
                    group_by(year, control.rule) %>%
                    summarise(mean.catch=mean(total.catch))
                    # # group_by(year, control.rule) %>%
                    # summarise(cum.catch.quant=quantile(total.catch, c(0.025, 0.5, 0.975))) %>%
                    # mutate(percentile=c("2.5%", "50%", "97.5%")) %>%
                    # pivot_wider(
                    #     names_from = percentile,
                    #     values_from = cum.catch.quant
                    # ) %>%
                    # group_by(year, control.rule) %>%
                    # summarise(cum.catch.median = cumsum(`50%`))

ggplot(cum.catch.df)+
    geom_point(aes(x=year, y=`50%`, color=control.rule)) +
    geom_line(aes(x=year, y=`50%`, group=control.rule, color=control.rule))
    #geom_errorbar(aes(x=year, ymin=`2.5%`, ymax=`97.5%`, color=control.rule))
