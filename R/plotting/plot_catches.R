library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
source(file=paste0(here::here("R/plotting/"), "plot_util_vals.R"))

nyr <- 25
sims <- c(197, 2255, 2386, 3709, 4716, 8388, 8904, 8634, 8935, 1094, 4288)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

catch.data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)
for(cr in hcr.names){
    cr.dat <- read.catch.data(cr, sims, nyr)
    catch.data <- catch.data %>% bind_rows(cr.dat)
}

catch.data <- catch.data %>%
                mutate(control.rule = recode_factor(control.rule, !!!hcr.levels))

low.regime.catch <- catch.data[catch.data$year %in% c(1, 2, 3, 4, 5, 21, 22, 23, 24, 25),]
high.regime.catch <- catch.data[catch.data$year %in% seq(6, 20, 1),]

# Total catch (all fisheries) by control rule
total.catch.df <- catch.data %>% na.omit() %>%
                    group_by(year, control.rule, sim) %>%
                    summarise(total.catch = sum(catch)) %>%
                    group_by(year, control.rule) %>%
                    median_qi(total.catch, .width=c(0.5, 0.95))

ggplot(total.catch.df) + 
    geom_col(
        data = total.catch.df %>% filter(.width == 0.5),
        aes(x=year, y=total.catch)
    ) +
    geom_pointrange(aes(x=year, y=total.catch, ymax=.upper, ymin=.lower))+
    facet_wrap(~control.rule)


# Catch by fishery and control rule
per.fishery.catch.df <- catch.data %>% na.omit() %>% 
                            group_by(year, fishery, control.rule, sim) %>%
                            summarise(total.catch = sum(catch)) %>%
                            group_by(year, control.rule, fishery) %>%
                            median_qi(total.catch, .width=c(0.5, 0.95))

ggplot(per.fishery.catch.df) + 
    geom_col(
        data = per.fishery.catch.df %>% filter(.width == 0.5),
        aes(x=year, y=total.catch, fill=fishery)
    ) +
    #geom_pointrange(data=total.catch.df, aes(x=year, y=total.catch, ymax=.upper, ymin=.lower)) +
    labs(x="Year", y="Catch (mt)", title="Fishery Catch", fill="Fishery")+
    facet_wrap(~control.rule)

ggplot(per.fishery.catch.df) + 
    geom_col(
        data = per.fishery.catch.df %>% filter(.width == 0.5),
        aes(x=year, y=total.catch, fill=fishery)
    ) +
    #geom_pointrange(aes(x=year, y=`50%`, ymax=`97.5%`, ymin=`2.5%`)) +
    labs(x="Year", y="Catch (mt)", title="Fishery Catch", fill="Fishery") +
    facet_grid(rows=vars(fishery), cols=vars(control.rule))

ggplot(per.fishery.catch.df) + 
    geom_col(
        data = per.fishery.catch.df %>% filter(.width == 0.5),
        aes(x=year, y=total.catch, fill=fishery)
    ) +
    #geom_pointrange(aes(x=year, y=`50%`, ymax=`97.5%`, ymin=`2.5%`)) +
    labs(x="Year", y="Catch (mt)", title="Fishery Catch", fill="Fishery") +
    facet_grid(rows=vars(control.rule), cols=vars(fishery))


# Cumulative catch by control rule
cum.catch.df <- catch.data %>% na.omit() %>%
    group_by(year, control.rule, sim) %>%
    summarise(
        total.catch = sum(catch)
    ) %>%
    group_by(control.rule, year) %>%
    summarise(
        total.catch = median(total.catch),
    ) %>%
    group_by(control.rule) %>%
    summarise(year=year, cum.catch = cumsum(total.catch)) %>%
    print(n=200)
    
# cum.catch.df <- data.frame(year=NA, control.rule=NA, sim=NA, cum.catch=NA)
# for(cr in hcr.names){
#     for(s in sims){
#         df <- data.frame(year=1:nyr, control.rule=cr, sim=s, cum.catch=cumsum(med.catch.df$total.catch[med.catch.df$control.rule == cr & med.catch.df$sim == s]))
#         cum.catch.df <- rbind(cum.catch.df, df)
#     }
# }
# cum.catch.df <- cum.catch.df %>% na.omit() %>%
#                 group_by(year, control.rule) %>%
#                 summarise(cum.catch=quantile(cum.catch, c(0.025, 0.5, 0.975))) %>%
#                 mutate(percentile=c("2.5%", "50%", "97.5%")) %>%
#                 pivot_wider(
#                     names_from = percentile,
#                     values_from = cum.catch
#                 )

cum.catch.plot <- ggplot(cum.catch.df)+
                    geom_point(aes(x=year, y=cum.catch, color=control.rule)) +
                    geom_line(aes(x=year, y=cum.catch, group=control.rule, color=control.rule), size=1.0)+
                    #geom_errorbar(aes(x=year, y=`50%`, ymin=`2.5%`, ymax=`97.5%`, color=control.rule))+
                    scale_x_continuous("Year", breaks=seq(0, nyr, by=2))+
                    scale_y_continuous("Catch (mt)", breaks=seq(0, 500000, 50000))+
                    scale_color_manual(values=as.vector(hcr.colors))+
                    labs(x="Year", y="Catch (mt)", title="Cumulative Catch", color="Control Rule")+
                    theme(panel.grid.minor = element_blank())
cum.catch.plot

library(patchwork)

catch.data.for.years <- function(years){
    return(
        catch.data %>% na.omit() %>%
            filter(year %in% years) %>%
            mutate(control.rule = recode_factor(control.rule, !!!hcr.levels)) %>%
            group_by(control.rule, sim) %>%
            summarise(
                tot.catch=sum(catch),
                ann.catch=tot.catch/length(years)
            ) %>%
            group_by(control.rule) %>%
            median_qi(tot.catch, ann.catch, .width=c(0.5, 0.75, 0.95))
    )
}

plot.catch.data <- function(data, var, title=NA){

    if(grepl("ann", var)){
        x.breaks <- seq(0, 40000, 4000)
        x.labs <- x.breaks/1000
        x.lab <- "Annual Catch"
    }else if(grepl("tot", var)){
        x.breaks <- seq(0, 800000, 100000)
        x.labs <- x.breaks/1000
        x.lab <- "Total Catch"
    }else if(grepl("biomass", var)){
        x.breaks <- seq(0, 200000, 10000)
        x.labs <- x.breaks/1000
        x.lab <- "Final Year Biomass"
    }else if(grepl("dep", var)){
        x.breaks <- seq(0, 4, 0.5)
        x.labs <- x.breaks
        x.lab <- "Final Year Depletion"
    }else if(grepl("aav", var)){
        x.breaks <- seq(0, 2, 0.25)
        x.labs <- x.breaks
        x.lab <- "Average Annual Catch Variation"
    }else{
        x.breaks <- seq(0, 1.0, 0.1)
        x.labs <- x.breaks
        x.lab <- "Proportion of Years Below Regulatory Threshold"
    }
    

    title <- ifelse(is.na(title), paste0(x.lab, " by Control Rule"), title)

    return(
        ggplot(data) + 
            geom_pointinterval(aes(
                x=.data[[var]], 
                y=control.rule, 
                xmin=.data[[paste0(var, ".lower")]], 
                xmax=.data[[paste0(var, ".upper")]], 
                color=control.rule
            )) +
            geom_vline(
                data = data %>% filter(control.rule == "Default" & .width == 0.5),
                aes(xintercept=.data[[var]])) + 
            scale_color_manual(values=as.vector(hcr.colors)) +
            scale_y_discrete(limits=rev) + 
            labs(x=paste0(x.lab), y="Control Rule", title=title) + 
            scale_x_continuous(breaks=x.breaks, labels=x.labs)+
            coord_cartesian(xlim=c(min(x.breaks), max(x.breaks)))+
            theme(
                legend.position = "none"
            )
    )
}

all.catch.data <- catch.data.for.years(1:nyr)
plot.catch.data(all.catch.data, "ann.catch")

lr.catch.data <- catch.data.for.years(c(1, 2, 3, 4, 5, 2, 21, 22, 23, 24, 25))
plot.catch.data(lr.catch.data, "ann.catch")

hr.catch.data <- catch.data.for.years(6:20)
plot.catch.data(hr.catch.data, "ann.catch")

library(patchwork)

(plot.catch.data(all.catch.data, "ann.catch") / 
    (plot.catch.data(lr.catch.data, "ann.catch", title="Low Recruitment") | plot.catch.data(hr.catch.data, "ann.catch", title="High Recruitment"))
)

plot.catch.data(catch.biomass.df, "tot_catch")
plot.catch.data(catch.biomass.df, "biomass")
plot.catch.data(catch.biomass.df, "depletion")
plot.catch.data(catch.biomass.df, "aav")
plot.catch.data(catch.biomass.df, "prob.below")    
