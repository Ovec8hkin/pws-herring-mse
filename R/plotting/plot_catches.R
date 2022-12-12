library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

fnames <- c("foodbait_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "seine_catch.csv")

total.sims <- 6

set.seed(1998)
sims <- sample(1:1e4, size=total.sims)
nyr <- 25
sims <- c(197, 2255, 2386, 3709, 4716, 8388, 8904, 8634, 8935, 1094, 4288)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

catch.data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)
for(cr in hcr.names){
    cr.dat <- read.catch.data(cr, sims, nyr)
    catch.data <- catch.data %>% bind_rows(cr.dat)
}

low.regime.catch <- catch.data[catch.data$year %in% c(1, 2, 3, 4, 5, 21, 22, 23, 24, 25),]
high.regime.catch <- catch.data[catch.data$year %in% seq(6, 20, 1),]

# Total catch (all fisheries) by control rule
total.catch.df <- catch.data %>% na.omit() %>%
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
per.fishery.catch.df <- catch.data %>% na.omit() %>% 
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
med.catch.df <- catch.data %>% na.omit() %>%
    group_by(year, control.rule, sim) %>%
    summarise(total.catch = sum(catch)) %>%
    # group_by(year, control.rule, sim) %>%
    # summarise(cum.catch = cumsum(total.catch)) %>%
    # group_by(year, control.rule) %>%
    # summarise(med.catch = median(total.catch)) %>%
    print(n=100)
    
cum.catch.df <- data.frame(year=NA, control.rule=NA, sim=NA, cum.catch=NA)
for(cr in hcr.names){
    for(s in sims){
        df <- data.frame(year=1:nyr, control.rule=cr, sim=s, cum.catch=cumsum(med.catch.df$total.catch[med.catch.df$control.rule == cr & med.catch.df$sim == s]))
        cum.catch.df <- rbind(cum.catch.df, df)
    }
}
cum.catch.df <- cum.catch.df %>% na.omit() %>%
                group_by(year, control.rule) %>%
                summarise(cum.catch=quantile(cum.catch, c(0.025, 0.5, 0.975))) %>%
                mutate(percentile=c("2.5%", "50%", "97.5%")) %>%
                pivot_wider(
                    names_from = percentile,
                    values_from = cum.catch
                )
cum.catch.df$control.rule <- factor(cum.catch.df$control.rule, 
                                levels=c("base", "high.harvest", "low.harvest", "low.biomass", "high.biomass", "evenness", "gradient", "three.step.thresh", "constant.f.00", "big.fish"),
                                labels=c("Default", "High F", "Low F", "Low Threshold", "High Threshold", "Evenness", "Gradient", "Three-Step Threshold", "No Fishing", "Big Fish Only"))

cum.catch.plot <- ggplot(cum.catch.df)+
                    geom_point(aes(x=year, y=`50%`, color=control.rule)) +
                    geom_line(aes(x=year, y=`50%`, group=control.rule, color=control.rule), size=1.0)+
                    #geom_errorbar(aes(x=year, y=`50%`, ymin=`2.5%`, ymax=`97.5%`, color=control.rule))+
                    scale_x_continuous("Year", breaks=seq(0, nyr, by=2))+
                    scale_y_continuous("Catch (mt)", breaks=seq(0, 500000, 50000))+
                    scale_color_manual(values=c("black", "red", "blue", "#00A600", "#530d7e", "#E97902", "#AF0092", "#A6A6A6", "#B8DD4F", "#31aef1"))+
                    labs(x="Year", y="Catch (mt)", title="Cumulative Catch", color="Control Rule")+
                    theme(panel.grid.minor = element_blank())
cum.catch.plot
    
average.catch.df <- catch.data %>% na.omit() %>%
                        group_by(control.rule, sim) %>%
                        summarise(
                            tot.catch=sum(catch),
                            ann.catch=tot.catch/nyr
                        ) %>%
                        print(n=10)

average.catch.df$control.rule <- set.factor.levels(average.catch.df)

lr.average.catch.df <- low.regime.catch %>% na.omit() %>%
                        group_by(control.rule, sim) %>%
                        summarise(
                            tot.catch=sum(catch),
                            ann.catch=tot.catch/10
                        ) %>%
                        print(n=10)

lr.average.catch.df$control.rule <- set.factor.levels(lr.average.catch.df)

hr.average.catch.df <- high.regime.catch %>% na.omit() %>%
                        group_by(control.rule, sim) %>%
                        summarise(
                            tot.catch=sum(catch),
                            ann.catch=tot.catch/15
                        ) %>%
                        print(n=10)

hr.average.catch.df$control.rule <- set.factor.levels(hr.average.catch.df)

avg.catch.boxplot <- ggplot(average.catch.df) + 
                        geom_boxplot(aes(x=control.rule, y=ann.catch, fill=control.rule)) +
                        coord_flip() + 
                        geom_vline(aes(xintercept=1000))+
                        scale_y_continuous("Average Catch", limits=c(0, 20000), expand=c(0, 0)) +
                        scale_x_discrete("Control Rule")+
                        scale_fill_manual(values=c("black", "red", "blue", "#00A600", "#530d7e", "#E97902", "#AF0092", "#A6A6A6", "#B8DD4F", "#31aef1"))+
                        ggtitle("Average Annual Catch")+
                        theme(legend.position="none")
avg.catch.boxplot

library(patchwork)

(cum.catch.plot / avg.catch.boxplot)

set.factor.levels <- function(data){
    return(
        factor(data$control.rule, 
                levels=c("base", "high.harvest", "low.harvest", "low.biomass", "high.biomass", "evenness", "gradient", "three.step.thresh", "constant.f.00", "big.fish"),
                labels=c("Default", "High F", "Low F", "Low Threshold", "High Threshold", "Evenness", "Gradient", "Three-Step Threshold", "No Fishing", "Big Fish Only"))
    )
}

lr.plot <- ggplot(lr.average.catch.df) + 
        geom_boxplot(aes(x=control.rule, y=ann.catch, fill=control.rule)) +
        coord_flip() + 
        geom_vline(aes(xintercept=1000))+
        scale_y_continuous("Average Catch", limits=c(0, 20000), expand=c(0, 0)) +
        scale_x_discrete("Control Rule")+
        scale_fill_manual(values=c("black", "red", "blue", "#00A600", "#530d7e", "#E97902", "#AF0092", "#A6A6A6", "#B8DD4F", "#31aef1"))+
        ggtitle("Average Annual Catch (Low Recruitment)")+
        theme(legend.position="none")
lr.plot

hr.plot <- ggplot(hr.average.catch.df) + 
        geom_boxplot(aes(x=control.rule, y=ann.catch, fill=control.rule)) +
        coord_flip() + 
        geom_vline(aes(xintercept=1000))+
        scale_y_continuous("Average Catch", limits=c(0, 20000), expand=c(0, 0)) +
        scale_x_discrete("Control Rule")+
        scale_fill_manual(values=c("black", "red", "blue", "#00A600", "#530d7e", "#E97902", "#AF0092", "#A6A6A6", "#B8DD4F", "#31aef1"))+
        ggtitle("Average Annual Catch (High Recruitment)")+
        theme(legend.position="none")
hr.plot

avg.catch.boxplot

(avg.catch.boxplot / (lr.plot | hr.plot))
(avg.catch.boxplot | (lr.plot / hr.plot))
