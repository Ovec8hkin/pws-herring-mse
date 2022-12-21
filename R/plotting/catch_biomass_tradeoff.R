library(tidyverse)
library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_util_vals.r")))
source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

fnames <- c("foodbait_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "seine_catch.csv")
nyr <- 25
years <- seq(1980, 1980+42+nyr-1)

sims <- c(197, 2255, 2386, 3709, 4716, 8388, 8904, 8634, 8935, 1094, 4288)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

model.dir <- here::here("results/base/sim_197/year_25/model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)

bio.traj.df <- data.frame(year=NA, biomass=NA, control.rule=NA, sim=NA)
catch.data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)
for(cr in hcr.names){
    print(cr)
    cr.dat <- read.catch.data(cr, sims, nyr)
    catch.data <- catch.data %>% bind_rows(cr.dat)
    biomass.dat <- read.biomass.data(cr, sims, 25)
    bio.traj.df <- bio.traj.df %>% bind_rows(biomass.dat)
    # for(s in sims){
    #     #model.dir <- paste0(here::here("results/"), cr, "/sim_", s, "/year_25/model/")
    #     #biomass.df <- compute.biomass.traj(model.dir, length(years), years)
    #     # biomass.df <- read.biomass.data(cr)
    #     # tmp <- data.frame(year=biomass.df$year, biomass=biomass.df$biomass, sim=s, cr=cr)
    #     # bio.traj.df <- rbind(bio.traj.df, tmp)
    #     # rm(tmp)
    # }
}

aav <- function(catches){
    total.catch <- sum(catches)
    catch.diffs <- abs(diff(catches))
    aav <- mean((sum(catch.diffs)/total.catch))
    return(ifelse(is.nan(aav), 0, aav))
}

aav.df <- catch.data %>% na.omit() %>% 
            group_by(control.rule, sim, year) %>% 
            summarise(
                total.catch = sum(catch)
            ) %>%
            group_by(control.rule, sim) %>%
            summarise(
                aav = aav(total.catch)
            ) %>%
            print(n=10)

prob.threshold.df <- as_tibble(bio.traj.df) %>% na.omit() %>%
                    filter(year > 2021) %>%
                    group_by(control.rule, sim) %>%
                    summarise(
                        n=n(),
                        n.below=sum(biomass < 19958),
                        prob.below=n.below/n
                    ) %>%
                    print(n=10)

biomass.df <- as_tibble(bio.traj.df) %>% na.omit() %>%
                filter(year == max(bio.traj.df$year, na.rm=TRUE)) %>%
                mutate(depletion=biomass/40000) %>%
                select(-c(year)) %>%
                relocate(c(control.rule, sim, biomass, depletion)) %>%
                print(n=10)

average.catch.df <- catch.data %>% na.omit() %>%
                        group_by(control.rule, sim) %>%
                        summarise(
                            tot_catch=sum(catch),
                            ann_catch=tot_catch/nyr
                        ) %>%
                        print(n=10)

catch.biomass.df <- biomass.df %>% 
                    inner_join(average.catch.df, by=c("control.rule", "sim")) %>%
                    inner_join(aav.df, by=c("control.rule", "sim")) %>%
                    inner_join(prob.threshold.df, by=c("control.rule", "sim")) %>%
                    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels))


catch.biomass.df <- as_tibble(catch.biomass.df) %>%
                        select(-c(n, n.below)) %>%
                        group_by(control.rule) %>%
                        median_qi(
                            tot_catch, biomass, depletion, aav, prob.below,
                            .width=c(0.5, 0.95)
                        )
catch.biomass.df$cr <- factor(catch.biomass.df$cr, 
                                levels=c("base", "high.harvest", "low.harvest", "lower.b0", "low.biomass", "higher.b0", "high.biomass", "evenness", "gradient", "three.step.thresh", "constant.f.00", "big.fish"),
                                labels=c("Default", "High F", "Low F", "Low B0", "Low Threshold", "High B0", "High Threshold", "Evenness", "Gradient", "Three-Step Threshold", "No Fishing", "Big Fish Only"))

colnames(catch.biomass.df)[1] <- "control.rule"

tradeoff.df.long <- reformat.metric.df("tot_catch") %>%
                        bind_rows(
                            reformat.metric.df("biomass"),
                            reformat.metric.df("depletion"),
                            reformat.metric.df("aav"),
                            reformat.metric.df("prob.below")
                        )

fit.tradeoff.line <- function(formula, xrange, name){
    lm.model <- lm(formula, data=catch.biomass.df[catch.biomass.df$.width == 0.5,])
    newdata <- data.frame(x=xrange)
    colnames(newdata) <- c(name)
    preds <- predict(lm.model, newdata=newdata, se=TRUE)
    y = preds$fit
    ci <- preds$se.fit * qt(0.95 / 2 + .5, preds$df)
    ymin = y - ci
    ymax = y + ci
    pred.catch.biomass <- data.frame(x=xrange, y, ymin, ymax, se=preds$se.fit)

    return(list(
        model=lm.model,
        preds=pred.catch.biomass
    ))

}

generate.tradeoff.plot <- function(vars){

    d1 <- tradeoff.df.long %>% filter(metric == vars[1])
    d2 <- tradeoff.df.long %>% filter(metric == vars[2]) %>% select(c(control.rule, median, lower, upper, .width))

    d <- d1 %>% inner_join(d2, by=c("control.rule", ".width"))

    axis.breaks <- list(
        tot_catch = seq(0, 500000, 50000),
        biomass = seq(0, 150000, 10000),
        depletion = seq(0, 5, 0.5),
        aav = seq(1, 0, -0.1),
        prob.below = seq(0.5, 0, -0.05)
    )

    axis.labels <- list(
        tot_catch = seq(0, 500, 50),
        biomass = seq(0, 150, 10),
        depletion = seq(0, 5, 0.5),
        aav = seq(1, 0, -0.1),
        prob.below = seq(0.5, 0, -0.05)
    )

    metric.names <- list(
        tot_catch = "Total Catch (10000 mt)",
        biomass = "Final Year Biomass (10000 mt)",
        depletion = "Final Year Depletion Level",
        aav = "Average Annual Catch Variation",
        prob.below = "Proportion of Years Biomass Below Lower Regulatory Threshold"
    )

    metric.names.short <- list(
        tot_catch = "Total Catch",
        biomass = "Final Year Biomass",
        depletion = "Final Year Depletion",
        aav = "Catch Variation",
        prob.below = "Proportion of Years Below Threshold"

    )

    plot.title <- paste0(metric.names.short[[vars[1]]], " - ", metric.names.short[[vars[2]]], " Tradeoff Plot")

    x.breaks <- axis.breaks[[vars[1]]]
    x.labels <- axis.labels[[vars[1]]]

    y.breaks <- axis.breaks[[vars[2]]]
    y.labels <- axis.labels[[vars[2]]]

    x.axis.name <- metric.names[[vars[1]]]
    y.axis.name <- metric.names[[vars[2]]]

    plot <- ggplot(d)+
            geom_pointinterval(aes(x=median.x, xmin=lower.x, xmax=upper.x,
                                y=median.y,
                                color=control.rule))+
            geom_pointinterval(aes(x=median.x,
                                y=median.y, ymin=lower.y, ymax=upper.y,
                                color=control.rule))+
            #geom_smooth(aes_auto(model$preds), data=model$preds, stat="identity", color="black")+
            scale_x_continuous(x.axis.name, breaks=x.breaks, labels=x.labels)+
            scale_y_reverse(y.axis.name, breaks=y.breaks, labels=y.labels)+
            coord_cartesian(
                xlim=c(x.breaks[1], x.breaks[length(x.breaks)]), 
                ylim=c(y.breaks[1], y.breaks[length(y.breaks)]), 
                expand=c(0.1, 0.1)
            )+
            ggtitle(plot.title) + 
            theme(
                panel.grid.minor = element_blank()
            )

    show(plot)

    return(plot)
}

generate.tradeoff.plot(c("biomass", "tot_catch"))
generate.tradeoff.plot(c("biomass", "aav"))
generate.tradeoff.plot(c("tot_catch", "aav"))
generate.tradeoff.plot(c("tot_catch", "prob.below"))

library(patchwork)
(cvb.toff.plot+guides(color=FALSE) | aavvb.toff.plot+guides(color=FALSE) | aavvc.toff.plot)


cvb.lm.model   <- lm(tot_catch ~ biomass,     data=catch.biomass.df %>% filter(.width == 0.5) %>% select(tot_catch, biomass))
aavvb.lm.model <- lm(aav ~ biomass,           data=catch.biomass.df %>% filter(.width == 0.5) %>% select(aav, biomass))
aavvc.lm.model <- lm(aav ~ tot_catch,           data=catch.biomass.df %>% filter(.width == 0.5) %>% select(aav, tot_catch))
cvprob.lm.model <- lm(tot_catch ~ prob.below,   data=catch.biomass.df %>% filter(.width == 0.5) %>% select(tot_catch, prob.below) %>% mutate(prob.below = 100*prob.below))

intercepts  <- c(cvb.lm.model$coefficients[1],      aavvb.lm.model$coefficients[1],     aavvc.lm.model$coefficients[1],     cvprob.lm.model$coefficients[1])
slopes      <- c(cvb.lm.model$coefficients[2],      aavvb.lm.model$coefficients[2],     aavvc.lm.model$coefficients[2],     cvprob.lm.model$coefficients[2])
rsq         <- c(summary(cvb.lm.model)$r.squared,   summary(aavvb.lm.model)$r.squared , summary(aavvc.lm.model)$r.squared,  summary(cvprob.lm.model)$r.squared)

data.frame(name=c("Catch v Biomass", "AAV v Biomass", "AAV v Catch", "Catch v Prob"), intercept=intercepts, slope=slopes, r.squared=rsq)



test.df <- reformat.metric.df("tot_catch") %>%
    bind_rows(
        reformat.metric.df("biomass"),
        reformat.metric.df("depletion"),
        reformat.metric.df("aav"),
        reformat.metric.df("prob.below")
    )


generate.tradeoff.plot <- function(vars){

    d1 <- test.df %>% filter(metric == vars[1])
    d2 <- test.df %>% filter(metric == vars[2]) %>% select(c(control.rule, median, lower, upper, .width))

    d <- d1 %>% inner_join(d2, by=c("control.rule", ".width"))

    axis.breaks <- list(
        tot_catch = seq(0, 500000, 50000),
        biomass = seq(0, 150000, 10000),
        depletion = seq(0, 5, 0.5),
        aav = seq(1, 0, -0.1),
        prob.below = seq(0.5, 0, -0.05)
    )

    axis.labels <- list(
        tot_catch = seq(0, 500, 50),
        biomass = seq(0, 150, 10),
        depletion = seq(0, 5, 0.5),
        aav = seq(1, 0, -0.1),
        prob.below = seq(0.5, 0, -0.05)
    )

    metric.names <- list(
        tot_catch = "Total Catch (10000 mt)",
        biomass = "Final Year Biomass (10000 mt)",
        depletion = "Final Year Depletion Level",
        aav = "Average Annual Catch Variation",
        prob.below = "Probability Biomass Below Lower Regulatory Threshold"
    )

    metric.names.short <- list(
        tot_catch = "Total Catch",
        biomass = "Final Year Biomass",
        depletion = "Final Year Depletion",
        aav = "Catch Variation",
        prob.below = "Probability Biomass Below Threshold"

    )

    plot.title <- paste0(metric.names.short[[vars[1]]], " - ", metric.names.short[[vars[2]]], " Tradeoff Plot")

    x.breaks <- axis.breaks[[vars[1]]]
    x.labels <- axis.labels[[vars[1]]]

    y.breaks <- axis.breaks[[vars[2]]]
    y.labels <- axis.labels[[vars[2]]]

    x.axis.name <- metric.names[[vars[1]]]
    y.axis.name <- metric.names[[vars[2]]]

    plot <- ggplot(d)+
            geom_pointinterval(aes(x=median.x, xmin=lower.x, xmax=upper.x,
                                y=median.y,
                                color=control.rule))+
            geom_pointinterval(aes(x=median.x,
                                y=median.y, ymin=lower.y, ymax=upper.y,
                                color=control.rule))+
            #geom_smooth(aes_auto(model$preds), data=model$preds, stat="identity", color="black")+
            scale_x_continuous(x.axis.name, breaks=x.breaks, labels=x.labels)+
            scale_y_reverse(y.axis.name, breaks=y.breaks, labels=y.labels)+
            coord_cartesian(
                xlim=c(x.breaks[1], x.breaks[length(x.breaks)]), 
                ylim=c(y.breaks[1], y.breaks[length(y.breaks)]), 
                expand=c(0.1, 0.1)
            )+
            ggtitle(plot.title) + 
            theme(
                panel.grid.minor = element_blank()
            )

    show(plot)

    return(plot)
}



reformat.metric.df <- function(metric.name){
    return(
        catch.biomass.df %>% 
            select(c("control.rule", starts_with(metric.name), ".width", ".point", ".interval")) %>%
            rename_with(~ c("median", "lower", "upper"), c(metric.name, paste0(metric.name, ".lower"), paste0(metric.name, ".upper"))) %>%
            mutate(metric=metric.name) %>%
            relocate(c("control.rule", "metric", "median", "lower", "upper", ".width", ".point", ".interval"))

    )
}

