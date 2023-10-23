library(tidyverse)
library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)
library(doParallel)
library(patchwork)

source(file=paste0(here::here(), "/R/calc_utility.R"))
source(file=paste0(here::here("R/plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("R/plotting/", "plot_util_vals.r")))
source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

reformat.metric.df <- function(metric.name){
    return(
        catch.biomass.df %>% 
            select(c("control.rule", starts_with(metric.name), ".width", ".point", ".interval")) %>%
            rename_with(~ c("median", "lower", "upper"), c(metric.name, paste0(metric.name, ".lower"), paste0(metric.name, ".upper"))) %>%
            mutate(metric=metric.name) %>%
            relocate(c("control.rule", "metric", "median", "lower", "upper", ".width", ".point", ".interval"))

    )
}

nyr <- 30
years <- seq(1980, 1980+42+nyr-1)

#sims <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 3169, 3709, 4288, 4716, 4775, 6460, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)
set.seed(1120)
sims <- sample(1e4, 150)
hcr.names <- c("base", "low.harvest", "high.harvest", "low.biomass", "high.biomass", "constant.f.00", "evenness", "gradient", "three.step.thresh", "big.fish")

# cores <- parallel::detectCores()
# cl <- makeCluster(min(cores[1]-1, length(hcr.names)), outfile="")
# registerDoParallel(cl)

# bio.traj.df <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
#     source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
#     biomass.dat <- read.true.biomass.data(cr, seeds, nyr)
# }, seeds=sims, nyr=nyr, cl=cl)
# bio.traj.df <- bind_rows(bio.traj.df)

# catch.data <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
#     source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
#     biomass.dat <- read.catch.data(cr, seeds, nyr)
# }, seeds=sims, nyr=nyr, cl=cl)
# catch.data <- bind_rows(catch.data)

# #unregister_dopar()
# stopCluster(cl)

good.sims <- get.good.sims()

catch.data <- read_csv(file.path(here::here(), "results", "om_catch.csv"), col_names=TRUE, show_col_types=FALSE) %>%
    select(year, catch, fishery, control.rule, sim) %>%
    filter(sim %in% good.sims)  %>%
    print(n=10)
bio.traj.df <- read_csv(file.path(here::here(), "results", "om_biomass.csv"), col_names = TRUE, show_col_types = FALSE) %>%
    select(year, biomass, control.rule, sim) %>%
    filter(sim %in% good.sims)  %>%
    print(n=10)

aav <- function(data){
    total <- mean(data)
    diffs <- abs(diff(data))
    aav <- sum(diffs/total)/(length(data)-1)
    return(ifelse(is.nan(aav), 0, aav)) # If all data is 0, return 0 rather than NA
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
                left_join(
                    as_tibble(bio.traj.df) %>% na.omit() %>%                          
                      filter(year == max(bio.traj.df$year, na.rm=TRUE) & control.rule == "constant.f.00") %>%
                      rename(no.fish.biomass = biomass) %>%
                      select(c(sim, no.fish.biomass)),
                    by = "sim"
                ) %>%
                mutate(
                    depletion = biomass/40000,
                    dyn.b0 = biomass/no.fish.biomass
                ) %>%
                select(-c(year)) %>%
                relocate(c(control.rule, sim, biomass, depletion, dyn.b0)) %>%
                print(n=10)

avg.biomass.df <- as_tibble(bio.traj.df) %>% na.omit() %>%                      # Remove NAs
                    left_join(
                      as_tibble(bio.traj.df) %>% na.omit() %>%                          
                        filter(control.rule == "constant.f.00") %>%
                        rename(no.fish.biomass = biomass) %>%
                        select(c(sim, no.fish.biomass)),
                      by = c("sim")
                    ) %>%
                    group_by(control.rule, sim) %>%                             # Group
                    summarise(
                        avg.bio = median(biomass),                              # Median biomass across entire sim
                        avg.dep = avg.bio/40000,                                # Median depletion across entire sim
                        avg.db0 = median(biomass/no.fish.biomass)
                    ) %>%
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
                    inner_join(avg.biomass.df, by=c("control.rule", "sim")) %>%
                    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels))


catch.biomass.df <- as_tibble(catch.biomass.df) %>%
                        select(-c(n, n.below)) %>%
                        group_by(control.rule) %>%
                        median_qi(
                            tot_catch, ann_catch, biomass, depletion, aav, prob.below, dyn.b0, avg.db0,
                            .width=c(0.5, 0.80)
                        )
catch.biomass.df$cr <- factor(catch.biomass.df$control.rule, 
                                levels=c("base", "high.harvest", "low.harvest", "lower.b0", "low.biomass", "higher.b0", "high.biomass", "evenness", "gradient", "three.step.thresh", "constant.f.00", "big.fish"),
                                labels=c("Default", "High F", "Low F", "Low B0", "Low Threshold", "High B0", "High Threshold", "Evenness", "Gradient", "Three-Step Threshold", "No Fishing", "Big Fish Only"))

colnames(catch.biomass.df)[1] <- "control.rule"

tradeoff.df.long <- reformat.metric.df("ann_catch") %>%
                        bind_rows(
                            reformat.metric.df("biomass"),
                            reformat.metric.df("depletion"),
                            reformat.metric.df("aav"),
                            reformat.metric.df("prob.below"),
                            reformat.metric.df("dyn.b0"),
                            reformat.metric.df("avg.db0")
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
        tot_catch = seq(-50000, 550000, 50000),
        ann_catch = seq(-5000, 20000, 5000),
        biomass = seq(-10000, 160000, 10000),
        depletion = seq(-1, 5, 1),
        aav = round(seq(1.2, -0.2, -0.2), 2),
        prob.below = seq(0.5, 0, -0.05),
        dyn.b0 = seq(0.2, 1.2, 0.2),
        avg.db0 = seq(0.2, 1.2, 0.2)
    )

    axis.labels <- list(
        tot_catch = seq(-50, 550, 50),
        ann_catch = seq(-5, 20, 5),
        biomass = seq(-10, 160, 10),
        depletion = seq(-1, 5, 1),
        aav = round(seq(1.2, -0.2, -0.2), 2),
        prob.below = seq(0.5, 0, -0.05),
        dyn.b0 = seq(0.2, 1.2, 0.2),
        avg.db0 = seq(0.2, 1.2, 0.2)
    )

    metric.names <- list(
        tot_catch = "Total Catch (1000 mt)",
        ann_catch = "Annual Catch (1000 mt)",
        biomass = "Final Year Biomass (1000 mt)",
        depletion = "Final Year Depletion Level",
        aav = "Average Annual Catch Variation",
        prob.below = "Proportion of Years Biomass Below Lower Regulatory Threshold",
        dyn.b0 = "Final Year Biomass Relative to Unfished Conditions"
    )

    metric.names.short <- list(
        tot_catch = "Total Catch",
        ann_catch = "Annual Catch",
        biomass = "Final Year Biomass",
        depletion = "Final Year Depletion",
        aav = "Catch Variation",
        prob.below = "Proportion of Years Below Threshold",
        dyn.b0 = "Relative Biomass",
        avg.db0 = "Relative Biomass"

    )

    plot.title <- paste0(metric.names.short[[vars[1]]], " - ", metric.names.short[[vars[2]]], " Tradeoff Plot")

    x.breaks <- axis.breaks[[vars[1]]]
    x.labels <- axis.labels[[vars[1]]]

    y.breaks <- axis.breaks[[vars[2]]]
    y.labels <- axis.labels[[vars[2]]]

    x.axis.name <- metric.names.short[[vars[1]]]
    y.axis.name <- metric.names.short[[vars[2]]]

    print(x.breaks)

    model <- fit.tradeoff.line(as.formula(paste0(vars[2], "~", vars[1])), seq(min(x.breaks), max(x.breaks), length.out=100), name=vars[1])

    plot <- ggplot(d)+
            geom_smooth(aes(x=x, y=y, ymin=ymin, ymax=ymax), data=model$preds, stat="identity", color="black")+
            geom_point(aes(x=median.x, y=median.y, color=control.rule), size=4)+
            geom_pointinterval(aes(x=median.x, xmin=lower.x, xmax=upper.x,
                                y=median.y,
                                color=control.rule), alpha=0.33)+
            geom_pointinterval(aes(x=median.x,
                                y=median.y, ymin=lower.y, ymax=upper.y,
                                color=control.rule), alpha=0.33)+
            scale_color_manual(values=hcr.colors.named, name="Control Rule") + 
            scale_x_continuous(x.axis.name, breaks=x.breaks, labels=x.labels)+
            scale_y_reverse(y.axis.name, breaks=y.breaks, labels=y.labels)+
            coord_cartesian(
                xlim=c(x.breaks[2], x.breaks[length(x.breaks)-1]), 
                ylim=c(y.breaks[2], y.breaks[length(y.breaks)-1]), 
                expand=0
            )+
            ggtitle(plot.title) +
            theme_minimal()+ 
            theme(
                axis.line.x = element_line(),
                axis.line.y = element_line(),
                axis.text = element_text(size=15),
                axis.title = element_text(size=22, face="bold"),
                plot.title = element_blank(),
                panel.grid.minor = element_blank()
            )

    #show(plot)

    return(plot)
}

generate.utility.plot <- function(vars){

    yvar <- vars[1]
    xvar <- vars[2]

    d1 <- tradeoff.df.long %>% filter(metric == xvar)
    d2 <- tradeoff.df.long %>% filter(metric == yvar) %>% select(c(control.rule, median, lower, upper, .width))

    d <- d1 %>% inner_join(d2, by=c("control.rule", ".width"))

    bounds <- list(
        ann_catch = c(0, 15000),
        depletion = c(0, 4),
        dyn.b0 = c(0.4, 1),
        avg.db0 = c(0.4, 1),
        aav = c(0, 1)
    )

    x.bounds <- bounds[[xvar]]
    y.bounds <- bounds[[yvar]]

    x.seq <- seq(x.bounds[1], x.bounds[2], length.out=1000)
    y.seq <- seq(y.bounds[1], y.bounds[2], length.out=1000)

    u.x <- apply(as.matrix(x.seq), 1, function(x) calc.utility(x, xvar))
    u.y <- apply(as.matrix(y.seq), 1, function(x) calc.utility(x, yvar))

    cb.util.df <- expand.grid(y.seq, x.seq)
    cb.total.utility <- apply(expand.grid(u.y, u.x), 1, total.utility)
    cb.util.df$total.util <- cb.total.utility
    colnames(cb.util.df) <- c("y", "x", "total.util")

    plot <- ggplot(cb.util.df)+
                geom_raster(aes(x=x, y=y, fill=total.util, z=total.util))+
                geom_contour(aes(x=x, y=y, fill=total.util, z=total.util), breaks=c(0.25, 0.50, 0.75, 1.0), color="black", size=1)+
                geom_label_contour(aes(x=x, y=y, fill=total.util, z=total.util), breaks=c(0.25, 0.50, 0.75, 1.0), skip=0, label.placer=label_placer_fraction(0.5))+
                geom_point(data=d, aes(x=median.x, y=median.y, color=control.rule), size=4)+
                scale_color_manual(values=hcr.colors.named, name="Control Rule") +
                scale_fill_gradient("Utility", low="white", high="red", limits=c(0, 1))+
                coord_cartesian(expand=0)+
                guides(color="none")+
                theme(
                    axis.line = element_line()
                )

    if(xvar == "aav"){
        plot <- plot + scale_x_reverse()
    }

    return(plot)

}

strip.gg.axes <- function(plot, x=TRUE, y=TRUE){
    if(x){
        plot <- plot +
            theme(
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank()
            )
    }

    if(y){
        plot <- plot +
            theme(
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank()
            )
    }

    return(plot)
}

generate.fake.plot <- function(bounds, name, strip.x=TRUE, strip.y=TRUE, rev="none"){
    data <- data.frame(
        x=seq(bounds[1], bounds[2], length=3),
        y=seq(bounds[1], bounds[2], length=3)
    )

    plot <- ggplot(data)+
        geom_point(aes(x=x, y=y), color="white", alpha=0.0)+
        labs(x=name, y=name)+
        theme(
            panel.background = element_blank(),
            axis.text = element_text(size=15),
            axis.title = element_text(size=22, face="bold")
        )

    if(strip.x){
        plot <- plot+theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )
    }

    if(strip.y){
        plot <- plot+theme(
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
        )
    }

    if(rev == "x"){
        plot <- plot+scale_x_reverse()
    }

    return(plot)

}

library(metR)

dvc.trade.plot <- generate.tradeoff.plot(c("ann_catch", "avg.db0"))
dvv.trade.plot <- generate.tradeoff.plot(c("avg.db0", "aav"))
cvv.trade.plot <- generate.tradeoff.plot(c("ann_catch", "aav"))

cvd.util.plot <- generate.utility.plot(c("ann_catch", "avg.db0"))
cvv.util.plot <- generate.utility.plot(c("ann_catch", "aav"))
dvv.util.plot <- generate.utility.plot(c("avg.db0", "aav"))

fake.1 <- generate.fake.plot(c(0, 15), name="Annual Catch", strip.y=FALSE)
fake.2 <- generate.fake.plot(c(0, 4), name="Relative Biomass")
fake.3 <- generate.fake.plot(c(1, 0), name="Catch Variation", strip.x=FALSE, rev="x")

(fake.1                                 + strip.gg.axes(cvd.util.plot)              + strip.gg.axes(cvv.util.plot)) /
(strip.gg.axes(dvc.trade.plot, y=FALSE) + fake.2                                    + strip.gg.axes(dvv.util.plot)) /
(cvv.trade.plot                         + strip.gg.axes(dvv.trade.plot, x=FALSE)    + fake.3) +
plot_layout(guides="collect")#+plot_annotation(tag_levels = "A")

ggsave(file.path(here::here(), "figures", "tradeoff.png"), dpi=300, width=14, height=12, units="in")
ggsave(file.path(here::here(), "figures", "present", "tradeoff.png"), dpi=300, width=14, height=12, units="in")

(dvc.trade.plot| dvv.trade.plot| cvv.trade.plot) + plot_layout(guides="collect") & theme(legend.position="top-right")

###################

cvb.lm.model   <- lm(ann_catch ~ dyn.b0,     data=catch.biomass.df %>% filter(.width == 0.5) %>% select(ann_catch, dyn.b0))
aavvb.lm.model <- lm(aav ~ dyn.b0,           data=catch.biomass.df %>% filter(.width == 0.5) %>% select(aav, dyn.b0))
aavvc.lm.model <- lm(aav ~ ann_catch,           data=catch.biomass.df %>% filter(.width == 0.5) %>% select(aav, ann_catch))
#cvprob.lm.model <- lm(tot_catch ~ prob.below,   data=catch.biomass.df %>% filter(.width == 0.5) %>% select(tot_catch, prob.below) %>% mutate(prob.below = 100*prob.below))

intercepts  <- c(cvb.lm.model$coefficients[1],      aavvb.lm.model$coefficients[1],     aavvc.lm.model$coefficients[1])
slopes      <- c(cvb.lm.model$coefficients[2],      aavvb.lm.model$coefficients[2],     aavvc.lm.model$coefficients[2])
rsq         <- c(summary(cvb.lm.model)$r.squared,   summary(aavvb.lm.model)$r.squared , summary(aavvc.lm.model)$r.squared)

data.frame(name=c("Catch v Biomass", "AAV v Biomass", "AAV v Catch"), intercept=intercepts, slope=slopes, r.squared=rsq)
