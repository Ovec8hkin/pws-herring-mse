library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)
library(doParallel)

source(paste0(here::here("R/utils"), "/fun_read_dat.R"))
source(paste0(here::here("R/plotting"), "/plot_util_vals.R"))

nyr <- 30
years <- seq(1980, 1980+42+nyr-1)
#sims <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4050, 4288, 4716, 4775, 6460, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)

set.seed(1120)
sims <- sample(1:1e4, size=150)

# set.seed(1)
# sims <- c(sims, sample(1:1e4, size=40))
#sims <- c(1017, 4775, 9725, 8462, 8789, 8522, 1799, 8229, 1129, 878, 7845, 5922, 6526, 5071, 4650, 2159, 3476, 2580, 1530, 7289, 4633, 4344, 1222, 2858, 5400, 526, 1069)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

year.0.fname <- paste0(here::here("results/base/sim_197/year_0/model/mcmc_out/"), "PFRBiomass.csv")
curr.biomass.data <- read_csv(year.0.fname, col_names=as.character(years)) %>% pivot_longer(everything(), names_to="year", values_to="biomass")
curr <- data.frame(
    year=as.numeric(rep(pull(select(curr.biomass.data, "year")), length(hcr.names))), 
    biomass=rep(pull(select(curr.biomass.data, "biomass")), length(hcr.names)), 
    control.rule=rep(hcr.names, each=length(pull(curr.biomass.data, "biomass"))),
    sim=0
)

# cores <- parallel::detectCores()
# cl <- makeCluster(min(cores[1]-1, length(hcr.names)), outfile="")
# registerDoParallel(cl)

# biomass.traj.raw <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
#     source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
#     biomass.dat <- read.true.biomass.data(cr, seeds, nyr)
# }, seeds=sims, nyr=nyr, cl=cl)
# biomass.traj.raw <- bind_rows(biomass.traj.raw)

# catch.data <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
#     source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
#     biomass.dat <- read.catch.data(cr, seeds, nyr)
# }, seeds=sims, nyr=nyr, cl=cl)
# catch.data <- bind_rows(catch.data)

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

ann.catch.traj <- catch.data %>% na.omit() %>%
                    mutate(
                        year = 2021+year
                    ) %>%
                    group_by(year, control.rule, sim) %>%
                    summarise(
                        total.catch = sum(catch)
                    )

sample.sims <- sample(sims, size=5)

sample.trajectories <- bio.traj.df %>% 
    filter(sim %in% sample.sims) %>% 
    left_join(
        ann.catch.traj %>% filter(sim %in% sample.sims), 
        by=c("year", "control.rule", "sim")
    ) %>%
    replace_na(list(year=0, biomass=0, control.rule=0, sim=0, total.catch=0)) %>%
    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels))

traj <- bio.traj.df %>% na.omit() %>%
                    left_join(ann.catch.traj, by=c("year", "control.rule", "sim")) %>%
                    replace_na(list(year=0, biomass=0, control.rule=0, sim=0, total.catch=0)) %>%
                    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels)) %>%
                    group_by(year, control.rule) %>% 
                    median_qi(biomass, total.catch, .width = c(.80))

bio.plot.1 <- ggplot(traj %>% filter(control.rule %in% c("Default", "Low Harvest", "High Harvest", "Low Threshold", "High Threshold"))) +
    geom_lineribbon(aes(x=year, y=biomass/1000, ymin=biomass.lower/1000, ymax=biomass.upper/1000, color=control.rule))+
    geom_line(data=sample.trajectories %>% filter(control.rule %in% c("Default", "Low Harvest", "High Harvest", "Low Threshold", "High Threshold")), aes(x=year, y=biomass/1000, group=sim), alpha=0.3)+
    geom_hline(yintercept=20000/1000, linetype="dashed")+
    geom_hline(yintercept=40000/1000, linetype="dashed")+
    scale_color_manual(values=hcr.colors.named)+
    scale_fill_grey(start=0.9, end=0.7)+
    coord_cartesian(ylim=c(0, 125000)/1000, expand=0.1)+
    labs(x="Year", y="Biomass\n(1000 mt)")+
    facet_wrap(~control.rule, nrow=1) +
    theme(
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=alpha("white", 0)),
        panel.spacing = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10)
    )

bio.plot.1 <- tag_facet(bio.plot.1, tag_pool=toupper(letters), size=3.5)

catch.plot.1 <- ggplot(traj %>% filter(control.rule %in% c("Default", "Low Harvest", "High Harvest", "Low Threshold", "High Threshold"))) +
    geom_lineribbon(aes(x=year, y=total.catch/1000, ymin=total.catch.lower/1000, ymax=total.catch.upper/1000, color=control.rule))+
    geom_line(data=sample.trajectories %>% filter(control.rule %in% c("Default", "Low Harvest", "High Harvest", "Low Threshold", "High Threshold")), aes(x=year, y=total.catch/1000, group=sim), alpha=0.3)+
    scale_color_manual(values=hcr.colors.named)+
    scale_fill_grey(start=0.9, end=0.7)+
    scale_y_continuous(breaks=seq(0, 45000, 15000)/1000)+
    coord_cartesian(ylim=c(0, 50000)/1000, expand=0)+
    labs(x="Year", y="Catch\n(1000 mt)")+
    facet_wrap(~control.rule, nrow=1) +
    theme(
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(fill=alpha("white", 0)),
        panel.spacing = unit(0, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10)
    )

#catch.plot.1 <- tag_facet(catch.plot.1, tag_pool=letters[6:26])

bio.plot.2 <- ggplot(traj %>% filter(control.rule %in% c("Evenness", "Gradient", "Three Step", "Big Fish", "No Fishing"))) +
    geom_lineribbon(aes(x=year, y=biomass/1000, ymin=biomass.lower/1000, ymax=biomass.upper/1000, color=control.rule))+
    geom_line(data=sample.trajectories %>% filter(control.rule %in% c("Evenness", "Gradient", "Three Step", "Big Fish", "No Fishing")), aes(x=year, y=biomass/1000, group=sim), alpha=0.3)+
    geom_hline(yintercept=20000/1000, linetype="dashed")+
    geom_hline(yintercept=40000/1000, linetype="dashed")+
    scale_color_manual(values=hcr.colors.named)+
    scale_fill_grey(start=0.9, end=0.7)+
    coord_cartesian(ylim=c(0, 125000)/1000, expand=0.1)+
    labs(x="Year", y="Biomass\n(1000 mt)")+
    facet_wrap(~control.rule, nrow=1) +
    theme(
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=alpha("white", 0)),
        panel.spacing = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10)
    )

bio.plot.2 <- tag_facet(bio.plot.2, tag_pool=toupper(letters[6:26]), size=3.5)

catch.plot.2 <- ggplot(traj %>% filter(control.rule %in% c("Evenness", "Gradient", "Three Step", "Big Fish", "No Fishing"))) +
    geom_lineribbon(aes(x=year, y=total.catch/1000, ymin=total.catch.lower/1000, ymax=total.catch.upper/1000, color=control.rule))+
    geom_line(data=sample.trajectories %>% filter(control.rule %in% c("Evenness", "Gradient", "Three Step", "Big Fish", "No Fishing")), aes(x=year, y=total.catch/1000, group=sim), alpha=0.3)+
    scale_color_manual(values=hcr.colors.named)+
    scale_fill_grey(start=0.9, end=0.7)+
    scale_y_continuous(breaks=seq(0, 45000, 15000)/1000)+
    coord_cartesian(ylim=c(0, 50000)/1000, expand=0)+
    labs(x="Year", y="Catch\n(1000 mt)")+
    facet_wrap(~control.rule, nrow=1) +
    theme(
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(fill=alpha("white", 0)),
        panel.spacing = unit(0, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10)
    )

#catch.plot.2 <- tag_facet(catch.plot.2, tag_pool=letters[16:26])

library(patchwork)
library(cowplot)

patch1 <- (bio.plot.1 / plot_spacer() / catch.plot.1) + plot_layout(guide="collect", heights=c(4, -0.5, 4)) & 
    theme(
        legend.position = 'none',
    )

patch2 <- (bio.plot.2 / plot_spacer() / catch.plot.2) + plot_layout(guide="collect", heights=c(4, -0.5, 4)) & 
    theme(
        legend.position = 'none',
    )

plot_grid(patch1+theme(legend.position="none"), patch2+theme(legend.position="none"), nrow=2)

#ggsave(file.path(here::here(), "figures", "bio_catch.png"), dpi=300, width=11, height=8.5, units="in")

#ggsave(file.path(here::here(), "figures", "publication", "Fig6_biomass_catch.jpg"), dpi=300, width=170, height=132, units="mm")

ggsave(file.path(here::here(), "figures", "publication", "Fig6_biomass_catch.pdf"), dpi=300, width=170, height=132, units="mm")
