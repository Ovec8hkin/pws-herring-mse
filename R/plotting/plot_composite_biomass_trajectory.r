library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(paste0(here::here("R/utils"), "/fun_read_dat.R"))
source(paste0(here::here("R/plotting"), "/plot_util_vals.R"))
source(paste0(here::here("R/utils/other", "get_good_sims.R")))

nyr <- 30
years <- seq(1980, 1980+42+nyr-1)

set.seed(1120)
sims <- sample(1:1e4, size=150)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

year.0.fname <- paste0(here::here("results/admb123/sim_358/year_0/model/mcmc_out/"), "PFRBiomass.csv")
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

# stopCluster(cl)

good.sims <- get.good.sims()

biomass.traj.raw <- read_csv(file.path(here::here(), "results", "om_biomass.csv"), col_names = TRUE, show_col_types = FALSE) %>%
  select(year, biomass, control.rule, sim) %>%
  filter(sim %in% good.sims)

biomass.traj <- biomass.traj.raw %>% na.omit() %>% 
                    bind_rows(curr) %>%
                    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels)) %>%
                    group_by(year, control.rule) %>% 
                    median_qi(biomass, .width = c(.50, .95))    # Compute confidence intervals around biomass

# Side-by-side trajectories 
p1 <- ggplot(biomass.traj, aes(x=year, y=biomass, color=control.rule, group=control.rule))+
        geom_line(size=1.0)+
        scale_color_manual("Control Rule", values=as.vector(hcr.colors))+
        geom_lineribbon(
            data = biomass.traj %>% filter(control.rule == "Default" & year < 2022), 
            aes(ymin=.lower, ymax=.upper), 
            color="black", size=0.75
        )+
        geom_vline(xintercept=2022,    linetype="longdash", color="gray50")+
        geom_vline(xintercept=2022+15, linetype="longdash", color="gray50")+
        geom_hline(yintercept = 20000, linetype="longdash", color="gray50")+
        geom_hline(yintercept = 40000, linetype="longdash", color="gray50")+
        scale_fill_grey(start=0.8, end=0.6)+
        scale_x_continuous("Year", breaks=seq(1980, 2022+nyr, by=10), expand=c(0,0))+
        scale_y_continuous("Biomass (1,000 mt)", breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), labels=c(0, 20, 40, 50, 100, 150, 200), expand=c(0,0))+
        coord_cartesian(ylim=c(0, 200000))+
        ggtitle("Median Biomass Trajectories under Different Control Rules")+
        theme(
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.title = element_text(size=12),
            plot.title = element_blank()
        )+
        guides(fill="none", color="none")
p1

ggsave(file.path(here::here(), "figures", "present", "biomass_traj.png"), dpi=300, width=8, height=4, units="in")

constant.f.biomass <- as_tibble(biomass.traj.raw) %>% na.omit() %>%
                        bind_rows(curr) %>%
                        filter(control.rule == "constant.f.00") %>% 
                        select(year, biomass, sim) %>%
                        group_by(year, sim) %>%
                        summarise(biomass = median(biomass)) %>%
                        rename(biomass.unfished = biomass) %>%
                        print(n=100)

dynamic.b0.raw <- as_tibble(biomass.traj.raw) %>% na.omit() %>%
                      bind_rows(curr) %>%
                      left_join(
                        constant.f.biomass,
                        by = c("year", "sim")
                      ) %>%
                      group_by(year, sim, control.rule) %>%
                      summarise(
                        biomass = median(biomass),
                        biomass.unfished = median(biomass.unfished),
                        dyn.b0 = biomass/biomass.unfished
                      ) %>%
                      filter(dyn.b0 != 0) %>%
                      mutate(control.rule=recode_factor(control.rule, !!!hcr.levels)) %>%
                      group_by(year, control.rule) %>% 
                      median_qi(dyn.b0, .width = c(.50, .80))

p2 <- ggplot(dynamic.b0.raw %>% filter(year >= 2022), aes(x=year, y=dyn.b0, color=control.rule, group=control.rule))+
        geom_line(size=1.0)+
        geom_line(
          data = dynamic.b0.raw %>% filter(control.rule == "Default" & year <= 2022), 
          color="black", size=0.75
        )+
        geom_vline(xintercept=2022, linetype="longdash", color="gray50")+
        geom_vline(xintercept=2022+15, linetype="longdash", color="gray50")+
        scale_color_manual("Control Rule", values=as.vector(hcr.colors))+
        scale_x_continuous("Year", breaks=seq(1980, 2022+nyr, by=10), expand=c(0,0))+
        scale_y_continuous("Biomass Relative to Unfished", breaks=seq(0.0, 1.0, 0.2), expand=c(0,0))+
        coord_cartesian(ylim=c(0, 1.05))+
        ggtitle("Dynamic Biomass Trajectories under Different Control Rules")+
        theme(
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size=12),
          plot.title = element_blank()
        )+
        guides(fill="none", color=guide_legend(nrow=3,byrow=TRUE))
p2

ggsave(file.path(here::here(), "figures", "present", "rel_biomass_traj.png"), dpi=300, width=8, height=4, units="in")

library(patchwork)

(p1/p2)+
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides="collect") & 
  theme(
        legend.text = element_text(size=10),
        legend.position = 'bottom', 
        legend.direction = "horizontal"
  )

#ggsave(file.path(here::here(), "figures", "biomass_trajectories.png"), dpi=300, width=8, height=8, units="in")

#ggsave(file.path(here::here(), "figures", "publication", "Fig5_biomass_trajectories.jpg"), dpi=300, width=170, height=190, units="mm")
ggsave(file.path(here::here(), "figures", "publication", "Fig5_biomass_trajectories.pdf"), dpi=300, width=170, height=190, units="mm")




# Panelled trajectories by control rule
ggplot(biomass.traj, aes(x=year, y=biomass, ymin=.lower, ymax=.upper, group=1)) +
    geom_lineribbon(size=0.75)+
    geom_point(size=1.2)+
    geom_vline(xintercept=2022)+
    geom_hline(yintercept = 20000, linetype="longdash")+
    geom_hline(yintercept = 40000, linetype="longdash")+
    scale_fill_grey(start=0.8, end=0.6)+
    scale_x_continuous("Year", breaks=seq(1980, 2052, by=5))+
    scale_y_continuous("Pre-Fishery Biomass", breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), expand=c(0,0))+
    coord_cartesian(ylim=c(0, 220000))+
    facet_wrap(~control.rule, drop=TRUE, ncol=2)+
  theme_minimal()+
    theme(panel.grid.minor = element_blank())
