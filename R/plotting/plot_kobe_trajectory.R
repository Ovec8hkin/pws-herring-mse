library(ggplot2)
library(ggdist)
library(tidyverse)

source(paste0(here::here("R/utils/"), "fun_read_dat.R"))

model.dir <- here::here("results/save/base/sim_100/year_15/model/")

b.star <- 40000
f.star <- 0.20

biomass.estimates <- read.biomass.estimates(model.dir)
exploit.rate.estimates <- read.exploit.rates(model.dir)

biomass.est.rel <- biomass.estimates/b.star
exploit.est.rel <- exploit.rate.estimates/f.star

exploit.rate.df <- as_tibble(exploit.est.rel) %>%
                    pivot_longer(everything(), names_to="year", values_to="exploit") %>%
                    group_by(year) %>%
                    median_qi(exploit, .width=c(0.95)) %>%
                    print(n=10)

biomass.df <- as_tibble(biomass.est.rel) %>%
                pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                group_by(year) %>%
                median_qi(biomass, .width=c(0.95)) %>%
                print(n=10)

biomass.exploit.df <- biomass.df %>%
                        left_join(exploit.rate.df, by="year", suffix=c("bio", "exp"))

rect_dat <- data.frame(panel = c("bottom_left", "top_right",
                                    "bottom_right", "top_left"),
                          x_min = c(-Inf, 1, 1, -Inf),
                          x_max = c(1, Inf, Inf, 1), y_min = c(-Inf,1, -Inf, 1),
                          y_max = c(1, Inf, 1, Inf))



ggplot(biomass.exploit.df)+
    geom_rect(data=rect_dat, aes(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max, fill = panel))+
    geom_point(aes(x=biomass, y=exploit), size=3)+
    geom_label_repel(data=biomass.exploit.df[(as.numeric(biomass.exploit.df$year) %% 5)== 0,], aes(x=biomass, y=exploit, label=year), box.padding = 0.35, point.padding = 0.5, max.overlaps=5)+
    geom_errorbar(data=tail(biomass.exploit.df, n=1), aes(x=biomass, y=exploit, ymin=.lowerexp, ymax=.upperexp), width=0.1, color="white")+
    geom_errorbarh(data=tail(biomass.exploit.df, n=1), aes(x=biomass, y=exploit, xmin=.lowerbio, xmax=.upperbio), height=0.1, color="white")+
    #geom_pointinterval(aes(x=biomass, y=exploit, ymin=.lowerexp, ymax=.upperexp, xmin=.lowerbio, xmax=.upperbio))+
    geom_segment(aes(x=biomass, y=exploit, xend=c(tail(biomass, n=-1), NA), yend=c(tail(exploit, n=-1), NA)), arrow=arrow(length=unit(0.5, "cm")))+
    geom_hline(yintercept=1, size=1)+
    geom_vline(xintercept=1, size=1)+
    scale_fill_manual(values = c("#ffc400","limegreen", "#FF4500", "#ffc400"))+
    coord_cartesian(xlim=c(0, 4), ylim=c(0, 2), expand=FALSE)   