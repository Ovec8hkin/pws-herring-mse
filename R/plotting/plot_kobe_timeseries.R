library(ggplot2)
library(ggdist)
library(tidyverse)

source(paste0(here::here("R/utils/"), "fun_read_dat.R"))


b.star <- 40000
f.star <- 0.20

set.seed(1998)
seeds <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4050, 4288, 4716, 4775, 6460, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)

#hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "constant.f.00")
nyr <- 25
sims <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4288, 4716, 4775, 7251, 7915, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)

control.rules <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

data <- data.frame(year=NA, biomass=NA, exploit=NA, sim=NA, cr=NA)
for(c in control.rules){
    for(s in sims){
        print(s)
        model.dir <- paste0(here::here("results/"), c, "/sim_", s, "/year_", nyr, "/model/")

        biomass.estimates <- read.true.biomass.data(model.dir)
        exploit.rate.estimates <- read.exploit.rates(model.dir)

        biomass.est.rel <- as_tibble(biomass.estimates/b.star) %>% pivot_longer(everything(), "year", values_to="biomass")
        exploit.est.rel <- as_tibble(exploit.rate.estimates/f.star) %>% pivot_longer(everything(), "year", values_to="exploit")

        d <- biomass.est.rel 
        d$exploit <- exploit.est.rel$exploit
        d$sim <- s
        d$cr <- c
        
        data <- rbind(data, d)
    }
}

cores <- parallel::detectCores()
cl <- makeCluster(min(cores[1]-1, length(hcr.names)), outfile="")
registerDoParallel(cl)

exploit.data <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
    source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
    e = read.true.exploit.rates(cr, seeds, nyr)
}, seeds=seeds, nyr=nyr, cl=cl)
exploit.data <- bind_rows(exploit.data)

unregister_dopar()
stopCluster(cl)

data$year <- as.numeric(data$year)

data.df <- exploit.data %>% na.omit() %>%
            group_by(year, control.rule) %>%
            summarise(
                biomass = median(biomass),
                exploit = median(exploit)
            ) %>%
            print(n=10)

rect_dat <- data.frame(panel = c("bottom_left", "top_right",
                                    "bottom_right", "top_left"),
                          x_min = c(-Inf, 1, 1, -Inf),
                          x_max = c(1, Inf, Inf, 1), y_min = c(-Inf,1, -Inf, 1),
                          y_max = c(1, Inf, 1, Inf))



kobe.plot <- ggplot(data.df[data.df$year >= 2022, ])+
    geom_rect(data=rect_dat, aes(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max, fill = panel))+
    geom_point(aes(x=biomass, y=exploit, color=cr), size=3)+
    #geom_label_repel(data=data.df[(as.numeric(data.df$year) %% 5)== 0,], aes(x=biomass, y=exploit, label=year), box.padding = 0.35, point.padding = 0.5, max.overlaps=5)+
    geom_hline(yintercept=1, size=1)+
    geom_vline(xintercept=1, size=1)+
    scale_fill_manual(values = c("orange","limegreen", "#FF4000", "#fcd526"))+
    coord_cartesian(xlim=c(0, 4), ylim=c(0, 2), expand=FALSE)

for(c in control.rules){
    kobe.plot <- kobe.plot + geom_segment(data=data.df[data.df$year >= 2022 & data.df$control.rule==c,], aes(x=biomass, y=exploit, color=control.rule, xend=c(tail(biomass, n=-1), NA), yend=c(tail(exploit, n=-1), NA), group=control.rule), arrow=arrow(length=unit(0.5, "cm")))
}

n.tot <- 600*total.sims

kobe.plot

kobe.df <- exploit.data %>% 
                mutate(
                    biomass = biomass/b.star,
                    exploit = exploit/f.star,
                    # kobe.color = ifelse(
                    #     exploit.data$biomass < 1 & exploit.data$exploit > 1,
                    #     "red",
                    #     ifelse(
                    #         exploit.data$biomass > 1 & exploit.data$exploit < 1,
                    #         "green",
                    #         ifelse(
                    #             exploit.data$biomass < 1 & exploit.data$exploit < 1,
                    #             "orange",
                    #             "yellow"
                    #         )
                    #     )
                    # )
                ) %>% 
                na.omit() %>%
                group_by(year, control.rule, kobe.color) %>%
                summarise(
                    n=n()
                ) %>%
                mutate(freq=n/n.tot) %>%
                mutate(across(kobe.color, factor, levels=c("red", "orange", "yellow", "green")))


ggplot(kobe.df) +
    geom_col(aes(x=year, y=freq, fill=kobe.color))+
    geom_vline(xintercept=2022, color="black", size=1, linetype="dashed")+
    scale_fill_manual(values=c("red", "orange", "#fedd1f", "limegreen")) +
    coord_cartesian(expand=FALSE)+
    scale_x_continuous(breaks=as.integer(seq(1980, 2022+nyr+1, by=10)))+
    labs(x="Year", y="Proportion of Outcomes", title="Kobe Timeseries")+
    facet_wrap(~cr, nrow=2)+
    theme(
        legend.position = "bottom",
        panel.spacing.x = unit(0.4, "in")
    )               

data.df.2 <- data %>% na.omit() %>%
                group_by(year, cr) %>%
                summarise(
                    biomass.quants=quantile(biomass, c(0.025, 0.5, 0.975)),
                    exploit.quants=quantile(exploit, c(0.025, 0.5, 0.975))
                ) %>%
                mutate(perc=c("2.5%", "50%", "97.5%")) %>%
                pivot_wider(
                    names_from = perc,
                    values_from = c(biomass.quants, exploit.quants)
                ) %>%
                filter(year == max(data$year, na.rm=TRUE)) %>%
                print(n=10)


ggplot(data.df.2)+
    geom_rect(data=rect_dat, aes(xmin = x_min, ymin = y_min, xmax = x_max, ymax = y_max, fill = panel))+
    geom_point(aes(x=`biomass.quants_50%`, y=`exploit.quants_50%`), size=3, color="white")+
    geom_errorbar(aes(x=`biomass.quants_50%`, y=`exploit.quants_50%`, ymin=`exploit.quants_2.5%`, ymax=`exploit.quants_97.5%`), width=0.5, size=1, color="white")+
    geom_errorbarh(aes(x=`biomass.quants_50%`, y=`exploit.quants_50%`, xmin=`biomass.quants_2.5%`, xmax=`biomass.quants_97.5%`), height=0.1, size=1, color="white")+
    gghighlight(use_direct_label = FALSE)+
    geom_hline(yintercept=1, size=1)+
    geom_vline(xintercept=1, size=1)+
    scale_fill_manual(values = c("orange","limegreen", "#FF4000", "#fcd526"))+
    coord_cartesian(xlim=c(0, 5), ylim=c(0, 2), expand=FALSE)+
    facet_wrap(~cr, nrow=2)+
    theme(legend.position = "bottom") 
