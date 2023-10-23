library(ggplot2)
library(ggdist)
library(here)
library(dplyr)
library(magrittr)
library(tidyverse)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
source(file=paste0(here::here("R/plotting/"), "plot_util_vals.R"))

get.historical.catches <- function(){

    a<-read.table(file.path(here::here(), "supp_data", "historical_FoodBaitCatch.csv"), header=TRUE, sep=",")
    b<-read.table(file.path(here::here(), "supp_data", "historical_PoundCatch.csv"), header=TRUE, sep=",")
    c<-read.table(file.path(here::here(), "supp_data", "historical_GillnetCatch.csv"), header=TRUE, sep=",")
    foodbait<-rowSums(a[,-1])
    pound<-rowSums(b[,-1])
    gillnet<-rowSums(c[,-1])

    # in metric tonnes
    seineYield<-read.table(file.path(here::here(), "supp_data", "historical_SeineYieldCatch.csv"), header=TRUE, sep=",")
    seineYield<-seineYield[,-1]
    seineYield<-replace(seineYield, seineYield == 0, NA)
    # Years of the model
    years<-seq(1980, 2012)

    # Weight-at-age
    weight<-read.table(file.path(here::here(), "supp_data", "historical_Weight.csv"), header=TRUE, sep=",")
    weight<-weight[,-1]
    a<-a[,-1]
    abiomass<-weight*a
    b<-b[,-1]
    bbiomass<-weight*b
    c<-c[,-1]
    cbiomass<-weight*c
    aSum<-rowSums(abiomass)
    aSum<-replace(aSum, aSum == 0, NA)
    bSum<-rowSums(bbiomass)
    bSum<-replace(bSum, bSum == 0, NA)
    cSum<-rowSums(cbiomass)
    cSum<-replace(cSum, cSum == 0, NA)
    # Matrix of catches by gear type in mt, used below
    Tot<-cbind(aSum,bSum,cSum, seineYield)
    # Total catches by year in mt, used below
    Tot1<-Tot
    Tot1[is.na(Tot1)]<-0 
    TotalCatchMass<-rowSums(Tot1)

    archDataYrs<-as.numeric(1900:2021)
    archDataYrs<-1900:2021

    archData<-rep(0,121) # so that plugging numbers into indices out of order is easier to handle
    archData[14:27] <- c(100, 200, 200, 300, 4000, 5000, 9000, 8500, 20000, 16000, 8500, 12000, 5500, 5000) 
    archData[28:34]<-c(7000, 6500, 14800, 7000, 10900, 14800, 22000)
    archData[35:43] <- c(40000, 41000, 42000, 51000, 48000, 28000, 31000, 1500, 6000)
    # 1938 should be 51000
    archData[44:51] <- c(21000, 16000, 18000, 2000, 18000, 250, 22000, 20000)
    archData[52:72]<-c( 3000, 350, 9000, 16000, 7000, 10000, 4000, 0, rep(50, 7), 100, 100, 300, 50, 800, 1900)
    # the zero should be 1959, and 1969 should be 300
    archData[73:79]<-c(6400, 5600, 5400, 2300, 2200, 2400, 4500)
    archData[80:112] <- TotalCatchMass
    return(archData)
}

nyr <- 25
seeds <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4050, 4288, 4716, 4775, 6460, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)

hcr.names <- c("base", "high.harvest", "low.harvest", "high.biomass", "low.biomass", "evenness", "gradient", "three.step.thresh", "big.fish", "constant.f.00")

cores <- parallel::detectCores()
cl <- makeCluster(min(cores[1]-1, length(hcr.names)), outfile="")
registerDoParallel(cl)

catch.data <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
    source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
    biomass.dat <- read.catch.data(cr, seeds, nyr)
}, seeds=seeds, nyr=nyr, cl=cl)
catch.data <- bind_rows(catch.data)

unregister_dopar()
stopCluster(cl)

catch.data <- catch.data %>%
                mutate(control.rule = recode_factor(control.rule, !!!hcr.levels))

low.regime.catch <- catch.data[catch.data$year %in% c(1, 2, 3, 4, 5, 21, 22, 23, 24, 25),]
high.regime.catch <- catch.data[catch.data$year %in% seq(6, 20, 1),]

hist.catch.seq <- get.historical.catches()
hist.catches <- data.frame(
    year=rep(1900:2020, length(hcr.names)*length(seeds)),
    control.rule = rep(factor(hcr.names, labels=hcr.levels), each=length(1900:2020)*length(seeds)),
    sim=rep(seeds, length(1900:2020)*length(hcr.names)),
    total.catch=rep(hist.catch.seq, length(hcr.names)*length(seeds))
)

# Total catch (all fisheries) by control rule
total.catch.df <- catch.data %>% na.omit() %>%
                    mutate(year = 2021+year) %>%
                    group_by(year, control.rule, sim) %>%
                    summarise(
                        total.catch = sum(catch)
                    ) %>%
                    bind_rows(hist.catches) %>%
                    group_by(year, control.rule) %>%
                    median_qi(total.catch, .width=c(0.5, 0.95)) %>%
                    print(n=100)

ggplot(total.catch.df) + 
    geom_col(
        data = total.catch.df %>% filter(.width == 0.5),
        aes(x=year, y=total.catch)
    ) +
    geom_pointrange(aes(x=year, y=total.catch, ymax=.upper, ymin=.lower))+
    facet_wrap(~control.rule)

ggplot(total.catch.df, aes(x=year, y=total.catch, group=control.rule, color=control.rule))+
    geom_line(
        data = total.catch.df %>% filter(.width == 0.5),
        size=0.9
    )+
    geom_line(
        data = total.catch.df %>% filter(.width == 0.5 & control.rule == "Default"),
        color="black",
        size=0.9
    )+
    scale_color_manual(values=as.vector(hcr.colors))+
    scale_x_continuous("Year", breaks=seq(1900, 2050, 20))+
    scale_y_continuous("Annual Catch (mt)", breaks=seq(0, 60000, 10000), labels = scales::comma)+
    ggtitle("Annual Catch of PWS Pacific Herring")+
    coord_cartesian(expand=0, ylim=c(0, 60000))+
    theme(
        panel.background = element_blank(),
        axis.line = element_line()
    )


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
