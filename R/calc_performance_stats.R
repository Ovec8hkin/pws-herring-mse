source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

total.sims <- 6
cr <- "base"
set.seed(1998)
sims <- sample(1:1e4, 6)
nyr <- 10
control.rules <- c("base", "low.harvest", "high.harvest", "low.biomass", "high.biomass", "lower.b0", "higher.b0", "constant.f.00")


catch.data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)
for(cr in control.rules){
    cr.dat <- read.catch.data(cr, sims, nyr)
    catch.data <- catch.data %>% bind_rows(cr.dat)
}

total.catch <- catch.data %>% na.omit() %>%
                group_by(year, fishery, control.rule) %>%
                summarise(catch = median(catch)) %>%
                group_by(control.rule) %>%
                summarise(tot.catch = sum(catch)) %>%
                print(n=10)

## Is average annual catch just the same thing as total catch
## divided by the number of years of the simulation?
avg.ann.catch <- catch.data %>% na.omit() %>%
                    group_by(year, fishery, control.rule) %>%
                    summarise(catch = median(catch)) %>%
                    group_by(year, control.rule) %>%
                    summarise(catch = sum(catch)) %>%
                    group_by(control.rule) %>%
                    summarise(avg.catch = mean(catch)) %>%
                    print(n=10)

avg.ann.catch.var <- catch.data %>% na.omit() %>%
                        group_by(year, fishery, control.rule) %>%
                        summarise(catch = median(catch)) %>%
                        group_by(year, control.rule) %>%
                        summarise(catch = sum(catch)) %>%
                        group_by(control.rule) %>%
                        summarise(aav = aav(catch)) %>%
                        print(n=10)


catch.stats <- total.catch %>% 
                left_join(avg.ann.catch) %>% 
                left_join(avg.ann.catch.var) %>% 
                print(n=10)


aav <- function(catches){
    total.catch <- sum(catches)
    catch.diffs <- abs(diff(catches))
    aav <- mean(catch.diffs/total.catch)*100
    return(ifelse(is.nan(aav), 0, aav))
}

#############################

read.biomass.data <- function(cr, sims, nyr){
    data <- data.frame(year=NA, biomass=NA, control.rule=NA, sim.num=NA) 
    for(s in sims){
        for(i in 1:nyr){
            fname <- paste0(here::here("results/"), cr, "/sim_", s, "/year_", i, "/model/mcmc_out/PFRBiomass.csv")
            biomass.data <- read_csv(fname, col_names=FALSE, show_col_types = FALSE) %>% select(last_col())
            bio.df <- data.frame(year=as.character(2022+i), biomass=pull(biomass.data), control.rule=cr, sim.num=s)
            data <- data %>% bind_rows(bio.df) %>% na.omit()
        }
    }
    return(data)
}

biomass.data <- data.frame(year=NA, biomass=NA, control.rule=NA, sim.num=NA)
for(cr in control.rules){
    cr.dat <- read.biomass.data(cr, sims, nyr)
    biomass.data <- biomass.data %>% bind_rows(cr.dat)
}

ref.point <- 40000

avg.depletion <- biomass.data %>% na.omit() %>%
                    mutate(depletion=biomass/ref.point) %>%
                    group_by(year, control.rule) %>%
                    summarise(depletion=median(depletion)) %>%
                    group_by(control.rule) %>%
                    summarise(avg.depletion=mean(depletion)) %>%
                    print(n=10)

final.depletion <- biomass.data %>% na.omit() %>%
                        filter(year == max(biomass.data$year, na.rm=TRUE)) %>%
                        mutate(depletion=biomass/ref.point) %>%
                        group_by(year, control.rule) %>%
                        summarise(depletion=median(depletion)) %>%
                        group_by(control.rule) %>%
                        summarise(fin.depletion=mean(depletion)) %>%
                        print(n=10)

biomass.threshold <- 20000
prob.below.thresh <- biomass.data %>% na.omit() %>%
                        filter(year == max(biomass.data$year, na.rm=TRUE)) %>%
                        group_by(control.rule) %>%
                        summarise(
                            below = sum(biomass < biomass.threshold),
                            n=n()
                        ) %>%
                        mutate(below.thresh=below/n) %>%
                        select(-c(below, n)) %>%
                        print(n=100)

biomass.stats <- avg.depletion %>% 
                    left_join(final.depletion) %>% 
                    left_join(prob.below.thresh) %>%
                    print(n=10)


perf.stats <- catch.stats %>% left_join(biomass.stats) %>% print(n=10)
