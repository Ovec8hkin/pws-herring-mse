library(tidyverse)
library(ggplot2)
library(ggdist)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
source(file=paste0(here::here("R/plotting/"), "plot_util_vals.r"))
source(file=paste0(here::here("R/plotting/"), "ggplot_facet_scaling.R"))

## Calculates Average Annual Variation in a metric
## -----------------------------------------------
## data - annual data for which to compute AAV
##        (assumed to be in order by year)
## ----------------------------------------------- 
aav <- function(data){
    total <- sum(data)
    diffs <- abs(diff(data))
    aav <- mean((sum(diffs)/total))
    return(ifelse(is.nan(aav), 0, aav)) # If all data is 0, return 0 rather than NA
}

## Converts performance metric data from wide to long format
## ---------------------------------------------------------
## metric.name - name of performance metric to reformat
## ---------------------------------------------------------
reformat.metric.df <- function(metric.name){
    return(
        catch.biomass.df %>% 
            select(
                c("control.rule", starts_with(metric.name), ".width", ".point", ".interval")
            ) %>%
            rename_with(
                ~ c("median", "lower", "upper"), 
                c(metric.name, paste0(metric.name, ".lower"), paste0(metric.name, ".upper"))
            ) %>%
            mutate(metric=metric.name) %>%
            relocate(
                c("control.rule", "metric", "median", "lower", "upper", ".width", ".point", ".interval")
            )

    )
}

set.seed(1)
seeds <- sample(1e4, 40)#c(1017, 4775, 9725, 8462, 8789, 8522, 1799, 8229, 1129, 878, 7845, 5922, 6526, 5071, 4650, 2159, 3476, 2580, 1530, 7289, 4633, 4344, 1222, 2858, 5400, 526, 1069)
nyr <- 30
hcr.names <- c("base", "low.harvest", "high.harvest", "low.biomass", "high.biomass", "constant.f.00", "evenness", "gradient", "three.step.thresh", "big.fish")


## NOTE THAT IF SOME SIMULATIONS FAILED AT SOME STAGE
## BIOMASS AND CATCH DATA FROM THAT SIMULATION SEED
## IS NOT INCLUDED IN THE DATA. THUS, SOME CONTROL RULES
## WILL HAVE MORE DATA BASED ON SIMULATION STABILITY.
#bio.traj.df <- data.frame(year=NA, biomass=NA, control.rule=NA, sim=NA)
#catch.data <- data.frame(year=NA, catch=NA, fishery=NA, control.rule=NA, sim=NA)

cores <- parallel::detectCores()
cl <- makeCluster(min(cores[1]-1, length(hcr.names)), outfile="")
registerDoParallel(cl)

bio.traj.df <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
    source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
    biomass.dat <- read.true.biomass.data(cr, seeds, nyr)
}, seeds=seeds, nyr=nyr, cl=cl)
bio.traj.df <- bind_rows(bio.traj.df)

catch.data <- pbapply::pblapply(hcr.names, function(cr, seeds, nyr){
    source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))
    biomass.dat <- read.catch.data(cr, seeds, nyr)
}, seeds=seeds, nyr=nyr, cl=cl)
catch.data <- bind_rows(catch.data)

unregister_dopar()
stopCluster(cl)

# for(cr in hcr.names){
#     print(cr)
#     cr.dat <- read.catch.data(cr, seeds, nyr)
#     catch.data <- catch.data %>% bind_rows(cr.dat)
#     #biomass.dat <- read.biomass.data(cr, seeds, 25)
#     biomass.dat <- read.true.biomass.data(cr, seeds, nyr)
#     bio.traj.df <- bio.traj.df %>% bind_rows(biomass.dat)
# }

# Compute AAV for all control rules and simulations
aav.df <- catch.data %>% na.omit() %>%                                          # Remove NAs
            group_by(control.rule, sim, year) %>%                               # Group
            summarise(
                total.catch = sum(catch)                                        # sum catches across fisheries for each year
            ) %>%
            group_by(control.rule, sim) %>%                                     # Group
            summarise(
                aav = aav(total.catch)                                          # Compute AAV
            ) %>%
            print(n=10)

# Compute probability of biomass being below the lower
# regulatory threshold for each control rule and simulation
limit.thresholds <- list(
    "base" = 19958,
    "low.harvest" = 19958,
    "high.harvest" = 19958,
    "low.threshold" = 10000,
    "high.threshold" = 30000,
    "evenness" = 19958,
    "gradient" = 19958,
    "three.step.thresh" = 19958,
    "big.fish" = 19958,
    "constant.f.00" = 1e10
)
prob.threshold.df <- as_tibble(bio.traj.df) %>% na.omit() %>%                   # Remove NAs
                    filter(year > 2021) %>%                                     # Only use years after 2021
                    group_by(control.rule, sim) %>%                             # Group
                    summarise(
                        n=n(),                                                  # Total number of samples
                        n.below=sum(biomass < 19958),                           # Number of samples where biomass below threshold
                        prob.below=n.below/n,                                   # Proportion of samples below threshold
                        n.closed=sum(biomass < limit.thresholds[control.rule][[1]]), # Number of samples where fishery is closed
                        prop.closed=n.closed/n                                  # Proportion of years where fishery is closed
                    ) %>%
                    print(n=10)

# Compute biomass and depletion level in final year of simulation
biomass.df <- as_tibble(bio.traj.df) %>% na.omit() %>%                          # Remove NAs
                filter(year == max(bio.traj.df$year, na.rm=TRUE)) %>%           # Only use final year of sim
                left_join(
                    as_tibble(bio.traj.df) %>% na.omit() %>%                          
                      filter(year == max(bio.traj.df$year, na.rm=TRUE) & control.rule == "constant.f.00") %>%
                      rename(no.fish.biomass = biomass) %>%
                      select(c(sim, no.fish.biomass)),
                    by = "sim"
                ) %>%
                mutate(
                  depletion = biomass/40000,                                      # Compute depletion from biomass
                  dyn.b0 = biomass/no.fish.biomass
                ) %>%                             
                select(-c(year)) %>%                                            # Drop year column
                relocate(c(control.rule, sim, biomass, depletion)) %>%          # Reorder columns
                print(n=10)

# Compute lowest depletion level for each control rule and
# simulation (at least 5 years into the sim).
lowest.depletion.df <- as_tibble(bio.traj.df) %>% na.omit() %>%                 # Reomve NAs
                        filter(year > 2021+5 & biomass > 0) %>%                 # Only use years 5 years into sim
                        left_join(
                          as_tibble(bio.traj.df) %>% na.omit() %>%                          
                            filter(year > 2021+5 & biomass > 0 & control.rule == "constant.f.00") %>%
                            rename(no.fish.biomass = biomass) %>%
                            select(c(year, sim, no.fish.biomass)),
                          by = c("year", "sim")
                        ) %>%
                        mutate(
                          depletion=biomass/40000,                             # Compute depletion from biomass
                          dyn.b0 = biomass/no.fish.biomass
                        ) %>%                     
                        group_by(control.rule, sim, year) %>%                   # Group
                        summarise(
                            depletion = median(depletion),                      # Median depletion by year
                            dyn.b0    = median(dyn.b0)
                        ) %>%
                        group_by(control.rule, sim) %>%                         # Group
                        summarise(
                            low.dep = min(depletion),                           # Lowest single year depletion level
                            low.year = as.numeric(year[depletion == low.dep]),  # Year in which low depletion occurred
                            low.dynb0 = min(dyn.b0),
                            low.year.db0 = as.numeric(year[dyn.b0 == low.dynb0])
                        ) %>%
                        print(n=100)

# Compute average biomass/depletion level for each control rule
# and simulation.
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

# Compute "biomass stability" (eg. biomass AAV) for each control
# rule and simulation.
biomass.stab <- as_tibble(bio.traj.df) %>% na.omit() %>%                        # Remove NAs
                    group_by(control.rule, sim) %>%                             # Group
                    summarise(
                        stab=aav(biomass)                                       # Biomass AAV
                    ) %>% print(n=10)

# Compute total and annual catch for each control rule and 
# simulation.
catch.df <- catch.data %>% na.omit() %>%                                        # Remove NAs
                group_by(control.rule, sim) %>%                                 # Group
                summarise(
                    tot_catch=sum(catch),                                       # Total catch by sim
                    ann_catch=tot_catch/nyr                                     # Annual catch by sim
                ) %>%
                print(n=10)

# Compute realized harvest rate for each control rule and 
# simulation.
harvest.rate.df <- as_tibble(catch.data) %>% na.omit() %>%                      # Remove NAs
                    group_by(control.rule, sim, year) %>%                       # Group
                    summarise(
                        ann_catch = sum(catch)                                  # Total catch by year
                    ) %>%
                    inner_join(                                                 # Join biomass data
                        as_tibble(bio.traj.df) %>% na.omit() %>%                # Remove NAs
                            mutate(year = as.numeric(year)-2022) %>%            # Convert calendar years to sim years
                            group_by(control.rule, sim, year) %>%               # Group
                            summarise(
                                biomass=median(biomass)                         # Mediam biomass by year
                            ),
                        by = c("control.rule", "sim", "year")
                    ) %>%
                    mutate(
                        harvest.rate = ann_catch/biomass                        # Compute annual harvest rate
                    ) %>%
                    group_by(control.rule, sim) %>%                             # Group
                    summarise(
                        harvest.rate = median(harvest.rate)                     # Median harvest rate
                    ) %>%
                    print(n=10)

# Join all performance metrics by control.rule and simulation
catch.biomass.df <- biomass.df %>% 
                    inner_join(catch.df,            by=c("control.rule", "sim")) %>%
                    inner_join(aav.df,              by=c("control.rule", "sim")) %>%
                    inner_join(prob.threshold.df,   by=c("control.rule", "sim")) %>%
                    inner_join(biomass.stab,        by=c("control.rule", "sim")) %>%
                    inner_join(avg.biomass.df,      by=c("control.rule", "sim")) %>%
                    inner_join(lowest.depletion.df, by=c("control.rule", "sim")) %>%
                    inner_join(harvest.rate.df,     by=c("control.rule", "sim")) %>%
                    mutate(control.rule=recode_factor(control.rule, !!!hcr.levels)) %>%
                    filter(no.fish.biomass != 0)

# Compute median, 50% (25-75), and 95% (2.5-97.5) confidence intervals
# for all performane metrics.
catch.biomass.df <- as_tibble(catch.biomass.df) %>%
                        select(-c(n, n.below, n.closed)) %>%
                        group_by(control.rule) %>%
                        median_qi(
                            ann_catch, dyn.b0, avg.db0, low.dynb0, aav, prob.below, prop.closed, harvest.rate, #, avg.dep, depletion, low.dep,
                            .width=c(0.50, 0.95)
                        )

# Reformat perforance data into long format for ease of plotting
perf.data <- 
    bind_rows(
        reformat.metric.df("ann_catch"),
        reformat.metric.df("aav"),
        reformat.metric.df("harvest.rate"),
        reformat.metric.df("dyn.b0"),
        reformat.metric.df("avg.db0"),
        reformat.metric.df("low.dynb0"),
        #reformat.metric.df("depletion"),
        #reformat.metric.df("avg.dep"),
        #reformat.metric.df("stab"),
        #reformat.metric.df("low.dep"),
        reformat.metric.df("prob.below"),
        reformat.metric.df("prop.closed")
        
    ) %>% 
    mutate(
        metric.long =
            recode_factor(
                metric, 
                !!!c(
                    tot_catch = "Total Catch",
                    ann_catch = "Annual Catch",
                    aav = "Average Annual Catch Variation",
                    harvest.rate = "Realized Harvest Rate",
                    biomass = "Final Year Biomass",
                    depletion = "Final Year Depletion Level",
                    avg.bio = "Average Biomass",
                    avg.dep = "Average Depletion Level",
                    dyn.b0 = "Final Year Biomass Relative to Unfished",
                    avg.db0 = "Average Biomass Relative to Unfished",
                    low.dynb0 = "Lowest Biomass Relative to Unfished",
                    stab = "Annual Biomass Variability",
                    low.dep = "Lowest Depletion Level",
                    prob.below = "Proportion of Years Below Threshold",
                    prop.closed = "Proportion of Year Fishery is Closed"
                )
            )
    )
#     print(n=75)

# perf.data %>% filter(.width == 0.95) %>%
#     group_by(control.rule) %>%
#     summarise(
#         metric = metric,
#         median = median,
#         lower = lower,
#         upper = upper
#     ) %>%
#     write_csv(file=paste0(here::here("results/"), "performance_summary.csv"))

ggplot(perf.data) +
        geom_pointinterval(aes(
            x = median,
            y = control.rule,
            xmin = lower,
            xmax = upper,
            color = control.rule
        ), point_size=6) +
        geom_vline(
            data = perf.data %>% filter(control.rule == "Default" & .width == 0.5),
            aes(xintercept=median)) +
        scale_color_manual(values=as.vector(hcr.colors)) +
        scale_y_discrete(limits=rev, labels=function(x) str_wrap(x, width=15)) +
        facet_wrap_custom(~metric.long, ncol=3, scale="free_x", shrink=TRUE, scale_overrides = list(
            #scale_override(1, scale_x_continuous(breaks=seq(0, 1000000, 100000),   labels=seq(0, 1000, 100),  limits = c(0, 1000000))),
            scale_override(1, scale_x_continuous(breaks=seq(0, 50000,  5000),      labels=seq(0, 50, 5),      limits = c(0, 50000))),
            scale_override(4, scale_x_continuous(breaks=seq(0, 1, 0.2),           labels=seq(0, 1, 0.2),     limits = c(0, 1))),
            scale_override(5, scale_x_continuous(breaks=seq(0, 1, 0.2),           labels=seq(0, 1, 0.2),     limits = c(0, 1))),
            scale_override(6, scale_x_continuous(breaks=seq(0, 1, 0.2),           labels=seq(0, 1, 0.2),     limits = c(0, 1))),
            scale_override(2, scale_x_continuous(breaks=seq(0, 2.0, 0.25),        labels=seq(0, 2, 0.25),   limits = c(0, 2))),
            #scale_override(4, scale_x_continuous(breaks=seq(0, 500000, 100000),   labels=seq(0, 500, 100),  limits = c(0, 500000))),
            #scale_override(4, scale_x_continuous(breaks=seq(0, 15, 1),            labels=seq(0, 15, 1),     limits = c(0, 15))),
            #scale_override(5, scale_x_continuous(breaks=seq(0, 7, 1),    labels=seq(0, 7, 1),   limits = c(0, 7))),
            #scale_override(5, scale_x_continuous(breaks=seq(0, 200000, 25000),    labels=seq(0, 200, 25),   limits = c(0, 200000))),
            #scale_override(6, scale_x_continuous(breaks=seq(0, 0.5, 0.1),         labels=seq(0, 0.5, 0.1),  limits = c(0, 0.5))),
            #scale_override(7, scale_x_continuous(breaks=seq(0, 2.0, 0.1),         labels=seq(0, 2.0, 0.1),  limits = c(0, 2.0))),
            scale_override(7, scale_x_continuous(breaks=seq(0, 1.0, 0.2),         labels=seq(0, 1.0, 0.2),  limits = c(0, 1.0))),
            scale_override(8, scale_x_continuous(breaks=seq(0, 1.0, 0.2),         labels=seq(0, 1.0, 0.2),  limits = c(0, 1.0))),
            scale_override(3, scale_x_continuous(breaks=seq(0, 1.0, 0.1),         labels=seq(0, 1.0, 0.1),  limits = c(0, 1.0)))
        ))+
        labs(x="", y="", title="Performance Metric Summaries")+ 
        theme(
            panel.grid.minor.x = element_blank(),
            legend.position = "none",
            plot.title = element_blank(),
            strip.text = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14)
        )

ggsave("/Users/jzahner/Desktop/plot.png")

## -------------------------------------
## Bentley et al. 2003 Objectve Function
## -------------------------------------

ms <- list(
    #tot_catch = 150000,
    ann_catch = 6000,
    aav = 1.0,
    dyn.b0 = 0.50,
    avg.db0 = 0.50,
    low.db0 = 0.35,
    #biomass = 20000,
    #avg.bio = 20000,
    #avg.dep = 0.35,
    #depletion = 0.35,
    #low.dep = 0.2,
    #stab = 0.5,
    prob.below = 0.2,
    prop.closed = 0.33,
    harvest.rate = 0.0
)

ls <- list(
    #tot_catch  = 500000,  # max annual catch 1970-1990 * nyr
    ann_catch  = 20000,   # approximate max annual catch 1970-1990
    aav        = 0.0,     # a constant catch/F rule has AAV 0
    dyn.b0     = 1.0,
    avg.db0    = 1.0,
    low.db0    = 1.5,     
    #biomass    = 80000,   # 50% biomass peak in 1990
    #avg.bio    = 80000,   # 50% of biomass peak
    #avg.dep    = 2.0,
    #depletion  = 2.0,     # matches avg_bio
    #low.dep    = 2.0,
    #stab       = 0.3,     # approximate median 
    prob.below = 0.0,      # we want the probability to be tiny,
    prop.closed = 0.0,
    harvest.rate = 0.3  
)

calc.utility <- function(value, metric){
    l <- ls[[metric]]
    m <- ms[[metric]]

    value <- as.numeric(value)

    if(metric %in% c("aav", "stab", "prob.below", "prop.closed")){
        value <- 1-value
        m <- 1-m
        l <- 1-l
    }
    
    print(metric)
    if(value < m) return(0)
    if(value > l) return(1)

    return((value - m)/(l-m))

}

total.utility <- function(utilities){
    n <- length(utilities)
    print(utilities)
    return(prod(as.numeric(utilities))^(1/n))
}

perf.matrix <- perf.data %>% select(control.rule, metric, median, lower, upper) %>% filter(metric %in% c("ann_catch", "dyn.b0", "aav", "prop.closed")) %>% as.matrix

utility.matrix <- perf.matrix
median.utilities <- apply(perf.matrix, 1, function(x) calc.utility(x[["median"]], x[["metric"]]))
lower.utilities <- apply(perf.matrix, 1, function(x) calc.utility(x[["lower"]], x[["metric"]]))
upper.utilities <- apply(perf.matrix, 1, function(x) calc.utility(x[["upper"]], x[["metric"]]))

utility.matrix[,3] <- as.numeric(median.utilities)
utility.matrix[,4] <- lower.utilities
utility.matrix[,5] <- upper.utilities

utility.matrix[,3:5] <- as.numeric(utility.matrix[,3:5])

utility.df <- as_tibble(utility.matrix) %>%
                group_by(control.rule) %>%
                summarise(
                    median.total.utility = total.utility(median),
                    lower.total.utility  = total.utility(lower),
                    upper.total.utility  = total.utility(upper)
                ) %>%
                mutate(
                  rel.util = median.total.utility/max(median.total.utility)
                ) %>%
                arrange(desc(rel.util)) %>%
                print(n=10)

