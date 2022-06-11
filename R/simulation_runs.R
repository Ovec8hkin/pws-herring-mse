library(doParallel)
library(tictoc)
library(ggplot2)
library(gridExtra)
library(tidyr)

source(file=paste0(here::here("R/"), "mse_loop.R"))
source(file=paste0(here::here("R/utils/", "fun_read_dat.R")))
source(file=paste0(here::here("R/utils/", "performance.R")))

clean.results.dir <- function(){
    dirs <- list.dirs(here::here("results"))
    for(d in dirs){
        unlink(d, recursive = TRUE, force=TRUE)
    }
}

#clean.results.dir()

nyr.sim <- 10
total.sims <- 2

#seeds <- rep(0, total.sims)
set.seed(1120)
seeds <- sample(1:1e4, size=total.sims)
control.rules <- list(
    base            = list(type="hcr.hockey.stick", lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    high.harvest    = list(type="hcr.hockey.stick", lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.30),
    low.harvest     = list(type="hcr.hockey.stick", lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.15),
    # lower.b0        = list(type="hcr.hockey.stick", lower.threshold=10000, upper.threshold=30000, min.harvest = 0.0, max.harvest=0.20),
    # low.biomass     = list(type="hcr.hockey.stick", lower.threshold=10000, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    # higher.b0       = list(type="hcr.hockey.stick", lower.threshold=30000, upper.threshold=50000, min.harvest = 0.0, max.harvest=0.20),
    # high.biomass    = list(type="hcr.hockey.stick", lower.threshold=30000, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    constant.f.00   = list(type="hcr.constant.f",   f.rate=0.0)
    ##constant.escape = list(type="hcr.constant.escapement", threshold=30000, proportion=1.0)
)
hcr.names <- names(control.rules)

# cores <- parallel::detectCores()
# cl <- makeCluster(cores[1]-1, type="FORK") #not to overload your computer
# registerDoParallel(cl)

for(n in 1:length(hcr.names)){
    tic()
    # foreach(s=1:total.sims) %dopar% {
    #     sim.dir <- paste0(here::here("results"), "/", hcr.names[n], "/sim_", seeds[s], "/")
    #     if(!dir.exists(sim.dir)){ dir.create(sim.dir, recursive = TRUE) }
    #     run.simulation(control.rules[[n]], nyr.sim, sim.seed=seeds[s], write=sim.dir)
    #     return()
    # }
    for(s in 1:total.sims){
        sim.dir <- paste0(here::here("results"), "/", hcr.names[n], "/sim_", seeds[s], "/")
        # TODO: include some code to not repeat runs when HCR and sim seed are repeated
        if(!dir.exists(sim.dir)){ dir.create(sim.dir, recursive = TRUE) }
        run.simulation(control.rules[[n]], nyr.sim, sim.seed=seeds[s], write=sim.dir, start.year=9)
    }
    
    #return()
    toc()
}
# stopCluster(cl)

# Accumulate catch data across all simulations for each of the four 
# major herring fisheries. Sum together to get a total catch for the
# entire fishery. Output is (nyr.sim * n.ages * total.sims)


get.harvest.rate.sim.results <- function(nyr.sim, total.sims, seeds, trial){
    harvest.rate <- accumulate.results.data(nyr.sim, total.sims, seeds, trial, c("harvest.csv"))
    harvest.median <- as.vector(apply(harvest.rate[[1]], c(1, 2), median))
    return(harvest.median)
}

get.harvest.results.all.trials <- function(nyr.sim, total.sims, seeds, trial){
    sim.results <- data.frame(Year = 1:nyr.sim)

    trials <- list.files(here::here("results"))
    for(t in trials){
        trial.results <- get.harvest.rate.sim.results(nyr.sim, total.sims, seeds, t)
        sim.results[,t] <- trial.results
    }

    return(sim.results)

}

harvest.rate <- get.harvest.results.all.trials(nyr.sim, total.sims, seeds, "")
harvest.rate
get.biomass.sim.results <- function(nyr.sim, total.sims, seeds, trial){

    prefish.spawn.biomass <- accumulate.results.data(nyr.sim, total.sims, seeds, trial, c("prefish_spawn_biomass.csv"), byage=TRUE)
    prefish.biomass.quants <- apply(apply(prefish.spawn.biomass[[1]], c(1, 3), sum), 1, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

    assessment.biomass <- accumulate.assessment.posteriors(nyr.sim, total.sims, seeds, trial)
    assessment.biomass.quants <- apply(assessment.biomass[[1]], 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    annual.assessment.biomass <- apply(assessment.biomass[[1]], c(2,3), median)

    prob.below.threshold <- compute.prob.below.threshold(annual.assessment.biomass, 19958)

    data <- cbind(prefish.biomass.quants, assessment.biomass.quants)

    sim.results <- data.frame(year=2021:(2021+nyr.sim-1), rep(trial, nyr.sim), rep(c("det", "ass"), each=nyr.sim), t(data), thresh.prob=prob.below.threshold)
    names(sim.results) <- c("Year", "Trial", "Type", "Biomass2.5", "Biomass25", "Biomass50", "Biomass75", "Biomass97.5", "thresh.prob")

    return(sim.results)
}

get.biomass.results.all.trials <- function(nyr.sim, total.sims, seeds){
    trials <- list.files(here::here("results"))
    sim.results <- data.frame()
    for(tr in trials){
        print(tr)
        trial.results <- get.biomass.sim.results(nyr.sim, total.sims, seeds, tr)
        sim.results <- rbind(sim.results, trial.results)
    }
    return(sim.results)
}

sim.results <- get.biomass.results.all.trials(nyr.sim, total.sims, seeds) #get.biomass.sim.results(nyr.sim, total.sims, seeds, "test")

#base.sim.results <- get.biomass.sim.results(nyr.sim, total.sims, seeds, "base")
#no.harvest.sim.results <- get.biomass.sim.results(nyr.sim, total.sims, seeds, "no.harvest")

sim.results$Trial <- factor(sim.results$Trial, 
                            levels=c("base",    "high.harvest", "low.harvest",  "lower.b0", "low.biomass",      "higher.b0",    "high.biomass",     "constant.f.00"), 
                            labels=c("Current", "High F",       "Low F",        "Low B0",   "Low Threshold",    "High B0",      "High Threshold",   "No Fishing"))

ggplot(sim.results)+
    geom_line(aes(x=Year, y=20000), color="blue", alpha=0.5)+
    geom_line(aes(x=Year, y=40000), color="blue", alpha=0.5)+
    geom_point(aes(x=Year, y=Biomass50, linetype=Type, color=Type))+
    geom_line(aes(x=Year, y=Biomass50, linetype=Type, color=Type))+
    geom_ribbon(aes(x=Year, ymin=Biomass25, ymax=Biomass75, linetype=Type, fill=Type), alpha=0.25)+
    geom_ribbon(aes(x=Year, ymin=Biomass2.5, ymax=Biomass97.5, linetype=Type, fill=Type), alpha=0.125)+
    geom_point(aes(x=Year, y=thresh.prob*60000), color="purple")+
    geom_line(aes(x=Year, y=thresh.prob*60000), color="purple")+
    scale_color_manual(values=c("red", "blue"))+
    scale_y_continuous(
        name = "Spawning Biomass (metric tons)",
        limits = c(0, 120000),
        sec.axis = sec_axis(trans=~.*1/80000, name="Probability below 20k metric tons")
    )+xlab("Year")+ggtitle("Spawning Biomass Predictions")+facet_wrap(vars(Trial))

ssb.traj.plot <- ggplot(sim.results)+
                    #geom_ribbon(aes(x=Year, ymin=Biomass25, ymax=Biomass75, fill=Trial), alpha=0.25)+
                    geom_point(aes(x=Year, y=Biomass50, color=Trial, linetype=Type), size=3)+
                    geom_line(aes(x=Year, y=Biomass50, color=Trial, linetype=Type), size=1.5)+
                    #geom_ribbon(aes(x=Year, ymin=Biomass2.5, ymax=Biomass97.5, fill=Trial), alpha=0.125)+
                    ylim(0, 75000)+
                    scale_color_manual(values=c("black", "red", "blue", "green", "purple", "orange", "darkslategray3", "coral4"))+
                    scale_fill_manual(values=c("black", "red", "blue", "green", "purple", "orange", "darkslategray3", "coral4"))+
                    labs(title="SSB Trajectories Under Candidate HCRs", 
                         x="Year of Simulation", 
                         y="Spawning Stock Biomass (mt)", 
                         color="HCR")
                    # theme(axis.text=element_text(size=18),
                    #       axis.title=element_text(size=22),
                    #       axis.title.y=element_text(margin=margin(r=20)),
                    #       axis.title.x=element_text(margin=margin(t=20)),
                    #       plot.title=element_text(size=32),
                    #       legend.text=element_text(size=18),
                    #       legend.key.size=unit(1.0, 'cm'),
                    #       legend.title=element_text(size=22))
ssb.traj.plot

# Performance Metrics
# compute.total.catch(total.catch)
# compute.average.catch(total.catch)
# compute.catch.aav(total.catch)
# compute.prob.threshold.final(annual.spawning.biomass, 19958)
# compute.number.of.closures(annual.spawning.biomass, 19958)
# compute.final.biomass(prefish.spawn.biomass)






compute.catch.performance.metrics <- function(nyr.sim, total.sims, seeds, trial, qs=c(0.025, 0.50, 0.975)){
    catches <- accumulate.results.data(nyr.sim, total.sims, seeds, trial,
                                   c("seine_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "foodbait_catch.csv"), 
                                   byage=TRUE)
    total.catch <- Reduce("+", catches)

    return(list(
        total.catch     = quantile(compute.total.catch(total.catch), prob=qs),
        average.catch   = quantile(compute.average.catch(total.catch), prob=qs),
        catch.aav       = quantile(compute.catch.aav(total.catch),  prob=qs)
    ))

}

catch.performance.metrics.all.trials <- function(nyr.sim, total.sims, seeds, qs=c(0.025, 0.50, 0.975)){
    trials <- list.files(here::here("results"))
    sim.results <- data.frame()

    metric.names <- c("total.catch", "average.catch", "catch.aav")

    for(m in metric.names){
        for(t in trials){
            perf.metric <- compute.catch.performance.metrics(nyr.sim, total.sims, seeds, t, qs)[[m]]
            data <- c(m, t, perf.metric)
            sim.results <- rbind(sim.results, data)
        }
    }
    colnames(sim.results) <- c("Metric", "Trial", "2.5%", "50%", "97.5%")

    return(sim.results)
    
}
compute.catch.performance.metrics(nyr.sim, total.sims, seeds, "base")
catch.performance.metrics.all.trials(nyr.sim, total.sims, seeds)
















cumulative.catch <- function(catch.quantiles){
    cum.catch <- array(dim=ncol(catch.quantiles))
    cum.catch[1] <- 0
    for(i in 2:ncol(catch.quantiles)){
        cum.catch[i] <- cum.catch[i-1]+catch.quantiles["50%", i]
    }
    return(cum.catch)
}

cum.catch.df <- data.frame(Year=1:25)
for(h in hcr.names){
    catches <- accumulate.results.data(nyr.sim, total.sims, seeds, h,
                                   c("seine_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "foodbait_catch.csv"), 
                                   byage=TRUE)
    total.catch <- Reduce("+", catches)
    catch.quantiles <- apply(apply(total.catch, c(1, 3), sum), 1, quantile, prob=c(0.025, 0.50, 0.975)) 
    cum.catch.df[h] <- cumulative.catch(catch.quantiles)
}

# catches <- accumulate.results.data(nyr.sim, total.sims, seeds, "constant.f.00",
#                                    c("seine_catch.csv", "gillnet_catch.csv", "pound_catch.csv", "foodbait_catch.csv"), 
#                                    byage=TRUE)
# total.catch <- Reduce("+", catches)
# catch.quantiles <- apply(apply(total.catch, c(1, 3), sum), 1, quantile, prob=c(0.025, 0.50, 0.975)) 

# low.biomass.cum.catch <- cumulative.catch(catch.quantiles)
# lower.b0.cum.catch <- cumulative.catch(catch.quantiles)
# low.harvest.cum.catch <- cumulative.catch(catch.quantiles)
# high.harvest.cum.catch <- cumulative.catch(catch.quantiles)
# base.cum.catch <- cumulative.catch(catch.quantiles)
# constant.f.00.cum.catch <- cumulative.catch(catch.quantiles)

# cum.catch.df <- data.frame(base=t(t(base.cum.catch)), 
#                            low.harvest=t(t(low.harvest.cum.catch)),
#                            high.harvest=t(t(high.harvest.cum.catch)),
#                            constant.f.00=t(t(constant.f.00.cum.catch)),
#                            low.b0=t(t(lower.b0.cum.catch)),
#                            low.biomass=t(t(low.biomass.cum.catch)),
#                            Year=1:25) %>% pivot_longer(!Year, names_to="control.rule", values_to="cum.catch")

cum.catch.df$control.rule <- factor(cum.catch.df$control.rule, levels=c("base", "high.harvest", "low.harvest", "low.b0", "low.biomass", "constant.f.00"), labels=c("Current", "High F", "Low F", "Low B0", "Low Threshold", "No Fishing"))
cum.catch.df <- as.data.frame(cum.catch.df)

cum.catch.plot <- ggplot(cum.catch.df)+
                    geom_point(aes(x=Year, y=cum.catch, color=control.rule), size=2)+
                    geom_line(aes(x=Year, y=cum.catch, color=control.rule), size=1.5)+
                    ylim(0, 350)+
                    scale_color_manual(values=colors)+
                    labs(title="Cumulative Catch Under Candidate HCRs", 
                         x="Year of Simulation", 
                         y="Cumulative Catch (mt)", 
                         color="HCR")

grid.arrange(ssb.traj.plot, cum.catch.plot, ncol=2)
