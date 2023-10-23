library(doParallel)
library(tictoc)
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

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

### THIS ACTUALLY CALLS THE MSE SIMULATION
run.mse <- function(cr, cr.name, seed, start.year=1, hindcast=FALSE){
    source(file=paste0(here::here("R/"), "mse_loop.R")) # source this here for parallel support
    
    sim.dir <- paste0(here::here("results"), "/", cr.name, "/sim_", seed, "/")
    # TODO: include some code to not repeat runs when HCR and sim seed are repeated
    if(!dir.exists(sim.dir)){ dir.create(sim.dir, recursive = TRUE) }
    res <- run.simulation(cr, nyr.sim, sim.seed=seed, 
                            write=sim.dir, start.year=start.year, hindcast=hindcast, cr.name=cr.name)
    if(!res$success){
        res <- run.simulation(cr, nyr.sim, sim.seed=seed, 
                        write=sim.dir, start.year=res$last.year-1, hindcast=hindcast, cr.name=cr.name)
    }
    return(res)
}

#clean.results.dir()

nyr.sim <- 30
total.sims <- 1

set.seed(99)
#seeds <- sample(1:1e4, size=total.sims) # 5552 7457 2294 2922 7102 3200  358 8973 9188 7071
#seeds <- c(2922)

seeds <- c(197, 649, 1017, 1094, 1144, 1787, 1998)#, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 2512, 3169, 3709, 4288, 4716, 4775, 6460, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)
control.rules <- list(
    base                = list(type="hcr.threshold.linear",      lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    high.harvest        = list(type="hcr.threshold.linear",      lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.40),
    low.harvest         = list(type="hcr.threshold.linear",      lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.10),
    low.biomass         = list(type="hcr.threshold.linear",      lower.threshold=10000, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    high.biomass        = list(type="hcr.threshold.linear",      lower.threshold=30000, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    evenness            = list(type="hcr.agestructure",          lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    gradient            = list(type="hcr.gradient",              lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20, p=0.5),
    constant.f.00       = list(type="hcr.constant.f",            f.rate=0.0),
    three.step.thresh   = list(type="hcr.threshold.multi", lower.threshold=19958, middle.threshold=38555, upper.threshold=60000, min.harvest=0.0, mid.harvest=0.20, max.harvest=0.5),
    big.fish            = list(type="hcr.threshold.weight",                lower.threshold=19958, upper.threshold=38555, min.harvest = 0.0, max.harvest=0.20),
    none = list()
)
hcr.names <- names(control.rules)

# Create all combinations of control rules and simulation seeds
cr.seed.combs <- expand.grid(seeds, hcr.names)
cr.seed.combs <- cr.seed.combs %>% filter(Var2 != "none")
total.sims <- nrow(cr.seed.combs)

# Set up parallel processing with doParallel
cores <- parallel::detectCores()
cl <- makeCluster(min(cores[1]-1, as.integer(total.sims/1)), outfile="")
registerDoParallel(cl)

tic()
foreach(i=1:total.sims) %dopar% {
    seed <- cr.seed.combs[i, 1]
    cr.name <- cr.seed.combs[i, 2]
    cr <- control.rules[[cr.name]]
    
    mse <- run.mse(cr, cr.name, seed, start.year=1)
    
}
toc()

unregister_dopar()
stopCluster(cl)

## Check for which control rule simulations failed
get.year <- function(f){
    return(as.numeric(str_split(f, "_")[[1]][2]))
}

results.dir <- "/Users/jzahner/Desktop/Projects/pwsh_mse/results/"

failed.simulations <- list()
for(cr in hcr.names){
    if(cr == "none") next;
    dir.name <- paste0(results.dir, cr)
    if(dir.exists(dir.name)){
        fs <- list.files(dir.name)
        for(f in fs){
            sim.name <- paste0(dir.name, "/", f)
            years <- list.files(sim.name)
            ys <- sort(apply(as.matrix(years), 1, get.year))
            final.year <- ys[length(ys)]
            if(final.year < nyr.sim){
                failed.simulations[[length(failed.simulations)+1]] <- list(cr=cr, sim=as.numeric(str_split(f, "_")[[1]][2]), fail.year=final.year)
            }
        }
    }
}
print(failed.simulations)

reruns <- list()
for(f in failed.simulations){
    if(f$sim %in% seeds){
        print(paste0(f$cr, ": ", f$sim, " (", f$fail.year, ")"))
        reruns[[length(reruns)+1]] <- f
    } 
}

## Rerun failed simulations a second time to see if some will finish this time
total.sims <- length(reruns)
cores <- parallel::detectCores()
cl <- makeCluster(min(cores[1]-1, as.integer(total.sims/1)), outfile="")
registerDoParallel(cl)
foreach(i=1:total.sims) %dopar% {
    fs <- reruns[i][[1]]
    cr.name <- fs$cr
    cr <- control.rules[[cr.name]]
    sim <- fs$sim
    yr <- fs$fail.year
    
    mse <- run.mse(cr, cr.name, sim, start.year=yr-1)
}
unregister_dopar()
stopCluster(cl)


##################
# total.sims <- length(seeds)
# for(n in 1:length(hcr.names)){
#     tic()
#     foreach(s=1:total.sims) %dopar% {

#         source(file=paste0(here::here("R/"), "mse_loop.R"))

#         sim.dir <- paste0(here::here("results"), "/", hcr.names[n], "/sim_", seeds[s], "/")
#         # TODO: include some code to not repeat runs when HCR and sim seed are repeated
#         if(!dir.exists(sim.dir)){ dir.create(sim.dir, recursive = TRUE) }
#         res <- run.simulation(control.rules[[n]], nyr.sim, sim.seed=seeds[s], 
#                               write=sim.dir, start.year=1, hindcast = FALSE, cr.name=hcr.names[n])
#         if(!res$success){
#             run.simulation(control.rules[[n]], nyr.sim, sim.seed=seeds[s], 
#                            write=sim.dir, start.year=res$last.year-1, hindcast=FALSE, cr.name=hcr.names[n])
#         }
#     }
#     toc()
# }