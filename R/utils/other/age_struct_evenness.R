library(ggplot2)
library(tidyverse)

source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

shanon.weiner.evenness <- function(age.struct){
    if(sum(age.struct) != 1.0){
        age.struct <- age.struct/sum(age.struct)
    }

    as.log <- log(age.struct)
    as.log[is.infinite(as.log)] <- 0

    H = -sum(age.struct*as.log)
    return(H/log(length(age.struct)))
}

generate.post.pred <- function(seine.age.comp, spawn.age.comp, seine.ess, spawn.ess, ncol=10){
    pp.seine.age.comp<-matrix(0, nrow(seine.age.comp), ncol(seine.age.comp)) # CHANGE
    pp.spawn.age.comp<-matrix(0, nrow(spawn.age.comp), ncol(spawn.age.comp))
    # sample using the MCMC draws 
    #set.seed(100)
    for(i in 1:nrow(seine.age.comp)){ # Loop through the MCMC draws
        for(j in 1:length(seine.ess)){ # Loop through each year
            pp.seine.age.comp[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1, size=as.integer(seine.ess[j]),seine.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)])) / seine.ess[j]*100
            if(all(spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)] == 0)){
                pp.spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)] <- rep(0, times=ncol)
            }else{
                pp.spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1, spawn.ess[j], spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)])) / spawn.ess[j] *100
            }
        }
    }
    return(list(
        seine.pp = pp.seine.age.comp,
        spawn.pp = pp.spawn.age.comp
    ))
}

model.dir <- paste0(here::here("results/base/sim_197/year_0/model/"))
years <- 1980:2021
raw.data <- read.data.files(model.dir)
model.data <- list(
    ess.seine = raw.data$PWS_ASA_ESS.ctl$seine_ess,
    ess.spawn = raw.data$PWS_ASA_ESS.ctl$spawn_ess
)

model.data$ess.seine[model.data$ess.seine == -9] <- 0
model.data$ess.spawn[model.data$ess.spawn == -9] <- 0

nyr=42
nburn=1

seine.age.comp<-read.csv(paste0(model.dir, "mcmc_out/SeAC.csv"), header = FALSE, dec=".") 
seine.age.comp<-seine.age.comp[-c(1:nburn), 1:(nyr*10)]*100
spawn.age.comp<-read.csv(paste0(model.dir, "mcmc_out/SpAC.csv"), header = FALSE, dec=".") 
spawn.age.comp<-spawn.age.comp[-c(1:nburn), 1:(nyr*10)]*100

pps <- generate.post.pred(seine.age.comp, spawn.age.comp, model.data$ess.seine, model.data$ess.spawn)

seine.age.comp.quants <- apply(pps$seine.pp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)
spawn.age.comp.quants <- apply(pps$spawn.pp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)

cnames <- rep(NA, ncol(seine.age.comp))
i=1
for(y in years){
    for(a in 0:9){
        cnames[i] <- paste0(y, "_", a)
        i = i+1
    }
}
colnames(seine.age.comp.quants) <- cnames
colnames(spawn.age.comp.quants) <- cnames

seac.df <- as_tibble(seine.age.comp.quants) %>%
            pivot_longer(everything(), names_to="year_age", values_to="val") %>%
            separate(year_age, c("year", "age"), sep="_") %>%
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) %>%
            pivot_wider(names_from=percentile, values_from=val) %>%
            filter(age > 2) %>%
            mutate(
                type="seine"
            ) %>%
            print(n=30)

spac.df <- as_tibble(spawn.age.comp.quants) %>%
            pivot_longer(everything(), names_to="year_age", values_to="val") %>%
            separate(year_age, c("year", "age"), sep="_") %>%
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) %>%
            pivot_wider(names_from=percentile, values_from=val) %>%
            filter(age > 2) %>%
            mutate(
                type="spawn"
            ) %>%
            print(n=30)

spac.25 <- spac.df %>% select(year, age, `2.5%`) %>%
            filter(!is.na(`2.5%`)) %>%
            pivot_wider(names_from=age, values_from=`2.5%`) %>%
            print(n=10)
spac.25 <- as.matrix(spac.25 %>% select(c(`3`, `4`, `5`, `6`, `7`, `8`, `9`)))
rownames(spac.25) <- 1982:2021

spac.50 <- spac.df %>% select(year, age, `50%`) %>%
            filter(!is.na(`50%`)) %>%
            pivot_wider(names_from=age, values_from=`50%`) %>%
            print(n=10)

spac.50 <- as.matrix(spac.50 %>% select(c(`3`, `4`, `5`, `6`, `7`, `8`, `9`)))
rownames(spac.50) <- 1982:2021

spac.975 <- spac.df %>% select(year, age, `97.5%`) %>%
            filter(!is.na(`97.5%`)) %>%
            pivot_wider(names_from=age, values_from=`97.5%`) %>%
            print(n=10)

spac.975 <- as.matrix(spac.975 %>% select(c(`3`, `4`, `5`, `6`, `7`, `8`, `9`)))
rownames(spac.975) <- 1982:2021

even.25 <- apply(spac.25, 1, shanon.weiner.evenness)
even.50 <- apply(spac.50, 1, shanon.weiner.evenness)
even.975 <- apply(spac.975, 1, shanon.weiner.evenness)

evenness.df <- data.frame(years=1982:2021, `2.5%`=even.25, `50%`=even.50, `97.5%`=even.975)

hist(even.50, breaks=seq(0, 1.0, 0.1), freq=TRUE, xlab="Shannon Evenness", main="PWS Herring AS Evenness")
