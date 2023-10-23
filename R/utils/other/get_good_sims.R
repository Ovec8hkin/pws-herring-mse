nyr <- 30

#sims <- c(197, 649, 1017, 1094, 1144, 1787, 1998, 2078, 2214, 2241, 2255, 2386, 3169, 3709, 4288, 4716, 4775, 6460, 7251, 7915, 8004, 8388, 8462, 8634, 8789, 8904, 8935, 9204, 9260, 9716, 9725)
set.seed(1120)
sims <- sample(1e4, 150)


get.good.sims <- function(){
    nyr <- 30
    hcr.names <- c("base", "low.harvest", "high.harvest", "low.biomass", "high.biomass", "constant.f.00", "evenness", "gradient", "three.step.thresh", "big.fish")

    bio.traj.df <- read_csv(file.path(here::here(), "results", "om_biomass.csv"), col_names = TRUE, show_col_types = FALSE) %>%
    select(year, biomass, control.rule, sim)

    final.year <- 1980+42+nyr-1

    sims <- unique(bio.traj.df$sim)

    d <- data.frame(matrix(0, nrow=length(sims), ncol=length(hcr.names)))
    rownames(d) <- as.character(sims)
    colnames(d) <- hcr.names

    for(s in sims){
        for(cr in hcr.names){
            b.final <- bio.traj.df %>% filter(sim == s & control.rule == cr & year==final.year) %>% slice_head(n=1) %>% pull(biomass) %>% as.integer
            d[as.character(s), cr] <- as.integer(b.final > 0)
        }
    }

    d$successful <- rowSums(d)
    nrow(d[d$successful == 10,])

    good.sims <- as.integer(rownames(d[d$successful == 10,]))

    return(good.sims)
}



bio.traj.df <- read_csv(file.path(here::here(), "results", "om_biomass_1.csv")) %>%
                    bind_rows(read_csv(file.path(here::here(), "results", "om_biomass_2.csv"))) %>%
                    write_csv(file.path(here::here(), "results", "om_biomass.csv"))

catch.traj.df <- read_csv(file.path(here::here(), "results", "om_catch_1.csv")) %>%
                    bind_rows(read_csv(file.path(here::here(), "results", "om_catch_2.csv"))) %>%
                    write_csv(file.path(here::here(), "results", "om_catch.csv"))
