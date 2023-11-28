library(tidyverse)
library(ggdist)

filename <- here::here("supp_data/asl_survey_data_1973-2021.csv")

asl.data <- read_csv(filename, show_col_types=FALSE) %>%
              select(Year, GI, LENGTH, WEIGHT, AGE) %>%
              filter(!is.na(GI) & GI < 7 & AGE < 16 & AGE > 2) %>%
              mutate(
                AGE = case_when(AGE >= 9 ~ 9, AGE < 9 ~ AGE)
              ) %>%
              mutate(
                AGE = as.factor(AGE)
              ) %>%
              print(n=10)

asl.summary.data <- asl.data %>% na.omit() %>%
                      group_by(AGE) %>%
                      summarise(
                        n=n(),
                        mean = mean(WEIGHT),
                        sd = sd(WEIGHT),
                        lower = quantile(WEIGHT, c(0.025)),
                        upper = quantile(WEIGHT, c(0.975)),
                        median = median(WEIGHT)
                      ) %>%
                      print(n=10)

weight.threshold <- 110
asl.bigfish <- asl.data %>% na.omit() %>%
                mutate(tot.prob = sum(WEIGHT > weight.threshold)/n()) %>%
                group_by(AGE) %>%
                summarise(
                  n=n(),
                  n.above = sum(WEIGHT > weight.threshold),
                  prob = round(n.above/n, 2),
                  tot.prob = mean(tot.prob)
                ) %>%
                print(n=10)

ggplot(asl.data) +
    stat_halfeye(aes(x=WEIGHT, y=AGE, fill=stat(abs(x) > weight.threshold))) +
    geom_vline(aes(xintercept=weight.threshold), size=0.5) + 
    geom_text(data = asl.bigfish, aes(x=315, y=AGE, label=prob), size=5)+
    scale_fill_manual(values=c("gray", "#008800")) +
    scale_x_continuous(limits = c(0, 325), breaks=seq(0, 300, 50))+
    labs(x="Weight (g)", y="Age", fill=NA, title="Age-Weight Distributions")+
    theme_classic()+
    theme(
      panel.grid.major.y = element_line(),
      legend.position="none",
      plot.title = element_text(size=16),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
    )

#ggsave(file.path(here::here(), "figures", "publication", "Fig3_weightdist.jpg"), dpi=300, width=170, height=105, units="mm")
ggsave(file.path(here::here(), "figures", "publication", "Fig3_weightdist.pdf"), dpi=300, width=170, height=105, units="mm")

  # ggplot(asl.data) +
  #   stat_slab(aes(x=WEIGHT, y=AGE, fill=after_stat(abs(x) > weight.threshold), thickness = after_stat(pdf)), scale = 0.7) +
  #   stat_dotsinterval(aes(x=WEIGHT, y=AGE, fill=stat(abs(x) > weight.threshold)), side="bottom", scale=0.4) +
  #   geom_vline(aes(xintercept=weight.threshold), size=0.5) + 
  #   geom_text(data = asl.bigfish, aes(x=315, y=AGE, label=prob), size=5)+
  #   scale_fill_manual(values=c("gray", "#008800")) +
  #   scale_x_continuous(limits = c(0, 325), breaks=seq(0, 300, 50))+
  #   labs(x="Weight (g)", y="Age", fill=NA, title="Age-Weight Distributions")+
  #   theme(
  #     legend.position="none"
  #   )
