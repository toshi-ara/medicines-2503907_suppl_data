library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(survminer)
library(ggplot2)
library(RColorBrewer)
library(readr)
source("common_functions.R")

## prepare data for survival analysis
## group1: raw data
## group2: simulation result1 (6 trials)
## group3: simulation result2 (8 trials)
load("rds/data_preparation.rds")         # dat (raw data)
load("rds/data_simulation_model2.rds")   # res1(6 trials), res2 (8 trials)

drug_name <- c("Pro", "Lid", "Mep", "Bup", "Lid+Adr")
cond <- c("Original", "Parameter1", "Parameter2")
condLegend <- c("Raw Data", "Parameter1", "Parameter2")


## x: logical vector [TRUE/FALSE] in each time point
## time: vector of time
## SEQ:
## return: time when the effect of local anesthesia finished
##         (otherwise Inf)
getDuration <- function(x, time, SEQ = 3) {
    n <- length(time)
    idx <- lapply(1:(n-SEQ+1),
                  function(i) x[i:(i+SEQ-1)]) |>
        sapply(all) |>
        match(x = TRUE)+SEQ-1
    return(ifelse(is.na(idx), Inf, time[idx]))
}

## x: data.frame(time, score)
## k: number of at least k successes
## SEQ: 
## return: time when the effect of local anesthesia finished
##         (otherwise Inf)
getMatchCondition <- function(x, k, SEQ = 3) {
    return(getDuration((x$score >= k), x$time, SEQ = SEQ))
}


########################################
## Kaplan-Meier method
########################################

## status: 0 (alive), 1 (dead)
## return is tibble
groupSurvData <- function(x, k) {
    res <- x |>
        group_nest(Drug, ID) |>
        mutate(time = map_dbl(data, ~ getMatchCondition(., k)),
               status = ifelse(time == Inf, 0, 1)) |>
        select(-data) |>
        mutate(ID = as.integer(ID))
    return(res)
}


## prepare data for survival analysis
## group1: raw data
## group2: parameter1
## group3: parameter2
dat0 <- groupSurvData(dat, k = 6)
dat1 <- groupSurvData(res[[1]], k = 6)
dat2 <- groupSurvData(res[[2]], k = 6)

## combine data in each drug
datD <- bind_rows( dat0, dat1, dat2, .id = "Condition") |>
    mutate(Drug = factor(Drug, levels = drug_name),
           Condition = factor(Condition,
                              levels = seq_len(length(cond)),
                              labels = cond))

## Kaplan-Meier method
resKM <- survfit(Surv(time, status) ~ Drug + Condition, data = datD)

## Save table
tbl1 <- surv_summary(resKM) |>
    mutate(Drug = factor(Drug, levels = drug_name),
           Condition = factor(Condition, labels = condLegend)) |>
    group_by(Drug, Condition) |>
    summarise(n = max(n.risk), events = sum(n.event), strata = unique(strata))
tbl2 <- surv_median(resKM)

tbl <- left_join(tbl1, tbl2, by = "strata") |>
    select(-strata)

write_csv(file = "result/result_KM.csv", tbl)
save(datD, resKM, tbl,
     file = "rds/survival_analysis.rds")


## Plot
mytheme <- list(
    labs(x = "Time (min)",
         y = "Rate that the effect of drug persists",
         color = "Condition"),
    geom_hline(yintercept = 0.5, linetype = "dotted"),
    scale_color_brewer(palette = "Dark2", labels = condLegend),
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
    )
)

p <- ggsurvplot_facet(
    resKM, datD,
    size = 0.7, alpha = 0.8,
    surv.scale = "percent",
    censor.shape = "|",
    short.panel.labs = TRUE,
    panel.labs.font = list(size = 16),
    panel.labs.background = list(color = "white", fill = "white"),
    facet.by = "Drug", nrow = 2,
    ggtheme = theme_bw())

p$data$Drug <- factor(p$data$Drug, levels = drug_name)
p <- p + mytheme


## Save plots (PDF)
cairo_pdf(filename = "result/fig5.pdf",
          width = 7, height = 5, onefile = TRUE)
print(p)
dev.off()

