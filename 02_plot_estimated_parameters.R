library(rstan)
library(dplyr)
library(tidyr)
library(patchwork)
library(readr)
library(xtable)
library(loo)
source("common_functions.R")

dir.create("result", showWarnings = FALSE)
dir.create("result/table", showWarnings = FALSE)

load(file = "rds/data_fit1.rds")
load(file = "rds/data_fit2.rds")
drug_name <- c("Pro", "Lid", "Mep", "Bup", "Lid+Adr")


## Model1
### traceplot
p_t1 <- stan_trace(fit1, pars = c('mu0', 'log_sigma0', 'adr'),
                   inc_warmup = TRUE, ncol = 4) +
    ggtitle("Model 1")
### density plot
p_d1 <- stan_dens(fit1, pars = c('mu0', 'log_sigma0', 'adr'),
              inc_warmup = FALSE, ncol = 4)
p1 <- p_t1 + p_d1

## Model2
### traceplot
p_t2 <- stan_trace(fit2, pars = c('mu0', 'log_sigma0', 'adr'),
                   inc_warmup = TRUE, ncol = 4) +
    ggtitle("Model 2")
### density plot
p_d2 <- stan_dens(fit2, pars = c('mu0', 'log_sigma0', 'adr'),
              inc_warmup = FALSE, ncol = 4)

p2 <- p_t2 + p_d2


## Save as PDF
## only parameters of fixed effect model
cairo_pdf(filename = "result/SFig2.pdf",
          width = 10, height = 5, onefile = TRUE)
print(p1)
print(p2)
dev.off()


## Get parameters from posterior distribution (mean)
N <- 51
pars <- c("mu0", "s_mu0", "log_sigma0", "log_s_sigma0")

## model1
tmp <- get_posterior_mean(fit1, pars = pars)[,"mean-all chains"] |>
    matrix(ncol = 4) |>
    data.frame()
colnames(tmp) <- pars

param1 <- tmp |>
    mutate(sigma0 = exp(log_sigma0),
           adr = get_posterior_mean(fit1, pars = "adr")[,"mean-all chains"],
           drug = factor(drug_name[1:4])) |>
    select(drug, mu0, s_mu0, sigma0, everything()) |>
    tibble()

## model2
d <- get_posterior_mean(fit2, pars = "d")[,"mean-all chains"]
tmp <- get_posterior_mean(fit2, pars = pars)[,"mean-all chains"] |>
    matrix(ncol = 4) |>
    data.frame()
colnames(tmp) <- pars

param2 <- tmp |>
    mutate(sigma0 = exp(log_sigma0),
           d_mean = mean(d), d_sd = sd(d),
           adr = get_posterior_mean(fit2, pars = "adr")[,"mean-all chains"],
           drug = factor(drug_name[1:4])) |>
    select(drug, mu0, s_mu0, sigma0, everything()) |>
    tibble()


## parameters of each individual (random effect model)
## model1
param1_individual <- tibble(
    drug = rep(drug_name[1:4], each = N),
    ID = rep(seq_len(N), 4),
    Mean = get_posterior_mean(fit1, pars = "mu")[,"mean-all chains"],
    Sigma = get_posterior_mean(fit1, pars = "sigma")[,"mean-all chains"]
)

## model2
param2_individual <- tibble(
    drug = rep(drug_name[1:4], each = N),
    ID = rep(seq_len(N), 4),
    Mean = get_posterior_mean(fit2, pars = "mu")[,"mean-all chains"],
    Sigma = get_posterior_mean(fit2, pars = "sigma")[,"mean-all chains"]
)


save(param1, param1_individual, param2, param2_individual,
     file = "rds/estimated_parameters.rds")

write_csv(file = "result/param1.csv", param1)
write_csv(file = "result/param1_individual.csv", param1_individual)

write_csv(file = "result/param2.csv", param2)
write_csv(file = "result/param2_individual.csv", param2_individual)


## WAIC / LOO-CV
## Model1
log_lik1 <- extract_log_lik(fit1)
waic1 <- waic(log_lik1)
loo1 <- loo::loo(log_lik1)

## Model2
log_lik2 <- extract_log_lik(fit2)
waic2 <- waic(log_lik2)
loo2 <- loo::loo(log_lik2)

## compare two models
loo_compare(waic1, waic2)
loo_compare(loo1, loo2)

