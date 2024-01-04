## Read data & Shaping
filename <- "local_anesteshia_data.csv"
drug_name <- c("Pro", "Lid", "Mep", "Bup", "Lid+Adr")

dat <- readr::read_csv(filename, show_col_types = FALSE) |>
    tidyr::pivot_longer(col = c(-ID, -Drug),
                        names_to = "time",
                        names_transform = list(time = as.numeric),
                        values_to = "score",
                        values_drop_na = TRUE) |>
    dplyr::mutate(Drug = factor(Drug, levels = drug_name),
                  drug1 = forcats::fct_collapse(Drug,
                              Pro = "Pro", Lid = c("Lid", "Lid+Adr")),
                  drug2 = ifelse(Drug == "Lid+Adr", 1L, 0L))

## Preparation of data for stan
dat_list <- list(
    N = dat |> nrow(),
    ND = dat$drug1 |> unique() |> length(),
    NG = dat$ID |> unique() |> length(),
    GROUP = dat$ID,
    DRUG1 = dat$drug1 |> as.integer(),
    DRUG2 = dat$drug2,
    TIME = dat$time,
    SCORE = dat$score
)

## Save raw and shaped data
dir.create("rds", showWarnings = FALSE)
save(dat, drug_name, dat_list, file = "rds/data_preparation.rds")



library(rstan)
options(mc.cores = parallel::detectCores())
load(file = "rds/data_preparation.rds")
dir.create("result", showWarnings = FALSE)


## model1
mod1 <- stan_model(file = "model/model1.stan")

## model2
## similar to model1
## Parameters of drugs (Mean) include the offset value for each individual
mod2 <- stan_model(file = "model/model2.stan")


## Condition of sampling
seed <- 1234
iter <- 10000
thin <- 10
warmup <- 2000
chains <- 4

## Estimate parametes
fit1 <- sampling(mod1, data = dat_list,
                 iter = iter, thin = thin, warmup = warmup, chains = chains,
                 seed = seed)

fit2 <- sampling(mod2, data = dat_list,
                 iter = iter, thin = thin, warmup = warmup, chains = chains,
                 seed = seed)

## Save models and stan results
save(mod1, fit1, file = "rds/data_fit1.rds")
save(mod2, fit2, file = "rds/data_fit2.rds")

