## Plot with fitted curves
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("result", showWarnings = FALSE)
dir.create("result/SVG", showWarnings = FALSE)

load(file = "rds/data_preparation.rds")
load(file = "rds/data_fit1.rds")
load(file = "rds/data_fit2.rds")
load(file = "rds/estimated_parameters.rds")
source("common_functions.R")
drug_name <- c("Pro", "Lid", "Mep", "Bup", "Lid+Adr")


## Function to get parameters from the result by stan
## x: raw data
## paramFix: estimated parameters (fix effect model)
## paramRnd: estimated parameters of each individual
set_param <- function(x, paramFix, paramRnd, drug_name) {
    param0 <- paramFix |>
        add_row(paramFix[2,]) |>
        select(mu0, sigma0, adr)
    param0$Drug <- factor(drug_name, levels = drug_name, labels = drug_name)
    param0$adr[1:4] <- 0

    param1 <- paramRnd |>
        rename(drug1 = drug, mu = Mean, sigma = Sigma)

    res <- x |>
        select(ID, Drug, drug1, drug2) |>
        unique() |>
        left_join(param0, by = "Drug") |>
        left_join(param1, by = c("drug1", "ID")) |>
        mutate(Drug = factor(Drug, levels = drug_name, labels = drug_name))
    return(res)
}

## Function to get predicted value from parameters
## x: raw data
## time: vector of time (minute)
## paramFix: estimated parameters (fix effect model)
## paramRnd: estimated parameters of each individual
get_pred <- function(x, time, paramFix, paramRnd, drug_name, score = 6) {
    pred <- set_param(x, paramFix, paramRnd, drug_name) |>
        expand_grid(time) |>
        mutate(pred0 = predScore(time, mu0, sigma0, adr, score),
               pred1 = predScore(time, mu, sigma, adr, score))
    return(pred)
}


## Functions for plot
mytheme <- list(
    theme(
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 16)
    )
)

plot_rawdata <- function(x) {
    p <- ggplot(x, aes(time, score)) +
            geom_line(color = "gray30") +
            geom_point(size = 1) +
            labs(x = "Time (min)", y = "Score") +
            facet_grid(Drug ~ ID) +
            theme_bw() +
            mytheme
    return(p)
}

plot_fitted <- function(x, pred) {
    p <- ggplot(x, aes(time, score)) +
            geom_line(color = "gray30") +
            geom_point(size = 1) +
            labs(x = "Time (min)", y = "Score") +
            geom_line(data = pred, aes(y = pred0), color = "red") +
            geom_line(data = pred, aes(y = pred1), color = "blue",
                      linetype = "twodash") +
            facet_grid(Drug ~ ID) +
            theme_bw() +
            mytheme
    return(p)
}


## Function to group for dividing figures
Grouping <- function(x) {
    cut(x, breaks = c((0:6) * 8 + 1, Inf), right = FALSE) |>
        as.integer()
}


## data for fitted curves
time <- seq(0, 100, by = 0.5)
pred1 <- get_pred(dat, time, param1, param1_individual, drug_name)
pred2 <- get_pred(dat, time, param2, param2_individual, drug_name)

datG <- dat |> mutate(group = Grouping(ID)) |> group_nest(group)
pred1G <- pred1 |> mutate(group = Grouping(ID)) |> group_nest(group)
pred2G <- pred2 |> mutate(group = Grouping(ID)) |> group_nest(group)


## Plot with fitted curves
n <- length(datG$data)
p_raw <- p1 <- p2 <- vector("list", n)

for (i in seq_len(n)) {
    p_raw[[i]] <- plot_rawdata(datG$data[[i]])
    p1[[i]] <- plot_fitted(datG$data[[i]], pred1G$data[[i]])
    p2[[i]] <- plot_fitted(datG$data[[i]], pred2G$data[[i]])
}


## Save plots (PDF)
cairo_pdf(filename = "result/fig1.pdf",
          width = 9, height = 6, onefile = TRUE)
for (i in seq_len(n)) {
    print(p_raw[[i]])
}
dev.off()

cairo_pdf(filename = "result/fig2.pdf",
          width = 9, height = 6, onefile = TRUE)
for (i in seq_len(n)) {
    print(p1[[i]])
}
dev.off()

cairo_pdf(filename = "result/fig3.pdf",
          width = 9, height = 6, onefile = TRUE)
for (i in seq_len(n)) {
    print(p2[[i]])
}
dev.off()

