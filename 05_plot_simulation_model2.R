library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
source("common_functions.R")


########################################
# Set parameters of each individual
########################################

drug_name <- c("Pro", "Lid", "Mep", "Bup", "Lid+Adr")
drug_param <- res <- p <- vector("list", 2)

## parameters
params <- list(
    ## estimated parameters
    data.frame(
        Drug = 1:4,
        mu0 = c(67.0, 60.9, 50.0, 29.1),
        s_mu0 = c(15.7, 5.1, 12.0, 20.3),
        log_sigma0 = c(2.20, 2.42, 2.41, 2.50),
        log_s_sigma0 = c(0.78, 0.51, 0.65, 0.83),
        adr = 0.665
    ),
    ## adjusted parameters
    data.frame(
        Drug = 1:4,
        mu0 = c(75, 67, 43, 30),
        s_mu0 = c(8, 5, 6, 10),
        log_sigma0 = c(2.2, 2.4, 2.4, 2.5),
        log_s_sigma0 = c(0.4, 0.4, 0.4, 0.5),
        adr = 0.7
    )
)

d_sigma <- 4
Time <- seq(0, 120, by = 5)
N <-100

## Offset value for each individual
set.seed(1234)
d <- rnorm(n = N, mean = 0, sd = d_sigma)

## drug_param
for (i in 1:2) {
    set.seed(123)
    ## generate parameters by random generator
    param <- by(params[[i]], params[[i]]$Drug, identity) |>
        lapply(function(x)
               data.frame(Drug = x$Drug,
                          ID = seq_len(N),
                          Mean = rnorm(N, x$mu0, x$s_mu0) + d,
                          Sigma = rlnorm(N, x$log_sigma0, x$log_s_sigma0),
                          adr = 0))
    ## Set parameter of Lid + Adr
    param_Lid <- param[[2]] |>
        mutate(Drug = 5, adr = params[[i]]$adr[1])
    ## Combine
    drug_param[[i]] <- param |>
        bind_rows() |>
        bind_rows(param_Lid) |>
        mutate(Drug = factor(Drug, labels = drug_name))
}


########################################
# Start simulation
########################################

## x: parameter of drugs (data.frame of Drug, ID, Mean, Sigma, adr)
## time: vector of time point
## n_stim: number of stimulation by needle
##
## return: data.frame(ID, Drug, time, score)
doSimulation <- function(x, time, n_stim) {
    n <- nrow(x)
    res_list <- vector("list", n)

    for (i in seq_len(n)) {
        Param <- x[i,]
        prob <- predProb(Time, Param$Mean, Param$Sigma, Param$adr)
        res_list[[i]] <- data.frame(
            ID = Param$ID,
            Drug = Param$Drug,
            time = time,
            score = rbinom(n = length(time), size = n_stim, prob = prob)
        )
    }
    res <- bind_rows(res_list) |>
        mutate(ID = factor(ID),
               Drug = factor(Drug))
    return(res)
}


mytheme <- list(
    labs(x = "Time (min)", y = "Score"),
    scale_x_continuous(breaks = seq(0, 120, by = 30)),
    theme_bw(),
    theme(
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 16)
    )
)

plotSim <- function(x) {
    p <- ggplot(x, aes(time, score)) +
        geom_line(color = "gray30") +
        geom_point(size = 1) +
        facet_grid(Drug ~ ID) +
        mytheme
    return(p)
}


# Start simulation
p <- pAll <- vector("list", 2)

for (i in 1:2) {
    set.seed(1234)
    res[[i]] <- doSimulation(drug_param[[i]], Time, n_stim = 6)
    p[[i]] <- plotSim(dplyr::filter(res[[i]], ID %in% 1:8))

    for (j in 0:floor(N/8)) {
        pAll[[i]][[j+1]] <- plotSim(dplyr::filter(res[[i]],
                                    ID %in% (8*j+1):(8*(j+1))))
    }
}

save(res, file = "rds/data_simulation_model2.rds")


## Save plots (PDF)
cairo_pdf(filename = "SFig4.pdf",
          width = 9, height = 6, onefile = TRUE)
print(pAll[[1]])
dev.off()

cairo_pdf(filename = "SFig5.pdf",
          width = 9, height = 6, onefile = TRUE)
print(pAll[[2]])
dev.off()

