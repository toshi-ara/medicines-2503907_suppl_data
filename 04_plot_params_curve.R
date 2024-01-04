library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
source("common_functions.R")

dir.create("result", showWarnings = FALSE)
dir.create("result/SVG", showWarnings = FALSE)


## Parameter: mu0, sigma0
mu <- c(20, 50, 80)
sigma <- c(5, 10, 15)

res_list <- vector("list", length(mu) * length(sigma))
k <- 1
for (i in mu)  {
    for (j in sigma) {
        time <- seq(0, 100, 0.5)
        prob <- predProb(time, i, j, 0)
        res_list[[k]] <- data.frame(mu = factor(i),
                                    sigma = factor(j), time, prob)
        k <- k + 1
    }
}
res <- bind_rows(res_list)

p1 <- ggplot(res, aes(time, prob,
                group = interaction(mu, sigma),
                color = sigma, linetype = sigma)) +
    geom_line(size = 1) +
    labs(color = expression(sigma), linetype = expression(sigma))


## Parameter: adr
mu <- 60
sigma <- 10
Adr <- seq(0, 1, 0.2)
time <- seq(0, 100, 0.5)

res_list <- vector("list", length(Adr))
for (i in seq_len(length(Adr)))  {
    prob <- predProb(time, mu, sigma, Adr[i])
    res_list[[i]] <- data.frame(adr = i, time, prob)
}
res <- dplyr::bind_rows(res_list) |>
    mutate(adr = factor(adr, labels = Adr))

p2 <- ggplot(res, aes(time, prob, group = adr,
                      color = adr, linetype = adr)) +
    geom_line(size = 1)


## Multipanel plot
p <- (p1 / p2) &
    # plot_annotation(tag_levels = "A") &
    labs(x = "Time (min)", y = "Probability") &
    scale_color_brewer(palette = "Dark2") &
    theme_bw() &
    theme(axis.text = element_text(size = 21, color = "black"),
          axis.title = element_text(size = 24),
          axis.ticks = element_line(size = 1),
          legend.position = "bottom",
          legend.box.background = element_rect(),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 21),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    )


WIDTH <- 6
HEIGHT <- 10

## Save plots (PDF)
cairo_pdf(filename = "result/fig2_params_curve.pdf",
          width = WIDTH, height = HEIGHT, onefile = TRUE)
print(p)
dev.off()

