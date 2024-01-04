## Return predict probability
predProb <- function(time, mu, sigma, adr) {
    conc <- 100 - (1 - adr) * time
    pnorm(conc, mu, sigma, lower.tail = FALSE)
}

## Return predict score
predScore <- function(time, mu, sigma, adr, score = 6) {
    predProb(time, mu, sigma, adr) * score
}

## add $ mark
addMathDollor <- function(x) {
    ifelse(is.na(x), NA, paste0("$", x, "$"))
}

