library(xtable)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcpp)
library(profvis)
library(microbenchmark)
library(parallel)
options(mc.cores=3, width=120)

source("functions/EM.R")
source("functions/misc.R")
source("functions/optim.R")
source("functions/preprocess.R")
source("functions/simulation.R")

## Model assumptions
beta <- 2
k <- 4
lambda <- 0.0058
pH <- 0.5

## Simulation parameters
p1 <- 0.2
alpha <- 4

## relative tolerance threshold in parameters to stop EM algorithm
EM_rel.tol <- 0.0005



