source("init.R")

set.seed(20161218)

################################################################
#### 3.4.8 Speed-up of grid-purging and memoization when evaluating log-likelihood on real data set

dat <- readRDS("data/01_dat_imputedparents.rds")

t1 <- system.time({ lik1 <- l(dat, p1=0.2, alpha=4, parallel=FALSE, gridpurging=FALSE, memoize=FALSE) })
t2 <- system.time({ lik2 <- l(dat, p1=0.2, alpha=4, parallel=FALSE, gridpurging=TRUE, memoize=TRUE) })

t2 / t1  # ca. 0.092

################################################################
#### 3.5.1 Simulation Study

#### Simulation setup (parameters)

p1 <- 0.2
alpha <- 4
pH <- 0.5

#### Simulate all datasets

Ns <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250)  # new
Ds <- c(15, 17, 19, 21, 23, 25)  # WARNING: This is manually set to be the same as the simulated families below.

## Simulate max-size families once. Then just subset for smaller N

all_raw_dat <- list()
set.seed(20161122)
all_raw_dat[["D5"]] <- simulate_dataset(n_familys=max(Ns), n_generations=2, p1, alpha, pH, n_lastgen=3)
all_raw_dat[["D9"]] <- simulate_dataset(n_familys=max(Ns), n_generations=3, p1, alpha, pH, n_lastgen=3)
all_raw_dat[["D15"]] <- simulate_dataset(n_familys=max(Ns), n_generations=4, p1, alpha, pH, n_lastgen=1)

raw_dat_5gen <- simulate_dataset(n_familys=max(Ns), n_generations=5, p1, alpha, pH, n_lastgen=1)

kill_founders <- 1:14  # simulate 17 persons by killing 14 founder off a 5 generation (31 persons) pedigree
raw_dat <- raw_dat_5gen
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D17"]] <- raw_dat
    
kill_founders <- 1:12
raw_dat <- raw_dat_5gen
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D19"]] <- raw_dat

kill_founders <- 1:10
raw_dat <- raw_dat_5gen
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D21"]] <- raw_dat

kill_founders <- 1:8
raw_dat <- raw_dat_5gen
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D23"]] <- raw_dat

kill_founders <- 1:6
raw_dat <- raw_dat_5gen
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D25"]] <- raw_dat

all_dat <- lapply(all_raw_dat, preprocess, impute=TRUE, pH=0.5, compute_PZL=TRUE)

#### Runtime improvement

results <- expand.grid(algo=c("optim", "EM"), D=Ds, N=Ns, t=NA)
theta_0 <- list(p1=.2, alpha=4)

for(D in Ds){
    print(paste0("* ", format(Sys.time()), ":: Starting D=", D))
    cat(paste0("* ", format(Sys.time()), ":: Starting D=", D, "\n"), file="log.txt", append=TRUE)
    for(N in Ns){  # WARNING: Do not mclapply this because it falsifies the elapsed time from system.time
        
        print(paste0("** ", format(Sys.time()), ":: Starting N=", N))
        cat(paste0("** ", format(Sys.time()), ":: Starting N=", N, "\n"), file="log.txt", append=TRUE)
        dat <- all_dat[[paste0("D", D)]]
        FamIDs <- unique(dat$FamID)
        dat <- subset(dat, FamID %in% FamIDs[1:N])
        attr(dat, "possible_Z_vectors_list") <- attr(all_dat[[paste0("D", D)]], "possible_Z_vectors_list")

        this_EM <- run_EM(dat, convergence_reltol=EM_rel.tol, theta_0=theta_0, log_every=100, SPA=TRUE)
        this_opt <- run_optim(dat, parallel=FALSE, theta_0=theta_0, gridpurging=TRUE)  # In diss, with grid purging!
        
        EM_idx  <- which(results$algo=="EM" & results$D==D & results$N==N)
        opt_idx <- which(results$algo=="optim" & results$D==D & results$N==N)
                
        results[EM_idx, "t"]  <- attr(this_EM,  "elapsed")[3]
        results[opt_idx, "t"] <- attr(this_opt, "elapsed")[3]
    }
    saveRDS(results, file=paste0("data/diss/runtime_comparison_results_D", D, ".rds"))
}

saveRDS(results, file="data/diss/runtime_comparison_results.rds")
results <- readRDS(file="data/diss/runtime_comparison_results.rds")
results$D <- as.factor(results$D)

## <add D=31> (Paper Revision)
D <- 31

results <- rbind(results, expand.grid(algo=c("optim","EM"), D=as.character(31), N=Ns, t=NA))
dat_5gen <- preprocess(raw_dat_5gen, impute=TRUE, pH=0.5, compute_PZL=FALSE)

print(paste0("* ", format(Sys.time()), ":: Starting D=", D))
cat(paste0("* ", format(Sys.time()), ":: Starting D=", D, "\n"), file="log.txt", append=TRUE)
for(N in Ns){  # WARNING: Do not mclapply this because it falsifies the elapsed time from system.time
    
    print(paste0("** ", format(Sys.time()), ":: Starting N=", N))
    cat(paste0("** ", format(Sys.time()), ":: Starting N=", N, "\n"), file="log.txt", append=TRUE)

    FamIDs <- unique(dat_5gen$FamID)
    dat <- subset(dat_5gen, FamID %in% FamIDs[1:N])
    
    this_EM <- run_EM(dat, convergence_reltol=EM_rel.tol, theta_0=theta_0, log_every=100, SPA=TRUE)
    EM_idx  <- which(results$algo=="EM" & results$D==D & results$N==N)
    results[EM_idx, "t"]  <- attr(this_EM,  "elapsed")[3]
}

saveRDS(results, file=paste0("data/diss/runtime_comparison_results_D", D, ".rds"))


## Table

xtab <- spread(results, key=algo, value=t) %>% transmute(D=D, N=N, fact=optim/EM) %>% spread(key=D, value=fact)
D <- ncol(xtab)-1
print(xtable(xtab, digits=c(0, 0, rep(2,D))), include.rownames=FALSE, booktabs=TRUE)

## Figure

legendtitle <- guide_legend("Family size (D)")

results$algo <- factor(ifelse(results$algo=="optim", "Nelder-Mead", "EM"), levels=c("Nelder-Mead", "EM"))

ggplot(results, aes(x=N, y=t, group=D, colour=D)) + geom_line() + geom_point() + facet_grid(~algo) + scale_y_log10() + labs(x="Number of families (N)", y="Runtime (s)") + guides(colour=legendtitle)

ggsave("data/diss/runtime_comparison_results.pdf", width=8, height=4)


#### Equivalence of EM and Nelder-Mead

options(mc.cores=10)

iterations <- 100
remove_20_percent <- FALSE  # you could loop over this \in {TRUE,FALSE}

EM_rel.tol.old <- EM_rel.tol
EM_rel.tol <- 5e-05  # For the Bland-Altman plot, be more strict with EM convergence

results <- as.data.frame(matrix(NA, nrow=iterations, ncol=4, dimnames=list(NULL, c("optim_p1","optim_alpha","optim_p1","optim_alpha"))))
set.seed(20161126)

fig4_res <- mclapply(1:iterations, function(iter){

    cat(paste0(Sys.time(), " :: Iteration ", iter, " of ", iterations, "\n"))
    
    raw_dat <- simulate_dataset(n_familys=500, n_generations=3, p1, alpha, pH, n_lastgen=3)
    if(remove_20_percent){
        raw_dat <- remove_persons(raw_dat, frac=0.2)
    }
    dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)

    theta_0 <- list(p1=runif(1,0.1,0.9), alpha=runif(1, 2, 20))
    
    EM_simdata  <- run_EM(dat, convergence_reltol=EM_rel.tol, theta_0=theta_0, log_every=100, SPA=FALSE)
    opt_simdata <- run_optim(dat, parallel=FALSE, theta_0=theta_0)

    c(optim_p1=invlogit(opt_simdata$par[1]),
      optim_alpha=opt_simdata$par[2],
      EM_p1=EM_simdata$p1[length(EM_simdata$p1)],
      EM_alpha=EM_simdata$alpha[length(EM_simdata$alpha)])
    
})

results <- Reduce(rbind, fig4_res)
rownames(results) <- NULL
results <- as.data.frame(results)

saveRDS(results, file="data/diss/fig4_results.rds")
results <- readRDS("data/diss/fig4_results.rds")

## Convergence Figure (Bland-Altman)

pdf(paste0("data/diss/convergence_comparison", ifelse(remove_20_percent, "_imputed", ""), ".pdf"))

layout(matrix(1:4, byrow=TRUE, nrow=2))

plot(results$optim_p1, results$EM_p1, main="(a) p1 convergence", xlab="Nelder-Mead", ylab="EM Algorithm", bty="n")
abline(h=0.2); abline(v=0.2)

plot(results$optim_alpha, results$EM_alpha, main="(b) alpha convergence", xlab="Nelder-Mead", ylab="EM Algorithm", bty="n")
abline(h=4); abline(v=4)

sdr <- 2* sd(results$optim_p1 - results$EM_p1)
plot((results$optim_p1 + results$EM_p1)/2, results$optim_p1 - results$EM_p1, main="(c) Bland-Altman plot of p1",
     xlab="Average", ylab="Difference", ylim=c(-1.25*sdr, 1.25*sdr))
abline(h=0)
abline(h=sdr , lty=2)
abline(h=-sdr, lty=2)

sdr <- 2* sd(results$optim_alpha - results$EM_alpha)
plot((results$optim_alpha + results$EM_alpha)/2, results$optim_alpha - results$EM_alpha, main="(d) Bland-Altman plot of alpha",
     xlab="Average", ylab="Difference", ylim=c(-1.25*sdr, 1.25*sdr))
abline(h=0)
abline(h=sdr , lty=2)
abline(h=-sdr, lty=2)

dev.off()

## Table

print(xtable(t(summary(results))), include.rownames=FALSE, booktabs=TRUE)
## Manually prettify this table.


## Convergence Figure (Bland-Altman) after removing and imputing 20%

remove_20_percent <- TRUE 

results <- as.data.frame(matrix(NA, nrow=iterations, ncol=4, dimnames=list(NULL, c("optim_p1","optim_alpha","optim_p1","optim_alpha"))))
set.seed(20161126)

fig4_res <- mclapply(1:iterations, function(iter){

    cat(paste0(Sys.time(), " :: Iteration ", iter, " of ", iterations, "\n"))
    
    raw_dat <- simulate_dataset(n_familys=500, n_generations=3, p1, alpha, pH, n_lastgen=3)
    if(remove_20_percent){
        raw_dat <- remove_persons(raw_dat, frac=0.2)
    }
    dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)

    theta_0 <- list(p1=runif(1,0.1,0.9), alpha=runif(1, 2, 20))
    
    EM_simdata  <- run_EM(dat, convergence_reltol=EM_rel.tol, theta_0=theta_0, log_every=100, SPA=FALSE)  # 16sec pro run
    opt_simdata <- run_optim(dat, parallel=FALSE, theta_0=theta_0)  # 10sec pro run

    c(optim_p1=invlogit(opt_simdata$par[1]),
      optim_alpha=opt_simdata$par[2],
      EM_p1=EM_simdata$p1[length(EM_simdata$p1)],
      EM_alpha=EM_simdata$alpha[length(EM_simdata$alpha)])
    
})

results <- Reduce(rbind, fig4_res)
rownames(results) <- NULL
results <- as.data.frame(results)

saveRDS(results, file="data/diss/fig4_results_imputed.rds")
results <- readRDS("data/diss/fig4_results_imputed.rds")

pdf(paste0("data/diss/convergence_comparison", ifelse(remove_20_percent, "_imputed", ""), ".pdf"))

layout(matrix(1:4, byrow=TRUE, nrow=2))

plot(results$optim_p1, results$EM_p1, main="(a) p1 convergence", xlab="Nelder-Mead", ylab="EM Algorithm", bty="n")
abline(h=0.2); abline(v=0.2)

plot(results$optim_alpha, results$EM_alpha, main="(b) alpha convergence", xlab="Nelder-Mead", ylab="EM Algorithm", bty="n")
abline(h=4); abline(v=4)

sdr <- 2* sd(results$optim_p1 - results$EM_p1)
plot((results$optim_p1 + results$EM_p1)/2, results$optim_p1 - results$EM_p1, main="(c) Bland-Altman plot of p1",
     xlab="Average", ylab="Difference", ylim=c(-1.25*sdr, 1.25*sdr))
abline(h=0)
abline(h=sdr , lty=2)
abline(h=-sdr, lty=2)

sdr <- 2* sd(results$optim_alpha - results$EM_alpha)
plot((results$optim_alpha + results$EM_alpha)/2, results$optim_alpha - results$EM_alpha, main="(d) Bland-Altman plot of alpha",
     xlab="Average", ylab="Difference", ylim=c(-1.25*sdr, 1.25*sdr))
abline(h=0)
abline(h=sdr , lty=2)
abline(h=-sdr, lty=2)

dev.off()

#### Estimation of P(riskfam) and ROC-Curve with AUC

dat <- preprocess(simulate_dataset(n_familys=1000, n_generations=3, p1=0.2, alpha=4, pH=0.5, n_lastgen=1, seed=20160921), compute_PZL=TRUE)

## compute: 1 - P(Z_I = 0 | X_I)
## only founders is enough. nonfounders come deterministically

opt <- run_optim(dat, parallel=FALSE)  # Same result than EM, but faster for family of size 9
saveRDS(opt, file="data/diss/opt_simdata.rds")
opt <- readRDS("data/diss/opt_simdata.rds")

theta_hat <- list(p1=invlogit(opt$par[1]), alpha=opt$par[2])
dat$T1hat <- E_SPA(dat, theta_t=theta_hat)  # compute T1 per sum-product algorithm

P_riskfam <- rep(NA, length(unique(dat$FamID)))
names(P_riskfam) <- as.character(unique(dat$FamID))

for(FamID in unique(dat$FamID)){
    fam <- dat[dat$FamID==FamID & dat$founder,]
    P_riskfam[as.character(FamID)] <- 1-prod(1-fam$T1hat)
}

## compare estimated probabilities against true status (obviously only works for simulated data where true status is known)
true_riskfam <- table(dat$FamID, dat$risk)[,2] != 0

## ROC curve:
library(ROCR)
pred <- prediction(P_riskfam, labels=true_riskfam)
perf <- performance(pred, measure="tpr", x.measure="fpr")
pdf("data/diss/ROC.pdf")
plot(perf, main=paste("ROC curve for P(riskfam)\nAUC:", round(performance(pred, measure="auc")@y.values[[1]], 2))); abline(0,1,lty=2)
dev.off()

################################################################
#### 3.5.2 Real data

EM_rel.tol <- EM_rel.tol.old  # switch back to old tolerance after Bland-Altman plots

#### Preprocessing

raw_dat <- load_family_study()
dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)

#### Parameter estimates
#### Runtime comparison

EM_realdata_1  <- run_EM(dat, convergence_reltol=EM_rel.tol, log_every=50, SPA=TRUE, SPA_cutoff=1)
EM_realdata_18  <- run_EM(dat, convergence_reltol=EM_rel.tol, log_every=50, SPA=TRUE, SPA_cutoff=18)

EM_realdata_20  <- run_EM(dat, convergence_reltol=EM_rel.tol, log_every=50, SPA=TRUE, SPA_cutoff=20)
opt_realdata <- run_optim(dat, parallel=FALSE)


EM_realdata_16  <- run_EM(dat, convergence_reltol=EM_rel.tol, log_every=50, SPA=TRUE, SPA_cutoff=16)
EM_realdata_21  <- run_EM(dat, convergence_reltol=EM_rel.tol, log_every=50, SPA=TRUE, SPA_cutoff=21)
save(EM_realdata_1, EM_realdata_16, EM_realdata_18, EM_realdata_20, EM_realdata_21, opt_realdata, file="data/diss/runtime_realdata.RData")

attr(EM_realdata_1, "elapsed")
attr(EM_realdata_16, "elapsed")
attr(EM_realdata_18, "elapsed")
attr(EM_realdata_20, "elapsed")
attr(EM_realdata_21, "elapsed")
attr(opt_realdata, "elapsed")

#### EM is slower. Try on largest family only:

largest_family <- as.numeric(names(sort(table(dat$FamID), decreasing=TRUE)[1]))
fam <- dat[dat$FamID==largest_family,]

opt_fam <- run_optim(fam, parallel=FALSE)
EM_fam  <- run_EM(fam, convergence_reltol=EM_rel.tol, log_every=50, SPA=TRUE)
attr(EM_fam, "elapsed") # 41.99sec
attr(opt_fam, "elapsed") # 172.347sec

################################################################
#### 3.5.3 Multiple starts of the EM algorithm are not necessary

#### Likelihood surfaces

p1s <- seq(0.05, 0.95, by=0.05)
alphas <- seq(2, 8, by=0.5)

likelihood_surface <- function(dat, p1s, alphas){
    surface <- expand.grid(p1=p1s, alpha=alphas)
    surface$l <- mclapply(1:nrow(surface), function(row){
        l(dat, p1=surface[row, "p1"], alpha=surface[row, "alpha"])
    }) %>% unlist()
    z <- matrix(surface$l, nrow=length(p1s), ncol=length(alphas))
    return(z)
}

contours <- list()

raw_dat <- simulate_dataset(n_familys=100, n_generations=3, p1, alpha, pH, n_lastgen=3, seed=20160921)
dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)
contours$z1 <- likelihood_surface(dat, p1s, alphas)

raw_dat <- simulate_dataset(n_familys=1000, n_generations=3, p1, alpha, pH, n_lastgen=3, seed=20160921)
dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)
contours$z2 <- likelihood_surface(dat, p1s, alphas)

raw_dat <- simulate_dataset(n_familys=100, n_generations=4, p1, alpha, pH, n_lastgen=3, seed=20160921)
dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)
contours$z3 <- likelihood_surface(dat, p1s, alphas)

dat <- readRDS("data/01_dat_imputedparents.rds")
contours$z4 <- likelihood_surface(dat, p1s, alphas)

saveRDS(contours, file="data/diss/contours.rds")

contours <- readRDS("data/diss/contours.rds")

pdf("data/diss/contours.pdf", width=10, height=5)
layout(matrix(1:4, byrow=T, nrow=2))
par(oma=c(0,0,0,0), mar=c(4,4,2,1))
idx <- c(1,2,3, 2^(2:floor(log2(length(pretty(range(contours$z1), 512))))))  # hack to show nicely spaced contour lines
contour(x=p1s, y=alphas, z=contours$z1, levels=rev(pretty(range(contours$z1), 512))[idx], main="(a) 100 families of 3 generations", xlab="p1", ylab="alpha")
idx <- c(1,2,3, 2^(2:floor(log2(length(pretty(range(contours$z2), 512))))))
contour(x=p1s, y=alphas, z=contours$z2, levels=rev(pretty(range(contours$z2), 1024))[idx], main="(b) 1000 families of 3 generations", xlab="p1", ylab="alpha")
idx <- c(1,2,3, 2^(2:floor(log2(length(pretty(range(contours$z3), 512))))))
contour(x=p1s, y=alphas, z=contours$z3, levels=rev(pretty(range(contours$z3), 512))[idx], main="(a) 100 families of 4 generations", xlab="p1", ylab="alpha")
idx <- c(1,2,3, 2^(2:floor(log2(length(pretty(range(contours$z4), 512))))))
contour(x=p1s, y=alphas, z=contours$z4, levels=rev(pretty(range(contours$z4), 512))[idx], main="(b) real data", xlab="p1", ylab="alpha")
dev.off()


################################################################
#### Bootstrapping for parameter estimation stability

## use real data!

B <- 100

dat_B <- bootstrap_families(dat, B=B, seed=20160624)

## saveRDS(dat_B, file=paste0("data/bootstrap_families_", B, ".rds"))  # 1GB and takes long

results <- mclapply(dat_B, function(dat_b){
    print(paste(format(Sys.time()), ":: Starting bootstrap analysis", attr(dat_b,"b"), "of", B))

    ## opt_res <- run_optim(dat_b)  # or NA if disregarding optim
    opt_res <- NA
    em_res  <- run_EM(dat_b, 100)

    return(list(optim=opt_res, EM=em_res))
})

saveRDS(results, file=paste0("data/diss/bootstrap_results_", B, ".rds"))

res <- t(sapply(results, function(x) sapply(x$EM[1:2], tail, 1)))
plotrix::std.error(res[,"p1"])
plotrix::std.error(res[,"alpha"])

## Plot bootstrap results

## par(mfrow=c(1,2))
results <- readRDS(paste0("data/diss/bootstrap_results_", B, ".rds"))
## plot(t(sapply(results, function(r) r$optim$par)), main="optim")
plot(t(sapply(results, function(x) sapply(x$EM[1:2], tail, 1))), main="EM")

