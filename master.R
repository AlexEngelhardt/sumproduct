source("init.R")

## Run one of these two lines (i.e. either simulate data or load the family study):
raw_dat <- simulate_dataset(n_familys=1000, n_generations=3, p1, alpha, pH, n_lastgen=3, seed=20160921)
raw_dat <- load_family_study()

#### Either preprocess now:
dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)
## saveRDS(dat, file="data/01_dat_imputedparents.rds")
## saveRDS(dat, file="data/simfam_100_3.rds")

#### Or load a previously preprocessed data set:
dat <- readRDS("data/01_dat_imputedparents.rds")  # real data
dat <- readRDS("data/simfam_100_3.rds")           # simulated data


################################################################
#### Results: Measure speed improvement of SPA over optim

#### One Family

dat <- readRDS("data/01_dat_imputedparents.rds")  # real data

largest_three <- as.numeric(names(sort(table(dat$FamID), dec=TRUE)[1:3]))  # family IDs with more than 20 members
fam <- dat[dat$FamID==largest_three[1],]  # largest family

attr(fam, "FG") <- make_factor_graph(fam)

theta_t <- list(p1=0.2, alpha=4)  # starting values

EM_largefam  <- run_EM(fam, convergence_reltol=EM_rel.tol, log_every=10, SPA=TRUE)
opt_largefam <- run_optim(fam, parallel=FALSE)
attr(EM_largefam, "elapsed")   # 56sec @fam104 with 23persons
attr(opt_largefam, "elapsed")  # 240sec

attr(EM_largefam, "elapsed") / attr(opt_largefam, "elapsed")

#### Whole real data set

dat <- readRDS("data/01_dat_imputedparents.rds")  # real data

EM_realdata  <- run_EM(dat, convergence_reltol=EM_rel.tol, log_every=1, SPA=TRUE)
attr(EM_realdata, "elapsed")   # 6855sec @ rel.tol 5e-04 (but SPA for all (i.e. also <10members) familys)
opt_realdata <- run_optim(dat, parallel=FALSE)
attr(opt_realdata, "elapsed")  # 505sec

## Estimated aprameters:
invlogit(opt_realdata$par[1])  # 0.9006657
opt_realdata$par[2]            # 5.723495

################################################################
#### Speed comparison within EM algo
#### Marginalization vs. SPA

fams <- data.frame(FamID=unique(dat$FamID), size=NA, margin.time=NA, SPA.time=NA)

for(row in 1:nrow(fams)){
    print(paste0(row, " of ", nrow(fams)))
    fam <- dat[dat$FamID==fams[row, "FamID"], ]
    fams[row, "size"] <- nrow(fam)
    fams[row, "margin.time"] <- system.time( E(fam, theta_t=list(p1=0.2, alpha=4)) )[3]
    fams[row, "SPA.time"] <- system.time( E_SPA(fam, theta_t=list(p1=0.2, alpha=4)) )[3]
}
fams$diff <- fams$margin.time-fams$SPA.time

saveRDS(fams, file="data/Estep_runtime_SPA_vs_marginalization.rds")
fams <- readRDS(file="data/Estep_runtime_SPA_vs_marginalization.rds")

plot(fams$size, fams$diff)  # Paper Figure this, or thesis only?
abline(v=19.5)  # cutoff from marginalization to EM makes sense at 20 persons
fams[567,]  # ok so linearify this summation in here. 13 kids!

################################################################
#### Bootstrapping for parameter estimation stability

B <- 3

dat_B <- bootstrap_families(dat, B=B, seed=20160624)

saveRDS(dat_B, file=paste0("data/results/bootstrap_families_", B, ".rds"))

results <- mclapply(dat_B, function(dat_b){
    print(paste(format(Sys.time()), ":: Starting bootstrap analysis", attr(dat_b,"b"), "of", B))

    opt_res <- run_optim(dat_b)  # or NA if disregarding optim 
    em_res  <- run_EM(dat_b, 100)

    return(list(optim=opt_res, EM=em_res))
})

saveRDS(results, file=paste0("data/results/bootstrap_results_", B, ".rds"))

## Plot bootstrap results

par(mfrow=c(2,2))
results <- readRDS(paste0("data/results/bootstrap_results_", "imputed", ".rds"))
plot(t(sapply(results, function(r) r$optim$par)), main="01, optim")
plot(t(sapply(results, function(x) sapply(x$EM[1:2], tail, 1))), main="01, EM")

results <- readRDS(paste0("data/results/bootstrap_results_", "censored80", ".rds"))
plot(t(sapply(results, function(r) r$optim$par)), main="02, optim")
plot(t(sapply(results, function(x) sapply(x$EM[1:2], tail, 1))), main="02, EM")


################################################################
################################################################
################################################################
#### Generate paper Figure 4: Convergence comparison of EM and optim with and without 20% removed/imputed

p1 <- 0.2
alpha <- 4
pH <- 0.5

options(mc.cores=7)

iterations <- 100
remove_20_percent <- TRUE  # TODO loop over this \in {TRUE,FALSE}

results <- as.data.frame(matrix(NA, nrow=iterations, ncol=4, dimnames=list(NULL, c("optim_p1","optim_alpha","optim_p1","optim_alpha"))))
set.seed(20161003)

fig4_res <- mclapply(1:iterations, function(iter){

    cat(paste0(Sys.time(), " :: Iteration ", iter, " of ", iterations, "\n"))
    
    raw_dat <- simulate_dataset(n_familys=500, n_generations=3, p1, alpha, pH, n_lastgen=3)
    if(remove_20_percent){
        raw_dat <- remove_persons(raw_dat, frac=0.2)
    }
    dat <-  preprocess(raw_dat, impute=TRUE, pH=pH, compute_PZL=TRUE)

    theta_0 <- list(p1=runif(1,0.1,0.9), alpha=runif(1, 2, 20))
    
    EM_simdata  <- run_EM(dat, convergence_reltol=EM_rel.tol, theta_0=theta_0, log_every=100, SPA=TRUE)  # 16sec pro run
    opt_simdata <- run_optim(dat, parallel=FALSE, theta_0=theta_0)  # 10sec pro run

    ## results[iter, "optim_p1"] <- invlogit(opt_simdata$par[1])
    ## results[iter, "optim_alpha"] <- opt_simdata$par[2]
    ## results[iter, "EM_p1"] <- EM_simdata$p1[length(EM_simdata$p1)]
    ## results[iter, "EM_alpha"] <- EM_simdata$alpha[length(EM_simdata$alpha)]

    c(optim_p1=invlogit(opt_simdata$par[1]),
      optim_alpha=opt_simdata$par[2],
      EM_p1=EM_simdata$p1[length(EM_simdata$p1)],
      EM_alpha=EM_simdata$alpha[length(EM_simdata$alpha)])
    
})
Sys.time()
results <- Reduce(rbind, fig4_res)
rownames(results) <- NULL
results <- as.data.frame(results)

saveRDS(results, file="data/fig4_results_20percentremoved.rds")
results <- readRDS("data/fig4_results_20percentremoved.rds")

relative_difference_p1 <- (results$optim_p1 - results$EM_p1) / (results$optim_p1)
relative_difference_alpha <- (results$optim_alpha - results$EM_alpha) / (results$optim_alpha)

pdf(paste0("data/convergence_comparison", ifelse(remove_20_percent, "_imputed", ""), ".pdf"))

layout(matrix(1:4, byrow=TRUE, nrow=2))
plot(results$optim_p1, results$optim_alpha, xlim=c(0,1), ylim=c(1,10), main="(a) Optim convergence")
plot(results$EM_p1, results$EM_alpha, xlim=c(0,1), ylim=c(1,10), main="(b) EM convergence")
boxplot(relative_difference_p1, main="(c) Relative difference of p1"); abline(h=0, lty=2)
boxplot(relative_difference_alpha, main="(d) Relative difference of alpha"); abline(h=0, lty=2)
dev.off()

################################################################
################################################################
################################################################
#### Generate paper Figure 5: Runtime comparison of EM/optim by increasing N and D

Ns <- seq(10, 100, by=10)

Ds <- c(9, 15, 20, 22, 24)  # WARNING: This is manually set to be the same as the simulated families below.
## Simulate max-size families once. Then just subset for smaller N

all_raw_dat <- list()

all_raw_dat[["D9"]] <- simulate_dataset(n_familys=max(Ns), n_generations=3, p1, alpha, pH, n_lastgen=3)
all_raw_dat[["D15"]] <- simulate_dataset(n_familys=max(Ns), n_generations=4, p1, alpha, pH, n_lastgen=1)

## make 4.5 generations hack:
purge_me_dat <- simulate_dataset(n_familys=max(Ns), n_generations=5, p1, alpha, pH, n_lastgen=2)

kill_founders <- 1:12  # 1:12 = 20 persons, 1:10 = 22 persons, 1:8 = 24 persons @ n_generations=5
raw_dat <- purge_me_dat
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D20"]] <- raw_dat

kill_founders <- 1:10  # 1:12 = 20 persons, 1:10 = 22 persons, 1:8 = 24 persons @ n_generations=5
raw_dat <- purge_me_dat
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D22"]] <- raw_dat

kill_founders <- 1:8  # 1:12 = 20 persons, 1:10 = 22 persons, 1:8 = 24 persons @ n_generations=5
raw_dat <- purge_me_dat
raw_dat <- raw_dat[-which(raw_dat$pos %in% kill_founders), ]
new_founders <- raw_dat$father_pos %in% kill_founders
raw_dat[new_founders, "father_pos"] <- NA
raw_dat[new_founders, "mother_pos"] <- NA
raw_dat[new_founders, "founder"] <- TRUE
all_raw_dat[["D24"]] <- raw_dat
## end 4.5 generations hack

## The preprocessing takes around 30secs now (was: 10 minutes because of slice imputing of data.frame instead of matrix)
all_dat <- lapply(all_raw_dat, preprocess, impute=TRUE, pH=0.5, compute_PZL=TRUE)


results <- expand.grid(algo=c("optim", "EM"), D=Ds, N=Ns, t=NA)
theta_0 <- list(p1=.2, alpha=4)

for(D in Ds){
    print(paste0("* ", format(Sys.time()), ":: Starting D=", D))
    for(N in Ns){  # WARNING: Do not mclapply this because it falsifies the elapsed time from system.time
        
        print(paste0("** ", format(Sys.time()), ":: Starting N=", N))
        dat <- all_dat[[paste0("D", D)]]
        FamIDs <- unique(dat$FamID)
        dat <- subset(dat, FamID %in% FamIDs[1:N])
        attr(dat, "possible_Z_vectors_list") <- attr(all_dat[[paste0("D", D)]], "possible_Z_vectors_list")

        this_EM <- run_EM(dat, convergence_reltol=EM_rel.tol, theta_0=theta_0, log_every=100, SPA=TRUE)
        this_opt <- run_optim(dat, parallel=FALSE, theta_0=theta_0)  # 10sec pro run
        
        EM_idx  <- which(results$algo=="EM" & results$D==D & results$N==N)
        opt_idx <- which(results$algo=="optim" & results$D==D & results$N==N)
                
        results[EM_idx, "t"]  <- attr(this_EM,  "elapsed")[3]
        results[opt_idx, "t"] <- attr(this_opt, "elapsed")[3]
    }
    
}

## saveRDS(results, file="data/runtime_comparison_results.rds")
results <- readRDS(file="data/runtime_comparison_results.rds")

ggplot(results, aes(x=N, y=t, group=D, colour=factor(D))) + geom_line() + facet_grid(~algo) + scale_y_log10()

ggsave("data/runtime_comparison_results.pdf")


################################################################
################################################################
################################################################
################################################################
#### Postprocessing


################################################################
#### Compute p_riskfam via E_SPA_onefam

## compute: 1 - P(Z_I = 0 | X_I)
## only founders is enough. nonfounders come deterministically

theta_hat <- list(p1=invlogit(opt_realdata$par[1]), alpha=opt_realdata$par[2])
dat$T1hat <- E_SPA(dat, theta_t=theta_hat)  # compute T1 per sum-product algorithm

P_riskfam <- rep(NA, length(unique(dat$FamID)))
names(P_riskfam) <- as.character(unique(dat$FamID))

for(FamID in unique(dat$FamID)){
    fam <- dat[dat$FamID==FamID & dat$founder,]
    P_riskfam[as.character(FamID)] <- 1-prod(1-fam$T1hat)
}

hist(P_riskfam)
quantile(P_riskfam, (1:100)/100)

## compare estimated probabilities against true status (obviously only works for simulated data where true status is known)
true_riskfam <- table(dat$FamID, dat$risk)[,2] != 0
boxplot(P_riskfam ~ true_riskfam, main="P(riskfam) vs. true riskfam")
## ROC curve:
library(ROCR)
pred <- prediction(P_riskfam, labels=true_riskfam)
perf <- performance(pred, measure="tpr", x.measure="fpr")
plot(perf, main=paste("ROC curve for P(riskfam)\nAUC:", round(performance(pred, measure="auc")@y.values[[1]], 2))); abline(0,1,lty=2)

## you could also use Z-scores instead:
## Z_riskfam <- scale(P_riskfam)
## names(Z_riskfam) <- as.character(unique(dat$FamID))

################################################################
#### Plot imputed data

library(kinship2)
dat$ID <- paste0(dat$FamID, "_", dat$position)

dat$dad_id <- NA
dat$mom_id <- NA
for(row in 1:nrow(dat)){
    if(!is.na(dat[row, "father_pos"])){
        dad_row <- which(dat$FamID==dat[row,"FamID"] & dat$position==dat[row,"father_pos"])
        dat[row, "dad_id"] <- dat[dad_row, "ID"]
    }
    if(!is.na(dat[row, "mother_pos"])){
        mom_row <- which(dat$FamID==dat[row,"FamID"] & dat$position==dat[row,"mother_pos"])
        dat[row, "mom_id"] <- dat[mom_row, "ID"]
    }
}

single_famIDs <- as.numeric(names(table(dat$FamID))[table(dat$FamID)==1])
dat <- dat[!(dat$FamID %in% single_famIDs), ]
pedAll <- pedigree(id=dat$ID,
                   dadid=dat$dad_id, momid=dat$mom_id,
                   sex=ifelse(dat$m==1,1,2), famid=dat$FamID)

## throws errors though if one or two parents are missing
for(fam in unique(dat$FamID)){
    ## a <- try( plot(pedAll[as.character(fam)]) )
    ## if(inherits(a, "try-error")) next
    png(paste0("plots/", fam, "_P_", round(P_riskfam[as.character(fam)], 4), ".png"))
    plot(pedAll[as.character(fam)],
         id=dat[dat$FamID==fam, "t"],
         affected=dat[dat$FamID==fam, "c"],
         status=(dat[dat$FamID==fam, "t"]==0))
    dev.off()
}
