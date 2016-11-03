################################################################
#### Load and prepare family study (i.e. real data example)

load_family_study <- function(filename="data/Datenpaket04neu.RData"){
    load(filename)  # mdat
    mdat <- mdat[-which(is.na(mdat$Mann)),]
    mdat <- mdat[-which(is.na(mdat$generation)),]

    ## Split up the largest family, 555, manually. They are three independent families, see a plot of the pedigree.
    maxID <- max(mdat$famID)
    fam_555_2 <- mdat$famID==555 & (mdat$position %in% c(18:22))
    fam_555_3 <- mdat$famID==555 & (mdat$position %in% c(14:17))
    mdat[fam_555_2, "famID"] <- maxID+1
    mdat[fam_555_3, "famID"] <- maxID+2

    dat <- transmute(mdat,
                     FamID=famID, generation=generation, position=position,
                     father_pos=positionVater, mother_pos=positionMutter,
                     founder=(is.na(mdat$positionVater) & is.na(mdat$positionMutter)),
                     m=Mann, c=CRC, t=AlterEvent
                     )
    dat <- dat[order(dat$FamID, dat$generation, dat$position),]
}


################################################################
################################################################
#### Sum-Product-Algorithm test with fixed family structure (9 people)
#### The E-step for families9.

## Testing the linear time summation of many children

E_SPA_onefam9 <- function(fam, theta_t){
    p1_t    <- theta_t$p1
    alpha_t <- theta_t$alpha

    f1 <- function(z1) f(fam$t[1], r=z1, fam$c[1], fam$m[1], alpha_t) * p1_t^z1 * (1-p1_t)^(1-z1)
    f2 <- function(z2) f(fam$t[2], r=z2, fam$c[2], fam$m[2], alpha_t) * p1_t^z2 * (1-p1_t)^(1-z2)
    f3 <- function(z3) f(fam$t[3], r=z3, fam$c[3], fam$m[3], alpha_t) * p1_t^z3 * (1-p1_t)^(1-z3)
    f4 <- function(z4) f(fam$t[4], r=z4, fam$c[4], fam$m[4], alpha_t) * p1_t^z4 * (1-p1_t)^(1-z4)

    f5 <- function(z1, z2, z5) f(fam$t[5], r=z5, fam$c[5], fam$m[5], alpha_t) * (pH*z1 + pH*z2 - pH*pH*z1*z2)^z5 * (1-(pH*z1 + pH*z2 - pH*pH*z1*z2))^(1-z5)
    f6 <- function(z3, z4, z6) f(fam$t[6], r=z6, fam$c[6], fam$m[6], alpha_t) * (pH*z3 + pH*z4 - pH*pH*z3*z4)^z6 * (1-(pH*z3 + pH*z4 - pH*pH*z3*z4))^(1-z6)

    f7 <- function(z5, z6, z7) f(fam$t[7], r=z7, fam$c[7], fam$m[7], alpha_t) * (pH*z5 + pH*z6 - pH*pH*z5*z6)^z7 * (1-(pH*z5 + pH*z6 - pH*pH*z5*z6))^(1-z7)
    f8 <- function(z5, z6, z8) f(fam$t[8], r=z8, fam$c[8], fam$m[8], alpha_t) * (pH*z5 + pH*z6 - pH*pH*z5*z6)^z8 * (1-(pH*z5 + pH*z6 - pH*pH*z5*z6))^(1-z8)
    f9 <- function(z5, z6, z9) f(fam$t[9], r=z9, fam$c[9], fam$m[9], alpha_t) * (pH*z5 + pH*z6 - pH*pH*z5*z6)^z9 * (1-(pH*z5 + pH*z6 - pH*pH*z5*z6))^(1-z9)

    f789 <- function(z5, z6, z7, z8, z9){
        f7(z5, z6, z7) * f8(z5, z6, z8) * f9(z5, z6, z9)
    }

    mu_f1_z1 <- f1
    mu_f2_z2 <- f2
    mu_f3_z3 <- f3
    mu_f4_z4 <- f4

    mu_z1_f5 <- mu_f1_z1
    mu_z2_f5 <- mu_f2_z2
    mu_z3_f6 <- mu_f3_z3
    mu_z4_f6 <- mu_f4_z4

    mu_f5_z5 <- function(z5){
        sapply(0:1, function(z1){
            sapply(0:1, function(z2){
                f5(z1, z2, z5) * mu_z1_f5(z1) * mu_z2_f5(z2)
            }) %>% sum()
        }) %>% sum()
    }

    mu_f6_z6 <- function(z6){
        sapply(0:1, function(z3){
            sapply(0:1, function(z4){
                f6(z3, z4, z6) * mu_z3_f6(z3) * mu_z4_f6(z4)
            }) %>% sum()
        }) %>% sum()
    }

    mu_z5_f789 <- mu_f5_z5
    mu_z6_f789 <- mu_f6_z6

    mu_z7_f789 <- function(z7) 1
    mu_z8_f789 <- function(z8) 1
    mu_z9_f789 <- function(z9) 1

    mu_f789_z5 <- function(z5) {  # FAST
        sapply(0:1, function(z6){
            mu6 <- mu_z6_f789(z6)
            f7mu7 <- f7(z5, z6, z7=0) * mu_z7_f789(z7=0) + f7(z5, z6, z7=1) * mu_z7_f789(z7=1)
            f8mu8 <- f8(z5, z6, z8=0) * mu_z8_f789(z8=0) + f8(z5, z6, z8=1) * mu_z8_f789(z8=1)
            f9mu9 <- f9(z5, z6, z9=0) * mu_z9_f789(z9=0) + f9(z5, z6, z9=1) * mu_z9_f789(z9=1)

            mu6 * f7mu7 * f8mu8 * f9mu9
        }) %>% sum()
    }

    mu_f789_z6 <- function(z6) {  # version 1 with intermediate variables
        sapply(0:1, function(z5){
            mu5 <- mu_z5_f789(z5)
            f7mu7 <- f7(z5, z6, z7=0) * mu_z7_f789(z7=0) + f7(z5, z6, z7=1) * mu_z7_f789(z7=1)
            f8mu8 <- f8(z5, z6, z8=0) * mu_z8_f789(z8=0) + f8(z5, z6, z8=1) * mu_z8_f789(z8=1)
            f9mu9 <- f9(z5, z6, z9=0) * mu_z9_f789(z9=0) + f9(z5, z6, z9=1) * mu_z9_f789(z9=1)

            mu5 * f7mu7 * f8mu8 * f9mu9
        }) %>% sum()
    }

    mu_f789_z7 <- function(z7) {  # version 2 without intermediate variables
        sapply(0:1, function(z5){
            sapply(0:1, function(z6){
                mu_z5_f789(z5) * mu_z6_f789(z6) * f7(z5, z6, z7) *
                    (f8(z5, z6, z8=0) * mu_z8_f789(z8=0) + f8(z5, z6, z8=1) * mu_z8_f789(z8=1)) *
                        (f9(z5, z6, z9=0) * mu_z9_f789(z9=0) + f9(z5, z6, z9=1) * mu_z9_f789(z9=1))
            }) %>% sum()
        }) %>% sum()
    }

    mu_f789_z8 <- function(z8) {
        sapply(0:1, function(z5){
            sapply(0:1, function(z6){
                mu_z5_f789(z5) * mu_z6_f789(z6) * f8(z5, z6, z8) *
                    (f7(z5, z6, z7=0) * mu_z7_f789(z7=0) + f7(z5, z6, z7=1) * mu_z7_f789(z7=1)) *
                        (f9(z5, z6, z9=0) * mu_z9_f789(z9=0) + f9(z5, z6, z9=1) * mu_z9_f789(z9=1))
            }) %>% sum()
        }) %>% sum()
    }

    mu_f789_z9 <- function(z9) {
        sapply(0:1, function(z5){
            sapply(0:1, function(z6){
                mu_z5_f789(z5) * mu_z6_f789(z6) * f9(z5, z6, z9) *
                    (f7(z5, z6, z7=0) * mu_z7_f789(z7=0) + f7(z5, z6, z7=1) * mu_z7_f789(z7=1)) *
                        (f8(z5, z6, z8=0) * mu_z8_f789(z8=0) + f8(z5, z6, z8=1) * mu_z8_f789(z8=1))
            }) %>% sum()
        }) %>% sum()
    }

    mu_z5_f5 <- mu_f789_z5
    mu_z6_f6 <- mu_f789_z6

    mu_f5_z1 <- function(z1){
        sapply(0:1, function(z2){
            sapply(0:1, function(z5){
                f5(z1, z2, z5) * mu_z2_f5(z2) * mu_z5_f5(z5)
            }) %>% sum()
        }) %>% sum()
    }
    mu_f5_z2 <- function(z2){
        sapply(0:1, function(z1){
            sapply(0:1, function(z5){
                f5(z1, z2, z5) * mu_z1_f5(z1) * mu_z5_f5(z5)
            }) %>% sum()
        }) %>% sum()
    }
    mu_f6_z3 <- function(z3){
        sapply(0:1, function(z4){
            sapply(0:1, function(z6){
                f6(z3, z4, z6) * mu_z4_f6(z4) * mu_z6_f6(z6)
            }) %>% sum()
        }) %>% sum()
    }
    mu_f6_z4 <- function(z4){
        sapply(0:1, function(z3){
            sapply(0:1, function(z6){
                f6(z3, z4, z6) * mu_z3_f6(z3) * mu_z6_f6(z6)
            }) %>% sum()
        }) %>% sum()
    }

    T1_SPA <- rep(NA, 9)

    f_z1 <- function(z1) mu_f1_z1(z1) * mu_f5_z1(z1)
    f_z2 <- function(z2) mu_f2_z2(z2) * mu_f5_z2(z2)
    f_z3 <- function(z3) mu_f3_z3(z3) * mu_f6_z3(z3)
    f_z4 <- function(z4) mu_f4_z4(z4) * mu_f6_z4(z4)
    f_z5 <- function(z5) mu_f5_z5(z5) * mu_f789_z5(z5)
    f_z6 <- function(z6) mu_f6_z6(z6) * mu_f789_z6(z6)
    f_z7 <- function(z7) mu_f789_z7(z7)
    f_z8 <- function(z8) mu_f789_z8(z8)
    f_z9 <- function(z9) mu_f789_z9(z9)

    T1_SPA[1] <- f_z1(1) / (f_z1(1) + f_z1(0))
    T1_SPA[2] <- f_z2(1) / (f_z2(1) + f_z2(0))
    T1_SPA[3] <- f_z3(1) / (f_z3(1) + f_z3(0))
    T1_SPA[4] <- f_z4(1) / (f_z4(1) + f_z4(0))
    T1_SPA[5] <- f_z5(1) / (f_z5(1) + f_z5(0))
    T1_SPA[6] <- f_z6(1) / (f_z6(1) + f_z6(0))
    T1_SPA[7] <- f_z7(1) / (f_z7(1) + f_z7(0))
    T1_SPA[8] <- f_z8(1) / (f_z8(1) + f_z8(0))
    T1_SPA[9] <- f_z9(1) / (f_z9(1) + f_z9(0))

    T1_SPA
}

E_SPA9 <- function(dat, theta_t){
    T1 <- rep(NA, nrow(dat))
    for(FamID in unique(dat$FamID)){
        idx <- dat$FamID==FamID
        fam <- dat[idx, ]
        T1_fam <- E_SPA_onefam9(fam, theta_t)
        T1[idx] <- T1_fam
    }
    T1
}



################################################################
## Bootstrapping

bootstrap_families <- function(dat, B=25, seed=NULL){
    ## Generate families
    ## Input: dat: Data set
    ##        B: number of bootstrap data sets
    ##        seed: If you want reproducible results, set a random seed here
    ## Output: A list of length B, with a bootstrapped (unpreprocessed!) data set per list element

    FamIDs <- unique(dat$FamID)

    if(!is.null(seed)){
        set.seed(seed)
    }

    dat_B <- list()
    for(b in 1:B){

        FamIDs_b <- sort(sample(FamIDs, length(FamIDs), replace=TRUE))

        dat_b <- dat[numeric(0),]
        possible_Z_vectors_list_b <- list()

        for(i in 1:length(FamIDs)){
            FamID_old <- FamIDs_b[i]
            fam <- dat[dat$FamID==FamID_old,]
            fam$FamID=i  # new FamID
            dat_b <- rbind(dat_b, fam)
            possible_Z_vectors_list_b[[i]] <- attr(dat, "possible_Z_vectors_list")[[as.character(FamID_old)]]
            names(possible_Z_vectors_list_b)[i] <- as.character(i)
        }
        attr(dat_b, "possible_Z_vectors_list") <- possible_Z_vectors_list_b
        attr(dat_b, "b") <- b

        dat_B[[b]] <- dat_b

    }

    return(dat_B)
}
