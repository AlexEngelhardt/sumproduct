################################################################
#### Simulation

invSurv <- function(u, risk, sex, alpha, lambda=0.0058, k=4, beta=2)
    ## inverse Weibull survival function. Used to simulate survival times from random uniformly distributed 'u'
    ((-log(u) / ((alpha^risk) * (beta^sex)))^(1/k))  /  lambda

simulate_single <- function(FamID, p1, alpha, pH){
    ## Function to simulate a pedigree of one person
    fam <- data.frame(ID=paste0(FamID, "_1"),
                      FamID=FamID,
                      generation=1,
                      position=1,
                      father_pos=NA,
                      mother_pos=NA,
                      founder=TRUE,
                      m=sample(0:1,1),
                      risk=runif(1) < p1,
                      censored_at=Inf
                      )
    fam$AlterCRC <- invSurv(runif(1), r, m, alpha) %>% round()
    return(fam)
}

simulate_family <- function(FamID, n_generations, p1, alpha, pH, n_lastgen=1){
    ## Function to simulate a whole family of arbitrary generations

    ## simulate a single person if wanted:
    if(n_generations == 1){
        return(simulate_single(FamID, p1, alpha, pH))
    }

    ## Simulate founder generation:
    n_founders <- 2^(n_generations - 1)

    fam <- data.frame(ID=paste0(FamID, "_", 1:n_founders),
                      FamID=FamID,
                      generation=1,
                      position=1:n_founders,
                      father_pos=NA,
                      mother_pos=NA,
                      founder=TRUE,
                      m=rep(c(1,0), length.out=n_founders),
                      risk=runif(n_founders) < p1,
                      censored_at=Inf
                      )


    ## Simulate nonfounder generations:
    if(n_generations >= 2){
        for(gen in 2:n_generations){

            n_thisgen <- 2^(n_generations-gen)

            if((gen == n_generations) && n_lastgen > 1){
                n_thisgen = n_lastgen
            }

            positions <- (max(fam$position)+1):(max(fam$position)+n_thisgen)

            fam_gen <- data.frame(ID=paste0(FamID, "_", positions),
                                  FamID=FamID,
                                  generation=gen,
                                  position=positions,
                                  father_pos=fam[fam$generation==(gen-1) & fam$m==1, "position"],
                                  mother_pos=fam[fam$generation==(gen-1) & fam$m==0, "position"],
                                  founder=FALSE,
                                  m=rep(c(1,0), length.out=n_thisgen),  # this way, the last child is always male. It should not matter in the simulation.
                                  risk=NA,  # will be filled in the next command
                                  censored_at=Inf
                                  )

            ## filling the risk property according to inheritance rules:
            risk_father <- fam[fam$generation == (gen-1) & fam$position %in% fam_gen$father_pos, "risk"]
            risk_mother <- fam[fam$generation == (gen-1) & fam$position %in% fam_gen$mother_pos, "risk"]
            p1_tilde <- pH * risk_father + pH * risk_mother + pH*pH*risk_father*risk_mother
            fam_gen$risk <- runif(n_thisgen) < p1_tilde

            fam <- rbind(fam, fam_gen)  # inefficient growing. Preallocate if time becomes an issue
        }
    }

    fam$ageCRC <- invSurv(runif(nrow(fam)), fam$risk, fam$m, alpha) %>% round()
    fam$ageCensoring <- round(rnorm(nrow(fam), mean=125, sd=10)-4*fam$m)  # 66% censored
    fam$c <- as.numeric(fam$ageCRC < fam$ageCensoring)
    fam$t <- pmin(fam$ageCRC, fam$ageCensoring)

    ## delete unnecessary columns:
    fam <- select(fam, -ID, -ageCRC, -ageCensoring)
    ## fam <- select(fam, -ID, -ageCRC, -ageCensoring, -risk)  # keep them (especially 'risk') if you want to see simulation details

    return(fam)
}

simulate_dataset <- function(n_familys, n_generations, p1, alpha, pH, n_lastgen=1, seed=NULL){
    ## Simulate an entire dataset of n families, each with specified number of generations
    if(!is.null(seed)){
        set.seed(seed)
    }
    dat <- lapply(1:n_familys, function(i) simulate_family(i, n_generations, p1, alpha, pH, n_lastgen)) %>% Reduce(rbind, .)
    return(dat)
}

kill_ancestors <- function(raw_fam, remove_pos){
    ## Helper function for remove_persons(). This takes a family and a position as input, and
    ##  removes that position and all of its ancestors from the data set.
    kill_pos <- integer(0)
    that_guy <- raw_fam[raw_fam$pos == remove_pos, ]

    if(that_guy$founder == TRUE){
        return(integer(0))  # watch out, when using this with negative indexing it removes everyone.
    }
    
    parents <- raw_fam[raw_fam$pos %in% c(that_guy$father_pos, that_guy$mother_pos), ]
    parents <- parents[!is.na(parents),]  # if a single parent is NA.
    
    kill_pos <- c(kill_pos, parents[, "position"], kill_ancestors(raw_fam, parents[1, "position"]))
    if(nrow(parents) == 2){  # only append second parent if it is not already NA
        kill_pos <- c(kill_pos, kill_ancestors(raw_fam, parents[2, "position"]))
    }
    
    return(kill_pos)
}

remove_persons <- function(raw_dat, frac){
    ## This function removes a fraction of 'frac' persons from the dataset raw_dat.
    ## It is used to create more realistic pedigrees where sometimes only one parent is available.
    ##  With this simulated data, we verify our imputing procedure.
    stopifnot(frac>0 && frac<1)
    
    nrow_before <- nrow(raw_dat)

    while(nrow(raw_dat) > (1-frac)*nrow_before){
        remove_idx <- sample(1:nrow(raw_dat), 1)
        FamID <- raw_dat[remove_idx, "FamID"]
        raw_fam <- raw_dat[raw_dat$FamID==FamID,]
        remove_pos <- raw_dat[remove_idx, "position"]

        removing_founder <- raw_dat[remove_idx, "founder"] == TRUE
        removing_terminal_node <- !(remove_pos %in% raw_fam$father_pos) && !(remove_pos %in% raw_fam$mother_pos)
        
        if( removing_founder ){
            ## set his/her children to father/mother=NA:
            NA_col <- ifelse(raw_dat[remove_idx, "m"]==1, "father_pos", "mother_pos")
            raw_dat[(raw_dat$FamID == FamID) & !is.na(raw_dat[[NA_col]]) & (raw_dat[[NA_col]] == remove_pos), NA_col] <- NA
        } else if ( !removing_terminal_node ) {
            ## set his/her children to father/mother=NA:
            NA_col <- ifelse(raw_dat[remove_idx, "m"]==1, "father_pos", "mother_pos")
            raw_dat[(raw_dat$FamID == FamID) & !is.na(raw_dat[[NA_col]]) & (raw_dat[[NA_col]] == remove_pos), NA_col] <- NA

            ## remove parents and all their parents too:
            kill_pos <- kill_ancestors(raw_fam, remove_pos)
            raw_dat <- raw_dat[!(raw_dat$FamID==FamID & raw_dat$position %in% kill_pos),]
            
        }  # else: removing terminal node, which would be fine, i.e. no additional processing necessary

        raw_dat <- raw_dat[!(raw_dat$FamID==FamID & raw_dat$position == remove_pos),]
        ## In case you removed both parents from some nonfounder in two different iterations, set them to founder=TRUE now:
        raw_dat[is.na(raw_dat$mother_pos) & is.na(raw_dat$father_pos), "founder"] <- TRUE

    }

    return(raw_dat)
}


