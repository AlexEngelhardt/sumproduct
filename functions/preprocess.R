################################################################
#### Preprocessing

censor_CRCs <- function(dat, censor_at=80, censor_IP_only=FALSE){
    ## If people got CRC after say 80 years of age, this function instead censors these observations at t=80 (with c=0)
    idx <- dat$c==1 & dat$t > censor_at
    if(censor_IP_only){
        idx <- idx & dat$position==1
    }
    dat[idx,"c"] <- 0
    dat[idx,"t"] <- censor_at
    dat$censored_at <- censor_at
    return(dat)
}

## Two helper functions for possible_Z_vectors:
get_parents <- function(fam, pos){
    dad_pos <- fam[pos, "father_pos"]
    mom_pos <- fam[pos, "mother_pos"]

    return(c(dad_pos, mom_pos))
}

get_related_founders <- function(fam, pos){
    ## This function only works with imputed second parents. But we assume this anyway, after preprocessing.
    if(fam[fam$pos==pos,"founder"]){  # end of recursion
        return(pos)
    }

    ## recurse over mother and father:
    return( c(get_related_founders(fam, get_parents(fam, pos)[1]), get_related_founders(fam, get_parents(fam, pos)[2])) )
}

possible_Z_vectors <- function(dat, famID, pH){
    fam <- dat[dat$FamID==famID,]

    if(nrow(fam)==1){  # treat singles separately
        candidates <- data.frame(a=0:1)
        colnames(candidates) <- paste0("Z", fam$position)
        return(candidates)
    }

    founder_possibilities <- as.matrix(expand.grid(lapply(fam$position[fam$founder], function(i) 0:1)))
    nonfounder_possibilities <- as.matrix(expand.grid(lapply(fam$position[!fam$founder], function(i) 0:1)))
    colnames(founder_possibilities) <- paste0("Z", fam$position[fam$founder])
    colnames(nonfounder_possibilities) <- paste0("Z", fam$position[!fam$founder])

    candidates <- list()

    for(row in 1:nrow(founder_possibilities)){
        founder_constellation <- founder_possibilities[row,]

        sub_candidates <- cbind( rep(1, nrow(nonfounder_possibilities)) %*% t(founder_constellation),
                                nonfounder_possibilities
                                )
        sub_candidates <- sub_candidates[,paste0("Z",fam$position)]  # column order is wrong if nonfounders are far back

        for(nonfounder_col in which(!fam$founder)){
            dad_pos <- fam$father_pos[nonfounder_col]
            mom_pos <- fam$mother_pos[nonfounder_col]
            dad_col <- which(fam$position==dad_pos)
            mom_col <- which(fam$position==mom_pos)

            ptilde_this_nonfounder <- pH * sub_candidates[,dad_col] + pH * sub_candidates[,mom_col] - pH*pH * sub_candidates[,dad_col] * sub_candidates[,mom_col]

            keep_these_rows <- !( (ptilde_this_nonfounder==0 & sub_candidates[,nonfounder_col]==1) | (ptilde_this_nonfounder==1 & sub_candidates[,nonfounder_col]==0) )

            sub_candidates <- sub_candidates[keep_these_rows,,drop=FALSE]
        }
        candidates[[row]] <- sub_candidates
    }

    nrows <- sapply(candidates, nrow) %>% sum()

    ## combine all possible Z vectors by preallocation instead of sequential rbind() to save time
    candidates_df <- matrix(NA, nrow=nrows, ncol=nrow(fam), dimnames=list(NULL, paste0("Z", fam$position)))
    rowctr <- 1
    for(i in 1:length(candidates)){
        cand_len <- nrow(candidates[[i]])
        candidates_df[rowctr:(rowctr+cand_len-1), ] <- candidates[[i]]
        rowctr <- rowctr + cand_len
    }

    return(as.data.frame(candidates_df))
}

PZL <- function(dat, pH=1){
    ## This function loops over each family in 'dat' and calls possible_Z_vectors_list() on each family
    possible_Z_vectors_list <- list()
    if(all(table(dat$FamID, dat$position) == 1)){  # for simulated data where each family has the exact same structure, compute only one PZL and copy to all families.
        ## cat("All families have the same amount of members and position IDs\n",
        ##     "I am assuming these are simulated data where each family \n",
        ##     "has the same structure and am computing possible_Z_vectors only once.\n")
        PZV <- possible_Z_vectors(dat, dat$FamID[1], pH)
        for(fam in unique(dat$FamID)){
            possible_Z_vectors_list[[as.character(fam)]] <- PZV
        }
    } else {
        for(fam in unique(dat$FamID)){  # if it's real data (i.e. familys have different structures)
            print(fam)
            possible_Z_vectors_list[[as.character(fam)]] <- possible_Z_vectors(dat, fam, pH)
        }
    }
    return(possible_Z_vectors_list)
}

make_dataset <- function(dat, censor_at, remove_older, pH){
    ## This function is a convenient wrapper that takes a preprocessed data set and censors it at some age.
    ##  The parameter 'remove_older' specifies whether families whose indexpatients (i.e. position==1) CRC
    ##  age is greater than the censoring age should be removed or censored
    if(remove_older==FALSE){  # d.h. if censor instead of remove
        dat_new <- censor_CRCs(dat, censor_at=censor_at, censor_IP_only=FALSE)

        attr(dat_new, "possible_Z_vectors_list") <- attr(dat, "possible_Z_vectors_list")
    } else {
        rm_idx <- dat$position==1 & dat$c==1 & dat$t > censor_at
        FamIDs_to_keep <- unique(dat$FamID[!rm_idx])
        dat_new <- dat[dat$FamID %in% FamIDs_to_keep, ]

        attr(dat_new, "possible_Z_vectors_list") <- attr(dat, "possible_Z_vectors_list")[as.character(FamIDs_to_keep)]
    }
    return(dat_new)
}

preprocess <- function(raw_dat, impute=TRUE, pH=0.5, compute_PZL=FALSE){

    ## Main preprocessing function that imputes missing parents and performs grid purging, i.e. computes the possible Z vectors.
    
    dat <- raw_dat

    if(impute){

################################################################
#### Special families

        ppl_per_generation <- (table(dat$FamID, dat$generation))
        singles_FamIDs <- as.numeric(names((rowSums(ppl_per_generation)==1)[(rowSums(ppl_per_generation)==1)]))

#### siblings where both parents are missing:

        persons_without_parents <-  is.na(dat$father_pos) & is.na(dat$mother_pos)

        persons_without_children <- rep(NA, nrow(dat))
        for(p_row in 1:nrow(dat)){
            hisfamily <- dat$FamID[p_row]
            has_children <- ((dat$position[p_row] %in% dat[dat$FamID==hisfamily, "father_pos"]) | (dat$position[p_row] %in% dat[dat$FamID==hisfamily, "mother_pos"]))
            persons_without_children[p_row] <- !has_children
        }

        dat_has_1or2_parents <- !persons_without_parents
        dat_has_child <- !persons_without_children

        weird_families <- unique(dat$FamID[persons_without_parents & persons_without_children])
        nonsingle_idx <- !(weird_families %in% singles_FamIDs)
        weird_families_without_singles <- weird_families[nonsingle_idx]

        for(weird_FamID in weird_families_without_singles){
            fam <- dat[dat$FamID==weird_FamID, ]

            ## If there is only one generation, everyone is sibling and we impute noninformative parents:
            if(length(unique(fam$generation)) == 1){
                imputed_parents <- data.frame(FamID=weird_FamID, generation=fam$generation[1]-1,
                                              position=max(fam$position)+(1:2), father_pos=NA, mother_pos=NA, founder=TRUE, m=c(1,0), censored_at=Inf, c=0, t=0)
                dat[dat$FamID==weird_FamID, "father_pos"] <- max(fam$position)+1
                dat[dat$FamID==weird_FamID, "mother_pos"] <- max(fam$position)+2
                dat[(dat$FamID==weird_FamID) & (dat$generation==min(fam$generation)), "founder"] <- FALSE
                dat <- rbind(dat, imputed_parents)
                next
            }

            ## If in the oldest generation, no one has a child together, they are all siblings, and we impute noninformative parents:
            fam_2ndgen <- fam[fam$generation==min(fam$generation)+1,]
            if(all(is.na(fam_2ndgen$father_pos) | is.na(fam_2ndgen$mother_pos))){  # TRUE if no one in the first generation has a child together
                imputed_parents <- data.frame(FamID=weird_FamID, generation=min(fam$generation)-1,
                                              position=max(fam$position)+(1:2), father_pos=NA, mother_pos=NA, founder=TRUE, m=c(1,0), censored_at=Inf, c=0, t=0)
                dat[(dat$FamID==weird_FamID) & (dat$generation==min(fam$generation)), "father_pos"] <- max(fam$position)+1
                dat[(dat$FamID==weird_FamID) & (dat$generation==min(fam$generation)), "mother_pos"] <- max(fam$position)+2
                dat[(dat$FamID==weird_FamID) & (dat$generation==min(fam$generation)), "founder"] <- FALSE
                dat <- rbind(dat, imputed_parents)
            }

        }
        dat <- dat[order(dat$FamID, dat$generation, dat$position),]

#### If only one parent is missing, impute them:

        for(FamID in unique(dat$FamID)){
            ## FamID=2  # <DEBUG>
            fam <- dat[dat$FamID==FamID,]
            ein_elternteil_fehlt <- xor( is.na(fam$father_pos) , is.na(fam$mother_pos) )

            while(sum(ein_elternteil_fehlt) > 0){
                halbwaisen_zeile <- which(ein_elternteil_fehlt)[1]
                halbwaisen_pos <- fam$pos[halbwaisen_zeile]
                halbwaisen_gen <- fam$generation[halbwaisen_zeile]
                mom_fehlt <- is.na(fam$mother_pos[halbwaisen_zeile])
                vorhandener_elter_pos <-  ifelse(mom_fehlt, fam$father_pos[halbwaisen_zeile], fam$mother_pos[halbwaisen_zeile])
                impute_geschlecht <-      ifelse(mom_fehlt, 0, 1)
                fehlender_elternteil <-   ifelse(mom_fehlt, "mother_pos", "father_pos")
                vorhandener_elternteil <- ifelse(mom_fehlt, "father_pos", "mother_pos")
                impute_pos <- max(fam$position)+1

                impute_row <- data.frame(FamID=FamID, generation=halbwaisen_gen-1, position=impute_pos,
                                         father_pos=NA, mother_pos=NA,
                                         founder=TRUE, m=impute_geschlecht, censored_at=Inf, c=0, t=0)  # impute: c=0, t=0

                dat <- rbind(dat, impute_row)

                ## find this persons siblings and give them all a position for the missing parent:
                parent_pos <- fam[halbwaisen_zeile, vorhandener_elternteil]
                sibling_rows <- which( (dat$FamID==FamID) & (dat[[vorhandener_elternteil]] == parent_pos) & is.na(dat[[fehlender_elternteil]]) )
                dat[sibling_rows , fehlender_elternteil] <- impute_pos

                ## recompute b/c you could have fixed multiple orphans in this iteration
                fam <- dat[dat$FamID==FamID,]
                ein_elternteil_fehlt <- xor( is.na(fam$father_pos) , is.na(fam$mother_pos) )
            }
        }

    }

    dat <- dat[order(dat$FamID, dat$generation, dat$position),]

    ## There are c=0, t=NA in the original data. Replace this by censored @ t=0
    idx <- is.na(dat$t)
    dat[idx, "t"] <- 0
    dat[idx, "c"] <- 0

    ## attr(dat, "factor_graphs") <- make_factor_graphs(dat)  # this results in a very large size RDS object on the hard disk. Compute on-the-fly instead

    dat$censored_at <- Inf

    if(compute_PZL==TRUE){
        attr(dat, "possible_Z_vectors_list") <- PZL(dat, pH=pH)
    } else {
        attr(dat, "possible_Z_vectors_list") <- lapply(unique(dat$FamID), function(i) NA)
        names(attr(dat, "possible_Z_vectors_list")) <- as.character(unique(dat$FamID))

    }

    return(dat)
}

