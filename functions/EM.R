################################################################
#### Run EM Algorithm

E <- function(dat, theta_t){
    ## This is the E step without the sum-product algorithm, i.e. with a marginalization within (see the expand.grid call)
    
    ## goes per family and computes T_{R,I} first, then marginalizes out the T_{Id,1}

    p1_t <- theta_t$p1
    alpha_t <- theta_t$alpha

    FamIDs <- unique(dat$FamID)

    T1 <- rep(NA, nrow(dat))

    for(FamID in FamIDs){
        fam_idx <- dat$FamID==FamID
        fam <- dat[fam_idx,]

        if(!is.null(attr(dat, "possible_Z_vectors_list"))){
            Z_vectors <- attr(dat, "possible_Z_vectors_list")[[as.character(FamID)]]
        } else {  # this runs if the data frame does not have possible_Z_vectors_list precomputed, i.e. preprocessing did not run completely. This step should be avoided since it takes time in every E step instead of only once.
            Z_vectors <- expand.grid(lapply(fam$position, function(i) 0:1))
            colnames(Z_vectors) <- paste0("Z", fam$position)
        }

        Z_vectors$zaehler <- f_XI.ZI(fam, Z_vectors, p1_t, alpha_t)

        nenner <- sum(Z_vectors$zaehler)
        Z_vectors$T <- Z_vectors$zaehler / nenner  # T = PZgivenX, but a whole Z_I vector for a family. We want T1, the marginalized prob:

        for(d in 1:nrow(fam)){
            pos <- fam[d,"position"]
            T1[fam_idx][d] <- sum(Z_vectors$T[Z_vectors[,d]==1])  # marginalized out the other Z_Id
        }

    }
    return(T1)
}

M <- function(dat, T1){
    ## The M step is straightforward
    p1 <- mean(T1[dat$founder])  
    alpha <- sum(T1 * dat$c) / sum(T1 * (dat$t * lambda)^k * beta^dat$m)

    return(list(p1=p1, alpha=alpha))
}

entropy_of_Z <- function(theta_t){
    ## Previously used to assess convergence (instead, we now use the relative difference of p1 and alpha)
    ## The expected_loglik() plus the entorpy_of_Z() would result in the function B() from Dellaert's EM tutorial, that converges to the true log-likelihood
    return(NA)
}

expected_loglik <- function(theta_t, theta_tplus1){
    ## Previously used to assess convergence (instead, we now use the relative difference of p1 and alpha)
    return(NA)
}

run_EM <- function(dat, convergence_reltol=NA, EM_iterations=NA, theta_0=list(p1=0.2, alpha=4), verbose=TRUE, log_every=25, SPA=TRUE){

    if(is.na(convergence_reltol) && is.na(EM_iterations)){
        stop("You must set one of the arguments convergence_reltol or EM_iterations to define a stopping criterion!")
    }
    if(!is.na(convergence_reltol) && !is.na(EM_iterations)){
        stop("You can only set one of the arguments convergence_reltol and EM_iterations to define a stopping criterion!")
    }

    theta_t <- theta_0

    if(!is.na(EM_iterations)){  # run the EM for a specified number of iterations
        theta_trace <- list(p1=rep(theta_t$p1, EM_iterations), alpha=rep(theta_t$alpha, EM_iterations),
                            optim_loglik=rep(NA, EM_iterations),
                            expected_loglik=rep(NA, EM_iterations),
                            entropy=rep(NA, EM_iterations),
                            B=rep(NA, EM_iterations)
                            )
        T1_trace <- matrix(NA, nrow=EM_iterations, ncol=nrow(dat))

        elapsed <- system.time({

            for(i in 2:EM_iterations){
                if(verbose && !(i%%log_every)) print(paste(date(), ":: Iteration", i))

                if(SPA == TRUE){
                    T1 <- E_SPA(dat, theta_t)
                } else {
                    T1 <- E(dat, theta_t)
                }

                theta_tplus1 <- M(dat, T1)
                ##
                ## if(!i%%10){
                ##     theta_trace$entropy[i] <- entropy_of_Z(theta_t)
                ##     theta_trace$expected_loglik[i] <- expected_loglik(theta_t, theta_tplus1)
                ##     theta_trace$B[i] <- theta_trace$expected_loglik[i] + theta_trace$entropy[i]
                ## }
                ##
                T1_trace[i,] <- T1
                theta_t <- theta_tplus1
                theta_trace$p1[i] <- theta_t$p1
                theta_trace$alpha[i] <- theta_t$alpha
            }

        })
    } else if(!is.na(convergence_reltol)){  # run the EM until the parameters don't change much anymore

        theta_trace <- list(p1=theta_t$p1, alpha=theta_t$alpha,
                            optim_loglik=NA,
                            expected_loglik=NA,
                            entropy=NA,
                            B=NA
                            )
        T1_trace <- matrix(NA, nrow=1, ncol=nrow(dat))


        reltol_this_iter <- convergence_reltol * 2
        i <- 1

        elapsed <- system.time({
            while(reltol_this_iter > convergence_reltol){
                i <- i+1

                if(SPA == TRUE){
                    T1 <- E_SPA(dat, theta_t)
                } else {
                    T1 <- E(dat, theta_t)
                }

                theta_tplus1 <- M(dat, T1)
                ##
                ## if(!i%%10){
                ##     theta_trace$entropy[i] <- entropy_of_Z(theta_t)
                ##     theta_trace$expected_loglik[i] <- expected_loglik(theta_t, theta_tplus1)
                ##     theta_trace$B[i] <- theta_trace$expected_loglik[i] + theta_trace$entropy[i]
                ## }
                ##
                T1_trace <- rbind(T1_trace, T1)
                theta_t <- theta_tplus1
                theta_trace$p1 <- c(theta_trace$p1, theta_t$p1)
                theta_trace$alpha <- c(theta_trace$alpha, theta_t$alpha)

                reltol_this_iter_p1 <- abs( (theta_trace$p1[i] - theta_trace$p1[i-1]) / theta_trace$p1[i-1] )
                reltol_this_iter_alpha <- abs( (theta_trace$alpha[i] - theta_trace$alpha[i-1]) / theta_trace$alpha[i-1] )
                reltol_this_iter <- max(reltol_this_iter_p1, reltol_this_iter_alpha)

                if(verbose && !(i%%log_every)) print(paste(date(), ":: Iteration", i, "rel.tol =", reltol_this_iter))
            }
        })
    }

    attr(theta_trace, "elapsed") <- elapsed

    return(theta_trace)
}

################################################################
#### EM Algo with the sum-product algorithm

make_factor_graph <- function(fam){
    ## creates a list in relational database structure containing a factor graph for a family
    FG <- list()

    ## Variable nodes
    dad_pos <- fam$position[match(fam$father_pos, fam$position)]
    mom_pos <- fam$position[match(fam$mother_pos, fam$position)]

    FG$variable_nodes <- data.frame(
        name=paste0("z_", fam$position),
        is_founder=fam$founder,
        is_male=fam$m,
        dad=ifelse(is.na(dad_pos), NA, paste0("z_", dad_pos)),
        mom=ifelse(is.na(mom_pos), NA, paste0("z_", mom_pos))
    )

    FG$variable_nodes$name <- as.character(FG$variable_nodes$name)
    FG$variable_nodes$dad <- as.character(FG$variable_nodes$dad)
    FG$variable_nodes$mom <- as.character(FG$variable_nodes$mom)

    ## Factor nodes
    founder_factors <- paste("f", fam[fam$founder, "position"], sep="_")
    parents_factors <- paste("f", fam$father_pos[!fam$founder], fam$mother_pos[!fam$founder], sep="_") %>% unique()

    FG$factor_nodes <- data.frame(
        name=c(founder_factors, parents_factors),
        is_founder=c(rep(TRUE, length(founder_factors)), rep(FALSE, length(parents_factors)))
    )
    FG$factor_nodes$name <- as.character(FG$factor_nodes$name)

    ## Edges
    FG$edges <- data.frame(variable_name=character(0), factor_name=character(0), var_is_parent=logical(0), z_pos=numeric(0))
    for(pos in fam$position){
        variable_node <- paste0("z_", pos)
        parents_pos <- fam[fam$position==pos, c("father_pos", "mother_pos")]
        factor_up <- ifelse(all(is.na(parents_pos)), paste0("f_", pos), paste("f", parents_pos[1], parents_pos[2], sep="_"))  # if founder bzw. nonfounder
        FG$edges <- rbind(FG$edges, data.frame(variable_name=variable_node, factor_name=factor_up, var_is_parent=FALSE, z_pos=pos))

        if(length(which(fam$father_pos == pos | fam$mother_pos == pos)) > 0){  # create factor for kids
            factors_down <- paste("f", fam[which(fam$father_pos == pos | fam$mother_pos == pos), "father_pos"], fam[which(fam$father_pos == pos | fam$mother_pos == pos), "mother_pos"], sep="_") %>% unique()

            FG$edges <- rbind(FG$edges, data.frame(variable_name=variable_node, factor_name=factors_down, var_is_parent=TRUE, z_pos=pos))
        }
    }
    FG$edges$variable_name <- as.character(FG$edges$variable_name)
    FG$edges$factor_name <- as.character(FG$edges$factor_name)

    ## Messages
    messages_from_variable <- data.frame(
        name = paste("mu", FG$edges$variable_name, FG$edges$factor_name, sep="."),
        from = FG$edges$variable_name,
        to = FG$edges$factor_name,
        z_pos = FG$edges$z_pos,
        from_factor = FALSE
    )
    messages_from_factor <- data.frame(
        name = paste("mu", FG$edges$factor_name, FG$edges$variable_name, sep="."),
        from = FG$edges$factor_name,
        to = FG$edges$variable_name,
        z_pos = FG$edges$z_pos,
        from_factor = TRUE
    )
    FG$messages <- rbind(messages_from_variable, messages_from_factor)
    FG$messages$status="waiting"
    FG$messages$zero=NA
    FG$messages$one=NA
    FG$messages$to <- as.character(FG$messages$to)
    FG$messages$from <- as.character(FG$messages$from)
    FG$messages$status <- as.character(FG$messages$status)
    FG$messages$name <- as.character(FG$messages$name)

    return(FG)
}

## make_factor_graphs <- function(dat){
##     FGs <- list()
##     for(famid in unique(dat$FamID)){
##         FGs[[as.character(famid)]] <- make_factor_graph(dat[dat$FamID==famid,])
##     }
##     return(FGs)
## }

make_factor <- function(fam, factor, theta_t){  # This is a closure, it returns a function. See Hadley Wickham, Advanced R

    founder_factors <- paste("f", fam[fam$founder, "position"], sep="_")
    parents_factors <- paste("f", fam$father_pos[!fam$founder], fam$mother_pos[!fam$founder], sep="_") %>% unique()

    FG <- attr(fam, "FG")
    if(is.null(FG)){
        stop("object 'fam' does not have the factor graph as an attribute named 'FG'!")
    }

    alpha_t <- theta_t$alpha
    p1_t <- theta_t$p1

    if(factor %in% founder_factors){

        position <- as.numeric(strsplit(factor, "_")[[1]][2])
        fam_row <- which(fam$position==position)

        return(
            function(z){
                fX <- f(t=fam[fam_row, "t"], r=z, c=fam[fam_row, "c"], m=fam[fam_row, "m"], alpha=alpha_t, indexpatient=as.numeric(fam[fam_row, "position"]==1), censor_age=fam[fam_row, "censored_at"])
                pZ <- p1_t^(z) * (1-p1_t)^(1-(z))  # pZ for founders

                fX*pZ
            }
        )

    } else {  # nonfounder factor (i.e. with parents)

        edges <- subset(FG$edges, factor_name==factor)
        edges <- left_join(edges, FG$variable_nodes, by=c("variable_name"="name"))

        dad_row <- which(edges$var_is_parent & edges$is_male)  # row within 'edges'!
        mom_row <- which(edges$var_is_parent & !edges$is_male) # row within 'edges'!
        kids_rows_z <- which(!edges$var_is_parent)
        kids_pos <- edges$z_pos[kids_rows_z]
        kids_rows <-   which(fam$position %in% kids_pos)       # row within 'fam'!


        return(
            function(z){
                stopifnot(length(z)==nrow(edges))
                fX <- f(t=fam[kids_rows, "t"], r=z[kids_rows_z], c=fam[kids_rows, "c"], m=fam[kids_rows, "m"], alpha=theta_t$alpha, indexpatient=as.numeric(fam[kids_rows, "position"]==1), censor_age=fam[kids_rows, "censored_at"])

                p1tilde <- pH*z[dad_row] + pH*z[mom_row] - pH*pH*z[dad_row]*z[mom_row]
                pZ <- p1tilde ^ z[kids_rows_z] * (1-p1tilde) ^ (1-z[kids_rows_z])

                prod(fX) * prod(pZ)
            }
        )
    }
}

make_factors <- function(fam, theta_t){
    founder_factors <- paste("f", fam[fam$founder, "position"], sep="_")
    parents_factors <- paste("f", fam$father_pos[!fam$founder], fam$mother_pos[!fam$founder], sep="_") %>% unique()

    factors <- list()

    for(factor in founder_factors){
        factors[[factor]] <- make_factor(fam, factor, theta_t)
    }
    for(factor in parents_factors){
        factors[[factor]] <- make_factor(fam, factor, theta_t)
    }
    return(factors)
}

compute_messages_onefam <- function(fam, theta_t){
    p1_t <- theta_t$p1
    alpha_t <- theta_t$alpha

    ## This function can be compartmentalized even further by making a list of closures:
    ##  make_message() -> msg[["f_5_6"]][["z_9"]](z9=0)

    fam_factors <- make_factors(fam, theta_t)

    messages <- attr(fam, "FG")$messages

    while(sum(messages$status=="waiting") > 0){  # for as many iterations until all messages are computed

        ## this for loop determines all messages that can be calculated in this iteration:
        for(row in which(messages$status=="waiting")){
            sender <- messages$from[row]
            destination <- messages$to[row]

            dependencies <- (messages$to == sender) & (messages$from != destination)

            if((sum(dependencies) == 0) || (all(messages$status[dependencies] == "done"))){
                messages[row, "status"] <- "calculating"
            }
        }

        ## this for loop computes all messages that can be calculated in this iteration:
        for(row in which(messages$status=="calculating")){
            sender <- messages$from[row]
            destination <- messages$to[row]
            from_variable <- substr(sender, 1, 1) == "z"

            if(from_variable){  # if message sender is variable node
                incoming_rows <- which( (messages$to == sender) & (messages$from != destination) )

                if(length(incoming_rows) != 0){
                    messages[row, "zero"] <- prod(messages[incoming_rows, "zero"])
                    messages[row, "one"] <- prod(messages[incoming_rows, "one"])
                } else {  # leaf node, i.e. no other incoming messages on variable node
                    messages[row, "zero"] <- 1
                    messages[row, "one"] <- 1
                }
            } else {  # if message sender is factor node
                all_connected_variables <- messages[messages$from==sender, "to"]
                connected_variables <- messages[messages$from==sender & messages$to!=destination, "to"]

                this_factor <- fam_factors[[sender]]

                if(length(connected_variables)==0){  # i.e. if this is a founder factor
                    messages[row, "zero"] <- this_factor(0)
                    messages[row, "one"] <- this_factor(1)
                } else {  # i.e. message leaves from a nonfounder factor
                    destination_pos <- as.numeric(strsplit(destination, "_")[[1]][2])  # position of the destination variable node (i.e. 5 if it's "z_5")
                    x_Id <- fam[fam$position==destination_pos, ]

                    ## find out if the destination is a child or parent *within the factor*:
                    var_is_parent <- subset(attr(fam, "FG")$edges, subset=(factor_name==sender & variable_name==destination), select="var_is_parent") %>% as.logical()

                    if(var_is_parent){  # if destination is a parent, "big" sum only over one parent

                        other_parent <- subset(attr(fam, "FG")$edges, subset=(factor_name==sender & var_is_parent==TRUE & variable_name!=destination), select="variable_name")[[1]]
                        other_parent_vals <- 0:1
                        child_variables <- subset(attr(fam, "FG")$edges, subset=factor_name==sender & var_is_parent==FALSE & variable_name!=destination,
                                                  select="variable_name")[[1]] # "& variable_name!=destination" is actually not necessary here because destination is always a parent in this loop. However, for further generalization of the code, I leave it in.

                        messages[row, "zero"] <- sapply(other_parent_vals, function(z_other){
                            messages[messages$from==other_parent & messages$to==sender, ifelse(z_other==0, "zero", "one")] *
                                sapply(child_variables, function(child_var){
                                    x_child <- fam[fam$position==as.numeric(strsplit(child_var, "_")[[1]][2]), ]
                                    sapply(0:1, function(z_kid){
                                        messages[messages$from==child_var & messages$to==sender, ifelse(z_kid==0, "zero", "one")] *
                                            f(x_child$t, z_kid, x_child$c, x_child$m, alpha_t, x_child$position==1, x_child$censored_at) *
                                                (pH*0 + pH*z_other - pH^2*0*z_other)^z_kid * (1-(pH*0 + pH*z_other - pH^2*0*z_other))^(1-z_kid)
                                    }) %>% sum()
                                }) %>% prod()
                        }) %>% sum()

                        messages[row, "one"] <- sapply(other_parent_vals, function(z_other){
                            messages[messages$from==other_parent & messages$to==sender, ifelse(z_other==0, "zero", "one")] *
                                sapply(child_variables, function(child_var){
                                    x_child <- fam[fam$position==as.numeric(strsplit(child_var, "_")[[1]][2]), ]
                                    sapply(0:1, function(z_kid){
                                        messages[messages$from==child_var & messages$to==sender, ifelse(z_kid==0, "zero", "one")] *
                                            f(x_child$t, z_kid, x_child$c, x_child$m, alpha_t, x_child$position==1, x_child$censored_at) *
                                                (pH*1 + pH*z_other - pH^2*1*z_other)^z_kid * (1-(pH*1 + pH*z_other - pH^2*1*z_other))^(1-z_kid)
                                    }) %>% sum()
                                }) %>% prod()
                        }) %>% sum()


                    } else {  # if destination is a child, "big" sum goes over both parents
                        mom_vals <- 0:1  # possible values for the mother's Z_Id
                        dad_vals <- 0:1  # possible values for the father's Z_Id

                        parent_variables <- subset(attr(fam, "FG")$edges, subset=(factor_name==sender & var_is_parent==TRUE), select="variable_name")[[1]]
                        dad_variable <- subset(attr(fam, "FG")$variable_nodes, name==destination)$dad
                        mom_variable <- subset(attr(fam, "FG")$variable_nodes, name==destination)$mom
                        child_variables <- subset(attr(fam, "FG")$edges, subset=(factor_name==sender & var_is_parent==FALSE & variable_name != destination), select="variable_name")[[1]]
                        messages[row, "zero"] <- sapply(dad_vals, function(z_dad){
                            sapply(mom_vals, function(z_mom){

                                if(length(child_variables) > 0){  # children other than the destination child
                                    child_prods <- sapply(child_variables, function(child_var){
                                        x_child <- fam[fam$position==as.numeric(strsplit(child_var, "_")[[1]][2]), ]
                                        sapply(0:1, function(z_kid){
                                            messages[messages$from==child_var & messages$to==sender, ifelse(z_kid==0, "zero", "one")] *
                                                f(x_child$t, r=z_kid, x_child$c, x_child$m, alpha_t, x_child$pos==1, x_child$censored_at) *
                                                    (1-(pH*z_dad + pH*z_mom - pH^2*z_dad*z_mom))^(1-z_kid) * (pH*z_dad + pH*z_mom - pH^2*z_dad*z_mom)^(z_kid)
                                        }) %>% sum()
                                    }) %>% prod()
                                } else {
                                    child_prods <- 1
                                }

                                prod(
                                    messages[messages$to == sender & messages$from == dad_variable, ifelse(z_dad==1, "one", "zero")],
                                    messages[messages$to == sender & messages$from == mom_variable, ifelse(z_mom==1, "one", "zero")]
                                ) * f(x_Id$t, r=0, x_Id$c, x_Id$m, alpha_t, destination_pos==1, x_Id$censored_at) * (1-(pH*z_dad + pH*z_mom - pH^2*z_dad*z_mom)) *
                                    child_prods
                            }) %>% sum()
                        }) %>% sum()

                        messages[row, "one"] <- sapply(dad_vals, function(z_dad){
                            sapply(mom_vals, function(z_mom){

                                if(length(child_variables) > 0){
                                    child_prods <- sapply(child_variables, function(child_var){
                                        x_child <- fam[fam$position==as.numeric(strsplit(child_var, "_")[[1]][2]), ]
                                        sapply(0:1, function(z_kid){
                                            messages[messages$from==child_var & messages$to==sender, ifelse(z_kid==0, "zero", "one")] *
                                                f(x_child$t, r=z_kid, x_child$c, x_child$m, alpha_t, x_child$pos==1, x_child$censored_at) *
                                                    (1-(pH*z_dad + pH*z_mom - pH^2*z_dad*z_mom))^(1-z_kid) * (pH*z_dad + pH*z_mom - pH^2*z_dad*z_mom)^(z_kid)
                                        }) %>% sum()
                                    }) %>% prod()
                                } else {
                                    child_prods <- 1
                                }

                                prod(
                                    messages[messages$to == sender & messages$from == dad_variable, ifelse(z_dad==1, "one", "zero")],
                                    messages[messages$to == sender & messages$from == mom_variable, ifelse(z_mom==1, "one", "zero")]
                                ) * f(x_Id$t, r=1, x_Id$c, x_Id$m, alpha_t, destination_pos==1, x_Id$censored_at) * (pH*z_dad + pH*z_mom - pH^2*z_dad*z_mom) *
                                    child_prods
                            }) %>% sum()
                        }) %>% sum()
                    }
                }
            }
        }

        ## set all statuses from "calculating" to "done"
        messages$status[messages$status=="calculating"] <- "done"
    }

    return(messages)
}

E_SPA_onefam <- function(fam, theta_t){
    T1 <- rep(NA, nrow(fam))
    msg <- compute_messages_onefam(fam, theta_t)
    for(row in 1:nrow(fam)){
        pos <- fam[row, "position"]
        this_msgs <- msg[msg$to == paste0("z_", pos), ]
        T1[row] <- prod(this_msgs$one) / (prod(this_msgs$zero) + prod(this_msgs$one))
    }

    return(T1)
}

E_SPA <- function(dat, theta_t){  # We run the E step either as a SPA or with a marginalization, depending on the family size (SPA_cutoff)
    
    SPA_cutoff <- 1  # The minimum family size at which we use the SPA instead of a marginalization
    
    T1 <- rep(NA, nrow(dat))
    for(FamID in unique(dat$FamID)){
        idx <- dat$FamID==FamID
        fam <- dat[idx, ]
        if(nrow(fam) > SPA_cutoff){
            attr(fam, "FG") <- make_factor_graph(fam)  # ideally, this would be done once, instead of in each E step. But saving closures somehow results in huge or infinite file sizes
            T1_fam <- E_SPA_onefam(fam, theta_t)
        } else {
            T1_fam <- E(fam, theta_t)
        }
        T1[idx] <- T1_fam
    }
    T1
}
