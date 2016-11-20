################################################################
#### f(), the Weibull densities

f <- function(t, r, c, m, alpha, indexpatient=rep(0, length(t)), censor_age=Inf){
    ## The numerator is the Weibull density.
    ## The denominator is the CDF F(), and is only in effect for indexpatients with CRC and turns the density into a conditional density,
    ##   conditioned on the fact that we know that CRC for indexpatients appeared before the age of 'censor_age'.
    ## If censor_age==Inf, then we don't use conditional densities:
    ##  The division by the CDF is deactivated by setting censor_age to Inf, since the term then reduces to 1
    (k * lambda^k * t^(k-1) * alpha^r * beta^m)^c * exp(-(t*lambda)^k * alpha^r * beta^m) / ((1-exp(-(censor_age*lambda)^k * alpha^r * beta^m))^(indexpatient*c))
}

## A C++ implementation of rowProds (like R's rowSums but with products) since this has to be fast in the Nelder-Mead optimization
cppFunction('
NumericVector rowProdsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 1;
    for (int j = 0; j < ncol; j++) {
      total *= x(i, j);
    }
    out[i] = total;
  }
  return out;
}
')

f_XI.ZI <- function(fam, Z_vectors, p1, alpha){
    ## f(X_I, Z_I | theta)
    ## Helper for L_fam() and E(). Computes for all Z possibilities: f(X_I, Z_I | theta),
    ##  then returns the vector of these f's.
    ## Optim will need the sum of this vector, the naive EM algorithm every single element for T_{R,I}^{(t)}
    Z_probs <- as.data.frame(matrix(0, nrow=nrow(Z_vectors), ncol=ncol(Z_vectors), dimnames=list(NULL, colnames(Z_vectors))))
    f_densities <- Z_probs

    ## precompute (memoize) densities instead of recomputing them in every iteration:
    fX_Id_givenZ <- cbind(f(fam$t, 0, fam$c, fam$m, alpha, indexpatient=as.numeric(fam$position==1), censor_age=fam$censored_at),
                          f(fam$t, 1, fam$c, fam$m, alpha, indexpatient=as.numeric(fam$position==1), censor_age=fam$censored_at))  

    for(row in 1:nrow(fam)){  # vector-value precomputing of P(Z_Id = 1 | other_Z_I, theta_t)
        if(fam$founder[row]){
            thisp1 <- p1
        } else {
            dad_row <- which(fam$position == fam[row, "father_pos"])
            mom_row <- which(fam$position == fam[row, "mother_pos"])
            thisp1 <- pH*Z_vectors[,dad_row] + pH*Z_vectors[,mom_row] - pH*Z_vectors[,dad_row]*pH*Z_vectors[,mom_row]
        }
        Z_probs[,row] <- thisp1^Z_vectors[,row] * (1-thisp1)^(1-Z_vectors[,row])  # geht evtl schneller wenn du 2 Werte precomputest und extractest
        ##
        f_densities[,row] <- fX_Id_givenZ[matrix(c(rep(row, nrow(f_densities)), Z_vectors[,row]+1), ncol=2)]
    }

    f_XI.ZI <- rowProdsC(as.matrix(Z_probs)) * rowProdsC(as.matrix(f_densities))

    return(f_XI.ZI)
}

################################################################
#### Grid search / optim full likelihood functions

L_fam <- function(fam, riskconstellations, p1, alpha, gridpurging=TRUE){  # un-logged L for one family
    stopifnot(diff(range(fam$FamID))==0)  # Throws an error if you supply more than one family. If so, then subset dat before.

    if(gridpurging==FALSE){
        riskconstellations <- expand.grid(lapply(fam$position, function(i) 0:1))
    }
        
    all_f_XI.ZI <- f_XI.ZI(fam, Z_vectors=riskconstellations, p1, alpha)


    return(sum(all_f_XI.ZI))
}

l <- function(dat, p1, alpha, parallel=FALSE, gridpurging=TRUE){  # logged l over all families
    FamIDs <- unique(dat$FamID)

    if(!parallel){
        my.apply <- lapply
    } else {
        my.apply <- mclapply
    }

    likes <- my.apply(FamIDs, function(fam){
        L_fam(dat[dat$FamID==fam,], attr(dat, "possible_Z_vectors_list")[[as.character(fam)]], p1, alpha, gridpurging=gridpurging)
    }) %>% unlist()
    loglikes <- log(likes)
    loglike <- sum(loglikes)
    return(loglike)
}

################################################################
#### Run Nelder-Mead optimization

logit <- function(p) log(p/(1-p))
invlogit <- function(x) exp(x) / (1+exp(x))

run_optim <- function(dat, parallel=FALSE, theta_0=list(p1=0.2, alpha=4)){
    trace_ <- 0  # level of verbose during optim run

    elapsed <- system.time({
        opt <- optim(par=c(logit(theta_0$p1), theta_0$alpha), function(pars) l(dat, invlogit(pars[1]), pars[2], parallel=parallel, gridpurging=FALSE),
                     control=list(REPORT=1, trace=trace_, fnscale=-1))
    })
    attr(opt, "elapsed") <- elapsed
    return(opt)
}
