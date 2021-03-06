## auxiliary functions for dense setting
# (01) scenario.2014CLX : two-means      test (4 means, 8 covariances)
# (02) scenario.2013CLX : two-covariance test (4 covariances)
# (03) gen.cov.2014CLX  : 8 covariance matrix generator (scenario.2014CLX)
# (04) gen.2013CLX.mine  : sparse covariance difference setting
# (05) gen.2013CLX.dense : dense  covariance difference setting



# (01) scenario.2014CLX ---------------------------------------------------
#' @keywords internal
#' @noRd
scenario.2014CLX <- function(n, p, no.mean=1, no.cov=1, htype=c("H0","H1")){
  # put all the relevant functions here for clarity
  scenario.2014CLX.mean1 = function(n, p){
    m    = max(1,floor(0.05*p))
    vals = sample(c(-1,1), m, replace=TRUE)*sqrt(log(p)/n)

    mu0  = rep(0,p)
    mu0[sample(1:p, m, replace=FALSE)] = vals
    return(mu0)
  }
  scenario.2014CLX.mean2 = function(n, p){
    m = max(1,floor(0.05*p))
    bound = sqrt(8*log(p)/n)
    vals  = runif(m, min=-bound, max=bound)

    mu0  = rep(0,p)
    mu0[sample(1:p, m, replace=FALSE)] = vals
    return(mu0)
  }
  scenario.2014CLX.mean3 = function(n, p){
    m = max(1,sqrt(p))
    vals = sample(c(-1,1), m, replace=TRUE)*sqrt(log(p)/n)

    mu0  = rep(0,p)
    mu0[sample(1:p, m, replace=FALSE)] = vals
    return(mu0)
  }
  scenario.2014CLX.mean4 = function(n, p){
    m = max(1,sqrt(p))
    bound = sqrt(8*log(p)/n)
    vals  = runif(m, min=-bound, max=bound)

    mu0  = rep(0,p)
    mu0[sample(1:p, m, replace=FALSE)] = vals
    return(mu0)
  }
  scenario.2014CLX.cov1 = function(n, p){
    Sig = array(0,c(p,p))
    for (k in 1:floor(p/2)){
      lowk = as.integer(2*(k-1)+1)
      uppk = as.integer(2*k)
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          Sig[i,j] = 0.8
        }
      }
    }
    diag(Sig) = 1
    return(Sig)
  }
  scenario.2014CLX.cov2 = function(n, p){
    Sig = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        Sig[i,j] = 0.6^(abs(i-j))
      }
    }
    return(Sig)
  }
  scenario.2014CLX.cov3 = function(n, p){
    Omega = array(0,c(p,p))
    for (i in 1:(p-1)){
      Omega[i,i+1] = 0.8
    }
    for (i in 1:(p-2)){
      Omega[i,i+2] = 0.4
    }
    for (i in 1:(p-3)){
      Omega[i,i+3] = 0.2
    }
    Omega = Omega + t(Omega)
    diag(Omega) = 2
    Sig = solve(Omega, diag(p))
    return(Sig)
  }
  scenario.2014CLX.cov4 = function(n, p){
    Omega = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = 0.6^(abs(i-j))
        Omega[i,j] = theval
        Omega[j,i] = theval
      }
    }
    diag(Omega) = 1
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Sig   = Dhalf%*%solve(Omega,Dhalf)
    return(Sig)
  }
  scenario.2014CLX.cov5 = function(n, p){
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Ohalf = array(0,c(p,p))
    for (k in 1:floor(p/2)){
      lowk = as.integer(2*(k-1)+1)
      uppk = as.integer(2*k)
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          Ohalf[i,j] = 0.8
        }
      }
    }
    diag(Ohalf) = 1

    Omega = Dhalf%*%Ohalf%*%Ohalf%*%Dhalf
    Sig   = solve(Omega, diag(p))
  }
  scenario.2014CLX.cov6 = function(n, p){
    SigStar = array(0,c(p,p))
    for (k in 1:floor(p/2)){
      lowk = as.integer(2*(k-1)+1)
      uppk = as.integer(2*k)
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          SigStar[i,j] = 0.8
        }
      }
    }
    diag(SigStar) = 1
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))

    E = array(0,c(p,p))
    noff = sum(1:(p-1))
    E[lower.tri(E)] = rbinom(noff, 1, 0.3)*runif(noff, min=-0.2, max=0.2)
    E = E + t(E)

    triplet = Dhalf%*%SigStar%*%Dhalf
    delta   = abs(min(eigen(triplet + E)$values))

    Sig = triplet + E + delta*diag(p)
    return(Sig)
  }
  scenario.2014CLX.cov7 = function(n, p){
    SigStar = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = ((abs(i-j))^{-5})/2
        SigStar[i,j] = theval
        SigStar[j,i] = theval
      }
    }
    diag(SigStar) = 1
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Sig   = Dhalf%*%SigStar%*%Dhalf
    return(Sig)
  }
  scenario.2014CLX.cov8 = function(n, p){
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Fmat  = diag(p)
    for (i in 1:(p-1)){
      j = i+1
      Fmat[i,j] = 0.5
      Fmat[j,i] = 0.5
    }
    UmatOriginal = matrix(rnorm(p*p), nrow=p)
    UmatOriginal = UmatOriginal%*%t(UmatOriginal)
    U3 = qr.Q(qr(UmatOriginal))[,1:3]

    Sig = Dhalf%*%(Fmat + (U3%*%t(U3)))%*%Dhalf
    return(Sig)
  }

  # 1. generate Ha's mean
  no.mean = as.integer(no.mean)
  if (!(no.mean %in% c(1,2,3,4))){
    stop("* scenario.2014CLX : 'no.mean' is from 1 to 4.")
  }
  mu.alpha = switch(no.mean,
                    "1" = scenario.2014CLX.mean1(n,p),
                    "2" = scenario.2014CLX.mean2(n,p),
                    "3" = scenario.2014CLX.mean3(n,p),
                    "4" = scenario.2014CLX.mean4(n,p))
  # 2. generate common covariance
  no.cov = as.integer(no.cov)
  if (!(no.cov %in% 1:8)){
    stop("* scenario.2014CLX : 'no.cov' is from 1 to 8.")
  }
  sig.alpha = switch(no.cov,
                     "1" = scenario.2014CLX.cov1(n,p),
                     "2" = scenario.2014CLX.cov2(n,p),
                     "3" = scenario.2014CLX.cov3(n,p),
                     "4" = scenario.2014CLX.cov4(n,p),
                     "5" = scenario.2014CLX.cov5(n,p),
                     "6" = scenario.2014CLX.cov6(n,p),
                     "7" = scenario.2014CLX.cov7(n,p),
                     "8" = scenario.2014CLX.cov8(n,p))
  # generation and return

  htype=match.arg(htype)
  output = list()
  if (htype=="H0"){
    output$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig.alpha)
    output$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig.alpha)
  } else if (htype=="H1"){
    output$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig.alpha)
    output$Y = mvtnorm::rmvnorm(n, mean=mu.alpha, sigma=sig.alpha)
  }
  return(output)
}

# (02) scenario.2013CLX ---------------------------------------------------
#' @keywords internal
#' @noRd
scenario.2013CLX <- function(n, p, no.cov=1, htype=c("H0","H1")){
  # covariance modeling 1
  scenario.2013CLX.cov1 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    SStar = array(0,c(p,p))
    for (k in 1:floor(p/5)){
      lowk = 5*(k-1)+1
      uppk = 5*k
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          SStar[i,j] = 0.5
        }
      }
    }
    diag(SStar) = 1.0
    mSig = Dhalf%*%SStar%*%Dhalf



    return(mSig)
  }
  scenario.2013CLX.cov2 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        sstar[i,j] = (0.5^abs(i-j))
      }
    }
    msig = Dhalf%*%sstar%*%Dhalf
    return(msig)
  }
  scenario.2013CLX.cov3 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = ifelse(rbinom(1, 1, 0.05) > 0.5, 0.5, 0)
        sstar[i,j] = theval
        sstar[j,i] = theval
      }
    }
    diag(sstar) = 1

    del  = abs(min(eigen(sstar)$values))+0.05
    ssig = (Dhalf%*%(sstar+del*diag(p))%*%Dhalf)/(1+del)
    return(ssig)
  }
  scenario.2013CLX.cov4 <- function(p){
    omat = diag(runif(p, min=1, max=5))
    tmat = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        tmat[i,j] = ((-1)^(i+j))*(0.4^((abs(i-j))^0.1))
      }
    }
    ssig = omat%*%tmat%*%omat
    return(ssig)
  }
  # 1. generate common covariance
  no.cov = as.integer(no.cov)
  if (!(no.cov %in% 1:4)){
    stop("* scenario.2013CLX : 'no.cov' is from 1 to 4.")
  }
  thesig = switch(no.cov,
                  "1" = scenario.2013CLX.cov1(p),
                  "2" = scenario.2013CLX.cov2(p),
                  "3" = scenario.2013CLX.cov3(p),
                  "4" = scenario.2013CLX.cov4(p))
  # 2. generate samples accordingly
  output = list()
  htype = match.arg(htype)
  if (htype=="H0"){
    output$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=thesig)
    output$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=thesig)
  } else {
    n.upptri = sum(1:(p-1))
    vec.utri = rep(0,n.upptri)
    vec.utri[sample(1:n.upptri, 4, replace=FALSE)] = max(diag(thesig))*runif(4, min=0, max=4)
    Umat = array(0,c(p,p))
    Umat[upper.tri(Umat)] = vec.utri
    Umat = (Umat + t(Umat))

    del  = abs(min(c(min(eigen(thesig+Umat)$values), min(eigen(thesig)$values))))+0.05
    sig1 = thesig + del*diag(p)
    sig2 = thesig + Umat + del*diag(p)

    output$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
    output$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig2)
  }
  return(output)
}

# (03) gen.cov.2014CLX ----------------------------------------------------
#' @keywords internal
#' @noRd
gen.cov.2014CLX <- function(p, no.cov=1){
  scenario.2014CLX.cov1 = function(n, p){
    Sig = array(0,c(p,p))
    for (k in 1:floor(p/2)){
      lowk = as.integer(2*(k-1)+1)
      uppk = as.integer(2*k)
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          Sig[i,j] = 0.8
        }
      }
    }
    diag(Sig) = 1
    return(Sig)
  }
  scenario.2014CLX.cov2 = function(n, p){
    Sig = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        Sig[i,j] = 0.6^(abs(i-j))
      }
    }
    return(Sig)
  }
  scenario.2014CLX.cov3 = function(n, p){
    Omega = array(0,c(p,p))
    for (i in 1:(p-1)){
      Omega[i,i+1] = 0.8
    }
    for (i in 1:(p-2)){
      Omega[i,i+2] = 0.4
    }
    for (i in 1:(p-3)){
      Omega[i,i+3] = 0.2
    }
    Omega = Omega + t(Omega)
    diag(Omega) = 2
    Sig = solve(Omega, diag(p))
    return(Sig)
  }
  scenario.2014CLX.cov4 = function(n, p){
    Omega = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = 0.6^(abs(i-j))
        Omega[i,j] = theval
        Omega[j,i] = theval
      }
    }
    diag(Omega) = 1
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Sig   = Dhalf%*%solve(Omega,Dhalf)
    return(Sig)
  }
  scenario.2014CLX.cov5 = function(n, p){
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Ohalf = array(0,c(p,p))
    for (k in 1:floor(p/2)){
      lowk = as.integer(2*(k-1)+1)
      uppk = as.integer(2*k)
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          Ohalf[i,j] = 0.8
        }
      }
    }
    diag(Ohalf) = 1

    Omega = Dhalf%*%Ohalf%*%Ohalf%*%Dhalf
    Sig   = solve(Omega, diag(p))
  }
  scenario.2014CLX.cov6 = function(n, p){
    SigStar = array(0,c(p,p))
    for (k in 1:floor(p/2)){
      lowk = as.integer(2*(k-1)+1)
      uppk = as.integer(2*k)
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          SigStar[i,j] = 0.8
        }
      }
    }
    diag(SigStar) = 1
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))

    E = array(0,c(p,p))
    noff = sum(1:(p-1))
    E[lower.tri(E)] = rbinom(noff, 1, 0.3)*runif(noff, min=-0.2, max=0.2)
    E = E + t(E)

    triplet = Dhalf%*%SigStar%*%Dhalf
    delta   = abs(min(eigen(triplet + E)$values))

    Sig = triplet + E + delta*diag(p)
    return(Sig)
  }
  scenario.2014CLX.cov7 = function(n, p){
    SigStar = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = ((abs(i-j))^{-5})/2
        SigStar[i,j] = theval
        SigStar[j,i] = theval
      }
    }
    diag(SigStar) = 1
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Sig   = Dhalf%*%SigStar%*%Dhalf
    return(Sig)
  }
  scenario.2014CLX.cov8 = function(n, p){
    Dhalf = diag(sqrt(runif(p, min=1, max=3)))
    Fmat  = diag(p)
    for (i in 1:(p-1)){
      j = i+1
      Fmat[i,j] = 0.5
      Fmat[j,i] = 0.5
    }
    UmatOriginal = matrix(rnorm(p*p), nrow=p)
    UmatOriginal = UmatOriginal%*%t(UmatOriginal)
    U3 = qr.Q(qr(UmatOriginal))[,1:3]

    Sig = Dhalf%*%(Fmat + (U3%*%t(U3)))%*%Dhalf
    return(Sig)
  }
  sig.alpha = switch(no.cov,
                     "1" = scenario.2014CLX.cov1(n,p),
                     "2" = scenario.2014CLX.cov2(n,p),
                     "3" = scenario.2014CLX.cov3(n,p),
                     "4" = scenario.2014CLX.cov4(n,p),
                     "5" = scenario.2014CLX.cov5(n,p),
                     "6" = scenario.2014CLX.cov6(n,p),
                     "7" = scenario.2014CLX.cov7(n,p),
                     "8" = scenario.2014CLX.cov8(n,p))
  return(sig.alpha)
}

# (04) gen.2013CLX.mine ---------------------------------------------------
#' @keywords internal
#' @noRd
gen.2013CLX.mine <- function(n, p, no.cov=1, m, rho){
  # m   : number of non-zero elements in U
  # rho : non-zero signal strength

  # covariance modeling 1
  scenario.2013CLX.cov1 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    SStar = array(0,c(p,p))
    for (k in 1:floor(p/5)){
      lowk = 5*(k-1)+1
      uppk = 5*k
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          SStar[i,j] = 0.5
        }
      }
    }
    diag(SStar) = 1.0
    mSig = Dhalf%*%SStar%*%Dhalf



    return(mSig)
  }
  scenario.2013CLX.cov2 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        sstar[i,j] = (0.5^abs(i-j))
      }
    }
    msig = Dhalf%*%sstar%*%Dhalf
    return(msig)
  }
  scenario.2013CLX.cov3 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = ifelse(rbinom(1, 1, 0.05) > 0.5, 0.5, 0)
        sstar[i,j] = theval
        sstar[j,i] = theval
      }
    }
    diag(sstar) = 1

    del  = abs(min(eigen(sstar)$values))+0.05
    ssig = (Dhalf%*%(sstar+del*diag(p))%*%Dhalf)/(1+del)
    return(ssig)
  }
  scenario.2013CLX.cov4 <- function(p){
    omat = diag(runif(p, min=1, max=5))
    tmat = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        tmat[i,j] = ((-1)^(i+j))*(0.4^((abs(i-j))^0.1))
      }
    }
    ssig = omat%*%tmat%*%omat
    return(ssig)
  }
  # 1. generate common covariance
  no.cov = as.integer(no.cov)
  if (!(no.cov %in% 1:4)){
    stop("* scenario.2013CLX : 'no.cov' is from 1 to 4.")
  }
  thesig = switch(no.cov,
                  "1" = scenario.2013CLX.cov1(p),
                  "2" = scenario.2013CLX.cov2(p),
                  "3" = scenario.2013CLX.cov3(p),
                  "4" = scenario.2013CLX.cov4(p))
  # 2. generate samples accordingly
  output = list()
  #   2-1. Let's adjust for the signal strength
  n.upptri = sum(1:(p-1))
  vec.utri = rep(0,n.upptri)
  vec.utri[sample(1:n.upptri, m, replace=FALSE)] = rho # fixed strength
  Umat = array(0,c(p,p))
  Umat[upper.tri(Umat)] = vec.utri
  Umat = (Umat + t(Umat))

  del  = abs(min(c(min(eigen(thesig+Umat)$values), min(eigen(thesig)$values))))+0.05
  sig1 = thesig + del*diag(p)
  sig2 = thesig + Umat + del*diag(p)
  #   2-2. generate data for H0
  tmpH0 = list()
  tmpH0$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  tmpH0$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  #   2-3. generate data for H1
  tmpH1 = list()
  tmpH1$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  tmpH1$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig2)
  #   2-4. wrap-up and return
  output$H0 = tmpH0
  output$H1 = tmpH1
  return(output)
}

# (05) gen.2013CLX.dense --------------------------------------------------
#' @keywords internal
#' @noRd
gen.2013CLX.dense <- function(n, p, no.cov=1, rho){
  # rho : non-zero signal strength

  # covariance modeling 1
  scenario.2013CLX.cov1 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    SStar = array(0,c(p,p))
    for (k in 1:floor(p/5)){
      lowk = 5*(k-1)+1
      uppk = 5*k
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          SStar[i,j] = 0.5
        }
      }
    }
    diag(SStar) = 1.0
    mSig = Dhalf%*%SStar%*%Dhalf



    return(mSig)
  }
  scenario.2013CLX.cov2 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        sstar[i,j] = (0.5^abs(i-j))
      }
    }
    msig = Dhalf%*%sstar%*%Dhalf
    return(msig)
  }
  scenario.2013CLX.cov3 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = ifelse(rbinom(1, 1, 0.05) > 0.5, 0.5, 0)
        sstar[i,j] = theval
        sstar[j,i] = theval
      }
    }
    diag(sstar) = 1

    del  = abs(min(eigen(sstar)$values))+0.05
    ssig = (Dhalf%*%(sstar+del*diag(p))%*%Dhalf)/(1+del)
    return(ssig)
  }
  scenario.2013CLX.cov4 <- function(p){
    omat = diag(runif(p, min=1, max=5))
    tmat = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        tmat[i,j] = ((-1)^(i+j))*(0.4^((abs(i-j))^0.1))
      }
    }
    ssig = omat%*%tmat%*%omat
    return(ssig)
  }
  # 1. generate common covariance
  no.cov = as.integer(no.cov)
  if (!(no.cov %in% 1:4)){
    stop("* scenario.2013CLX : 'no.cov' is from 1 to 4.")
  }
  thesig = switch(no.cov,
                  "1" = scenario.2013CLX.cov1(p),
                  "2" = scenario.2013CLX.cov2(p),
                  "3" = scenario.2013CLX.cov3(p),
                  "4" = scenario.2013CLX.cov4(p))
  # 2. generate samples accordingly
  output = list()
  #   2-1. Let's adjust for the signal strength by Kyoungjae's method
  vvec = rnorm(p)
  vvec = vvec/sqrt(sum(vvec^2))
  sig1 = thesig
  sig2 = thesig + outer(vvec,vvec)*rho
  #   2-2. generate data for H0
  tmpH0 = list()
  tmpH0$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  tmpH0$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  #   2-3. generate data for H1
  tmpH1 = list()
  tmpH1$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  tmpH1$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig2)
  #   2-4. wrap-up and return
  output$H0 = tmpH0
  output$H1 = tmpH1
  return(output)
}
#' @keywords internal
#' @noRd
gen.2013CLX.dense.ones <- function(n, p, no.cov=1, rho){
  # rho : non-zero signal strength

  # covariance modeling 1
  scenario.2013CLX.cov1 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    SStar = array(0,c(p,p))
    for (k in 1:floor(p/5)){
      lowk = 5*(k-1)+1
      uppk = 5*k
      for (i in lowk:uppk){
        for (j in lowk:uppk){
          SStar[i,j] = 0.5
        }
      }
    }
    diag(SStar) = 1.0
    mSig = Dhalf%*%SStar%*%Dhalf



    return(mSig)
  }
  scenario.2013CLX.cov2 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        sstar[i,j] = (0.5^abs(i-j))
      }
    }
    msig = Dhalf%*%sstar%*%Dhalf
    return(msig)
  }
  scenario.2013CLX.cov3 <- function(p){
    Dhalf = diag(sqrt(runif(p, min=0.5, max=2.5)))
    sstar = array(0,c(p,p))
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        theval = ifelse(rbinom(1, 1, 0.05) > 0.5, 0.5, 0)
        sstar[i,j] = theval
        sstar[j,i] = theval
      }
    }
    diag(sstar) = 1

    del  = abs(min(eigen(sstar)$values))+0.05
    ssig = (Dhalf%*%(sstar+del*diag(p))%*%Dhalf)/(1+del)
    return(ssig)
  }
  scenario.2013CLX.cov4 <- function(p){
    omat = diag(runif(p, min=1, max=5))
    tmat = array(0,c(p,p))
    for (i in 1:p){
      for (j in 1:p){
        tmat[i,j] = ((-1)^(i+j))*(0.4^((abs(i-j))^0.1))
      }
    }
    ssig = omat%*%tmat%*%omat
    return(ssig)
  }
  # 1. generate common covariance
  no.cov = as.integer(no.cov)
  if (!(no.cov %in% 1:4)){
    stop("* scenario.2013CLX : 'no.cov' is from 1 to 4.")
  }
  thesig = switch(no.cov,
                  "1" = scenario.2013CLX.cov1(p),
                  "2" = scenario.2013CLX.cov2(p),
                  "3" = scenario.2013CLX.cov3(p),
                  "4" = scenario.2013CLX.cov4(p))
  # 2. generate samples accordingly
  output = list()
  #   2-1. Let's adjust for the signal strength by Kyoungjae's method
  vvec = rep(1,p)
  sig1 = thesig
  sig2 = thesig + outer(vvec,vvec)*rho
  #   2-2. generate data for H0
  tmpH0 = list()
  tmpH0$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  tmpH0$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  #   2-3. generate data for H1
  tmpH1 = list()
  tmpH1$X = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig1)
  tmpH1$Y = mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sig2)
  #   2-4. wrap-up and return
  output$H0 = tmpH0
  output$H1 = tmpH1
  return(output)
}
