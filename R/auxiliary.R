# auxiliary functions -----------------------------------------------------
# (01) gen_table  : generating table
# (02) ROC.pvalue : TPR and FPR for ROC curve
# (03) ROC.BF     : for BF, either vector or matrix is fine
# (04) simple_auc : compute AUC value
# (05) gen.data   : two mean data generation

# (01) gen_table ----------------------------------------------------------
#' @keywords internal
#' @noRd
gen_table = function(true.label, est.label){
  output = array(0,c(2,2))
  output[1,1] = sum((true.label<0.5)*(est.label<0.5))
  output[1,2] = sum((true.label>0.5)*(est.label<0.5))
  output[2,1] = sum((true.label<0.5)*(est.label>0.5))
  output[2,2] = sum((true.label>0.5)*(est.label>0.5))
  return(output)
}

# (02) ROC.pvalue ---------------------------------------------------------
#' @keywords internal
#' @noRd
ROC.pvalue <- function(true.label, p.values, nsep=100){
  if (length(true.label)!=length(p.values)){
    stop("* ROC.pvalue : something wrong.")
  }
  nlabel  = length(true.label)
  thrvec  = seq(from=min(p.values)*1.000001, to=0.999999*max(p.values), length.out=nsep)

  vec.tpr = rep(0, nsep)
  vec.fpr = rep(0, nsep)
  for (i in 1:nsep){
    thrtgt    = thrvec[i]                              # target significance
    est.label = ifelse(p.values < thrtgt, 1, 0)        # if rejected, label as 1

    conf.table = gen_table(true.label, est.label)
    vec.fpr[i] = as.double(conf.table[2,1]/sum(conf.table[,1]))
    vec.tpr[i] = as.double(conf.table[2,2]/sum(conf.table[,2]))
  }
  output = list()
  output$TPR = vec.tpr
  output$FPR = vec.fpr
  return(output)
}


# (03) ROC.BF -------------------------------------------------------------
#' @keywords internal
#' @noRd
ROC.BF <- function(true.label, BF, nsep=100){
  if (length(true.label)!=length(BF)){
    stop("* ROC.BF : something wrong.")
  }

  minBF = min(BF)
  maxBF = max(BF)
  nlabel  = length(true.label)
  if (any(BF < 0)){
    minrange = 1.000001*minBF
  } else {
    minrange = 0.999999*minBF
  }
  if (all(BF <0)){
    maxrange = 0.999999*maxBF
  } else {
    maxrange = 1.000001*maxBF
  }
  thrvec  = seq(from=minrange, to=maxrange, length.out=nsep)

  vec.tpr = rep(0, nsep)
  vec.fpr = rep(0, nsep)
  for (i in 1:nsep){
    thrtgt    = thrvec[i]                              # target threshold
    est.label = ifelse(BF > thrtgt, 1, 0)           # if rejected, label as 1

    conf.table = gen_table(true.label, est.label)
    vec.fpr[i] = as.double(conf.table[2,1]/sum(conf.table[,1]))
    vec.tpr[i] = as.double(conf.table[2,2]/sum(conf.table[,2]))
  }
  output = list()
  output$TPR = vec.tpr
  output$FPR = vec.fpr
  return(output)
}


# (04) simple_auc ---------------------------------------------------------
#' @keywords internal
#' @noRd
simple_auc <- function(FPR, TPR){
  # inputs already sorted, best scores first
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  return(sum(TPR * dFPR) + sum(dTPR * dFPR)/2)
}

# (05) gen.data -----------------------------------------------------------
#' @keywords internal
#' @noRd
gen.data <- function(n, mu0, mu1, sig, htype=c("H0","H1")){
  # 1-1. check the dimensionality
  p = length(mu0)
  if (nrow(sig)!=p){
    stop("gen.data : param dimension unmatch.")
  }
  # 1-2. real generation
  htype=match.arg(htype)
  output = list()
  if (htype=="H0"){
    output$X = mvtnorm::rmvnorm(n, mean=mu0, sigma=sig)
    output$Y = mvtnorm::rmvnorm(n, mean=mu0, sigma=sig)
  } else if (htype=="H1"){
    output$X = mvtnorm::rmvnorm(n, mean=mu0, sigma=sig)
    output$Y = mvtnorm::rmvnorm(n, mean=mu1, sigma=sig)
  }
  return(output)
}
