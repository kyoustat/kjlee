#' Sparse Covariance Difference Test : covid [1,4]
#'
#'  estA = Schott k-sample    (2007)
#'  estB = Li and Chen        (2012)
#'  estC = Cai, Liu, and Xia  (2013)
#'  estX = mxPBF
#'
#'  @examples
#'  ## run
#'  out = sp_covs(p=5, n=25, mreplicate=100, covid=1, ndev=1)
#'
#' @export
sp_covs <- function(p, n, mreplicate, covid, ndev=1, rhos=c(seq(0.01, 2, length.out=20), seq(2.2, 5, length.out=10))){
  ## parameters
  myn    = round(n)
  myp    = round(p)
  myrhos = rhos
  nrhos  = length(myrhos)
  aucA = list()
  aucB = list()
  aucC = list()
  aucX = list() # X stands for mxPBF

  ## main computation
  for (i in 1:nrhos){ # iterate over 'nrhos' values
    # 1. data generation ----------------------------------- ********** VERY IMPORTANT
    tgtrho = myrhos[i]
    true.label = c(rep(0,mreplicate), rep(1,mreplicate))
    dat.H0 = list()
    dat.H1 = list()

    for (k in 1:mreplicate){
      gen.mine = gen.2013CLX.mine(myn, myp, no.cov=covid, m=ndev, rho=tgtrho)
      dat.H0[[k]] = gen.mine$H0
      dat.H1[[k]] = gen.mine$H1
    }

    # 2. apply tests
    estA = rep(0,2*mreplicate) # label recorder
    estB = rep(0,2*mreplicate)
    estC = rep(0,2*mreplicate)
    estX = rep(0,2*mreplicate)
    for (k in 1:(2*mreplicate)){
      if (k<=mreplicate){
        tgtdata = dat.H0[[k]]
      } else {
        tgtdata = dat.H1[[k-mreplicate]]
      }
      tgtX = tgtdata$X
      tgtY = tgtdata$Y

      estA[k] = SHT::useknd(tgtX, tgtY, "covk.2007Schott")$p.value
      estB[k] = SHT::cov2.2012LC(tgtX, tgtY)$p.value
      estC[k] = SHT::cov2.2013CLX(tgtX, tgtY)$p.value
      estX[k] = max(mxPBF::testcov2(tgtX,tgtY)$log.BF.mat)
    }
    #   2-3. compute ROC curve and AUC value
    ROCA = ROC.pvalue(true.label, estA, nsep=1000)
    ROCB = ROC.pvalue(true.label, estB, nsep=1000)
    ROCC = ROC.pvalue(true.label, estC, nsep=1000)
    ROCX = ROC.BF(true.label, estX, nsep=1000)

    aucA[[i]] = cbind(ROCA$FPR, ROCA$TPR)
    aucB[[i]] = cbind(ROCB$FPR, ROCB$TPR)
    aucC[[i]] = cbind(ROCC$FPR, ROCC$TPR)
    aucX[[i]] = cbind(ROCX$FPR, ROCX$TPR)
  }

  ## report output
  output = list(rhos=myrhos, A=aucA, B=aucB, C=aucC, X=aucX)
  return(output)
}
