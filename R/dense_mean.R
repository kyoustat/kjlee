#' Dense Mean Difference Test : covid [1,8]
#'
#'  estA = Bai and Saranadasa (1996)
#'  estB = Srivastava and Du  (2008)
#'  estC = Cai, Liu, and Xia  (2014)
#'  estX = mxPBF
#'
#'
#' @examples
#' # run
#' out = dense_mean(p=5, n=25, mreplicate=100, covid=1)
#'
#' @export
dense_mean <- function(p, n, mreplicate, covid, rhos=c(seq(0.01, 2, length.out=20), seq(2.2, 5, length.out=10))){
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
    # 1. data generation
    tgtrho = myrhos[i]
    true.label = c(rep(0,mreplicate), rep(1,mreplicate))
    dat.H0 = list()
    dat.H1 = list()

    my.mu0 = rep(0,myp)                            # mu0
    my.mu1 = rnorm(myp)                            # mu1
    my.mu1 = (my.mu1/sqrt(sum(my.mu1^2)))*tgtrho
    my.sig = gen.cov.2014CLX(myp, no.cov=covid)    # common covariance

    for (k in 1:mreplicate){
      dat.H0[[k]] = gen.data(myn, my.mu0, my.mu1, my.sig, htype="H0")
    }
    for (k in 1:mreplicate){
      dat.H1[[k]] = gen.data(myn, my.mu0, my.mu1, my.sig, htype="H1")
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

      estA[k] = highmean::apval_Bai1996(tgtX, tgtY)$pval
      estB[k] = SHT::mean2.2008SD(tgtX, tgtY)$p.value
      estC[k] = highmean::apval_Cai2014(tgtX, tgtY)$pval
      estX[k] = max((mxPBF::testmean2(tgtX, tgtY)$log.BF.vec))
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
