readf.binned = function(fin,res) {
    dfs    = read.delim(fin, header=F, sep='\t')
    #    nout   = round(max(dfs[,1])/res)+1
    #    bmatch = dfs[,1]/res
    #    bmatch.int =findInterval(bmatch, 1:nout, all.inside=T)
    # note that one change 1 -> 0 from Joshs code
    nout   = round(max(dfs[,1])/res)+1
    bmatch = dfs[,1]/res
    bmatch.int =findInterval(bmatch, 0:nout, all.inside=T)

    means  = rep(0, nout)
    means[bmatch.int]=dfs[,2]                # np  (i.e. 1000 * .5)
   
    counts = rep(0, nout)
    counts[bmatch.int]= dfs[,2] + dfs[,3]      # n   (i.e. 1000)
   
    p =rep(0,nout)
    p[bmatch.int] = dfs[,2]/(dfs[,2]+dfs[,3])   # p   (i.e. .5)
   
    variances =  counts * p *( 1-p)  # np(1-p)  (i.e  1000*.5(1-.5)
    variances[variances==0]=NA
    # equal to y, y_var, and d 
    # means = y
    # variances = y_var 
    # counts = d 
    return(list(means=means, variances=variances, counts=counts))
}

# same as above, but instead of reading a file, it wants to be fed the counts and positions
# first column: SNP position on the chromosome
# 2nd column: ref counts
# 3rd column: alt counts
feedf.binned = function(finObject, res) {
    dfs    = finObject
    #    nout   = round(max(dfs[,1])/res)+1
    #    bmatch = dfs[,1]/res
    #    bmatch.int =findInterval(bmatch, 1:nout, all.inside=T)
    # note that one change 1 -> 0 from Joshs code
    nout   = round(max(dfs[,1])/res)+1
    bmatch = dfs[,1]/res
    bmatch.int =findInterval(bmatch, 0:nout, all.inside=T)
    
    means  = rep(0, nout)
    means[bmatch.int]=dfs[,2]                # np  (i.e. 1000 * .5)
    
    counts = rep(0, nout)
    counts[bmatch.int]= dfs[,2] + dfs[,3]      # n   (i.e. 1000)
    
    p =rep(0,nout)
    p[bmatch.int] = dfs[,2]/(dfs[,2]+dfs[,3])   # p   (i.e. .5)
    
    variances =  counts * p *( 1-p)  # np(1-p)  (i.e  1000*.5(1-.5)
    variances[variances==0]=NA
    # equal to y, y_var, and d
    # means = y
    # variances = y_var
    # counts = d
    return(list(means=means, variances=variances, counts=counts))
}

kalman=function(y, y_var, d, TT, N, p) {
    mu = rep(0,TT)#numpy.zeros(T)
    V = rep(0,TT)#numpy.zeros(T)
    P = rep(0,TT)# numpy.zeros(T)

    V_pstr = rep(0,TT)# numpy.zeros(T)
    mu_pstr = rep(0,TT) #numpy.zeros(T)

    cc = rep(1,TT)

    mu_initial = 0.5*N # Initial parameters, assumed given (binomial distribution)
    V_initial = 0.25*N

    A = (1.0 - 2.0*p)
    C = 1.0 * d / N
    S = p*(1.0-p)*N
    
    explode=C[1]/(C[1]^2.0*V_initial + ifelse(is.na(y_var[1]),0, y_var[1]))
    
    K = V_initial* ifelse(is.na(explode),0, explode)
    # C[1]/(C[1]^2.0*V_initial + ifelse(is.na(y_var[1]),0, y_var[1]))
    mu[1] = mu_initial + K*(y[1] - C[1]*mu_initial)
    V[1] = (1.0-K*C[1])*V_initial
    if(is.na(y_var[1])) { cc[1]=1 } else {
        cc[1]= dnorm(y[1], C[1]*mu_initial, sqrt(C[1]^2.0*V_initial + y_var[1])) }
    
    # Forward pass:
    for (i in 2:(TT)) {
        if (i == 2) { 
            P[i-1] = A^2.0*V_initial + S }
        else { 
            P[i-1] = A^2.0*V[i-1] + S 
        }

        if ( is.na(y_var[i] )) {
            K = 0
            cc[i] = 1.0
         }
        else{
            K = P[i-1]*C[i]/(C[i]^2.0*P[i-1]+y_var[i])
            cc[i] = dnorm(y[i], C[i]*(A*mu[i-1]+p*N), sqrt(C[i]^2.0*P[i-1] + y_var[i]))
        }
        mu[i] = A * mu[i-1] + N*p + K * (y[i] - C[i]*(A*mu[i-1] + N*p))
        V[i] = (1.0-K*C[i])*P[i-1]
     }

    V_pstr[length(V_pstr)] = V[length(V)]
    mu_pstr[length(V_pstr)] = mu[length(V)]
    logLik = sum(log(cc))

    # Backwards pass:
    for (i in (TT-1):1) {
        J = V[i]*A/P[i]
        mu_pstr[i] = mu[i] + J * (mu_pstr[i+1] - A*(mu[i]) - N*p)
        V_pstr[i] = V[i] + J^2.0 * (V_pstr[i+1] - P[i])
    }
    
    return( list(mu_pstr=mu_pstr, V_pstr=V_pstr, logLik=logLik))

}

#k1=kalman(y,y_var,d, TT, N,p)
#k2=kalman(y2,y_var2,d2, TT, N,p)
  

lognormpdf=function(x, mu, sigma){-0.5*log(2*pi) - log(sigma) + (-(x-mu)^2.0/2.0/sigma^2.0) }

calcLOD = function(mu_pstr_vec, v_pstr_vec, TT, N) {
   # mu_pstr_vec =  k1$mu_pstr
   # v_pstr_vec  = k1$V_pstr
    LOD = rep(0,TT) #numpy.zeros(T)
    mu_MLE = rep(0,TT) #numpy.zeros(T)

    mu_initial = 0.5*N
    V_initial = 0.25*N
    delta = 0.0025
    x=seq(delta, 1-delta+delta/2, delta)
    p_alt=x

    #row-wise
    #log_p_precomp=(sapply(x, function(x) dnorm(N*x, N*p_alt, sqrt(p_alt*(1.0-p_alt)*N),log=T)))
    p_precomp=(sapply(x, function(x) dnorm(N*x, N*p_alt, sqrt(p_alt*(1.0-p_alt)*N))))
    logreweighter = lognormpdf(N*x, mu_initial, sqrt(V_initial))

    if(is.null(dim(mu_pstr_vec))) {
        for (i in 1:TT) {
            logallsums = rep(0, length(x))
            logtemp = lognormpdf(N*x, mu_pstr_vec[i], sqrt(v_pstr_vec[i])) - logreweighter
            scaler= max(logtemp)
            # for the logic of this step see: 
            # http://www.johndcook.com/blog/2012/07/26/avoiding-underflow-in-bayesian-computations/
            logallsums = logallsums + scaler+ log(p_precomp %*% exp(logtemp - scaler))
            #logallsums = logallsums + scaler+ log(jbMult(log_p_precomp,logtemp-scaler  )) 
            p_alt=x[which.max(logallsums)]*N
            mu_MLE[i] = p_alt
            LOD[i] = log10(N) + log10(x[2]-x[1])+ max(logallsums) / log(10.0)
        }
    } else {
        for (i in 1:TT) {
            logallsums = rep(0, length(x))
            logtemp = lognormpdf(N*x, mu_pstr_vec[i,1], sqrt(v_pstr_vec[i,1])) - logreweighter
            scaler  = max(logtemp)
            toadd1= scaler+log(p_precomp %*% exp(logtemp - scaler))
            # log(jbMult(log_p_precomp,logtemp-scaler  ))
            logallsums = logallsums + toadd1
            logtemp = lognormpdf(N*x, mu_pstr_vec[i,2], sqrt(v_pstr_vec[i,2])) - logreweighter
            scaler  = max(logtemp)
            toadd2= scaler+ log(p_precomp %*% exp(logtemp - scaler))
            #log(jbMult(log_p_precomp,logtemp-scaler  ))
            logallsums = logallsums + toadd2
            # inf correction
            # see the discussion on top on why this can be problematic
            # (when the two x distributions do still overlap (small, but non-inf sum), overwriting the flanks with the higher values from one of the distributions is problematic
            # so, hack this by only doing the correction if the joint distribution is always inf?
            # this will still make problems when we go from no AF difference to a very large AF difference
            # the correction will then suddenly kick in, and pull down the LOD a lot!
            # but for such extreme cases, the LOD at the peak will presumably be much too high to be affected?
            # still, we probably just moved the problem to a different region of the parameter space
            # the jumps in correction still occur, just elsewhere
            # ideally, we would not have to deal with these INFs at all, or replace them with something very small, or something else that is guaranteed to stay consistent along the chromosome
            if (length(logallsums[is.finite(logallsums)]) == 0){
								# grey out this hack:
                #inf1.notinf2=!is.finite(toadd1) & is.finite(toadd2)
                #inf2.notinf1=is.finite(toadd1)  & !is.finite(toadd2)
                #logallsums[inf1.notinf2]=toadd2[inf1.notinf2]
                #logallsums[inf2.notinf1]=toadd1[inf2.notinf1]
								
								# and try this instead:
								# set all of logallsums to log (not log10!) of smallest possible value
								logallsums <- rep(log(.Machine$double.xmin), length(logallsums))
								# set the point between the peaks of toadd1 and toadd2 to twice this minimum value
								logallsums[round(mean(c(which.max(toadd1), which.max(toadd2))))] <- log(.Machine$double.xmin * 2)
								# still no good - now the correction results in the LOD jumping up at the peak
								
								# OK, one more idea
								# toadd1 and toadd2 are in log space and seem to form exponentical-like looking relationships
								# can I numerically approximate these and find their intersect
								# all the while remaining in log space?
								# trying splinefun
								# method "fmm" is a cubic fit that sometimes comes bouncing back up at the fringes, not good
								# "natural" would just use a stright line beyond the ends, also not sure about that
								#f1 <- splinefun(x[is.finite(toadd1)], toadd1[is.finite(toadd1)], method="fmm")
								#f2 <- splinefun(x[is.finite(toadd2)], toadd2[is.finite(toadd2)], method="fmm")
								#fSum <- f1(x) + f2(x)
								
								# try a quadratic fit?
								# looking good for row 315
								toFit1y = toadd1[is.finite(toadd1)]
								toFit1x = x[is.finite(toadd1)]
								f1 <- lm(toFit1y ~ poly(toFit1x, 2, raw=TRUE))
								toFit2y = toadd2[is.finite(toadd2)]
								toFit2x = x[is.finite(toadd2)]
								f2 <- lm(toFit2y ~ poly(toFit2x, 2, raw=TRUE))
								fSum <- predict(f1, data.frame(toFit1x=x)) + predict(f2, data.frame(toFit2x=x))
												
								logallsums <- fSum
								
            }
						
            p_alt=x[which.max(logallsums)]*N
            mu_MLE[i] = p_alt
            LOD[i] = log10(N) + log10(x[2]-x[1])+  max(logallsums) / log(10.0)

        }
    }

    return(list(LOD=LOD, mu_MLE=mu_MLE))
}

								# debugging within the function
								# pre-fill variables
								#mu_pstr_vec = cbind(k1$mu_pstr,k2$mu_pstr)
								#v_pstr_vec = cbind(k1$V_pstr,k2$V_pstr)
								
								#plot stuff
#pdf("multiPoolExplainer_OneSiteGap.pdf", width=6, height=4)
                                #plot(x, toadd1, ylim=c(min(c(toadd1[is.finite(toadd1)], toadd2[is.finite(toadd2)])), max(c(toadd1[is.finite(toadd1)], toadd2[is.finite(toadd2)]))), type="l", xlab="Allele frequency", ylab="LOD")
#plot(x, toadd1, ylim=c(min(c(toadd1[is.finite(toadd1)], toadd2[is.finite(toadd2)], -1000)), max(c(toadd1[is.finite(toadd1)], toadd2[is.finite(toadd2)]))), type="l", xlab="Allele frequency", ylab="LOD")
                                #plot(x, toadd1, ylim=c(min(c(toadd1[is.finite(toadd1)], toadd2[is.finite(toadd2)], fSum)), max(c(toadd1[is.finite(toadd1)], toadd2[is.finite(toadd2)], fSum))), type="l")
#								points(x, toadd2, type="l", col="blue")
								#points(x, logallsums, pch=19, col="red")
#								points(x, f1(x))
#								points(x, f2(x), col="blue")
#                               points(x, predict(f1, data.frame(toFit1x=x)))
#								points(x, predict(f2, data.frame(toFit2x=x)), col="blue")
#points(x, fSum, type="l", col="green", lwd=2)
#dev.off()
								
#jbMult=function(a,b) {apply(a, 1, function(x) { sum(exp((x+b))) }) } 

# this wants to calculate 'temp', the fraction of the sum of the 10^LOD at each point
# this gets too large in my simulations!
# see my replacement hack
# note that the resulting CIs seem to be somehwat smaller than those from the python version
# is this due to the hack, that seems to set many values far away from the peak to zero?
# also note that I am reporting the mean and mode that implemented in the python version
# seems like the mean is a weighted peak estimate (eg if the peak leans to one side)?
getCIs = function(LOD, TT, confidenceArea=0.9){
    # from original MP:
    #    temp = 10^LOD / sum(10^LOD)
    # how about using this trick to address 'log-sum-exp'
    #    https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/
    #    http://machineintelligence.tumblr.com/post/4998477107/the-log-sum-exp-trick
    # (this looks like the same trick Josh found above)
    temp = 10^(LOD - (max(LOD) + log10(sum(10 ^ (LOD - max(LOD))))))
    left = NA
    right = NA
    maxLOD = max(LOD)
    maxPos = which(LOD == maxLOD)[ceiling(length(which(LOD == maxLOD)) / 2)]
    # we can replace Matts python explicit loop using cumsum in R
    credibleIndeces = which(cumsum(temp) >= 0.5 - confidenceArea/2 & cumsum(temp) <= 0.5 + confidenceArea/2)
    # this can happen for very extreme simulations: nearly all the mass is at one point
    if(max(temp) > confidenceArea){credibleIndeces = which(temp == max(temp))[ceiling(length(which(temp == max(temp))) / 2)]}
    left = credibleIndeces[1]
    right = credibleIndeces[length(credibleIndeces)]
    if (is.na(left)){left = 1}
    if (is.na(right)){right = TT}
    CImean = sum(temp * (1:length(temp) - 1))
    CImode = which(temp == max(temp))[ceiling(length(which(temp == max(temp))) / 2)]
    return(list(left=left, right=right, maxPos = maxPos, maxLOD = maxLOD, CImean=CImean, CImode=CImode))
}

getLODCIs = function(LOD, TT, confidenceDrop = 1){
    maxLOD = max(LOD)
    maxPos = which(LOD == maxLOD)[ceiling(length(which(LOD == maxLOD)) / 2)]
    left = maxPos
    right = maxPos
    thisPos = maxPos
    while (thisPos > 0 & LOD[thisPos] > maxLOD - confidenceDrop){
        if(thisPos <= 1){break()}
        thisPos = thisPos - 1
    }
    left = thisPos
    thisPos = maxPos
    while (thisPos <= TT & LOD[thisPos] > maxLOD - confidenceDrop){
        if(thisPos >= length(LOD)){break()}
        thisPos = thisPos + 1
    }
    right = thisPos
    if (is.na(left)){left = 1}
    if (is.na(right)){right = TT}
    return(list(left=left, right=right, maxPos = maxPos, maxLOD = maxLOD))
}

#L1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
#L2=calcLOD( k2$mu_pstr, k2$V_pstr, TT, N)
#L3= calcLOD( cbind(k1$mu_pstr,k2$mu_pstr), cbind(k1$V_pstr,k2$V_pstr), TT, N)
#L=(L1$LOD + L2$LOD) - L3$LOD
##L=(L1$LOD + L2$LOD)

#getLODCIs(L, TT, 1)
#getCIs(L, TT, 0.9)

# for external calling
# make this as seamless as possible with the python version
# this leads to some inelegance now!
# e.g. could feed the counts in directly instead of running through these temp files
doMultiPoolAsIfItWerePython <- function(hiFile, loFile, N, res, cM, sFile = "multipoolSummaryOutDefault.txt", LODFile = "multipoolOutDefault.txt"){

    p= res/100.0/cM
    
    hin = readf.binned(hiFile, res)
    lin = readf.binned(loFile, res)

    y=lin$means
    y_var=lin$variances
    d=lin$counts
    
    y2=hin$means
    y_var2=hin$variances
    d2=hin$counts
    
    TT=length(hin$means)

    k1=kalman(y,y_var,d, TT, N,p)
    k2=kalman(y2,y_var2,d2, TT, N,p)

    L1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
    L2=calcLOD( k2$mu_pstr, k2$V_pstr, TT, N)
    L3= calcLOD( cbind(k1$mu_pstr,k2$mu_pstr), cbind(k1$V_pstr,k2$V_pstr), TT, N)
    L=(L1$LOD + L2$LOD) - L3$LOD

    LODCIs = getLODCIs(L, TT, 1)
    CIs = getCIs(L, TT, 0.9)

# make files that mimic the python output
# note that we need to shift all the indeces down by 1
# this is because R vectors are 1 based (such that the first bin would be 100)
# but the chromosome starts at zero!
# in python, this is not an issue b/c there the vectors are zero based
    write.table(t(c(LODCIs$maxLOD, (LODCIs$maxPos-1)*res, (LODCIs$left-1)*res, (LODCIs$right-1)*res, (CIs$left-1)*res, (CIs$right-1)*res, (CIs$CImean-1)*res, (CIs$CImode-1)*res)), file=sFile, quote = FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    write.table(format(cbind((1:length(L) - 1) * res, round(L, 2)), scientific=FALSE, trim=TRUE), file=LODFile, quote = FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

#doMultiPoolAsIfItWerePython(hiFile = 'row_3518iteration8N1e+06depth1e+06effect1HGmultiPoolIn.txt', loFile = 'row_3518iteration8N1e+06depth1e+06effect1LGmultiPoolIn.txt', N = 1e6, res=100, cM=2200)

# this takes the counts directly without writing
# and spits out the output "files" as lists
doMultiPoolFromWithinR <- function(hiCounts, loCounts, N=1000, res=100, cM=2200){
    
    p= res/100.0/cM
    hin = feedf.binned(hiCounts, res)
    lin = feedf.binned(loCounts, res)
    
    y=lin$means
    y_var=lin$variances
    d=lin$counts
    
    y2=hin$means
    y_var2=hin$variances
    d2=hin$counts
    
    TT=length(hin$means)
    
    k1=kalman(y,y_var,d, TT, N,p)
    k2=kalman(y2,y_var2,d2, TT, N,p)
    
    L1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
    L2=calcLOD( k2$mu_pstr, k2$V_pstr, TT, N)
    L3= calcLOD( cbind(k1$mu_pstr,k2$mu_pstr), cbind(k1$V_pstr,k2$V_pstr), TT, N)
    L=(L1$LOD + L2$LOD) - L3$LOD
    
    LODCIs = getLODCIs(L, TT, 1)
    CIs = getCIs(L, TT, 0.9)
    
    # make files that mimic the python output
    # note that we need to shift all the indeces down by 1
    # this is because R vectors are 1 based (such that the first bin would be 100)
    # but the chromosome starts at zero!
    # in python, this is not an issue b/c there the vectors are zero based
    peakOut <- data.frame(LODCIs$maxLOD, (LODCIs$maxPos-1)*res, (LODCIs$left-1)*res, (LODCIs$right-1)*res, (CIs$left-1)*res, (CIs$right-1)*res, (CIs$CImean-1)*res, (CIs$CImode-1)*res)
    names(peakOut) <- c("maxLOD", "maxLODPosition", "LOD_CI_left", "LOD_CI_right", "CI_left", "CI_right", "CI_mean", "CI_mode")
    LODOut <- cbind((1:length(L) - 1) * res, round(L, 2))
    return(list(peakOut, LODOut))
}



