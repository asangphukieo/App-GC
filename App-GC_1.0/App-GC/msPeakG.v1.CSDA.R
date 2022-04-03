############################################
############################################
# The R package msPeakG for Peak detection #
############################################
############################################

########################################
# Author: Seongho Kim                  #  
# Email: biostatistician.kim@gmail.com #
# Version: v.1                         #
# Date: June. 19, 2015                 #
########################################

########################################
# Author: Apiwat Sangphukieo           #  
# Email: apiwat.bif@mal.kmutt.ac.th    #
# Version: v.1.2                       #
# Comment: modify for using in         #
#         auto_GCxGC pipeline          #
########################################

######################
# Required libraries #
######################
library(msm)
library(pracma)
library(snowfall)
#install.packages(c("msm","pracma","snowfall"))

####################
# From NormalGamma #
####################

dnormgam <- function (par, x = NULL, N0 = 65536, plot = TRUE, log = FALSE, 
    tail.cor = TRUE, cor = 1e-15, mu = par[1], sigma = par[2], 
    k = par[3], theta = par[4]) 
{
    if (sum(par != c(mu, sigma, k, theta))) 
        stop("Inconsistent definition of the parameters")
    if (sigma < 0) 
        stop("negative argument 'sigma'")
    if (k < 0) 
        stop("negative argument 'k'")
    if (theta < 0) 
        stop("negative argument 'theta'")
    xout <- x
    if (!is.null(xout)) {
        Trange = diff(range(xout))
        Tmin = min(xout) - 0.1 * Trange
        Tmax = max(xout) + 0.1 * Trange
    }
    else {
        Tmin = mu - 5 * sigma
        Tmax = mu + 5 * sigma + qgamma(0.99999, shape = k, scale = theta)
    }
    A = pi * N0/(Tmax - Tmin)
    mu = mu - Tmin
    L = -A + 2 * A * (0:(N0 - 1))/N0
    a = 1 - (0+1i) * L * theta
    Cs = a^(-k)
    b = -sigma^2 * L^2/2 + (0+1i) * mu * L
    Cb = exp(b)
    Y = Cs * Cb
    T = pi * (0:(N0 - 1))/A
    fT = Re(A/N0/pi * exp((0+1i) * A * T) * fft(Y))
    fT[fT < 0] = 0
    if (tail.cor) {
        ind = which((T > mu + k * theta) & (abs(fT) < cor))
    }
    else {
        ind = NULL
    }
    if (length(ind) > 0) {
        ind = ind[1]:length(fT)
        badfT = fT[ind]
        fT[ind] = dgamma(T[ind] - mu, shape = k, scale = theta)
    }
    else badfT = NULL
    if (log) {
        fT = log(fT)
        if (length(ind) > 0) 
            fT[ind] = dgamma(T[ind] - mu, shape = k, scale = theta, 
                log = TRUE)
    }
    T = T + Tmin
    jnd = which(T < Tmax)
    T = T[jnd]
    fT = fT[jnd]
    if (!is.null(xout)) {
        fT = approx(x = T, y = fT, xout = xout, rule = 2)$y
        T = xout
    }
    if (plot) 
        plot(T, fT, type = "l")
    list(dout = fT, xout = T)
}

normgam.signal <- function (x, par, tail.cor = TRUE, cor = 1e-15, gshift = FALSE, 
    mu = par[1], sigma = par[2], k = par[3], theta = par[4]) 
{
    xout <- x
    Sh = xout
    tri = sort(xout, index.return = TRUE)
    xout2 = tri$x
    Iout = tri$ix
    D1 = dnormgam(par = c(par[1], par[2], par[3] + 1, par[4]), 
        x = xout2, tail.cor = tail.cor, cor = cor, plot = FALSE)
    D2 = dnormgam(par = par, x = xout2, tail.cor = tail.cor, 
        cor = cor, plot = FALSE)
    R = D1$dout/D2$dout
    I0 = which.min(R)
    if (I0 > 1) 
        R[1:(I0 - 1)] = approx(c(0, I0), c(0, R[I0]), xout = 1:(I0 - 
            1), rule = 2)$y
    Shat = approx(D1$xout, k * theta * R, xout = xout2, rule = 2)$y
    if ((gshift) && (k < 1)) {
        jnd = which(Shat < k * theta/10)
        B = hist(Shat[jnd], breaks = 5000, plot = FALSE)
        shift = B$mids[which.max(B$counts)]
        Shat <- Shat - shift
        Shat[Shat < 0] = 0
    }
    Sh[Iout] = Shat
    return(Sh)
}

###################################
###################################
###################################
# Normal-Gamma-Bernoulli model
# Gamma(alpha,beta)
# alpha: shape
# beta: scale
# m = alpha*beta
# v = alpha*beta^2
###################################
###################################
###################################

quadinf2 <- function (f, xa, xb, tol = 1e-12, ...) 
{
    stopifnot(is.numeric(xa), length(xa) == 1, is.numeric(xb), 
        length(xb) == 1, is.numeric(tol), length(tol) == 1)
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    if (xa == xb) {
        return(list(Q = 0, reltol = 0))
    }
    else if (xa > xb) {
        stop("For the integration limits 'xa <= xb' is required.")
    }
    u <- -1
    v <- 1
    if (is.finite(xa) && is.finite(xb)) {
        duv <- (xb - xa)/(v - u)
        ff <- function(y) duv * f(xa + duv * (y - u))
    }
    else if (is.finite(xa) && is.infinite(xb)) {
        ff <- function(y) {
            duv <- (v - u)/(1 - y)^2
            z <- (y - u)/(v - u)
            duv * f(xa + z/(1 - z))
        }
    }
    else if (is.infinite(xa) && is.finite(xb)) {
        xa <- -xb
        ff <- function(y) {
            duv = (v - u)/(1 - y)^2
            z <- (y - u)/(v - u)
            duv * f(-(xa + z/(1 - z)))
        }
    }
    else if (is.infinite(xa) && is.infinite(xb)) {
        ff <- function(y) {
            z <- pi * (y - u)/(v - u) - 0.5
            duv <- (1/cos(z)^2) * pi/(v - u)
            duv * f(tan(z))
        }
    }
    else {
        stop("Other integration domains will be treated later.")
    }
    cc <- .quadinf_pre()
    xx <- cc$nodes
    ww <- cc$weights
    h <- 1/2
    x <- xx[[1]]
    w <- ww[[1]]
    s <- w[7] * ff(x[7])
    for (j in 1:6) {
        s <- s + w[j] * (ff(x[j]) + ff(-x[j]))
    }
    Q <- s * h
    for (k in 2:7) {
        x <- xx[[k]]
        w <- ww[[k]]
        s <- 0
        for (j in 1:length(w)) {
            s <- s + w[j] * (ff(x[j]) + ff(-x[j]))
        }
        h <- h/2
        newQ <- s * h + Q/2
        delta <- abs(newQ - Q)
        Q <- newQ
        if (delta < tol) 
            break
    }
    return(list(Q = Q, relerr = delta, niter = k))
}
environment(quadinf2) = environment(quadinf)


true.test <- function(par,x){
	mu = par[1]
	s2 = par[2]
	phi = par[3]
	tp = mu+s2/phi
	td = x-tp + sqrt(s2)*dnorm((x-tp)/sqrt(s2))/pnorm((x-tp)/sqrt(s2))
	td
}


true.test2 <- function(par,x){
	n3 = exp(lne.test(x,par[1],par[2],par[3]))
	post.d <- function(val){
		n1 = dnorm(x,(par[1]+val),sqrt(par[2]),log=FALSE)	
		n2 = dexp(val,1/par[3],log=FALSE)
		val*n1*n2/n3
	}
	quadinf(post.d,0,Inf)$Q
}


lne.test2 <- function(x,mu,sigma2,phi){
	sigma = sqrt(sigma2)
	p1 = -(x-mu-.5*sigma2/phi)/phi - log(phi)
	p2 = pnorm((x-mu-sigma2/phi)/sigma,log.p=TRUE)
	p = p1+p2
	p
}

lne.test <- function(x,mu,sigma2,phi){
	sigma = sqrt(sigma2)
	ne <- function(val){
		n1 = dnorm(x,(mu+val),sigma,log=FALSE)	
		n2 = dexp(val,1/phi,log=FALSE)
		n1*n2
	}
	log(quadinf(ne,0,Inf)$Q)
}

dng <- function(par,x,log=T,cpus=2){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	lterm1 = -(x-mu-sigma2/(2*beta))
	lterm2 = log(1-pnorm(-(x-mu-sigma2/beta)/sigma))
	lterm3 = log(gamma(alpha))+alpha*log(beta)
	lterm4 = (alpha-1)*log(x-mu-sigma2/beta
			+(dnorm(-(x-mu-sigma2/beta)/sigma)*sigma)/(1-pnorm(-(x-mu-sigma2/beta)/sigma)))

	logval = lterm1+lterm2+lterm4-lterm3

	if(log){
		logval
	}else{
		exp(logval)
	}
}

dng.mode <- function(par,x,log=T,cpus=2){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	tmean = x-mu-sigma2/beta
	lesszero = which(tmean<0)
	tmode = tmean
	tmode[lesszero] = 0

	lterm1 = -(x-mu-sigma2/(2*beta))/beta
	lterm2 = log(1-pnorm(-(x-mu-sigma2/beta)/sigma))
	lterm3 = log(gamma(alpha))+alpha*log(beta)

	lterm4 = (alpha-1)*log(tmode)

	if(length(which(is.nan(lterm4)))>0){
		stop()
	}

	logval = lterm1+lterm2+lterm4-lterm3

	if(log){
		logval
	}else{
		exp(logval)
	}
}

dng.mean <- function(par,x,log=T,cpus=2){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	tmean = x-mu-sigma2/beta
	lesszero = which(tmean<0)
	tmode = tmean
	tmode[lesszero] = 0

	lterm1 = -(x-mu-sigma2/(2*beta))/beta
	lterm2 = log(1-pnorm(-(x-mu-sigma2/beta)/sigma))
	lterm3 = log(gamma(alpha))+alpha*log(beta)
	lterm4 = (alpha-1)*log(tmode+(dnorm(-tmode/sigma)*sigma)/(1-pnorm(-tmode/sigma)))

	logval = lterm1+lterm2+lterm4-lterm3

	if(log){
		logval
	}else{
		exp(logval)
	}
}
dng.delta1 <- function(part,x,log=TRUE){
	mu = part[1]
	sigma = part[2]
	alpha = part[3]
	beta = part[4]
	sigma2 = (sigma)^2

	tmean = x-mu-sigma2/beta
	lesszero = which(tmean<0)
	tmode = tmean
	tmode[lesszero] = 0

	lterm1 = -(x-mu-sigma2/(2*beta))/beta
	lterm2 = log(1-pnorm(-(x-mu-sigma2/beta)/sigma))
	lterm3 = log(gamma(alpha))+alpha*log(beta)

	term40 = (dnorm(-tmode/sigma))/(1-pnorm(-tmode/sigma))
	term40[which(is.nan(term40))] = 1
	lterm4 = (alpha-1)*log(tmode+term40*sigma)

	logval = lterm1+lterm2+lterm4-lterm3

	if(log){
		dout = logval
	}else{
		dout = exp(logval)
	}
	dout
}

dng.delta2 <- function(part,x,log=TRUE){
	mu = part[1]
	sigma = part[2]
	alpha = part[3]
	beta = part[4]
	sigma2 = (sigma)^2

	tmean = x-mu-sigma2/beta
	lesszero = which(tmean<0)
	tmode = tmean
	tmode[lesszero] = 0

	lterm1 = -(x-mu-sigma2/(2*beta))/beta
	lterm2 = log(1-pnorm(-(x-mu-sigma2/beta)/sigma))
	lterm3 = log(gamma(alpha))+alpha*log(beta)

	term40 = (dnorm(-tmode/sigma))/(1-pnorm(-tmode/sigma))
	term40[which(is.nan(term40))] = 1
	lterm41 = (alpha-1-2)*log(tmode+term40*sigma)
	term42 = sigma2*(1-(tmode/sigma)*term40-(term40)^2)
	term43 = (tmode+term40*sigma)^2 + (alpha-1)*(alpha-2)*term42/2
	lterm4 = lterm41 + log(term43)

	logval = lterm1+lterm2+lterm4-lterm3

	if(log){
		dout = logval
	}else{
		dout = exp(logval)
	}
	dout
}


dngN <- function(par,x,log=T,cpus=2){
	ng <- function(val,xv){
		n1 = dnorm(xv,(par[1]+val),sqrt(par[2]),log=FALSE)	
		n2 = dgamma(val,par[3],1/par[4],log=FALSE)
		n1*n2
	}
	len = length(x)

	wrapper <- function(x){
		quadinf(ng,0,Inf,xv=x)$Q
	}

	if(len > 1){
		result <- sfLapply(x, wrapper)
		val =  unlist( rbind( result ) )
	}else{
		val = quadinf(ng,0,Inf,xv=x)$Q
	}

	if(log){
		log(val)
	}else{
		val
	}
}	

dng0 <- function(par,x,log=T){
	ng <- function(val,xv){
		n1 = dnorm(xv,(par[1]+val),sqrt(par[2]),log=FALSE)	
		n2 = dgamma(val,par[3],1/par[4],log=FALSE)
		n1*n2
	}
	len = length(x)
	val = rep(0,len)
	for(i in 1:len){
		val[i] = quadinf(ng,0,Inf,xv=x[i])$Q
	}

	if(log){
		log(val)
	}else{
		val
	}
}

dng2 <- function(par,x,log=T){
	par[2] = sqrt(par[2])
	lenx = length(x)
	x = c(x,0)
	tmp = dnormgam(par=par,x=x,log=log,plot=FALSE)$dout
	as.numeric(tmp[1:lenx])
}

conv.ng.signal <- function(x,par){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	lterm1 = (alpha)*log(x-mu-sigma2/beta
			+(dnorm(-(x-mu-sigma2/beta)/sigma)*sigma)/(1-pnorm(-(x-mu-sigma2/beta)/sigma)))
	lterm2 = (alpha-1)*log(x-mu-sigma2/beta
			+(dnorm(-(x-mu-sigma2/beta)/sigma)*sigma)/(1-pnorm(-(x-mu-sigma2/beta)/sigma)))
	
	exp(lterm1-lterm2)
}

conv.ng.signal.delta1 <- function(x,par){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	tmode = x-mu-sigma2/beta

	term40 = (dnorm(-tmode/sigma))/(1-pnorm(-tmode/sigma))
	term40[which(is.nan(term40))] = 1

	(tmode+term40*sigma)
}

conv.ng.signal.delta2 <- function(x,par){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	tmode = x-mu-sigma2/beta

	term40 = (dnorm(-tmode/sigma))/(1-pnorm(-tmode/sigma))
	term40[which(is.nan(term40))] = 1

	term41 = (tmode+term40*sigma)

	term42 = sigma2*(1-(tmode/sigma)*term40-(term40)^2)
	term43 = (tmode+term40*sigma)^2 + (alpha)*(alpha-1)*term42/2
	term44 = (tmode+term40*sigma)^2 + (alpha-1)*(alpha-2)*term42/2

	term41*(term43/term44)
}

conv.ng.signal.mode <- function(x,par){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	tmean = x-mu-sigma2/beta
	if(tmean < 0)
		tmode = 0
	else
		tmode = tmean			

	tmode
}

conv.ng.signal.mean <- function(x,par){
	mu = par[1]
	sigma2 = par[2]
	alpha = par[3]
	beta = par[4]
	sigma = sqrt(sigma2)

	tmean = x-mu-sigma2/beta
	if(tmean < 0)
		tmode = 0
	else
		tmode = tmean			

	tmode+(dnorm(-tmode/sigma)*sigma)/(1-pnorm(-tmode/sigma))
}

conv.ng.signalN <- function(x,par){
	n3 = dng(par=par,x=x,log=FALSE)
	post.d <- function(val){
		n1 = dnorm(x,(par[1]+val),sqrt(par[2]),log=FALSE)	
		n2 = dgamma(val,par[3],1/par[4],log=FALSE)
		val*n1*n2/n3
	}
	quadinf(post.d,0,Inf)$Q
}

conv.ng.signal2 <- function(x,par){
	par[2] = sqrt(par[2])
	normgam.signal(x=x,par=par)
}

###################################
# Sub functions of peak detection #
###################################

lng.p0.x <- function(x,mu,sigma2){
	sigma = sqrt(sigma2)
	p = -((x-mu)^2)/(2*sigma2) - log(sigma*sqrt(2*pi))
	p
}

lng.pa.x <- function(x,mu,sigma2,alpha,beta){
	p = dng(c(mu,sigma2,alpha,beta),x,log=TRUE)
	p
}

lng.pa.x.mean <- function(x,mu,sigma2,alpha,beta){
	p = dng.mean(c(mu,sigma2,alpha,beta),x,log=TRUE)
	p
}

lng.pa.x.mode <- function(x,mu,sigma2,alpha,beta){
	p = dng.mode(c(mu,sigma2,alpha,beta),x,log=TRUE)
	p
}

lng.pa.x.delta1 <- function(x,mu,sigma2,alpha,beta){
	p = dng.delta1(c(mu,sigma2,alpha,beta),x,log=TRUE)
	p
}

lng.pa.x.delta2 <- function(x,mu,sigma2,alpha,beta){
	p = dng.delta2(c(mu,sigma2,alpha,beta),x,log=TRUE)
	p
}

lng.zk <- function(x,mu,sigma2,alpha,beta,p){
	p0 = lng.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = log(dnormgam(c(mu,sqrt(sigma2),alpha,beta),x,plot=FALSE)$dout)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	posz = which(z==0)
	z
}

lng.lc <- function(mu,sigma2,alpha,beta,p,x.v,z.v){
	lpa = log(dnormgam(c(mu,sqrt(sigma2),alpha,beta),x.v,plot=FALSE)$dout)
	lv1 = z.v*(lpa-lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}

lng.zk.delta1 <- function(x,mu,sigma2,alpha,beta,p){
	p0 = lng.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = lng.pa.x.delta1(x=x,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	posz = which(z==0)
	z
}

lng.lc.delta1 <- function(mu,sigma2,alpha,beta,p,x.v,z.v){
	lpa = lng.pa.x.delta1(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	lv1 = z.v*(lpa-lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}


lng.zk.delta2 <- function(x,mu,sigma2,alpha,beta,p){
	p0 = lng.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = lng.pa.x.delta2(x=x,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	posz = which(z==0)
	z
}

lng.lc.delta2 <- function(mu,sigma2,alpha,beta,p,x.v,z.v){
	lpa = lng.pa.x.delta2(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	lv1 = z.v*(lpa-lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}

lng.zk.mean <- function(x,mu,sigma2,alpha,beta,p){
	p0 = lng.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = lng.pa.x.mean(x=x,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	posz = which(z==0)
	z
}

lng.lc.mean <- function(mu,sigma2,alpha,beta,p,x.v,z.v){
	lpa = lng.pa.x.mean(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	lv1 = z.v*(lpa-lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}


lng.zk.mode <- function(x,mu,sigma2,alpha,beta,p){
	p0 = lng.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = lng.pa.x.mode(x=x,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	posz = which(z==0)
	z
}

lng.lc.mode <- function(mu,sigma2,alpha,beta,p,x.v,z.v){
	lpa = lng.pa.x.mode(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	lv1 = z.v*(lpa-lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}

lng.zk0 <- function(x,mu,sigma2,alpha,beta,p){
	p0 = lng.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = lng.pa.x(x=x,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	posz = which(z==0)
	z
}

lng.lc0 <- function(mu,sigma2,alpha,beta,p,x.v,z.v){
	lv1 = z.v*(lng.pa.x(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta)-lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lng.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}

lng.em.log <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		alpha = est[3]
		beta = est[4]
		p = est[5]
		
		z.v = lng.zk(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				alpha = exp(par[3])
				beta = exp(par[4])
				val = lng.lc(mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p,x.v=x.v,z.v=z.v)
				if(is.nan(val)){
					val = -10000000000000
				}
				-val
			}
		fit = nlminb(c(log(est[1:4])),ft)
		est[1:4] = exp(c(fit$par))
		est[5] = p
		rlt = rbind(rlt,est)
	}
	
	list(z=z.v,estimate=rlt,R=R)
}

lng.em.log.delta1 <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		alpha = est[3]
		beta = est[4]
		p = est[5]
		
		z.v = lng.zk.delta1(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				alpha = exp(par[3])
				beta = exp(par[4])
				val = lng.lc.delta1(mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p,x.v=x.v,z.v=z.v)
				if(is.nan(val)){
					val = -10000000000000
				}
				-val
			}
		fit = nlminb(c(log(est[1:4])),ft)
		est[1:4] = exp(c(fit$par))
		est[5] = p
		rlt = rbind(rlt,est)
	}
	
	list(z=z.v,estimate=rlt,R=R)
}

lng.em.log.delta2 <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		alpha = est[3]
		beta = est[4]
		p = est[5]
		
		z.v = lng.zk.delta2(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				alpha = exp(par[3])
				beta = exp(par[4])
				val = lng.lc.delta2(mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p,x.v=x.v,z.v=z.v)
				if(is.nan(val)){
					val = -1000000000000
				}
				-val
			}
		fit = nlminb(c(log(est[1:4])),ft)
		est[1:4] = exp(c(fit$par))
		est[5] = p
		rlt = rbind(rlt,est)
	}
	
	list(z=z.v,estimate=rlt,R=R)
}

lng.em.log.mode <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		alpha = est[3]
		beta = est[4]
		p = est[5]
		
		z.v = lng.zk.mode(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				alpha = exp(par[3])
				beta = exp(par[4])
				val = lng.lc.mode(mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p,x.v=x.v,z.v=z.v)
				if(is.nan(val)){
					val = -10000000000000
				}
				-val
			}
		fit = nlminb(c(log(est[1:4])),ft)
		est[1:4] = exp(c(fit$par))
		est[5] = p
		rlt = rbind(rlt,est)
	}
	
	list(z=z.v,estimate=rlt,R=R)
}

lng.em.log.mean <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		alpha = est[3]
		beta = est[4]
		p = est[5]
		
		z.v = lng.zk.mean(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				alpha = exp(par[3])
				beta = exp(par[4])
				val = lng.lc.mean(mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p,x.v=x.v,z.v=z.v)
				-val
			}
		fit = nlminb(c(log(est[1:4])),ft)
		est[1:4] = exp(c(fit$par))
		est[5] = p
		rlt = rbind(rlt,est)
	}
	
	list(z=z.v,estimate=rlt,R=R)
}

lng.em.log0 <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		alpha = est[3]
		beta = est[4]
		p = est[5]
		
		z.v = lng.zk(x=x.v,mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				alpha = exp(par[3])
				beta = exp(par[4])
				val = lng.lc(mu=mu,sigma2=sigma2,alpha=alpha,beta=beta,p=p,x.v=x.v,z.v=z.v)
				-val
			}
		fit = nlminb(c(log(est[1:4])),ft)
		est[1:4] = exp(c(fit$par))
		est[5] = p
		rlt = rbind(rlt,est)
	}
	
	list(z=z.v,estimate=rlt,R=R)
}

ng.odds.cont.val <- function(data,par,size){
	mu = par[1]
	s2 = par[2]
	alpha = par[3]
	beta = par[4]
	p = par[5]	

     	gridx = 1:size[1]
     	gridy = 1:size[2]

      outmat = matrix(NA, size[1], size[2])
      for(i in 1:size[1]){
		for(j in 1:size[2]){
	      	x = data[(i+size[1]*(j-1))]
			p0 = exp(lng.p0.x(x=x,mu=mu,sigma2=s2))
			pA = exp(lng.pa.x(x=x,mu=mu,sigma2=s2,alpha=alpha,beta=beta))
            	outmat[i,j] = (pA/p0)*(p/(1-p))
		}
     }
     list(gridx=gridx, gridy=gridy, outmat=outmat)
}

ng.z.cont.val <- function(z){
	size = dim(z)
     	gridx = 1:size[1]
     	gridy = 1:size[2]

      outmat = matrix(NA, size[1], size[2])
      for(i in 1:size[1]){
		for(j in 1:size[2]){
	      	x = z[(i+size[1]*(j-1))]
            	outmat[i,j] = x
		}
     }
     list(gridx=gridx, gridy=gridy, outmat=outmat)
}

ini.normgam <- function(Regular){
    	estimation_param_normexp_RMA <- function(pm, n.pts = 2^14) {
        	set.seed(1234)
        	max.density <- function(x, n.pts) {
            	aux <- density(x, kernel = "epanechnikov", n = n.pts, na.rm = TRUE)
            	aux$x[order(-aux$y)[1]]
        	}
        	pmbg <- max.density(pm, n.pts)
        	bg.data <- pm[pm < pmbg]
        	pmbg <- max.density(bg.data, n.pts)
        	bg.data <- pm[pm < pmbg]
        	bg.data <- bg.data - pmbg
        	bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2)
        	sig.data <- pm[pm > pmbg]
        	sig.data <- sig.data - pmbg
        	alpha <- max.density(sig.data, n.pts)
        	mubg <- pmbg
        	list(mu = mubg, sigma = bgsd, alpha = alpha)
	}
      prma = estimation_param_normexp_RMA(Regular)
      mu0 = prma$mu
      sigma0 = prma$sigma
      theta0 = ((sd(Regular))^2 - sigma0^2)/(max(mean(Regular), 
            median(Regular)) - mu0)
      k0 = (max(mean(Regular), median(Regular)) - mu0)/theta0
      if (!((k0 > 0) & (theta0 > 0))) {
      	k0 = 1
            theta0 = max(mean(Regular) - mu0, mu0/10)
      }
	c(mu0,sigma0^2,k0,theta0)
}

sig.peak.ng.fft <- function(data,w=1,ini=NULL){
	d = dim(data)
	cat("# msPeakG-FFT #: Size of signal = ",d,"\n")
	ad = w * c(data)
	if(is.null(ini)){
		ini = c(ini.normgam(ad),.05)
	}
	cat("# msPeakG-FFT #: Initial = ",ini,"\n")
	ad.fit = lng.em.log(data=ad,initial=ini,R=50)
	cat("# msPeakG-FFT #: Estimate = ",ad.fit$est[50,],"\n")
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	ngv = ng.z.cont.val(ad.fit$z)
	list(fit=ad.fit,nodd=ngv)
}

sig.peak.ng.delta1 <- function(data,w=1,ini=NULL){
	d = dim(data)
	cat("# msPeakG-Delta1 #: Size of signal = ",d,"\n")
	ad = w * c(data)
	if(is.null(ini)){
		ini = c(ini.normgam(ad),.05)
	}
	cat("# msPeakG-Delta1 #: Initial = ",ini,"\n")
	ad.fit = lng.em.log.delta1(data=ad,initial=ini,R=50)
	cat("# msPeakG-Delta1 #: Estimate = ",ad.fit$est[50,],"\n")
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	ngv = ng.z.cont.val(ad.fit$z)
	list(fit=ad.fit,nodd=ngv)
}

sig.peak.ng.delta2 <- function(data,w=1,ini=NULL){
	d = dim(data)
	cat("# msPeakG-Delta2 #: Size of signal = ",d,"\n")
	ad = w * c(data)
	if(is.null(ini)){
		ini = c(ini.normgam(ad),.05)
	}
	cat("# msPeakG-Delta2 #: Initial = ",ini,"\n")
	ad.fit = lng.em.log.delta2(data=ad,initial=ini,R=50)
	cat("# msPeakG-Delta2 #: Estimate = ",ad.fit$est[50,],"\n")
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	ngv = ng.z.cont.val(ad.fit$z)
	list(fit=ad.fit,nodd=ngv)
}

sig.peak.ng <- function(data,w=1){
	d = dim(data)
	cat("# msPeakG #: Size of signal = ",d,"\n")
	ad = w * c(data)
	cpar = c(normgam.fit(ad)$par,.05)
	cpar[2] = cpar[2]^2
	cat("# msPeakG #: Initial = ",cpar,"\n")
	ad.fit = lng.em.log(data=ad,initial=cpar,R=50)
	cat("# msPeakG #: Estimate = ",ad.fit$est[50,],"\n")
	ngv = ng.odds.cont.val(ad,ad.fit$est[50,],d)
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	list(fit=ad.fit,nodd=ngv)
}

sig.peak.ng.mean <- function(data,w=1){
	d = dim(data)
	cat("# msPeakG #: Size of signal = ",d,"\n")
	ad = w * c(data)
	cpar = c(normgam.fit(ad)$par,.05)
	cpar[2] = cpar[2]^2
	cat("# msPeakG #: Initial = ",cpar,"\n")
	ad.fit = lng.em.log.mean(data=ad,initial=cpar,R=50)
	cat("# msPeakG #: Estimate = ",ad.fit$est[50,],"\n")
	ngv = ng.odds.cont.val(ad,ad.fit$est[50,],d)
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	list(fit=ad.fit,nodd=ngv)
}

sig.peak.ng.mode <- function(data,w=1){
	d = dim(data)
	cat("# msPeakG #: Size of signal = ",d,"\n")
	ad = w * c(data)
	cpar = c(normgam.fit(ad)$par,.05)
	cpar[2] = cpar[2]^2
	cat("# msPeakG #: Initial = ",cpar,"\n")
	ad.fit = lng.em.log.mode(data=ad,initial=cpar,R=50)
	cat("# msPeakG #: Estimate = ",ad.fit$est[50,],"\n")
	ngv = ng.odds.cont.val(ad,ad.fit$est[50,],d)
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	list(fit=ad.fit,nodd=ngv)
}
true.s.odd.ng.delta1 <- function(par,odd,cv=.5,data){
	odd = c(odd$outmat)

	len = length(odd)
	zpos = which(odd<cv)
	nd = conv.ng.signal.delta1(x=data,par=par)
	nd[zpos] = 0
	nd
}

true.s.odd.ng.delta2 <- function(par,odd,cv=.5,data){
	odd = c(odd$outmat)

	len = length(odd)
	zpos = which(odd<cv)
	nd = conv.ng.signal.delta2(x=data,par=par)
	nd[zpos] = 0
	nd
}

true.s.odd.ng.fft <- function(par,odd,cv=.5,data){
	odd = c(odd$outmat)

	len = length(odd)
	zpos = which(odd<cv)
	nd = normgam.signal(x=data,par=par[1:4])
	nd[zpos] = 0
	nd
}

true.s.odd.ng.mean <- function(par,odd,cv=1,data){
	odd = c(odd$outmat)
	len = length(odd)
	nd = c()
	for(i in 1:len){
		if(odd[i]<cv)
			nd = c(nd,0)
		else{
			td = conv.ng.signal.mean(x=data[i],par=par)
			nd = c(nd,td)
		}
	}
	nd
}

true.s.odd.ng.mode <- function(par,odd,cv=1,data){
	odd = c(odd$outmat)

	len = length(odd)
	nd = c()
	for(i in 1:len){
		if(odd[i]<cv)
			nd = c(nd,0)
		else{
			td = conv.ng.signal.mode(x=data[i],par=par)
			nd = c(nd,td)
		}
	}
	nd
}

true.s.odd.ng0 <- function(par,odd,cv=1,data){
	odd = c(odd$outmat)
	len = length(odd)
	nd = c()
	for(i in 1:len){
		if(odd[i]<cv)
			nd = c(nd,0)
		else{
			td = conv.ng.signal(x=data[i],par=par)
			nd = c(nd,td)
		}
	}
	nd
}

true.s.odd.ng <- function(par,odd,cv=1,data){
	odd = c(odd$outmat)
	len = length(odd)
	zpos = which(odd<cv)
	nd = normgam.signal(x=data,par=par[1:4])
	nd[zpos] = 0
	nd
}

true.s.odd.ng2 <- function(par,odd,cv=1,data){
	odd = c(odd$outmat)
	len = length(odd)
	nd = c()
	for(i in 1:len){
		if(odd[i]<cv)
			nd = c(nd,0)
		else{
			td = conv.ng.signal(x=data[i],par=par)
			nd = c(nd,td)
		}
	}
	nd
}

lne.p0.x <- function(x,mu,sigma2){
	sigma = sqrt(sigma2)
	p = -((x-mu)^2)/(2*sigma2) - log(sigma*sqrt(2*pi))
	p
}

lne.pa.x <- function(x,mu,sigma2,phi){
	sigma = sqrt(sigma2)
	p1 = -(x-mu-.5*sigma2/phi)/phi - log(phi)
	p2 = pnorm((x-mu-sigma2/phi)/sigma,log.p=TRUE)
	p = p1+p2
	p
}

lne.zk <- function(x,mu,sigma2,phi,p){
	p0 = lne.p0.x(x=x,mu=mu,sigma2=sigma2)
	pA = lne.pa.x(x=x,mu=mu,sigma2=sigma2,phi=phi)
	n.p = exp(p0-pA)*(1-p)/p
	z = 1/(1+n.p)
	z
}

lne.lc <- function(mu,sigma2,phi,p,x.v,z.v){
	lv1 = z.v*(lne.pa.x(x=x.v,mu=mu,sigma2=sigma2,phi=phi)-lne.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(p/(1-p)))
	lv2 = lne.p0.x(x=x.v,mu=mu,sigma2=sigma2)+log(1-p)
	lv = sum(lv1 + lv2)
	lv
}

lne.em.log <- function(data,initial,R){
	x.v = data
	est = initial

	rlt = c()
	
	for(i in 1:R){
		mu = est[1]
		sigma2 = est[2]
		phi = est[3]
		p = est[4]
		
		z.v = lne.zk(x=x.v,mu=mu,sigma2=sigma2,phi=phi,p=p)
		
		p = (2+sum(z.v))/(4+length(z.v))

		ft = function(par){
				mu = exp(par[1])
				sigma2 = exp(par[2])
				phi = exp(par[3])
				val = lne.lc(mu=mu,sigma2=sigma2,phi=phi,p=p,x.v=x.v,z.v=z.v)
				if(is.nan(val))
					val = 10000000000
				-val
			}
		fit = nlminb(c(log(est[1:3])),ft)
		est[1:3] = exp(c(fit$par))
		est[4] = p
		rlt = rbind(rlt,est)
	}
	list(z=z.v,estimate=rlt,R=R)
}

ne.odds.cont.val <- function(data,par,size){
	mu = par[1]
	s2 = par[2]
	phi = par[3]
	p = par[4]	

     	gridx = 1:size[1]
     	gridy = 1:size[2]

      outmat = matrix(NA, size[1], size[2])
      for(i in 1:size[1]){
		for(j in 1:size[2]){
	      	x = data[(i+size[1]*(j-1))]
			p0 = exp(lne.p0.x(x=x,mu=mu,sigma2=s2))
			pA = exp(lne.pa.x(x=x,mu=mu,sigma2=s2,phi=phi))
            	outmat[i,j] = (pA/p0)*(p/(1-p))
		}
     }
     list(gridx=gridx, gridy=gridy, outmat=outmat)
}

ne.z.cont.val <- function(z){
	size = dim(z)
     	gridx = 1:size[1]
     	gridy = 1:size[2]

      outmat = matrix(NA, size[1], size[2])
      for(i in 1:size[1]){
		for(j in 1:size[2]){
	      	x = z[(i+size[1]*(j-1))]
            	outmat[i,j] = x
		}
     }
     list(gridx=gridx, gridy=gridy, outmat=outmat)
}

ini.normexp <- function(pm, n.pts = 2^14) {
      set.seed(1234)
      max.density <- function(x, n.pts) {
         	aux <- density(x, kernel = "epanechnikov", n = n.pts, na.rm = TRUE)
         	aux$x[order(-aux$y)[1]]
      }
      pmbg <- max.density(pm, n.pts)
      bg.data <- pm[pm < pmbg]
      pmbg <- max.density(bg.data, n.pts)
      bg.data <- pm[pm < pmbg]
      bg.data <- bg.data - pmbg
      bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2)
      sig.data <- pm[pm > pmbg]
      sig.data <- sig.data - pmbg
      alpha <- max.density(sig.data, n.pts)
      mubg <- pmbg
      if(bgsd<=0 | alpha<0){
 	ad.q = quantile(pm,prob=c(0.95))
	mubg = mean(pm[pm<ad.q])
	bgsd = sqrt(var(pm[pm<ad.q]))
	alpha = mean(pm[pm>ad.q]) # shape
      }
	#list(mu = mubg, sigma = bgsd, alpha = alpha)
	c(mubg, bgsd^2, alpha)
}

sig.peak.ne <- function(data,w=1,ini=NULL){
	d = dim(data)
	ad = w * c(data)
	if(is.null(ini)){
		ini = c(ini.normexp(ad),.05)
	}
	cat("# msPeakG #: Initial = ",ini,"\n")
	ad.fit = lne.em.log(data=ad,initial=ini,R=50)
	cat("# msPeakG #: Estimate = ",ad.fit$est[50,],"\n")
	ad.fit$z = matrix(ad.fit$z,d[1],d[2])
	nev = ne.z.cont.val(ad.fit$z)
	list(fit=ad.fit,nodd=nev)
}

conv.ng.signal.ne <- function(x,par){
	mu = par[1]
	s2 = par[2]
	phi = par[3]
	tp = mu+s2/phi

	x-tp + sqrt(s2)*dnorm((x-tp)/sqrt(s2))/pnorm((x-tp)/sqrt(s2))
}

true.s.odd.ne <- function(par,odd,cv=.5,data){
	odd = c(odd$outmat)

	len = length(odd)
	zpos = which(odd<cv)
	nd = conv.ng.signal.ne(x=data,par=par)
	nd[zpos] = 0
	nd
}

hpd <- function (qfun, conf = 0.95, tol = 1e-08, ...) 
{

	conf <- min(conf, 1 - conf)
    	f <- function(x, qfun, conf, ...) {
        	qfun(1 - conf + x, ...) - qfun(x, ...)
    	}
    	out <- optimize(f, c(0, conf), qfun = qfun, conf = conf, tol = tol, ...)
    	return(c(qfun(out$minimum, ...), qfun(1 - conf + out$minimum, ...)))
}

pemg <- function (q, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) 
{
   	u <- (q - mu) * lambda
    	sigma1 <- lambda * sigma
    	p <- pnorm(u, 0, sigma1) - exp(-u + (sigma1 * sigma1)/2 + 
        pnorm(u, sigma1 * sigma1, sigma1, log.p = TRUE))
    	p <- if (lower.tail) {
        		p
    		}else {
        		1 - p
    		}
    	if (log.p) {
        	log(p)
    	}else {
        	p
    	}
}

demg <- function (x, mu = 0, sigma = 1, lambda = 1, log = FALSE){
	l <- max(length(x), length(mu), length(sigma), length(lambda))
    	x <- rep(x, times = ceiling(l/length(x)), length.out = l)
    	mu <- rep(mu, times = ceiling(l/length(mu)), length.out = l)
    	sigma <- rep(sigma, times = ceiling(l/length(sigma)), length.out = l)
    	lambda <- rep(lambda, times = ceiling(l/length(lambda)), length.out = l)
    	if (min(sigma) <= 0) {
        	stop("Sigma must be greater than zero")
    	}
    	if (min(lambda) <= 0) {
        	stop("Lambda must be greater than zero")
    	}
    	erfc <- pnorm((mu + lambda * sigma * sigma - x)/sigma, lower.tail = FALSE, log.p = log)
    	if (log) {
        	result <- lambda/2 * (2 * mu + lambda * sigma * sigma - 2 * x) + Re(erfc) + log(lambda)
    	}else {
        	result <- exp(lambda/2 * (2 * mu + lambda * sigma * sigma - 2 * x)) * Re(erfc) * lambda
    	}
    	result[is.nan(result)] <- 0
    	result
}

qemg <- function (p, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE){
    	if (length(p[p > 1]) > 1) {
        	stop("p must be equal to or less than 1")
    	}
    	if (length(p[p < 0]) > 1) {
        	stop("p must be equal to or greater than 0")
    	}
    	if (!lower.tail) {
        	p <- 1 - p
    	}
    	result <- rep(NaN, length(p))
    	for (i in 1:length(p)) {
        	result[i] <- optim(c(mu), fn = function(y) {
            						abs(p[i] - pemg(y, mu, sigma, lambda))
        						}, 
					lower = -Inf, upper = Inf, method = "L-BFGS-B")$par
    	}
    	if (log.p) {
        	log(result)
    	}else {
        	result
    	}
}

fsd.pd <- function(data){
	if(length(which(data==0))>0){
		stop("### ERROR!!! ZERO PEAKS IN NONZERO-SUB-REGION IN fad.pd \"fad.pd\" ###")
	}
	len = length(data)
	yder = rep(0,len)
	yarea = yder

	parea = 0 
	ppos = 0 
	for(i in 1:len){
		if(i==1){
			tyder = (data[i]-data[i+1])
			if(tyder>0){
				yder[i] = 1
				ppos = i
			}
			parea = parea + 1
		}else if(i==len){
			tyder = (data[i-1]-data[i])
			parea = parea + 1
			if(tyder<0){
				yder[i] = 1
				yarea[i] = parea
			}
			if(ppos!=0){
				yarea[ppos] = parea
			}
		}else{
			tyder1 = (data[i-1]-data[i])
			tyder2 = (data[i]-data[i+1])
			if(tyder1<=0 & tyder2>0){
				yder[i] = 1
				ppos = i
			}
			parea = parea + 1
			if(tyder1>0 & tyder2<=0){
				yder[i] = -1
				if(ppos==0){ # correct on 06072012
					stop("!!! ERROR for PPOS (ppos=",ppos,") in (i)=(",i,")\n")
				}
				yarea[ppos] = parea
				parea = 1
				ppos = 0
			}
		}
	}
	list(y=yder,area=yarea)
}

pois.pd <- function(data,fdata){
	data = data/sum(data)
	fpos = fdata$y
	farea = fdata$area

	len = length(data)
	ppos = (which(fpos>0))
	pnum = length(ppos)
	parnum = 1

	pnum.check = 0
	for(i in pnum:1){
		tmp.pnum = i
		num.point.4.fit = tmp.pnum*parnum + tmp.pnum - 1 + 1 # include sigma^2
		if(num.point.4.fit<len){
			pnum.check = 1
			break
		}
	}
			
	if(pnum.check==1){
		pnum = tmp.pnum	
		fun <- function(param,num.cluster=1){
			lamt = param[c(1:num.cluster)]
			sigma2 = exp(param[length(param)])
			if(length(which(is.nan(param)))>0){
				return(Inf)
			}
			if(num.cluster==1){
				wt = 1
			}else{
				wt = param[c((num.cluster*parnum+1):(length(param)-1))]
				if((1-sum(wt))<0){
					return(Inf)
				}
			}
			if(length(which(diff(lamt)<0))>0){
				return(Inf)
			}
			wt = c(wt,1-sum(wt))			

			fit = 0
			for(i in 1:len){
				yi = (data[i])
				fi = 0
				for(j in 1:num.cluster){
					tmp = wt[j]*dpois(i,lambda=lamt[j]) 
					fi = fi + tmp
				}
				fit = fit + (yi-fi)^2
			}
			fit/(2*sigma2) + len*log(sigma2)/2 + len*log(2*pi)/2
		}

		opt.pnum = 1
		opt.obj = Inf
		opt.fit = 0
		peak.order = rev(order(data[ppos]))
		for(i in 1:pnum){
			low = c(rep(c(1),i),rep(0,(i-1)),-Inf)
			up = c(rep(c(len),i),rep(1,(i-1)),Inf)
			ini.mt = ppos[peak.order[1:i]]
			ini.mt.od = order(ini.mt)
			ini.mt = ini.mt[ini.mt.od]
			ini = c(ini.mt,rep(1/i,(i-1)),-10)
			tmp.fit = nlminb(ini,fun,lower=low,upper=up,num.cluster=i)
			if(!is.nan(tmp.fit$objective) & opt.obj>tmp.fit$objective){
				opt.obj = tmp.fit$objective
				opt.fit = tmp.fit
				opt.pnum = i
			}
		}

		fpos = rep(0,length(fpos))
		farea = fpos

		tmp.lambda = floor(opt.fit$par[1:opt.pnum])
		tmp.ppos = pmax(tmp.lambda,1)
		tmp.ppos = pmin(tmp.ppos,len)
		tmp.area = tmp.ppos
		for(i in 1:opt.pnum){
			tmp.area[i] = diff(hpd(qpois,lambda=tmp.lambda[i]))
		}
		
		u.tmp.ppos = unique(tmp.ppos)
		if(length(u.tmp.ppos)<opt.pnum){
			n.pos = match(u.tmp.ppos,tmp.ppos)
			fpos[tmp.ppos[n.pos]] = 1
			farea[tmp.ppos[n.pos]] = tmp.area[n.pos]
		}else{
			fpos[tmp.ppos] = 1
			farea[tmp.ppos] = tmp.area
		}

		fdata = list(y=fpos,area=farea,fit=opt.fit,pnum=opt.pnum)
	}
	fdata
}

ab.01 <- function(x,a,b){
	(x-a)/(b-a)
}

i01.ab <- function(x,a,b){
	(b-a)*x+a
}

tgaus.pd <- function(data,fdata){
	data = data/sum(data)
	fpos = fdata$y
	farea = fdata$area

	len = length(data)
	ppos = (which(fpos>0))
	pnum = length(ppos)
	parnum = 2

	pnum.check = 0
	for(i in pnum:1){
		tmp.pnum = i
		num.point.4.fit = tmp.pnum*parnum + tmp.pnum - 1 + 1 # include sigma^2
		if(num.point.4.fit<len){
			pnum.check = 1
			break
		}
	}
			
	if(pnum.check==1){
		pnum = tmp.pnum	
		
		fun <- function(param,num.cluster=1){
			if(length(which(is.nan(param)))>0){
				return(Inf)
			}
			mt = param[c(1:num.cluster)]
			sdt = param[c((num.cluster+1):(num.cluster*parnum))]
			sigma2 = exp(param[length(param)])
			if(num.cluster==1){
				wt = 1
			}else{
				wt = param[c((num.cluster*parnum+1):(length(param)-1))]
				if(1-sum(wt)<0){
					return(Inf)
				}
			}
			if(length(which(diff(mt)<0))>0){
				return(Inf)
			}
			wt = c(wt,1-sum(wt))			

			fit = 0
			for(i in 1:len){
				yi = (data[i])
				fi = 0
				for(j in 1:num.cluster){
					tmp = wt[j]*dtnorm(i,mean=mt[j],sd=sdt[j],lower=1,upper=len) 
					fi = fi + tmp
				}
				fit = fit + (yi-fi)^2
			}
			fit/(2*sigma2) + len*log(sigma2)/2 + len*log(2*pi)/2
		}

		opt.pnum = 1
		opt.obj = Inf
		opt.fit = 0
		peak.order = rev(order(data[ppos]))
		for(i in 1:pnum){
			low = c(rep(c(1,0),each=i),rep(0,(i-1)),-Inf)
			up = c(rep(c(len,len),i),rep(1,(i-1)),Inf)
			ini.mt = ppos[peak.order[1:i]]
			ini.mt.od = order(ini.mt)
			ini.sdt = farea[ppos[peak.order[1:i]]]/2
			ini.mt = ini.mt[ini.mt.od]
			ini.sdt = ini.sdt[ini.mt.od]
			ini = c(ini.mt,ini.sdt,rep(1/i,(i-1)),-10)
			tmp.fit = nlminb(ini,fun,lower=low,upper=up,num.cluster=i)
			if(!is.nan(tmp.fit$objective) & opt.obj>tmp.fit$objective){
				opt.obj = tmp.fit$objective
				opt.fit = tmp.fit
				opt.pnum = i
			}
		}
		fpos = rep(0,length(fpos))
		farea = fpos

		tmp.ppos = round(opt.fit$par[1:opt.pnum])
		tmp.ppos = pmax(tmp.ppos,1)
		tmp.ppos = pmin(tmp.ppos,len)
		tmp.area = 2*qnorm(0.975)*opt.fit$par[(opt.pnum+1):(opt.pnum*2)]
		
		u.tmp.ppos = unique(tmp.ppos)
		if(length(u.tmp.ppos)<opt.pnum){
			n.pos = match(u.tmp.ppos,tmp.ppos)
			fpos[tmp.ppos[n.pos]] = 1
			farea[tmp.ppos[n.pos]] = tmp.area[n.pos]
		}else{
			fpos[tmp.ppos] = 1
			farea[tmp.ppos] = tmp.area
		}

		fdata = list(y=fpos,area=farea,fit=opt.fit,pnum=opt.pnum)
	}
	fdata
}

gaus.pd <- function(data,fdata){
	data = data/sum(data)
	fpos = fdata$y
	farea = fdata$area

	len = length(data)
	ppos = (which(fpos>0))
	pnum = length(ppos)
	parnum = 2

	pnum.check = 0
	for(i in pnum:1){
		tmp.pnum = i
		num.point.4.fit = tmp.pnum*parnum + tmp.pnum - 1 + 1 # include sigma^2
		if(num.point.4.fit<len){
			pnum.check = 1
			break
		}
	}
			
	if(pnum.check==1){
		pnum = tmp.pnum	
		
		fun <- function(param,num.cluster=1){
			if(length(which(is.nan(param)))>0){
				return(Inf)
			}
			mt = param[c(1:num.cluster)]
			sdt = param[c((num.cluster+1):(num.cluster*parnum))]
			sigma2 = exp(param[length(param)])
			if(num.cluster==1){
				wt = 1
			}else{
				wt = param[c((num.cluster*parnum+1):(length(param)-1))]
				if(1-sum(wt)<0){
					return(Inf)
				}
			}
			if(length(which(diff(mt)<0))>0){
				return(Inf)
			}
			wt = c(wt,1-sum(wt))			

			fit = 0
			for(i in 1:len){
				yi = (data[i])
				fi = 0
				for(j in 1:num.cluster){
					tmp = wt[j]*dnorm(i,mean=mt[j],sd=sdt[j]) 
					fi = fi + tmp
				}
				fit = fit + (yi-fi)^2
			}
			fit/(2*sigma2) + len*log(sigma2)/2 + len*log(2*pi)/2
		}

		opt.pnum = 1
		opt.obj = Inf
		opt.fit = 0
		peak.order = rev(order(data[ppos]))
		for(i in 1:pnum){
			low = c(rep(c(1,0),each=i),rep(0,(i-1)),-Inf)
			up = c(rep(c(len,len),i),rep(1,(i-1)),Inf)
			ini.mt = ppos[peak.order[1:i]]
			ini.mt.od = order(ini.mt)
			ini.sdt = farea[ppos[peak.order[1:i]]]/2
			ini.mt = ini.mt[ini.mt.od]
			ini.sdt = ini.sdt[ini.mt.od]
			ini = c(ini.mt,ini.sdt,rep(1/i,(i-1)),-10)
			tmp.fit = nlminb(ini,fun,lower=low,upper=up,num.cluster=i)			
			if(!is.nan(tmp.fit$objective) & opt.obj>tmp.fit$objective){
				opt.obj = tmp.fit$objective
				opt.fit = tmp.fit
				opt.pnum = i
			}
		}

		fpos = rep(0,length(fpos))
		farea = fpos

		tmp.ppos = round(opt.fit$par[1:opt.pnum])
		tmp.ppos = pmax(tmp.ppos,1)
		tmp.ppos = pmin(tmp.ppos,len)
		tmp.area = 2*qnorm(0.975)*opt.fit$par[(opt.pnum+1):(opt.pnum*2)]
		
		u.tmp.ppos = unique(tmp.ppos)
		if(length(u.tmp.ppos)<opt.pnum){
			n.pos = match(u.tmp.ppos,tmp.ppos)
			fpos[tmp.ppos[n.pos]] = 1
			farea[tmp.ppos[n.pos]] = tmp.area[n.pos]
		}else{
			fpos[tmp.ppos] = 1
			farea[tmp.ppos] = tmp.area
		}

		fdata = list(y=fpos,area=farea,fit=opt.fit,pnum=opt.pnum)
	}
	fdata
}

gamma.pd <- function(data,fdata){
	data = data/sum(data)
	fpos = fdata$y
	farea = fdata$area

	len = length(data)
	ppos = (which(fpos>0))
	pnum = length(ppos)
	parnum = 2

	pnum.check = 0
	for(i in pnum:1){
		tmp.pnum = i
		num.point.4.fit = tmp.pnum*parnum + tmp.pnum - 1 + 1 # include sigma^2
		if(num.point.4.fit<len){
			pnum.check = 1
			break
		}
	}
			
	if(pnum.check==1){
		pnum = tmp.pnum	
		
		fun <- function(param,num.cluster=1){
			if(length(which(is.nan(param)))>0){
				return(Inf)
			}
			mt = param[c(1:num.cluster)]
			sdt = param[c((num.cluster+1):(num.cluster*parnum))]
			sigma2 = exp(param[length(param)])
			if(num.cluster==1){
				wt = 1
			}else{
				wt = param[c((num.cluster*parnum+1):(length(param)-1))]
				if(1-sum(wt)<0){
					return(Inf)
				}
			}
			if(length(which(diff(mt)<0))>0){
				return(Inf)
			}
			wt = c(wt,1-sum(wt))			

			fit = 0
			for(i in 1:len){
				yi = (data[i])
				fi = 0
				for(j in 1:num.cluster){
					tmp = wt[j]*dgamma(i,shape=mt[j],rate=sdt[j]) 
					fi = fi + tmp
				}
				fit = fit + (yi-fi)^2
			}
			fit/(2*sigma2) + len*log(sigma2)/2 + len*log(2*pi)/2
		}

		opt.pnum = 1
		opt.obj = Inf
		opt.fit = 0
		peak.order = rev(order(data[ppos]))
		for(i in 1:pnum){
			low = c(rep(c(1,0),each=i),rep(0,(i-1)),-Inf)
			up = c(rep(c(len,len),i),rep(1,(i-1)),Inf)
			ini.mt = ppos[peak.order[1:i]]
			ini.mt.od = order(ini.mt)
			ini.mt = ini.mt[ini.mt.od]
			ini.sdt = rep(1,i) #ini.sdt[ini.mt.od]
			ini = c(ini.mt,ini.sdt,rep(1/i,(i-1)),-10)
			tmp.fit = nlminb(ini,fun,lower=low,upper=up,num.cluster=i)			
			if(!is.nan(tmp.fit$objective) & opt.obj>tmp.fit$objective){
				opt.obj = tmp.fit$objective
				opt.fit = tmp.fit
				opt.pnum = i
			}
		}

		fpos = rep(0,length(fpos))
		farea = fpos

		tmp.ppos = round(opt.fit$par[1:opt.pnum])
		tmp.ppos = pmax(tmp.ppos,1)
		tmp.ppos = pmin(tmp.ppos,len)
		tmp.scale = opt.fit$par[(opt.pnum+1):(opt.pnum*2)]
		tmp.area = tmp.ppos
		for(i in 1:opt.pnum){
			tmp.area[i] = diff(hpd(qgamma,shape=tmp.ppos[i],rate=tmp.scale[i]))
			fun.gamma <- function(param){
				fit = -dgamma(param,shape=tmp.ppos[i],rate=tmp.scale[i],log=T)
				fit
			}
			tmp.ppos[i] = round(nlminb(1,fun.gamma,lower=1,upper=len)$par)
		}
		
		u.tmp.ppos = unique(tmp.ppos)
		if(length(u.tmp.ppos)<opt.pnum){
			n.pos = match(u.tmp.ppos,tmp.ppos)
			fpos[tmp.ppos[n.pos]] = 1
			farea[tmp.ppos[n.pos]] = tmp.area[n.pos]
		}else{
			fpos[tmp.ppos] = 1
			farea[tmp.ppos] = tmp.area
		}

		fdata = list(y=fpos,area=farea,fit=opt.fit,pnum=opt.pnum)
	}
	fdata
}

emg.pd <- function(data,fdata){
	data = data/sum(data)
	fpos = fdata$y
	farea = fdata$area

	len = length(data)
	ppos = (which(fpos>0))
	pnum = length(ppos)
	parnum = 3

	pnum.check = 0
	for(i in pnum:1){
		tmp.pnum = i
		num.point.4.fit = tmp.pnum*parnum + tmp.pnum - 1 + 1 # include sigma^2
		if(num.point.4.fit<len){
			pnum.check = 1
			break
		}
	}
			
	if(pnum.check==1){
		pnum = tmp.pnum	
		
		fun <- function(param,num.cluster=1){
			if(length(which(is.nan(param)))>0){
				return(Inf)
			}
			mt = param[c(1:num.cluster)]
			sdt = param[c((num.cluster+1):(num.cluster*2))]
			lamt = param[c((num.cluster*2+1):(num.cluster*parnum))]
			sigma2 = exp(param[length(param)])
			if(length(which(sdt==0))>0){
				return(Inf)
			}
			if(length(which(lamt==0))>0){
				return(Inf)
			}
			if(num.cluster==1){
				wt = 1
			}else{
				wt = param[c((num.cluster*parnum+1):(length(param)-1))]
				if(1-sum(wt)<0){
					return(Inf)
				}
			}
			if(length(which(diff(mt)<0))>0){
				return(Inf)
			}
			wt = c(wt,1-sum(wt))			

			fit = 0
			for(i in 1:len){
				yi = (data[i])
				fi = 0
				for(j in 1:num.cluster){
					tmp = wt[j]*demg(i,mu=mt[j],sigma=sdt[j],lambda=lamt[j]) 
					fi = fi + tmp
				}
				fit = fit + (yi-fi)^2
			}
			fit/(2*sigma2) + len*log(sigma2)/2 + len*log(2*pi)/2
		}

		opt.pnum = 1
		opt.obj = Inf
		opt.fit = 0
		peak.order = rev(order(data[ppos]))
		for(i in 1:pnum){
			low = c(rep(c(1,0,0),each=i),rep(0,(i-1)),-Inf)
			up = c(rep(c(len,len,len),i),rep(1,(i-1)),Inf)
			ini.mt = ppos[peak.order[1:i]]
			ini.mt.od = order(ini.mt)
			ini.sdt = farea[ppos[peak.order[1:i]]]/2
			ini.mt = ini.mt[ini.mt.od]
			ini.sdt = ini.sdt[ini.mt.od]
			ini = c(ini.mt,ini.sdt,rep(0.1,i),rep(1/i,(i-1)),-10)
			tmp.fit = nlminb(ini,fun,lower=low,upper=up,num.cluster=i)			
			if(!is.nan(tmp.fit$objective) & opt.obj>tmp.fit$objective){
				opt.obj = tmp.fit$objective
				opt.fit = tmp.fit
				opt.pnum = i
			}
		}

		fpos = rep(0,length(fpos))
		farea = fpos

		tmp.sigma = opt.fit$par[(opt.pnum+1):(opt.pnum*2)]
		tmp.lambda = opt.fit$par[(opt.pnum*2+1):(opt.pnum*3)]
		tmp.mu = round(opt.fit$par[1:opt.pnum]) #+1/tmp.lambda)
		tmp.ppos = pmax(tmp.mu,1)
		tmp.ppos = pmin(tmp.ppos,len)
		tmp.area = tmp.ppos
		for(i in 1:opt.pnum){
			tmp.area[i] = diff(hpd(qemg,mu=tmp.ppos[i],sigma=tmp.sigma[i],lambda=tmp.lambda[i]))
			fun.emg <- function(param){
				fit = -demg(param,mu=tmp.ppos[i],sigma=tmp.sigma[i],lambda=tmp.lambda[i],log=T)
				fit
			}
			tmp.ppos[i] = round(nlminb(1,fun.emg,lower=1,upper=len)$par)
		}
		
		u.tmp.ppos = unique(tmp.ppos)
		if(length(u.tmp.ppos)<opt.pnum){
			n.pos = match(u.tmp.ppos,tmp.ppos)
			fpos[tmp.ppos[n.pos]] = 1
			farea[tmp.ppos[n.pos]] = tmp.area[n.pos]
		}else{
			fpos[tmp.ppos] = 1
			farea[tmp.ppos] = tmp.area
		}

		fdata = list(y=fpos,area=farea,fit=opt.fit,pnum=opt.pnum)
	}
	fdata
}

mb.pd <- function(data=ts.ad,method=1,type=1,sft=c(31,70),sst=c(401,500)){
	method.pd <- function(mdata,mfdata){
		if(method==2){
			return(pois.pd(data=mdata,fdata=mfdata))
		}else if(method==3){
			return(tgaus.pd(data=mdata,fdata=mfdata))
		}else if(method==4){
			return(gaus.pd(data=mdata,fdata=mfdata))
		}else if(method==5){
			return(gamma.pd(data=mdata,fdata=mfdata))
		}else if(method==6){
			return(emg.pd(data=mdata,fdata=mfdata))
		}
	}

	d = dim(data) # 1st rt by 2nd rt
	pos = c()
	gnum = 0
	measure = c()
	for(i in 1:d[1]){ # by 1st rt
		td = as.numeric(data[i,])
		gpos = which(td>0)
		glen = length(gpos)
		if(glen==0)
			next

		dpos = which(diff(gpos)>1)
		dlen = length(dpos)
		if(dlen==0){
			std = td[gpos]
			if(length(std)<2){
				next
			}
			tmp.pdata = fsd.pd(data=std)
			if(length(which(tmp.pdata$y>0))==0)
				next
			if(method!=1){
				tmp.pdata = method.pd(mdata=std,mfdata=tmp.pdata)
				tmp = tmp.pdata
				if(length(tmp.pdata$pnum)>0){
					tmp.len = length(std)
					if(type==1){
						tmp.sigma2 = exp(tmp$fit$par[length(tmp$fit$par)])
						tmp.mse = (tmp$fit$objective-tmp.len*log(tmp.sigma2)/2-tmp.len*log(2*pi)/2)*2*tmp.sigma2/tmp.len
						measure = c(measure,tmp.mse)
						rm(tmp.len)
						rm(tmp.sigma2)
						rm(tmp.mse)
					}else if(type==2)
						measure = c(measure,2*tmp$fit$objective+2*length(tmp$fit$par))
					else if(type==3)
						measure = c(measure,2*tmp$fit$objective+tmp.len*length(tmp$fit$par))
				}
				rm(tmp)
			}	

			pd = tmp.pdata$y
			parea = as.numeric(tmp.pdata$area)
			ppos = which(pd==1)
			plen = length(ppos)
			if(plen==0)
				next

			gnum = gnum+1

			trg = c(gnum,i,gpos[1],gpos[glen],(sft[1]-1+i),(sst[1]-1+gpos[1]),(sst[1]-1+gpos[glen]))
			for(j in 1:plen){
				top = ppos[j]+(trg[3]-1) # t2-based peak matrix position 
				pa = parea[ppos[j]] # peak area
				tpp = sst[1]-1+top # t2-based peak matrix position of the whole data
				if(tpp>=trg[6] & tpp<=trg[7]){
					pos = rbind(pos,c(trg,top,tpp,pa))
				}
			}
		}else{
			for(j in 1:(dlen+1)){
				gnum = gnum+1
				if(j==1){
					trg = c(gnum,i,gpos[1],gpos[dpos[j]],(sft[1]-1+i),(sst[1]-1+gpos[1]),(sst[1]-1+gpos[dpos[j]]))
				}else if(j<=dlen){
					trg = c(gnum,i,gpos[dpos[j-1]+1],gpos[dpos[j]],(sft[1]-1+i),(sst[1]-1+gpos[dpos[j-1]+1])
							,(sst[1]-1+gpos[dpos[j]]))
				}else{
					trg = c(gnum,i,gpos[dpos[j-1]+1],gpos[glen],(sft[1]-1+i),(sst[1]-1+gpos[dpos[j-1]+1])
							,(sst[1]-1+gpos[glen]))
				}

				std = td[(trg[3]):(trg[4])]
				if(length(std)<2){
					next
				}
				tmp.pdata = fsd.pd(data=std)
				if(length(which(tmp.pdata$y>0))==0)
					next
				if(method!=1){
					tmp.pdata = method.pd(mdata=std,mfdata=tmp.pdata)
					tmp = tmp.pdata
					if(length(tmp.pdata$pnum)>0){
						tmp.len = length(std)
						if(type==1){
							tmp.sigma2 = exp(tmp$fit$par[length(tmp$fit$par)])
							tmp.mse = (tmp$fit$objective-tmp.len*log(tmp.sigma2)/2-tmp.len*log(2*pi)/2)*2*tmp.sigma2/tmp.len
							measure = c(measure,tmp.mse)
							rm(tmp.len)
							rm(tmp.sigma2)
							rm(tmp.mse)
						}else if(type==2)
							measure = c(measure,2*tmp$fit$objective+2*length(tmp$fit$par))
						else if(type==3)
							measure = c(measure,2*tmp$fit$objective+tmp.len*length(tmp$fit$par))
					}
					rm(tmp)
				}
	
				pd = tmp.pdata$y
				parea = as.numeric(tmp.pdata$area)
				ppos = which(pd==1)
				plen = length(ppos)
				if(plen==0)
					next

				for(k in 1:plen){
					top = ppos[k]+trg[3]-1
					pa = parea[ppos[k]] # peak area
					tpp = sst[1]-1+top
					if(tpp>=trg[6] & tpp<=trg[7]){
						pos = rbind(pos,c(trg,top,tpp,pa))
					}
				}
			}
		}
	}
	if(is.null(pos)){
		cat("# msPeakG #: POS in mb.pd is NULL\n")
		list(pos=NULL,fit=NULL)
	}else{ 
		pos = cbind(c(1:dim(pos)[1]),pos)
		list(pos=pos,fit=mean(measure))
	}
}

mb.pd.opt <- function(data=ts.ad,type=1,sft=c(31,70),sst=c(401,500)){
	# type = 1: ll, 2: aic, 3: bic
	method.pd <- function(mdata,mfdata,vmethod=2){
		if(vmethod==2){
			return(pois.pd(data=mdata,fdata=mfdata))
		}else if(vmethod==3){
			return(tgaus.pd(data=mdata,fdata=mfdata))
		}else if(vmethod==4){
			return(gaus.pd(data=mdata,fdata=mfdata))
		}else if(vmethod==5){
			return(gamma.pd(data=mdata,fdata=mfdata))
		}else if(vmethod==6){
			return(emg.pd(data=mdata,fdata=mfdata))
		}
	}

	d = dim(data) # 1st rt by 2nd rt
	pos = c()
	gnum = 0

	measure = c()

	for(i in 1:d[1]){ # by 1st rt
		td = as.numeric(data[i,])
		gpos = which(td>0)
		glen = length(gpos)
		if(glen==0)
			next

		dpos = which(diff(gpos)>1)
		dlen = length(dpos)
		if(dlen==0){
			std = td[gpos]
			if(length(std)<2){
				next
			}
			tmp.pdata = fsd.pd(data=std)

			###################
			# Optimization
			###################

			selection = c()
			tsp = list()			
			for(i.method in 2:6){
				tmp = method.pd(mdata=std,mfdata=tmp.pdata,vmethod=i.method)
				tsp[[(i.method-1)]] = tmp
				if(length(tmp$pnum)>0){
					tmp.len = length(std)
					if(type==1){
						tmp.sigma2 = exp(tmp$fit$par[length(tmp$fit$par)])
						tmp.mse = (tmp$fit$objective-tmp.len*log(tmp.sigma2)/2-tmp.len*log(2*pi)/2)*2*tmp.sigma2/tmp.len
						selection = c(selection,tmp.mse)
						rm(tmp.len)
						rm(tmp.sigma2)
						rm(tmp.mse)
					}else if(type==2)
						selection = c(selection,2*tmp$fit$objective+2*length(tmp$fit$par))
					else if(type==3)
						selection = c(selection,2*tmp$fit$objective+tmp.len*length(tmp$fit$par))
				}else{
					selection = c(selection,Inf)
				}
				rm(tmp)
			}	
			if(length(which(selection==Inf))<5){
				best.pos = which(selection==min(selection)[1])[1]
				tmp.pdata = tsp[[best.pos]]
				measure = c(measure,selection[best.pos])
			}
			rm(tsp)
			rm(selection)		

			pd = tmp.pdata$y
			parea = as.numeric(tmp.pdata$area)
			ppos = which(pd==1)
			plen = length(ppos)
			if(plen==0)
				next

			gnum = gnum+1

			trg = c(gnum,i,gpos[1],gpos[glen],(sft[1]-1+i),(sst[1]-1+gpos[1]),(sst[1]-1+gpos[glen]))
			for(j in 1:plen){
				top = ppos[j]+(trg[3]-1) # t2-based peak matrix position 
				pa = parea[ppos[j]] # peak area
				tpp = sst[1]-1+top # t2-based peak matrix position of the whole data
				if(tpp>=trg[6] & tpp<=trg[7]){
					pos = rbind(pos,c(trg,top,tpp,pa))
				}
			}
		}else{
			for(j in 1:(dlen+1)){
				gnum = gnum+1
				if(j==1){
					trg = c(gnum,i,gpos[1],gpos[dpos[j]],(sft[1]-1+i),(sst[1]-1+gpos[1]),(sst[1]-1+gpos[dpos[j]]))
				}else if(j<=dlen){
					trg = c(gnum,i,gpos[dpos[j-1]+1],gpos[dpos[j]],(sft[1]-1+i),(sst[1]-1+gpos[dpos[j-1]+1])
							,(sst[1]-1+gpos[dpos[j]]))
				}else{
					trg = c(gnum,i,gpos[dpos[j-1]+1],gpos[glen],(sft[1]-1+i),(sst[1]-1+gpos[dpos[j-1]+1])
							,(sst[1]-1+gpos[glen]))
				}

				std = td[(trg[3]):(trg[4])]
				if(length(std)<2){
					next
				}
				tmp.pdata = fsd.pd(data=std)

				###################
				# Optimization
				###################

				selection = c()
				tsp = list()			
				for(i.method in 2:6){
					tmp = method.pd(mdata=std,mfdata=tmp.pdata,vmethod=i.method)
					tsp[[(i.method-1)]] = tmp
					if(length(tmp$pnum)>0){
						tmp.len = length(std)
						if(type==1){
							tmp.sigma2 = exp(tmp$fit$par[length(tmp$fit$par)])
							tmp.mse = (tmp$fit$objective-tmp.len*log(tmp.sigma2)/2-tmp.len*log(2*pi)/2)*2*tmp.sigma2/tmp.len
							selection = c(selection,tmp.mse)
							rm(tmp.len)
							rm(tmp.sigma2)
							rm(tmp.mse)
						}else if(type==2)
							selection = c(selection,2*tmp$fit$objective+2*length(tmp$fit$par))
						else if(type==3)
							selection = c(selection,2*tmp$fit$objective+tmp.len*length(tmp$fit$par))
					}else{
						selection = c(selection,Inf)
					}
					rm(tmp)
				}	
				if(length(which(selection==Inf))<5){
					best.pos = which(selection==min(selection)[1])[1]
					tmp.pdata = tsp[[best.pos]]
					measure = c(measure,selection[best.pos])
				}
				rm(tsp)
				rm(selection)	

				pd = tmp.pdata$y
				parea = as.numeric(tmp.pdata$area)
				ppos = which(pd==1)
				plen = length(ppos)
				if(plen==0)
					next

				for(k in 1:plen){
					top = ppos[k]+trg[3]-1
					pa = parea[ppos[k]] # peak area
					tpp = sst[1]-1+top
					if(tpp>=trg[6] & tpp<=trg[7]){
						pos = rbind(pos,c(trg,top,tpp,pa))
					}
				}
			}
		}
	}
	pos = cbind(c(1:dim(pos)[1]),pos)
	list(pos=pos,fit=mean(measure))
}

gp.rg <- function(data=plist){
	len = dim(data)[1]
	nd = data
	if(is.null(len)){
		nd = cbind(nd,rep(1,len))
		gp = 1
   	}else if(len==1){
		nd = cbind(nd,rep(1,len))
		gp = 1
   	}else{
	nd = cbind(data,rep(0,len))
	cls = list()
	for(i in 1:(len-1)){
		rd1 = c(nd[i,3:5])
		tc = c(i)
		for(j in (i+1):len){
			rd2 = c(nd[j,3:5])
			if(abs(rd1[1]-rd2[1])<=1){
				if((rd2[2]<=rd1[2] & rd1[3]<=rd2[3])
					| (rd1[2]<=rd2[2] & rd2[3]<=rd1[3])
					| ((rd2[2]-rd1[3])<=1 & (rd1[2]-rd2[3])<=1)
					| ((rd1[2]-rd2[3])<=1 & (rd2[2]-rd1[3])<=1)
				){
					tc = c(tc,j)
				}
			}
		}
		cls[[i]] = tc
	}
	cls[[len]] = len
	gp = list()
	gn = 0
	len = length(cls)
	for(i in 1:(len-1)){
		for(j in (i+1):(len)){
			tn = na.omit(match(cls[[i]],cls[[j]]))
			if(length(tn)>0){
				tg = unique(c(cls[[i]],cls[[j]]))
				if(gn>0){
					gflag = 0
					for(k in 1:gn){
						tn2 = na.omit(match(tg,gp[[k]]))
						if(length(tn2)>0){
							gp[[k]] = unique(c(gp[[k]],tg))
							gflag = gflag+1
						}
					}
					if(gflag==0){
						gn = gn+1
						gp[[gn]] = tg
					}
				}else if(gn==0){
					gn = gn+1
					gp[[gn]] = tg
				}
			}else{
				if(gn==0){
					gn = gn+1
					gp[[gn]] = cls[[i]]
					gn = gn+1
					gp[[gn]] = cls[[j]]
				}else{
					gflag1 = 0
					gflag2 = 0
					for(k in 1:gn){
						tn2 = na.omit(match(cls[[i]],gp[[k]]))
						tn3 = na.omit(match(cls[[j]],gp[[k]]))
						if(length(tn2)>0){
							gp[[k]] = unique(c(gp[[k]],cls[[i]]))
							gflag1 = gflag1+1
						}
						if(length(tn3)>0){
							gp[[k]] = unique(c(gp[[k]],cls[[j]]))
							gflag2 = gflag2+1
						}
					}
					if(gflag1==0){
						gn = gn+1
						gp[[gn]] = cls[[i]]
					}
					if(gflag2==0){
						gn = gn+1
						gp[[gn]] = cls[[j]]
					}
				}
			}
		}
	}
	# group assign
	for(i in 1:length(gp)){
		tg = gp[[i]]
		nd[tg,dim(nd)[2]]=i
	}	
   }

	list(g=gp,nd=nd)
}

get.ms <- function(ft=1,st=1,sic.name="sicdata"){
	fname = paste(sic.name,"_",ft,".txt",sep="")
	td = read.table(fname)
	if(td[1,1]==-333){
		if(td[1,2]!=ft){
			stop(paste("SIC File Name is not matched to Info:",ft," vs. ",td[1,2],sep=""))
		}
		start.st = td[1,3]
		st = st-start.st+1
		td = td[-1,]
	}
	if(st<=0){
		ms = c(t(td[1,]))
		ms =rep(0,length(ms))
	}else{
		ms = c(t(td[st,]))
	}
	ms
}

peak.mg <- function(gdata=g.ad$nd,data=ts.ad,cv=.95,sic.name="sicdata"){
	gdata = cbind(gdata,rep(0,dim(gdata)[1]))
	gdata = as.data.frame(gdata)
	g = unique(gdata[,12])
	len = length(g)
	for(i in 1:len){
		pos = gdata[gdata[,12]==g[i],c(3,9,6,10,11)]
		pos = data.frame(pos[,1],pos[,2],pos[,3],pos[,4],pos[,5])
		slen = dim(pos)[1]
		if(slen<=1)
			next
		sm = c()
		for(j in 1:slen){
			sm = cbind(sm,get.ms(pos[j,3],pos[j,4],sic.name=sic.name))
		}
		sms = cor(sm)
		sms[is.na(sms)] = 0
		for(j in 1:(slen-1)){
			for(k in (j+1):slen){
				if(as.numeric(sms[j,k])>=cv){
					pa.j = gdata[gdata[,3]==pos[j,1]
						& gdata[,9]==pos[j,2],13]
					pa.k = gdata[gdata[,3]==pos[k,1]
						& gdata[,9]==pos[k,2],13]

					if(pa.j==0)
						pa.j = pos[j,5]
					if(pa.k==0)
						pa.k = pos[k,5]

					tmp.pa = pa.j + pa.k
					if(data[pos[j,1],pos[j,2]]>=data[pos[k,1],pos[k,2]]){
						gdata[gdata[,3]==pos[k,1]
							& gdata[,9]==pos[k,2],c(12,13)]=c(0,0)						
						gdata[gdata[,3]==pos[j,1]
							& gdata[,9]==pos[j,2],13]=tmp.pa						
					}else{
						gdata[gdata[,3]==pos[j,1]
							& gdata[,9]==pos[j,2],c(12,13)]=c(0,0)
						gdata[gdata[,3]==pos[k,1]
							& gdata[,9]==pos[k,2],13]=tmp.pa						
					}
				}
			}
		}
	}
	gdata[gdata[,12]==0,13]=0
	gdata[gdata[,13]==0 & gdata[,12]!=0,13] = gdata[gdata[,13]==0 & gdata[,12]!=0,11]
	gdata	
}

ptable.wo.merg <- function(data=g.ad$nd,subt1=subrt1,subt2=subrt2,p1=6,p2=0.005,pmin1=300,pmin2=0.02,sic.name="sicdata"){
	len = dim(data)[1]
	rlt = c()
	for(i in 1:len){
		td = as.numeric(data[i,])
		tmp.ms = as.numeric(get.ms(td[6],td[10],sic.name=sic.name))
		tpos = c(td[3],td[9])
		tpos[1] = (subt1[1]-1+tpos[1])*p1+pmin1
		tpos[2] = (subt2[1]-1+tpos[2])*p2+pmin2
		rlt = rbind(rlt,c(i,tpos,td[3],td[9],td[11]))
	}
	rlt = data.frame(rlt)
	dimnames(rlt)[[2]] = c("id","t1","t2","tid1","tid2","area")
	rlt
}

ptable <- function(data=ms.ad.odd100.y,subt1=subrt1,subt2=subrt2,p1=6,p2=0.005,pmin1=300,pmin2=0.02,sic.name="sicdata"){
	len = dim(data)[1]
	rlt = c()
	for(i in 1:len){
		td = as.numeric(data[i,])
		tmp.ms = as.numeric(get.ms(td[6],td[10],sic.name=sic.name))
		tpos = c(td[3],td[9])
		tpos[1] = (subt1[1]-1+tpos[1])*p1+pmin1
		tpos[2] = (subt2[1]-1+tpos[2])*p2+pmin2
		rlt = rbind(rlt,c(i,tpos,td[3],td[9],td[12],td[11],td[13],tmp.ms))
	}
	rlt = data.frame(rlt)
	dimnames(rlt)[[2]][1:8] = c("id","t1","t2","tid1","tid2","pgroup","parea","area")
	rlt
}

#########
#########
#########
# msPeakG
#########

msPeakG <- function(
			x # output from "find.sub"
			,zcv=.5 # the cut-off of the z-prob for peak region detection
			,method=1 # 1=first derivative, 2=poisson, 3=tgaussian, 4=gaussian, 5=gamma, 6=emg
			,cv=0.95 # the cut-off of the mass spectrum similarity score for peak merging
			,out.fname=NA
			,tic.fname=NA # the filename of the TICs
			,merge=T # peak merging?
			,verbose=T
			,sic.name="sicdata"
			,ngb=T # normal-gamma-bernoulli = T; normal-exp-bernoulli = F
			,den="fft" # "fft","delta1","delta2"
			,imax=10 # maximum intensity
){
	if(is.na(out.fname)){
		out.fname = gsub("-","",Sys.Date())
	}
	
	if(is.na(tic.fname)){
		tic.fname = x$tic.fname
	}

	odds = zcv

	xarea = x

	subrt1 = x$idxt1 # the indices for the first dimension retention time
	subrt2 = x$idxt2 # the indices for the second dimension retention time

	r.t1 = x$range.t1
	r.t2 = x$range.t2

	# load the whole TIC map
	atic = as.matrix(read.table(tic.fname))
	if(verbose)
		cat("# msPeakG #: Loaded the TIC data set with the maximum intensity of ",imax," ...\n")
	
	# calculate the weight to magnify the whole signal; it will not affect the result
	if(is.null(imax)){
		imax = 0
		sad.weight = 1
	}else{
		sad.weight = imax/max(c(atic))
	}

	# magnify with the weight
	watic = atic * sad.weight

	if(verbose)
		cat("# msPeakG #: Selected the sub-data region between RT1 (",r.t1,") and RT2 (",r.t2,") ...\n")

	# select the interesting subregion
	s.ad =  t(watic)[subrt1,subrt2] # original TIC

	if(verbose)
		cat("# msPeakG #: Step I - Peak region finding with the z-cv ",odds," ...\n")

	# peak finding using hierarchical models
	if(ngb){
		if(den == "fft")
			fs.ad = sig.peak.ng.fft(s.ad) # peak finding # NGB model
		else if(den == "delta1")
			fs.ad = sig.peak.ng.delta1(s.ad) # peak finding # NGB model
		else
			fs.ad = sig.peak.ng.delta2(s.ad) # peak finding # NGB model

		ngblabel = paste("Gamma",".",den,sep="")
	}else{
		fs.ad = sig.peak.ne(s.ad) # NEB model
		ngblabel = "Exp"
	}

		# correct the data by doing baseline correction and denoising
		if(ngb){
			if(den == "fft")
				vts.ad.odd = true.s.odd.ng.fft(par=fs.ad$fit$est[fs.ad$fit$R,],fs.ad$nodd,cv=odds,data=c(s.ad)) # convoluted data
			else if(den == "delta1")
				vts.ad.odd = true.s.odd.ng.delta1(par=fs.ad$fit$est[fs.ad$fit$R,],fs.ad$nodd,cv=odds,data=c(s.ad)) # convoluted data
			else
				vts.ad.odd = true.s.odd.ng.delta2(par=fs.ad$fit$est[fs.ad$fit$R,],fs.ad$nodd,cv=odds,data=c(s.ad)) # convoluted data
		}else
			vts.ad.odd = true.s.odd.ne(par=fs.ad$fit$est[fs.ad$fit$R,],fs.ad$nodd,cv=odds,data=c(s.ad)) # convoluted data

		if(verbose)
			cat("# msPeakG #: Step II - Baseline removal and denoising ...\n")

		ts.ad.odd = matrix(vts.ad.odd,length(subrt1),length(subrt2)) # matrix of convoluted data

		method.lab = c("FDT","PMM","tGMM","GMM","GaMM","EGMM")

		if(verbose)
			cat("# msPeakG #: Step III - Peak peaking using ",method.lab[method]," and grouping ...\n")

		# peak finding and dgrouping
		gfs.ad.odd = mb.pd(ts.ad.odd,method,sft=range(subrt1),sst=range(subrt2))$pos
	
	########################
	########################
	########################

	gs.ad.odd = gp.rg(gfs.ad.odd)	

	if(merge){
		if(verbose)
			cat("# msPeakG #: Step IV - Peak merging with the cut-off ",cv," ...\n")

		# peak merging
		ms.ad.odd = peak.mg(gdata=gs.ad.odd$nd,data=ts.ad.odd,cv=cv,sic.name=sic.name)

		if(verbose)
			cat("# msPeakG #: Create the peak table ...\n")
		# peak table
		ptab.odd = ptable(ms.ad.odd,subt1=subrt1,subt2=subrt2,sic.name=sic.name)

		ptab1 = ptab.odd[,c(2,3,7)]
		ptab2 = ptab.odd[ptab.odd$pgroup!=0,c(2,3,8)]
	}else{
		if(verbose)
			cat("# msPeakG #: Create the peak table ...\n")
		# peak table
		ptab.odd = ptable.wo.merg(gs.ad.odd$nd,subt1=subrt1,subt2=subrt2,sic.name=sic.name)

		ptab1 = ptab.odd[,c(2,3,6)]
		ptab2 = ptab1
	}

	# output file name
		out.file = paste(out.fname,".rda",sep="")
		cat("# parameters for msPeakG\n")
		cat("# msPeakG.output.cv",cv,".zcv",odds,".m",method,".",out.fname,".",ngblabel,".max",imax,".rda\n",sep="")

	if(verbose)
		cat("# msPeakG #: Saving the result in the RData file \"",out.file,"\"\n")

	# save data for downstream analysis
	if(merge){
		save(
			sad.weight,subrt1,subrt2,cv,odds,method#,opt
			,tic.fname
			,s.ad,fs.ad
			,vts.ad.odd	,ts.ad.odd
			,gfs.ad.odd,gs.ad.odd,ms.ad.odd,ptab.odd
			,out.file,ptab1,ptab2
			,ngb
			,xarea
			,file=out.file
		)
	}else{
		save(
			sad.weight,subrt1,subrt2,cv,odds,method#,opt
			,tic.fname
			,s.ad,fs.ad
			,vts.ad.odd	,ts.ad.odd
			,gfs.ad.odd,gs.ad.odd,ptab.odd
			,out.file,ptab1,ptab2
			,ngb
			,xarea
			,file=out.file
		)
	}

	if(verbose)
		cat("# msPeakG #: =====> DONE!!! <=====\n")

	list(rda=out.file,bptable=ptab1,aptable=ptab2)
}

#########################
# Sub functions for plot #
#########################

# draw the raw data
sub.contour <- function(
				tic.fname = "ticmat_289.txt"
				,sub1=1:200
				,sub2=1:500
				,p1=6
				,p2=0.005
				,pmin1=300
				,pmin2=0.02
				,plot=T
				,col="grey"
){
	# TIC information: original: RT2 by RT1 so transpose necessary
	atic = t(as.matrix(read.table(tic.fname)))
	tic.size = dim(atic)

	r.rt1 = seq(pmin1,(pmin1+p1*(tic.size[1]-1)),by=p1)
	r.rt2 = seq(pmin2,(pmin2+p2*(tic.size[2]-1)),by=p2)

	sr.rt1 = r.rt1[sub1]
	sr.rt2 = r.rt2[sub2]
	satic = atic[sub1,sub2]

	if(plot){
		# plot the raw data
		plot(1:10,1:10,type="n"
				,xlab="First dimension retention time (s)",ylab="Second dimension reteintion time (s)"
			,xlim=range(sr.rt1),ylim=range(sr.rt2)
		)
		contour(sr.rt1,sr.rt2,satic,col=col,nlevels=15,add=T,labels="")
	}
}

# discrete boundary-like plot
peak.region <- function(data
				,add=F
				,col=2,cex=.1,pch=1,lty=2,lwd=2
				,ini
				,p1=6
				,p2=0.005
				,pmin1=300
				,pmin2=0.02
){
	d = dim(data)
	if(!add)
		plot(1:d[1],1:d[1],type="n",xlim=c(0,d[1]),ylim=c(0,d[2]),cex=cex,lwd=lwd)

	r1 <- function(x){
		(ini[1]-2+x)*p1+pmin1
	}

	r2 <- function(x){
		(ini[2]-2+x)*p2+pmin2
	}

	for(i in 1:d[1]){
		td = data[i,]
		nz.pos = which(td>0)
		len.nz = length(nz.pos)
		if(len.nz==0)
			next
		bnz.pos = which(abs(diff(nz.pos))>1)
		len.bnz = length(bnz.pos)
		if(len.bnz==0){
			segments(r1(i-0.5),r2(nz.pos[1]-0.5),r1(i+0.5),r2(nz.pos[1]-0.5),col=col,lty=lty,cex=cex,lwd=lwd)
			segments(r1(i-0.5),r2(nz.pos[len.nz]+0.5),r1(i+0.5),r2(nz.pos[len.nz]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
			segments(r1(i-0.5),r2(nz.pos[1]-0.5),r1(i-0.5),r2(nz.pos[len.nz]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
			segments(r1(i+0.5),r2(nz.pos[1]-0.5),r1(i+0.5),r2(nz.pos[len.nz]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
			next
		}		
		for(j in 1:len.bnz){
			if(j==1){
				segments(r1(i-0.5),r2(nz.pos[1]-0.5),r1(i+0.5),r2(nz.pos[1]-0.5),col=col,lty=lty,cex=cex,lwd=lwd)
				segments(r1(i-0.5),r2(nz.pos[bnz.pos[j]]+0.5),r1(i+0.5),r2(nz.pos[bnz.pos[j]]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
				segments(r1(i-0.5),r2(nz.pos[1]-0.5),r1(i-0.5),r2(nz.pos[bnz.pos[j]]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
				segments(r1(i+0.5),r2(nz.pos[1]-0.5),r1(i+0.5),r2(nz.pos[bnz.pos[j]]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
			}else{
				segments(r1(i-0.5),r2(nz.pos[(bnz.pos[j-1]+1)]-0.5),r1(i+0.5),r2(nz.pos[(bnz.pos[j-1]+1)]-0.5),col=col,lty=lty,cex=cex,lwd=lwd)
				segments(r1(i-0.5),r2(nz.pos[bnz.pos[j]]+0.5),r1(i+0.5),r2(nz.pos[bnz.pos[j]]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
				segments(r1(i-0.5),r2(nz.pos[(bnz.pos[j-1]+1)]-0.5),r1(i-0.5),r2(nz.pos[bnz.pos[j]]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
				segments(r1(i+0.5),r2(nz.pos[(bnz.pos[j-1]+1)]-0.5),r1(i+0.5),r2(nz.pos[bnz.pos[j]]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
			}
							
		}		
		segments(r1(i-0.5),r2(nz.pos[(bnz.pos[len.bnz]+1)]-0.5),r1(i+0.5),r2(nz.pos[(bnz.pos[len.bnz]+1)]-0.5),col=col,lty=lty,cex=cex,lwd=lwd)
		segments(r1(i-0.5),r2(nz.pos[len.nz]+0.5),r1(i+0.5),r2(nz.pos[len.nz]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
		segments(r1(i-0.5),r2(nz.pos[(bnz.pos[len.bnz]+1)]-0.5),r1(i-0.5),r2(nz.pos[len.nz]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
		segments(r1(i+0.5),r2(nz.pos[(bnz.pos[len.bnz]+1)]-0.5),r1(i+0.5),r2(nz.pos[len.nz]+0.5),col=col,lty=lty,cex=cex,lwd=lwd)
	}
}

# plot the detected peaks
pplotG <- function(x,merge=F,...){
	# data: the result from msPeak
	if(!is.list(x)){
		stop("x must be the outcome of msPeak\n")
	}
	if(merge){
		xd = x$aptable$t1
		yd = x$aptable$t2
	}else{
		xd = x$bptable$t1
		yd = x$bptable$t2
	}
	plot(xd,yd,xlab="First dimension retention time (s)"
				,ylab="Second dimension retention time (s)"
				,...
	)
}

# plot the detected peaks along with the selected sub-regions and the detected peak regions
tplotG <- function(x
			,merge=F # plot the merged peaks
			,col="purple" # the color of peaks
			,cex=1 # the size of the peak points
			,lwd=1 # the thickness of the lines
			,lpos="topleft" # the position of the legend
			,con.col="grey" # the color of the contour plot
){
	# data: the result from msPeak
	if(!is.list(x)){
		stop("x must be the outcome of msPeak\n")
	}

	cat("# CAUTION!!! # Make sure that ",x$rda," is in the same directory ...\n")

	load(x$rda)

	zcv = odds
	p1 = xarea$p1
	p2 = xarea$p2
	pmin1 = xarea$pmin1
	pmin2 = xarea$pmin2

	vmethod = c("fs","ps","tgs","gs","gas","egs")
	leg.met = c("FDT","PMM","tGMM","GMM","GaMM","EGMM","OPT-MSE","OPT-AIC","OPT-BIC")

	# var = 1:y, 2:z
	tdata = get("ts.ad.odd")
	sub.contour(tic.fname=tic.fname,sub1=subrt1,sub2=subrt2,col=con.col,p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2)
	tslab = paste("Peak region (cv=",zcv,")",sep="")
	peak.region(tdata,cex=cex,col=2,add=T,lwd=lwd,ini=c(subrt1[1],subrt2[1]),p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2)
	if(merge){
		mdata = get("ms.ad.odd")
		tp1 = (subrt1[1]-2+mdata[mdata[,12]!=0,3])*p1+pmin1
		tp2 = (subrt2[1]-2+mdata[mdata[,12]!=0,9])*p2+pmin2
		points(tp1,tp2,col=col,pch=18,cex=cex)
		pnum = dim(mdata[mdata[,12]!=0,])[1]
	}else{
		gdata = get("gs.ad.odd")
		gdata = gdata$nd
		tp1 = (subrt1[1]-2+gdata[,3])*p1+pmin1
		tp2 = (subrt2[1]-2+gdata[,9])*p2+pmin2
		points(tp1,tp2,col=col,pch=18,cex=cex)
		pnum = dim(gdata)[1]
	}
	vleg = c("TIC",tslab,paste(leg.met[method]," (",pnum,")",sep=""))
	legend(lpos,vleg
			,col=c(con.col,"red",col)
			,pch=c(-1,-1,18),lty=c(1,2,-1)
			,bty="n"
		)
}

###############
# extract info
###############
# estimate
estimate <- function(rda){
	load(rda)
	fit = fs.ad$fit$estimate[50,]
	fit
}

# compare histogram
ms.signal <- function(rda){
	load(rda)
	signal = s.ad
	tsignal = ts.ad.odd
	yr = range(c(c(signal),c(tsignal)))
	plot(c(signal),xlab="Scan",ylab="Intensity",cex=.5,ylim=yr)
	points(c(tsignal),cex=.5,col=2)
}

# compare histogram
ms2.signal <- function(rda1,rda2){
	load(rda1)
	signal = s.ad
	tsignal1 = ts.ad.odd
	load(rda2)
	#signal = s.ad
	tsignal2 = ts.ad.odd
	yr = range(c(c(signal),c(tsignal1),c(tsignal2)))
	par(mfrow=c(2,2))
	plot(c(signal),xlab="Scan",ylab="Intensity",cex=.5,ylim=yr)
	points(c(tsignal1),cex=.5,col=2)
	points(c(tsignal2),cex=.5,col=4)

	plot(c(signal),c(tsignal1),cex=.5,xlim=yr,ylim=yr,xlab="Signal",ylab="Converted")
	abline(0,1,col=2)
	plot(c(signal),c(tsignal2),cex=.5,xlim=yr,ylim=yr)
	abline(0,1,col=2)
	plot(c(tsignal1),c(tsignal2),cex=.5,xlim=yr,ylim=yr)
	abline(0,1,col=2)
}

# compare histogram
ms.hist <- function(rda,gam=TRUE){
	load(rda)
	signal = c(s.ad)
	tsignal = c(ts.ad.odd)
	
	tns = signal
	tns = tns[tsignal!=0]

	tnz = tsignal
	tnz = tnz[tsignal!=0]

	h1 = histogram(signal,type='irregular',verbose=FALSE,plot=FALSE)
	h2 = histogram(tns,type='irregular',verbose=FALSE,plot=FALSE)
	h3 = histogram(tsignal,type='irregular',verbose=FALSE,plot=FALSE)
	h4 = histogram(tnz,type='irregular',verbose=FALSE,plot=FALSE)
	
	xr = c(0,2)

	est = estimate(rda=rda)

	n = 10000
	if(gam){
		yy = rep(0,n)
		for(i in 1:n){
			tx = rgamma(1,est[3],1/est[4])
			ty = rnorm(1,(tx+est[1]),sqrt(est[2]))
			yy[i] = ty
		}
	}else{
		yy = rep(0,n)
		for(i in 1:n){
			tx = rexp(1,est[3])
			ty = rnorm(1,(tx+est[1]),sqrt(est[2]))
			yy[i] = ty
		}
	}		
	plot(h1,xlim=xr)
	points(density(yy),type="l",col=2)
	if(gam){
		F1 = dnormgam(est[1:4],plot=FALSE)
		lines(F1$xout,F1$dout,col="blue")
	}
}

# compare histogram
ms.hist2 <- function(rda,gam=TRUE){
	load(rda)
	signal = c(s.ad)
	tsignal = c(ts.ad.odd)
	par(mfrow=c(2,2))
	
	tns = signal
	tns = tns[tsignal!=0]

	tnz = tsignal
	tnz = tnz[tsignal!=0]

	h1 = histogram(signal,type='irregular',verbose=FALSE)
	h2 = histogram(tns,type='irregular',verbose=FALSE)
	h3 = histogram(tsignal,type='irregular',verbose=FALSE)
	h4 = histogram(tnz,type='irregular',verbose=FALSE)
	
	xr = c(0,2)

	plot(h1,xlim=xr)
	plot(h2,xlim=xr)

	plot(h4,xlim=xr)
	
	est = estimate(rda=rda)
	if(gam){
		ft <- function(x,est){
			dgamma(x,est[3],1/est[4])
		}
	}else{
		ft <- function(x,est){
			dexp(x,est[3])
		}
	}		
	
	xx = seq(xr[1],xr[2],length=100)
	yy = rep(0,length(xx))
	for(i in 1:length(xx)){
		yy[i] = ft(xx[i],est=est)
	}
	points(xx,yy,type="l",col=2)
	plot(density(tnz),xlim=xr)
	points(xx,yy,type="l",col=2)

}

gen.data <- function(par=c(0.5,0.0005,1,0.5),rn=1200,cn=289,p=.1){
	n = rn*cn
	y = rep(0,n)
	s = rep(0,n)
	for(i in 1:n){
		uv = runif(1)
		if(uv<p){
			tx = rgamma(1,par[3],1/par[4])
			s[i] = 1
		}else{
			tx = 0
		}
		ty = rnorm(1,(tx+par[1]),par[2])
		y[i] = ty
	}
	y = matrix(y,rn,cn)
	list(par=par,data=y,s=s)
}

gen2.data <- function(par=c(0.5,0.0005,1,0.5),rn=1200,cn=289,p=.1){
	n = rn*cn
	y = rep(0,n)
	s = rep(0,n)
	for(i in 1:n){
		uv = runif(1)
		ty = rnorm(1,par[1],par[2])
		if(uv<p){
			tx = rgamma(1,par[3],1/par[4])
			s[i] = 1
		}else{
			tx = 0
		}
		y[i] = ty + tx
	}
	y = matrix(y,rn,cn)
	list(par=par,data=y,s=s)
}

gen.tdata <- function(par=c(0.5,0.0002,.5,0.9),rn=1200,cn=289){
	n = rn*cn
	y = rep(0,n)
	for(i in 1:n){
		tx = rgamma(1,par[3],1/par[4])
		ty = rnorm(1,(tx+par[1]),par[2])
		y[i] = ty
	}
	list(par=par,data=y)
}

hist.eg <- function(fit,X,cv=.5,up=1){
	par = as.numeric(fit$fit$est[fit$fit$R,])
	par = par[-length(par)]
	z = fit$fit$z
	X = c(unlist(X))
	XS = X[z>cv]	
	NS = X[z<=cv]
	nx = seq(0,up,length=100)
	CS = normgam.signal(X,par)

	par(mfrow=c(2,2))
	Td = dnormgam(par, plot=FALSE)
	Fd = dnorm(nx,par[1],par[2], log=FALSE)
	H = histogram(X, type='irregular', verbose=FALSE, plot=FALSE)
	Hs = histogram(XS, type='irregular', verbose=FALSE, plot=FALSE)
	Ns = histogram(NS, type='irregular', verbose=FALSE, plot=FALSE)
	plot(H, xlim=c(0,up))	
	lines(Td$xout, Td$dout, col='red',lwd=.5)
	lines(nx, Fd, col='blue',lwd=1.5)
	plot(Hs, xlim=c(0,up))	
	lines(Td$xout, Td$dout, col='red',lwd=.5)
	plot(Ns, xlim=c(0,up))	
	lines(nx, Fd, col='red')
	plot(X,CS)
	abline(0,1,col=2)
}

sub.contourG <- function(
                                tic.fname = "ticmat_289.txt"
                                ,sub1=1:200
                                ,sub2=1:500
                                ,p1=6
                                ,p2=0.005
                                ,pmin1=300
                                ,pmin2=0.02
                                ,plot=T
                                ,col="grey"
					  ,vlevel=15
){
        # TIC information: original: RT2 by RT1 so transpose necessary
        atic = t(as.matrix(read.table(tic.fname)))
        tic.size = dim(atic)

        r.rt1 = seq(pmin1,(pmin1+p1*(tic.size[1]-1)),by=p1)
        r.rt2 = seq(pmin2,(pmin2+p2*(tic.size[2]-1)),by=p2)

        sr.rt1 = r.rt1[sub1]
        sr.rt2 = r.rt2[sub2]
        satic = atic[sub1,sub2]

        if(plot){
                # plot the raw data
                plot(1:10,1:10,type="n"
                                ,xlab="First dimension retention time (s)",ylab="Second dimension reteintion time (s)"
                        ,xlim=range(sr.rt1),ylim=range(sr.rt2)
                )
                contour(sr.rt1,sr.rt2,satic,col=col,nlevels=vlevel,add=T,labels="")
        }
}

find.subG <- function(
                                subrt1=31:70 # indices of the selected sub-region's 1st dimension rt
                                ,subrt2=401:500 # indices of the selected sub-region's 2nd dimension rt
                                ,tic.fname = "ticmat_289.txt" # TIC file name
                                ,col="grey"
                                ,p1=6 # modulation time of 1st RT
                                ,p2=0.005 # modulation time of 2nd RT
                                ,pmin1=300 # minimum of 1st RT
                                ,pmin2=0.02 # minimum of 2nd RT
                                ,plot=F
					  ,vlevel=15
){
        atic = t(as.matrix(read.table(tic.fname)))
        tic.size = dim(atic)

        # whole data
        xr = (c(range(subrt1))-1)*p1+pmin1
        yr = (c(range(subrt2))-1)*p2+pmin2
        if(plot){
                par(mfrow=c(2,1))
                sub.contourG(tic.fname=tic.fname,1:tic.size[1],1:tic.size[2],col=col,p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2,vlevel=vlevel)
                title("Entire Chromatogram")
                segments(xr[1],yr[1],xr[2],yr[1],col=3,lwd=2)
                segments(xr[1],yr[1],xr[1],yr[2],col=3,lwd=2)
                segments(xr[1],yr[2],xr[2],yr[2],col=3,lwd=2)
                segments(xr[2],yr[1],xr[2],yr[2],col=3,lwd=2)
                sub.contourG(tic.fname=tic.fname,subrt1,subrt2,col=col,p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2,vlevel=vlevel)
                title("Magnified Sub-chromatogram")
        }

        list(idxt1=subrt1,idxt2=subrt2,range.t1=xr,range.t2=yr,tic.fname=tic.fname,p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2)
}


