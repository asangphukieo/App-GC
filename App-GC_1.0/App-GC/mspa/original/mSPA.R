################################
# LOAD the ChromaToF DATA	 #
# 1: MIXTURE (STANDARD) (#=10) #
# 2: SPIKE-IN (DEBLANK) (#=5)  #
################################

# import the data
rt.import <- function(pattern="^STANDARD*"){
	fname = dir(pattern=pattern)
	rt = NULL
	lapply(1:length(fname),function(i) read.csv(fname[i]))
}

# extract the spec info
spec <- function(data=rt){
	len = length(data)
	nd = NULL
	nd = lapply(1:len,
				function(i){
					sp=data[[i]]$Spectra;
					strsplit(as.character(sp),split=" ");
				}
			)
	sd = NULL
	sd = lapply(1:len,
				function(i){
					td=nd[[i]];
					len2=length(td);
					lapply(1:len2,
							function(j){
								len3=length(td[[j]]);
								td2=t(matrix(as.numeric(unlist(strsplit(td[[j]],split=":"))),2,len3));
								td2[order(td2[,1]),];
							}
						)
				}
			)
	sd
}	

# load the ChromaToF data
load.data <- function(pattern){

	rt = rt.import(pattern=pattern) # real data rt[[1]][,]
	nrt = spec(rt) # spectrum nrt[[1]][[1]][,]
	
	# 1   2   3     4     5      6-
	# t1, t2, area, name, simil, spect...,

	trt = NULL
	for(i in 1:length(nrt)){
		tspec = NULL
		for(j in 1:length(nrt[[i]])){
			tspec = rbind(tspec,as.numeric(nrt[[i]][[j]][,2]))
		}
		trt[[i]] = data.frame(
				t1 = (rt[[i]]$X1)
				,t2 = (rt[[i]]$X2)
				,area = rt[[i]]$Area
				,name = rt[[i]]$Name
				,simil = rt[[i]]$Simil
				,tspec # spectrum # 6---
			)
	}
	for(i in 1:length(trt)){
		trt[[i]] = trt[[i]][match(unique(trt[[i]]$area),trt[[i]]$area),]
	}

	# multiple peak merging
	urt = list()
	for(i in 1:length(trt)){
		utmp = trt[[i]]
		total.name = as.character(unique(as.character(utmp$name)))
		n.len = length(total.name)
		for(j in 1:n.len){
			upos = match(as.character(utmp$name),as.character(total.name[j]))
			upos = which(!is.na(upos))
			if(length(upos)>1){
				sutmp = utmp[upos,]
				tpos = which(sutmp$area==max(sutmp$area))[1]
				dpos = upos[-tpos]
				utmp = utmp[-dpos,]
				utmp$name = factor(utmp$name)
			}
		}
		urt[[i]] = utmp
	}
	
	list(trt,urt)
}

# load the data
#pad = load.data(pattern="^Standard_*")
#pad = load.data(pattern="^Deblank_*")


###################
#  PEAL ALIGNMENT #
###################

# similarity product
sim.pd <- function(d1,d2,ms.pos,sm="dot"){
	d1 = d1[,ms.pos]
	d2 = d2[,ms.pos]
	if(sm=="dot"){
		dd1 = as.matrix(d1)/sqrt(apply((as.matrix(d1))^2,1,sum))
		dd2 = as.matrix(d2)/sqrt(apply((as.matrix(d2))^2,1,sum))
		sim = dd1 %*% t(dd2)
	}else if(sm=="pearson"){
		sim = cor(t(d1),t(d2))
	}
}

# distance from R
dist.r <- function(d1,d2,dm="canberra"){
	# dm = "euclidean", "maximum", "manhattan", "canberra"
	
	d1 = d1
	d2 = d2
	len1 = dim(d1)[1]
	len2 = dim(d2)[1]

	da = rbind(d1,d2)
	da = as.matrix(da)
	tdist = as.matrix(dist(da,method=dm,diag=TRUE,upper=TRUE))[(1:len1),((len1+1):(len1+len2))]
	tdist
}

# accuracy calculation
cscore <- function(d1,d2){
	d1.name = as.character(d1$name)
	d2.name = as.character(d2$name)
	total.pos = length(na.omit(match(d1.name,d2.name)))
	total.neg = length(d1.name)*length(d2.name) - total.pos
	
	d1 = d1[d1$nflag!=0,] # ref
	d2 = d2[d2$nflag!=0,] # tar
	d1 = d1[order(d1$nflag),]
	d2 = d2[order(d2$nflag),]
	d1$name = as.character(d1$name)
	d2$name = as.character(d2$name)
	total.match = dim(d1)[1]
	t.p = 0
	for(i in 1:total.match){
		if(d1$name[i] == d2$name[i]){
			t.p = t.p + 1
		}	
	}
	f.p = total.match - t.p
	f.n = total.pos - t.p
	t.n = total.neg - f.p
	t.p.r = t.p/total.pos
	p.p.v = t.p/(t.p+f.p)
	f1 = 2*t.p.r*p.p.v/(t.p.r+p.p.v)
		
	rlt=c(t.p.r,p.p.v,t.p,f.p,f.n,t.n,f1)
	rlt
}

# Peak alignment
mspa <- function(
		padata = pad
		,rid = 1
		,tid = 2
		,method="PMWO"
		,w = .5
		,k = 5
		,rho = .95
		,distance="canberra"
		,similarity="dot"
		,vn = c("t1","t2")
		,nseed = 3927
		,sp.pos=NA
		,anal=FALSE
		,plot=FALSE
	){
	# padata: peak alignment data; list
	# rid: index of the reference sample
	# tid: index of the target sample
	# distance (dm): "euclidean", "maximum", "manhattan", "canberra"
	# similarity (sm): "dot", "pearson"
	# method (pm): peak alignment method
	#   =1: Only Distance (PAD)
	#   =2: Only Similarity (PAS)
	#   =3: Distance window => Similarity (DW-PAS)
	#   =4: Similarity window => Distance (SW-PAD)
	#   =5: Mixture similarity (PAM)
	#   =6: Optimal version of Method 5 (OP-PAM)
	# anal: T if both data have the name information. i.e., identification information
	# sp.pos: positions of mass spectra in the dataset
	# vn: positions of RT1 and RT2
	# k (kdist): cutoff for the distance-based window (for pm=3)
	# rho (cv): cutoff for the similarity-based window (for pm=4)
	# w: weight for distance (for pm=5) and the initial weight for pm=6

	id1 = rid
	id2 = tid
	cv = rho
	kdist = k
	sm = similarity
	dm = distance
	if(method=="PAD"){
		pm = 1
	}else if(method=="PAS"){
		pm = 2
	}else if(method=="DW-PAS"){
		pm = 3
	}else if(method=="SW-PAD"){
		pm = 4
	}else if(method=="PAM"){
		pm = 5
	}else if(method=="OP-PAM"){
		pm = 6
	}
	w = as.numeric(w)

	set.seed(nseed)
	
	refd = padata[[id1]]
	tard = padata[[id2]]

	if(is.na(sp.pos)){
		sp.pos = sort(grep("X",colnames(refd)))
		sp.pos = sp.pos[1]:sp.pos[length(sp.pos)]
	}

	# index
	refd$ndx = c(1:dim(refd)[1])
	tard$ndx = c(1:dim(tard)[1])

	# flag
	refd$nflag = rep(0,dim(refd)[1])
	tard$nflag = rep(0,dim(tard)[1])

	# correlation
	refd$cor = rep(-2,dim(refd)[1])

	# distance
	refd$dist = refd$cor

	# underlying similarity (pm=1,4: distance; others: similarity)
	refd$sim = refd$cor

	# number of compounds
	lenr = dim(refd)[1]
	lent = dim(tard)[1]

	lntrans <- function(vdata){
		a = min(vdata)
		b = max(vdata)
		ndata = vdata/(b-a)-a/(b-a)
		ndata
	}

	# sort by RT1 and RT2
	refd = refd[order(refd[,vn[1]],refd[,vn[2]]),] # library
	tard = tard[order(tard[,vn[1]],tard[,vn[2]]),] # query

	# distance
	avec.dist = NULL
	avec.dist = dist.r(tard[,vn],refd[,vn],dm=dm)
	avec.dist = 1/(1+avec.dist)
	
	# spectra similarity
	amat.cor = NULL
	amat.cor = sim.pd(tard[,],refd[,],ms.pos=sp.pos,sm=sm)

	# search all the samples given a reference
	if(pm!=6){
		for(i in 1:lent){		

			sim2 = -2

			vec.dist = c(avec.dist[i,])
			mat.cor = c(amat.cor[i,])

			if(pm==1){ # distance
				k.dist.pos = 1:length(vec.dist)
				mat.cor.dist = vec.dist
				max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
				max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
				max.cor.k.pos = k.dist.pos[max.cor.pos]
			}else if(pm==2){ # similarity
				cv.cor.pos = 1:length(mat.cor)
				mat.cor.dist = mat.cor
				max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
				max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
				max.cor.k.pos = cv.cor.pos[max.cor.pos]
			}else if(pm==3){ # distance-window
				if(kdist<=-1){ # i.e.: pm=1
					k.dist.pos = 1:length(vec.dist)
				}else if(kdist<dim(tard)[1]){
					k.dist.pos = which(vec.dist>=sort(vec.dist,decreasing=T)[kdist])		
				}else{
					k.dist.pos = 1:length(vec.dist)
				}

				smat.cor = mat.cor[k.dist.pos]
				cv.cor.pos = 1:length(smat.cor)
				mat.cor.dist = smat.cor[cv.cor.pos]	
				max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
				max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
				max.cor.k.pos = k.dist.pos[cv.cor.pos[max.cor.pos]]
			}else if(pm==4){ # similarity-window
				if(cv<=-1){ # i.e.: pm=1
					cv.cor.pos = 1:length(mat.cor)
				}else{
					cv.cor.pos = which(mat.cor>=cv)
				}
				if(length(cv.cor.pos)==0)
					next

				svec.dist = vec.dist[cv.cor.pos]
				k.dist.pos = 1:length(svec.dist)
				mat.cor.dist = svec.dist[k.dist.pos]	
				max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
				max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
				max.cor.k.pos = cv.cor.pos[k.dist.pos[max.cor.pos]]
			}else if(pm==5){ # mixture
				k.dist.pos = 1:length(vec.dist)
				mat.cor.dist = w*vec.dist + (1-w)*mat.cor
				max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
				max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
				max.cor.k.pos = k.dist.pos[max.cor.pos]
			}

			sim2 = mat.cor.dist[max.cor.pos]

			if(refd$nflag[max.cor.k.pos]!=0){
				sim1 = refd$sim[max.cor.k.pos]
				if(sim2>sim1){
					tard$nflag[tard$ndx==refd$nflag[max.cor.k.pos]] = 0
					refd$nflag[max.cor.k.pos] = tard$ndx[i]
					refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
					refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
					refd$sim[max.cor.k.pos] = sim2
					tard$nflag[i] = tard$ndx[i]
				}
			}else{
				refd$nflag[max.cor.k.pos] = tard$ndx[i]
				refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
				refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
				refd$sim[max.cor.k.pos] = sim2
				tard$nflag[i] = tard$ndx[i]
			}				
		}
	}else{

		fun <- function(para,edm="canberra"){
			avec.dist = NULL
			avec.dist = dist.r(tard[,vn],refd[,vn],dm=edm)
			avec.dist = 1/(1+avec.dist)

			irefd = refd
			itard = tard
			
			for(i in 1:lent){
				sim2 = -2
	
				vec.dist = c(avec.dist[i,])
				mat.cor = c(amat.cor[i,])
	
				k.dist.pos = 1:length(vec.dist)
				mat.cor.dist = para*vec.dist + (1-para)*mat.cor
				max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
				max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
				max.cor.k.pos = k.dist.pos[max.cor.pos]
			
				sim2 = mat.cor.dist[max.cor.pos]

				if(refd$nflag[max.cor.k.pos]!=0){
					sim1 = refd$sim[max.cor.k.pos]
					if(sim2>sim1){
						itard$nflag[itard$ndx==irefd$nflag[max.cor.k.pos]] = 0
						irefd$nflag[max.cor.k.pos] = itard$ndx[i]
						irefd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
						irefd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
						irefd$sim[max.cor.k.pos] = sim2
						itard$nflag[i] = itard$ndx[i]
					}
				}else{
					irefd$nflag[max.cor.k.pos] = itard$ndx[i]
					irefd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
					irefd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
					irefd$sim[max.cor.k.pos] = sim2
					itard$nflag[i] = itard$ndx[i]
				}				
			}
			all.sim = irefd$sim[irefd$nflag!=0]
			all.cor = irefd$cor[irefd$nflag!=0]
			all.dist = irefd$dist[irefd$nflag!=0]

			s.sim = sum(all.sim)
			s.cor = sum(all.cor)
			s.dist = sum(all.dist)
			s.mix = s.sim + s.cor + s.dist
			ll = s.mix
			
			-ll
		}
		
		w1 = nlminb(w,fun,lower=0,upper=1,edm="euclidean")
		#cat("    *** euclidean = ",w1$par, "(",w1$obj,")\n")
		w2 = nlminb(w,fun,lower=0,upper=1,edm="maximum")
		#cat("    *** maximum = ",w2$par, "(",w2$obj,")\n")
		w3 = nlminb(w,fun,lower=0,upper=1,edm="manhattan")
		#cat("    *** manhattan = ",w3$par, "(",w3$obj,")\n")
		w4 = nlminb(w,fun,lower=0,upper=1,edm="canberra")
		#cat("    *** canberra = ",w4$par, "(",w4$obj,")\n")

		objs = c(w1$objective,w2$objective,w3$objective,w4$objective)
		ws = c(w1$par,w2$par,w3$par,w4$par)

		mim.pos = which(objs == min(objs))[1]
		w = ws[mim.pos]

		if(mim.pos==1){
			edm = "euclidean"
		}else if(mim.pos==2){
			edm = "maximum"
		}else if(mim.pos==3){
			edm = "manhattan"
		}else{
			edm = "canberra"
		}

		#cat("     *** ew=",w," at ",edm,"(",mim.pos,")\n")
		
		avec.dist = NULL
		avec.dist = dist.r(tard[,vn],refd[,vn],dm=edm)
		avec.dist = 1/(1+avec.dist)
		
		for(i in 1:lent){		

			sim2 = -2

			vec.dist = c(avec.dist[i,])
			mat.cor = c(amat.cor[i,])

			k.dist.pos = 1:length(vec.dist)
			mat.cor.dist = w*vec.dist + (1-w)*mat.cor
			max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
			max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
			max.cor.k.pos = k.dist.pos[max.cor.pos]

			sim2 = mat.cor.dist[max.cor.pos]

			if(refd$nflag[max.cor.k.pos]!=0){
				sim1 = refd$sim[max.cor.k.pos]
				if(sim2>sim1){
					tard$nflag[tard$ndx==refd$nflag[max.cor.k.pos]] = 0
					refd$nflag[max.cor.k.pos] = tard$ndx[i]
					refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
					refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
					refd$sim[max.cor.k.pos] = sim2
					tard$nflag[i] = tard$ndx[i]
				}
			}else{
				refd$nflag[max.cor.k.pos] = tard$ndx[i]
				refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
				refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
				refd$sim[max.cor.k.pos] = sim2
				tard$nflag[i] = tard$ndx[i]
			}				
		}			
	}
	
	all.sim = refd$sim[refd$nflag!=0]
	all.cor = refd$cor[refd$nflag!=0]
	all.dist = refd$dist[refd$nflag!=0]

	s.sim = sum(all.sim)
	s.cor = sum(all.cor)
	s.dist = sum(all.dist)
		
	s.mix = s.sim + s.cor + s.dist
	ll = -s.mix
			
	
	refr = refd[refd$nflag!=0,]
	tarr = tard[tard$nflag!=0,]
	refr = refr[order(refr$nflag),]
	tarr = tarr[order(tarr$nflag),]

	atable = data.frame(refr$ndx,tarr$ndx)
	dimnames(atable)[[2]] = paste("S",c(id1,id2),sep="")

	if(plot){
		plot(refd[refd$nflag==0,vn[1]],refd[refd$nflag==0,vn[2]],cex=.5,pch=19,col=3
			,xlim=c(min(refd[,vn[1]],tard[,vn[1]]),max(refd[,vn[1]],tard[,vn[1]]))
			,ylim=c(min(refd[,vn[2]],tard[,vn[2]]),max(refd[,vn[2]],tard[,vn[2]]))
			,xlab="RT1",ylab="RT2"
			)
		points(tard[tard$nflag==0,vn[1]],tard[tard$nflag==0,vn[2]],col=5,pch=19,cex=.5)

		points(refr[,vn[1]],refr[,vn[2]],cex=1,pch=3,col=1)
		points(tarr[,vn[1]],tarr[,vn[2]],col=2,pch=4,cex=1)
		
		for(i in 1:dim(refr)[1]){
			segments(refr[i,vn[1]],refr[i,vn[2]],tarr[i,vn[1]],tarr[i,vn[2]],col=4,lwd=.5)
		}
		legend("topleft",c("ref:nonmatch","tar:nonmatch","ref:match","tar:match","pair")
				,lty=c(-1,-1,-1,-1,1),col=c(3,5,1,2,4),pch=c(19,19,3,4,-1)
				,bty="n"
			)
	}

	if(anal){
		accur = cscore(refd,tard)
		summary = c(id1,id2,pm,dm,kdist,sm,cv,w,as.numeric(accur),ll)
		summary = data.frame(t(summary))
		dimnames(summary)[[2]] = c(	
					"id1","id2"
					,"method","distance","k","similarity","cv","w"
					,"tpr","ppv","tp","fp","fn","tn"
					,"f1","ll"
				)
	}else{
		summary = c(id1,id2,pm,dm,kdist,sm,cv,w,ll)
		summary = data.frame(t(summary))
		dimnames(summary)[[2]] = c(	
					"id1","id2"
					,"method","distance","k","similarity","cv","w"
					,"ll"
				)	}
	
	list(table=atable,summary=summary)
}

# alingment of all peak lists
align <- function(padata,alist,method="peak",plot=FALSE){
	# alist: it should be arranged.
	# e.g.: (1,2), (2,3), (3,4), (4,5), ...
	# Its order sould be consistent with that of padata
	# method: peak, name, area
	
	len = length(alist)
	wtable = alist[[1]]
	for(i in 2:len){
		ta = alist[[i]]
		if(colnames(ta)[1]!=colnames(alist[[i-1]])[2]){
			stop("!!! Please make sure the order of alist !!!\n")
		}
		t.pos = match(wtable[,i],ta[,1])
		wtable = data.frame(wtable[which(!is.na(t.pos)),],ta[na.omit(t.pos),2])
		dimnames(wtable)[[2]][i+1] = colnames(ta)[2]
	}
	wlen = dim(wtable)[2]
	if(method=="peak"){
		rtable = list()
		rtable$t1 = matrix(0,dim(wtable)[1],dim(wtable)[2])
		rtable$t2 = rtable$t1
		for(i in 1:wlen){
			tmp.wt = wtable[,i]
			rtable$t1[,i] = padata[[i]]$t1[tmp.wt]
			rtable$t2[,i] = padata[[i]]$t2[tmp.wt]
		}
	}else if(method=="name"){
		rtable = matrix("",dim(wtable)[1],dim(wtable)[2])
		for(i in 1:wlen){
			tmp.wt = wtable[,i]
			rtable[,i] = as.character(padata[[i]]$name[tmp.wt])
		}
	}else if(method=="area"){
		rtable = matrix(0,dim(wtable)[1],dim(wtable)[2])
		for(i in 1:wlen){
			tmp.wt = wtable[,i]
			rtable[,i] = padata[[i]]$area[tmp.wt]
		}
	}
	if(plot){
		if(method=="peak"){
			plot(1:10,1:10,type="n",xlim=range(rtable$t1),ylim=range(rtable$t2),xlab="RT1",ylab="RT2")
			for(i in 1:wlen){
				points(rtable$t1[,i],rtable$t2[,i],col=i,pch=i)
			}
		}else if(method=="area"){
			plot(1:10,1:10,type="n",xlim=c(1,dim(rtable)[2]),ylim=range(rtable),xlab="Peak",ylab="Peak Area")
			for(i in 1:dim(rtable)[1]){
				points(1:wlen,rtable[i,],col=i,pch=i,lty=i,type="o")
			}
		}
	}
		
			
	list(pos=wtable,info=rtable)
}

# PPV calculation based on the name table			
comp.name <- function(ntable){
	nt = factor(ntable)
	len = dim(ntable)
	total = len[1]
	correct = 0
	for(i in 1:len[1]){	
		flag = 1
		for(j in 2:len[2]){
			if(ntable[i,j-1]!=ntable[i,j]){
				flag=0
				break
			}
		}
		if(flag==1){
			correct = correct + 1
		}
	}
	c(TP=correct,P=total,PPV=correct/total)
}
