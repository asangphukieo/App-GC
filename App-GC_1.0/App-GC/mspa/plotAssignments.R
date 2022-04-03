
plotAssignments <- function(file, rid, tid) {
	pw <- read.table(file, header=T, sep="\t")
	d1 <- pw[,c(1,10,12,17)]
	names(d1) <- c("name", "t1", "t2", "sim")
	d2 <- pw[,c(1,11,13,17)]
	names(d2) <- c("name", "t1", "t2", "sim")

	cscore2(rid, tid, d1, d2)
}

isTP <- function(d1, d2, i, drt1, drt2) {
	v <- d2$t1[i] - d1$t1[i]
	w <- d2$t2[i] - d1$t2[i]
	s <- sum(drt1$y[c(tail(which(drt1$x<=v),n=1), head(which(drt1$x>v),n=1))])/sum(drt1$y)
	t <- sum(drt2$y[c(tail(which(drt2$x<=w),n=1), head(which(drt2$x>w),n=1))])/sum(drt2$y)
	#if (d1$name[i] == d2$name[i]) {
	#	print(paste(v,";",w,"->", s*t, sep=" "))
	#}
	s * t
}

# accuracy calculation
cscore2 <- function(rid,tid,d1,d2){
	cat("##################################","\n")
	cat("function cscore","\n")
	##convert levels to strings
	d1.name = as.character(d1$name)
	d2.name = as.character(d2$name)
	##count number of possible true matching names between d1 and d2
	total.pos = length(na.omit(match(d1.name,d2.name)))
	##all possible combinations - possible positives
	total.neg = length(d1.name)*length(d2.name) - total.pos

	cat("total.pos: ",total.pos,"\n")
	cat("total.neg: ",total.neg,"\n")
	cat("##################################","\n")	
	
	minx <- min(d1$t1, d2$t1)
	maxx <- max(d1$t1, d2$t1)
	miny <- min(d1$t2, d2$t2)
	maxy <- max(d1$t2, d2$t2)
        
	unmatchedRef = d1[d1$nflag==0,]
	unmatchedTar = d2[d2$nflag==0,]

	#d1 = d1[d1$nflag!=0,] # ref
	#d2 = d2[d2$nflag!=0,] # tar
	#d1 = d1[order(d1$nflag),]
	#d2 = d2[order(d2$nflag),]

	d1$name = as.character(d1$name)
	d2$name = as.character(d2$name)
	##number of matches between d1 and d2
	total.match = dim(d1)[1]
	cat("total.match: ",total.match,"\n")
	tpcols <- rep(2,total.match)

	diffs <- matrix(NA, total.match, 6)
	
	for(i in 1:total.match){
		diffs[i,5] <- d2$t1[i] - d1$t1[i]
                diffs[i,6] <- d2$t2[i] - d1$t2[i]
	}
	ds <- list()
	ds[[5]] <- density(na.omit(diffs[,5]), adjust=0.5)
	ds[[6]] <- density(na.omit(diffs[,6]), adjust=1.5)

	pdf(file=paste("rtMatchings_",rid,"-",tid,".pdf",sep=""), width=10, height=10)
	plot(d1$t1,d1$t2,type="n",xlab="RT1",ylab="RT2", ylim=c(miny,maxy), xlim=c(minx,maxx))
	pchcex <- 0.75

	t.p = 0
	ses <- matrix(NA, total.match, 1)
	for(i in 1:total.match){
		##which matches are true positive matches
		s <- isTP(d1,d2,i,ds[[5]], ds[[6]])
		ses[i,1] <- s
		if(s > 0.00005){
			t.p = t.p + 1
			tpcols[i] <- 3
			lines(c(d1$t1[i],d2$t1[i]),c(d1$t2[i],d2$t2[i]),col=3,lt=1)
			diffs[i,1] <- d2$t1[i] - d1$t1[i]
			diffs[i,2] <- d2$t2[i] - d1$t2[i]
		}else{
			#print(paste(d1$name[i], d2$name[i], sep=" - "))
			lines(c(d1$t1[i],d2$t1[i]),c(d1$t2[i],d2$t2[i]),col=2,lt=1)
			text(d2$t1[i], d2$t2[i], round(s, digits=4), cex=0.75)
			diffs[i,3] <- d2$t1[i] - d1$t1[i]
			diffs[i,4] <- d2$t2[i] - d1$t2[i]
		}
		#if (d1$sim[i] < 0.8) {
		#	text(d1$t1[i], d1$t2[i], paste(round(d1$sim[i], digits=2), d1$name[i], sep=" "), cex=0.75)
		#	if (d1$name[i] != d2$name[i]) {
		#		text(d2$t1[i], d2$t2[i], d2$name[i], cex=0.75)
		#	}
		#}
	}

	
	for (i in 1:4) {
		if (length(na.omit(diffs[,i])) > 0) {
			ds[[i]] <- density(na.omit(diffs[,i]), width=1, adjust=0.01)
		}
	}


	rs <- list()
	rs[[1]] <- range(ds[[1]]$x, ds[[3]]$x, ds[[5]]$x)
	rs[[2]] <- range(ds[[1]]$y, ds[[3]]$y, ds[[5]]$y)
	rs[[3]] <- range(ds[[2]]$x, ds[[4]]$x, ds[[6]]$x)
	rs[[4]] <- range(ds[[2]]$y, ds[[4]]$y, ds[[6]]$y)

	#transparent color
	#rgb(0,100,0,50,maxColorValue=255)
	
	points(d1$t1,d1$t2,col=1,pch=1:total.match%%25, cex=pchcex)
	points(d2$t1,d2$t2,col=tpcols,pch=1:total.match%%25, cex=pchcex)

	points(unmatchedRef$t1, unmatchedRef$t2, col=1, pch=20, cex=pchcex-0.3)
	points(unmatchedTar$t1, unmatchedTar$t2, col=4, pch=20, cex=pchcex-0.3)

	dev.off()
	
	pdf(file=paste("scores_",rid,"-",tid,".pdf",sep=""), width=10, height=10)
	plot(density(ses[,1], adjust=0.1))
	dev.off()
	
	pdf(file=paste("rt1Densities_",rid,"-",tid,".pdf",sep=""), width=20, height=10)
        plot(ds[[5]],xlab="RT1 Diff",ylab="Density", ylim=rs[[2]], xlim=rs[[1]])
	lines(ds[[1]], col=3)
	lines(ds[[3]], col=2)
	dev.off()

	pdf(file=paste("rt2Densities_",rid,"-",tid,".pdf",sep=""), width=20, height=10)
        plot(ds[[6]],xlab="RT2 Diff",ylab="Density", ylim=rs[[4]], xlim=rs[[3]])
	lines(ds[[2]], col=3)
	lines(ds[[4]], col=2)
	dev.off()	


	cat("##################################","\n")
	#all matched ones - true positive matches
	f.p = total.match - t.p
	#all non-matches that should be there
	f.n = total.pos - t.p
	#all true non-matches
	t.n = total.neg - f.p
	t.p.r = t.p/total.pos
	p.p.v = t.p/(t.p+f.p)
	f1 = 2*t.p.r*p.p.v/(t.p.r+p.p.v)
	cat("tp: ",t.p,"\n")
	cat("fp: ",f.p,"\n")
	cat("fn: ",f.n,"\n")
	cat("tn: ",t.n,"\n")
	cat("tpr: ",t.p.r,"\n")
	cat("ppv: ",p.p.v,"\n")
	cat("f1: ",f1,"\n")

	rlt=c(t.p.r,p.p.v,t.p,f.p,f.n,t.n,f1)
	rlt
}
