###################################################
## Author:	Jeffrey Blume
## Date: 	10/14/2017
## File: 	error-fdr.R
## Purpose:	Simulations of likelihood principle
###################################################


###################################################
## Parameter specifications
###################################################
n=14		## sample size
sd.both=1	## SD

mu.0=0		## null mean
mu.1=0.5	## alt mean, assumes mu.1 > mu.0

sims.0=1000000	## Null sims  ;## 200000
sims.1=1000000	## Alt Sims

p.0=sims.0/(sims.0+sims.1)	## P(H0)

rr=0.5		## Rejection region xbar>rr

###################################################
## Generate data
###################################################
#### X ~ N(mu, sd.both)

set.seed(464646)
x.0=matrix(rnorm(n*sims.0,mu.0,sd=sd.both),ncol=n,nrow=sims.0)
x.1=matrix(rnorm(n*sims.1,mu.1,sd=sd.both),ncol=n,nrow=sims.1)

#### Sample means
xbar.0=rowMeans(x.0)
xbar.1=rowMeans(x.1)

###################################################
## table sample means ; compute error rates & FDRs
###################################################
cutoff=cbind(sm=c(xbar.0,xbar.1),hyp=c(rep(0,sims.0),rep(1,sims.1)),reject=c(1*(xbar.0>rr),1*(xbar.1>rr)))

addmargins(table(hyp=cutoff[,"hyp"],reject=cutoff[,"reject"]))

mis.0=1-pnorm(rr,mean=mu.0,sd=sd.both/sqrt(n))
pow.1=1-pnorm(rr,mean=mu.1,sd=sd.both/sqrt(n))

fdr=1/(1+pow.1/mis.0*(p.0/(1-p.0)))
fcr=1/(1+(1-mis.0)/(1-pow.1)*((1-p.0)/p.0))

prop.table(table(hyp=cutoff[,"hyp"],reject=cutoff[,"reject"]),margin=1) ## Error rates
c(1-pow.1,mis.0)

prop.table(table(hyp=cutoff[,"hyp"],reject=cutoff[,"reject"]),margin=2) ## FDR rates
c(fcr,fdr)

###################################################
## Plots
###################################################

#### Plot two sampling distributions
hist(xbar.0,breaks=80,col=rgb(0,0,1,0.5),xlim=c(min(xbar.0,xbar.1),max(xbar.0,xbar.1)),
	xlab="sample means",main="Distribution of Sample Means",freq=FALSE)
hist(xbar.1,breaks=80,add=T,col=rgb(1,0,0,0.5),freq=FALSE)
	axis(side=1, at=0, labels=FALSE)
abline(v=rr,lty=1,lwd=3,col="black")

legend("topright",c("Null","Alternative","Rejection Line"),
		pch=c(15,15,NA),pt.cex=c(1.5,1.5,1),lwd=c(NA,NA,3),lty=c(NA,NA,1),
		col=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5),"black"),bty="n")

#### Plot distribuiton of sample means (mixture)
my.hist=hist(c(xbar.0,xbar.1), breaks=80, prob=TRUE,plot=F)
 
# Color vector:  'seagreen' rgb(0.2,0.8,0.5,0.5) & 'purple' rgb(0.2,0.2,0.2,0.2)
my.color= ifelse(my.hist$breaks < rr, 0 , 
			ifelse (my.hist$breaks >= rr, "purple", rgb(0.2,0.2,0.2,0.2) ))
 
# Final plot
plot(my.hist, col=my.color,freq=FALSE,xlab="sample means",main="Mixture Distribution of Sample Means")

den.mx=density(c(xbar.0,xbar.1),adjust=2)
den.0=density(c(xbar.0),adjust=2)
den.1=density(c(xbar.1),adjust=2)

lines(den.mx,lwd=2,col="black") 
lines(den.0$x,p.0*den.0$y,lwd=2,col="blue")  
lines(den.1$x,(1-p.0)*den.1$y,lwd=2,col="red")  

legend("topright",c("Mixture","Null","Alternative","Rejected"),
		lty=c(1,1,1,NA),lwd=c(2,2,2,NA),pch=c(NA,NA,NA,15),
		pt.cex=c(1,1,1,1.5),
		col=c("black","blue","red","purple"),bty="n")

mtext(paste("P(H0)= ",round(p.0,3),sep=""),1,line=4)

###
##
#
