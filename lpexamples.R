###################################################
## Author:	Jeffrey Blume
## Date: 	10/14/2017
## File: 	lpexamples.R
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
## Compute p-values
###################################################
#### X ~ N(mu, sd.both)

#### RR of sample mean > 0.5
alpha=1-pnorm(rr,mean=mu.0,sd=sd.both/sqrt(n))
power=1-pnorm(rr,mean=mu.1,sd=sd.both/sqrt(n))

#### observed sample mean = 0.6
pv.0=1-pnorm(0.6,mean=mu.0,sd=sd.both/sqrt(n))
pv.1=1-pnorm(0.6,mean=mu.1,sd=sd.both/sqrt(n))

###################################################
## Plot of likelihoods
###################################################
theta=seq(-0.25,1.5,0.01)

lf.pow=(1-pnorm(rr,mean=theta,sd=sd.both/sqrt(n)))/(1-pnorm(rr,mean=mu.0,sd=sd.both/sqrt(n)))
lf.pv=(1-pnorm(0.6,mean=theta,sd=sd.both/sqrt(n)))/(1-pnorm(0.6,mean=mu.0,sd=sd.both/sqrt(n)))
lf.nor=exp(n*(theta-mu.0)/(sd.both^2)*(0.6-(theta+mu.0)/2))/exp(n*(mu.0-mu.0)/(sd.both^2)*(0.6-(mu.0+mu.0)/2))
## check: lr.nor=dnorm(0.6,mean=theta,sd=1/sqrt(14))/dnorm(0.6,mean=0,sd=1/sqrt(14))

pv.bon=min(3*pv.0,1)
rr.bon=qnorm(1-pv.bon,mean=mu.0,sd=sd.both/sqrt(n))
lf.bon=(1-pnorm(rr.bon,mean=theta,sd=sd.both/sqrt(n)))/(1-pnorm(rr.bon,mean=mu.0,sd=sd.both/sqrt(n)))

plot(theta,lf.pv,type="n",ylab="Likelihood ratio L(H1)/L(H0)",xlab="Alternative Hypothesis (H1)")
lines(theta,lf.pow,col="forestgreen",lwd=2)
lines(theta,lf.pv,col="firebrick",lwd=2)
lines(theta,lf.nor,col="dodgerblue",lwd=2)
lines(theta,lf.bon,col="purple",lwd=2)
axis(side=1,at=0.6,label=0.6)

segments(0.5,-2.5,0.5,40,lty=2,lwd=1,col="black")

legend("topleft",c("LR | H0 rejected","LR | p = 0.0124","LR | xbar = 0.6","LR | p.adj = 0.037"),bty="n",
	col=c("forestgreen","firebrick","dodgerblue","purple"),lwd=2,lty=1)

###################################################
## Plot of FDRs
###################################################
plot(theta[lf.pv>=1],1/(1+lf.pv)[lf.pv>=1],type="n",ylab="False Discovery Rate",xlab="Alternative Hypothesis (H1)",ylim=c(0,0.5))
lines(theta[lf.pow>=1],1/(1+lf.pow)[lf.pow>=1],col="forestgreen",lwd=2)
lines(theta[lf.pv>=1],1/(1+lf.pv)[lf.pv>=1],col="firebrick",lwd=2)
lines(theta[lf.nor>=1],1/(1+lf.nor)[lf.nor>=1],col="dodgerblue",lwd=2)
lines(theta[lf.bon>=1],1/(1+lf.bon)[lf.bon>=1],col="purple",lwd=2)

axis(side=1,at=0.6,label=0.6)

segments(0.6,-0.02,0.6,0.2,lty=2,lwd=1,col="black")
text(0.6,0.21,"sample mean")

legend("top",c("FDR | H0 rejected","FDR | p = 0.0124","FDR | xbar = 0.6","FDR | p.adj = 0.037"),bty="n",
	col=c("forestgreen","firebrick","dodgerblue","purple"),lwd=2,lty=1)

###
##
#
