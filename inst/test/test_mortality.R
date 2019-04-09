#========================================
# Hierarchical Binomial Linear Regression
#========================================

library(twoe)

#== Generating data

# Constants
nobs <- 1000
nspecies <- 20
species <- c(1:nspecies,sample(c(1:nspecies),(nobs-nspecies),replace=TRUE))

# Covariates
X1 <- runif(n=nobs,min=-10,max=10)
X2 <- runif(n=nobs,min=-10,max=10)
X <- cbind(rep(1,nobs),X1,X2)
W <- X

# Target parameters
# beta
beta.target <- matrix(c(0.3,0.2,0.1),ncol=1)
# Vb
Vb.target <- c(0.5,0.05,0.05)
# b
b.target <- cbind(rnorm(nspecies,mean=0,sd=sqrt(Vb.target[1])),
				  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[2])),
				  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[3])))

# V 
V.target <- 0.5

# Response variable
e <- rnorm(nobs, 0, sd=sqrt(V.target))
logit_theta = X %*% beta.target + rowSums(W*b.target[IdentGroup,]) + e
theta_prim <- twoe::inv_logit(logit_theta)*Int_vect*SCell_vect
Y <- rbinom(nobs,size=1,theta_prim) 

# Data-set 
Data <- as.data.frame(cbind(Y,logit_theta,X1,X2,IdentGroup)) 
# Plot data
library(ggplot2)
ggplot(Data, aes(X1,logit_theta,color=as.factor(IdentGroup))) +
	geom_point()

#== Call to entry_mortality_gibbs
model <- entry_mortality_gibbs(
	fixed=Y~X1+X2, random=~X1+X2,
	group="IdentGroup",
	interval=1, data=Data, burnin=5000, mcmc=1000, thin=1,verbose=1,
	seed=NA, beta.start=0, sigma2.start=0.5,
	Vb.start=1, mubeta=0, Vbeta=1.0E6,
	r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001, FixOD=1)

#== MCMC analysis

# Graphics
pdf("Posteriors-entry_mortality_gibbs.pdf")
plot(model$mcmc)
dev.off()

# Summary
summary(model$mcmc)

# Predicted-Observed
plot(logit_theta, model$logit_theta_pred)
abline(a=0, b=1, col="red")
plot(twoe::inv_logit(logit_theta),
	 twoe::inv_logit(model$logit_theta_pred))
abline(a=0, b=1, col="red")

# End

