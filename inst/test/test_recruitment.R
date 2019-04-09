#========================================
# Hierarchical Poisson Linear Regression
#========================================

library(twoe)

#== Generating data

# Set seed
seed <- 1234
set.seed(seed)
# Ecological process (suitability)
# Constants 
nobs <- 1000 
ngroup <- 20 
np <- 3
nq <- 3
# Groups (at least one obervation per group)
IdentGroup <- c(1:ngroup,sample(c(1:ngroup),(nobs-ngroup),replace=TRUE))
# Time intervals (in years)
Int_vect <- rep(1, nobs) # runif(n=nobs,min=0.5,max=5)  
# Cell surface (in m2)
SCell_vect <- rep(1, nobs) # runif(n=nobs,min=0.5,max=10)
# Covariates 
X1 <- runif(n=nobs, min=-1, max=1) 
X2 <- runif(n=nobs, min=-1, max=1) 
X <- cbind(rep(1, nobs), X1, X2) 
W <- X
# Target parameters 
# beta 
beta.target <- matrix(c(0.3, 0.2, 01), ncol=1) 
# Vb 
Vb.target <- c(0.5, 0.1, 0.1) 
# b 
b.target <- cbind(rnorm(ngroup,mean=0,sd=sqrt(Vb.target[1])), 
				  rnorm(ngroup,mean=0,sd=sqrt(Vb.target[2])), 
				  rnorm(ngroup,mean=0,sd=sqrt(Vb.target[3]))) 
# V 
V.target <- 0.5

# Response variable
e <- rnorm(nobs, 0, sd=sqrt(V.target))
log_lambda = X %*% beta.target + rowSums(W*b.target[IdentGroup,]) + e
lambda_prim <- exp(log_lambda)*Int_vect*SCell_vect
Y <- rpois(nobs,lambda_prim) 

# Data-set 
Data <- as.data.frame(cbind(Y,log_lambda,X1,X2,IdentGroup)) 
# Plot data
library(ggplot2)
ggplot(Data, aes(X1,log_lambda,color=as.factor(IdentGroup))) +
	geom_point()

#== Call to entry_recruitment_gibbs
model <- entry_recruitment_gibbs(
	fixed=Y~X1+X2, random=~X1+X2, group="IdentGroup",
	interval=1, area=1, data=Data, burnin=5000, mcmc=1000, thin=1,verbose=1,
	seed=NA, beta.start=0, sigma2.start=0.5,
	Vb.start=1, mubeta=0, Vbeta=1.0E6,
	r=3, R=diag(c(0.1,0.1,0.1)), nu=0.001, delta=0.001, FixOD=1)

#== MCMC analysis

# Graphics
pdf("Posteriors-entry_recruitment_gibbs.pdf")
plot(model$mcmc)
dev.off()

# Summary
summary(model$mcmc)

# Predicted-Observed
plot(log_lambda, model$log_lambda_pred)
abline(a=0, b=1, col="red")
plot(exp(log_lambda), exp(model$log_lambda_pred))
abline(a=0, b=1, col="red")

# End