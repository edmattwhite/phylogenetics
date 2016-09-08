
# mcmc
# markov chain monte carlo
# author: Edward White
# date: 02/08/2016

# Purpose:
# Posterior inference under BDI model of gene cluster evolution. Assume that each cluster has birth rate drawn from a Gamma(theta, eps) distribution, and duplication rate drawn from a Beta(alpha, beta) distribution. Give these parameters uniform priors, and perform posterior inference by MCMC (Multi-try Metropolis algorithm).
# Start MCMC with values drawn from the sample so we have a sensible starting point.
# Under BDI we have cluster frequency following the Polya distribution, with mean
# LOTS OF DATA -> HIGH DIMENSIONS -> NEW MONTE CARLO METHOD -> MULTI PROPOSAL METROPOLIS

#### ------------------------------------------------------------------

# Number of iterations for MCMC
reps = 10000

# Number of clusters to consider
cluster.setup <- 1000

#### Attach dataframe for analysis
#Here the table is kept in the wide format, more logical for the MCMC step

cluster.frequency.wide.unfiltered <- read.table(paste(getwd(),"/comparativeAnalysis_table_clusters_counts.txt",sep=""), header = T)
attach(cluster.frequency.wide.unfiltered)
# Selecting the number of clusters to analyse
cluster.frequency.wide <- cluster.frequency.wide.unfiltered[1:cluster.setup,]
detach(cluster.frequency.wide.unfiltered)
species <- nspecies <- (dim(cluster.frequency.wide)[2]-1)
n <- dim(cluster.frequency.wide)[1]
frequency.matrix <- data.matrix(cluster.frequency.wide[,2:(dim(cluster.frequency.wide)[2]-1)])

#### MCMC ------------------------------------------------------------------
# Initialize vectors of parameters. Setting lambda = 0.5 gives simple expression for mean from data.

theta <- rep(5,reps)
eps <- rep(5,reps)
alpha <- rep(10,reps)
beta <- rep(3,reps)
r.vec <- rep(10,n)
lambda.vec <- rep(0.75,n)
loglike <- rep(-Inf,n)

# Begin MCMC loop

for (rep in 2:reps) {

	# Spit out current iteration to the console periodically, just so we can see it is working.
	if (floor(rep/1000)==(rep/1000)) {
		print(paste("mcmc completed iterations:", rep, "of", reps, sep = " " ))
		flush.console()
	}

# Move forward one step (new values are equal to draws from the previous iteration)

theta[rep] <- theta[rep-1]
eps[rep] <- eps[rep-1]
alpha[rep] <- alpha[rep-1]
beta[rep] <- beta[rep-1]

#### Re-sample r values. This is essentially a Metropolis update step, applied to the entire vector of r values simultaneously to save on computational cost.
# Propose new r values. These must be positive.

prop_r <- abs(rnorm(n,mean=r.vec,sd=1))

# Calculate log-likelihood ratio of proposed values compared with previous values

loglikeRatio <- rowSums(dnbinom(frequency.matrix,size=outer(prop_r/lambda.vec,rep(1,species)),prob=1-outer(lambda.vec,rep(1,species)),log=TRUE))+dgamma(prop_r,shape=theta[rep],rate=eps[rep],log=TRUE)-rowSums(dnbinom(frequency.matrix,size=outer(r.vec/lambda.vec,rep(1,species)),prob=1-outer(lambda.vec,rep(1,species)),log=TRUE))-dgamma(r.vec,shape=theta[rep],rate=eps[rep],log=TRUE)

# Update values if they pass Metropolis test

MetropolisTest <- runif(n)<exp(loglikeRatio)
r.vec[MetropolisTest] <- prop_r[MetropolisTest]

#### Re-sample lambda.vec values
# Propose new lambda.vec values. These must be between 0 and 1.

prop_lambda <- rnorm(n,mean=lambda.vec,sd=0.1)

# Shifting lambda.vec to the correct range

while(length(prop_lambda[prop_lambda<0])>0 | length(prop_lambda[prop_lambda>1])>0) {
	prop_lambda[prop_lambda<0] <- -prop_lambda[prop_lambda<0]
	prop_lambda[prop_lambda>1] <- 2-prop_lambda[prop_lambda>1]
}

# Calculate log-likelihood ratio of proposed values compared with previous values

loglikeRatio <- rowSums(dnbinom(frequency.matrix,size=outer(r.vec/prop_lambda,rep(1,species)),prob=1-outer(prop_lambda,rep(1,species)),log=TRUE))+dbeta(prop_lambda,alpha[rep],beta[rep],log=TRUE)-rowSums(dnbinom(frequency.matrix,size=outer(r.vec/lambda.vec,rep(1,species)),prob=1-outer(lambda.vec,rep(1,species)),log=TRUE))-dbeta(lambda.vec,alpha[rep],beta[rep],log=TRUE)

# Update values if they pass Metropolis test

MetropolisTest <- runif(n)<exp(loglikeRatio)
lambda.vec[MetropolisTest] <- prop_lambda[MetropolisTest]


#### Re-sample theta and eps values
# I have found that the MCMC runs better if we draw new values of the mean and variance of the Gamma distribution, rather than updating the values of theta and eps directly. The proposed mean and variance is a normal draw around the old mean and variance.

prop_gammaMean <- abs(rnorm(1,mean=theta[rep]/eps[rep],sd=0.05))
prop_gammaVar <- abs(rnorm(1,mean=theta[rep]/eps[rep]^2,sd=0.05))

prop_theta <- prop_gammaMean^2/prop_gammaVar
prop_eps <- prop_gammaMean/prop_gammaVar

# Resample theta

# Calculate log-likelihood ratio of proposed value compared with previous value
loglikeRatio <- sum(dgamma(r.vec,shape=prop_theta,rate=eps[rep],log=TRUE)-dgamma(r.vec,shape=theta[rep],rate=eps[rep],log=TRUE))

# Update value if it passes Metropolis test
if (runif(1)<exp(loglikeRatio)){
	theta[rep] <- prop_theta
	}

# Re-sample eps

# Calculate log-likelihood ratio of proposed value compared with previous value
loglikeRatio <- sum(dgamma(r.vec,shape=theta[rep],rate=prop_eps,log=TRUE)-dgamma(r.vec,shape=theta[rep],rate=eps[rep],log=TRUE))

# Update value if it passes Metropolis test
if (runif(1)<exp(loglikeRatio)){
	eps[rep] <- prop_eps
	}


#### Re-sample alpha and beta values
# As above, I have found that the MCMC runs better if we draw new values of the mean and variance of the Beta distribution, rather than updating the values of alpha and beta directly. The proposed mean must be between 0 and 1, and the proposed variance must be between 0 and mean*(1-mean).

prop_betaMean <- -1
prop_betaVar <- -1
while(prop_betaMean<0 | prop_betaMean>1 | prop_betaVar<0 | prop_betaVar>(prop_betaMean*(1-prop_betaMean))) {
	prop_betaMean <- rnorm(1,mean=alpha[rep]/(alpha[rep]+beta[rep]),sd=0.005)
	prop_betaVar <- rnorm(1,mean=alpha[rep]*beta[rep]/(alpha[rep]+beta[rep])^2/(alpha[rep]+beta[rep]+1),sd=0.005)
}
prop_beta <- prop_betaMean*(1-prop_betaMean)^2/prop_betaVar - 1 + prop_betaMean
prop_alpha <- prop_betaMean*prop_beta/(1-prop_betaMean)

#resample alpha

# Calculate log-likelihood ratio of proposed value compared with previous value
loglikeRatio <- sum(dbeta(lambda.vec,prop_alpha,beta[rep],log=TRUE)-dbeta(lambda.vec,alpha[rep],beta[rep],log=TRUE))

# Update value if it passes Metropolis test
if (runif(1)<exp(loglikeRatio)) {
	alpha[rep] <- prop_alpha
	}


# re-sample beta

# Calculate log-likelihood ratio of proposed value compared with previous value
loglikeRatio <- sum(dbeta(lambda.vec,alpha[rep],prop_beta,log=TRUE)-dbeta(lambda.vec,alpha[rep],beta[rep],log=TRUE))

# Update value if it passes Metropolis test
if (runif(1)<exp(loglikeRatio)) {
	beta[rep] <- prop_beta
	}


# Calculate the log-probability of all free parameters and the data. This may be useful later on (e.g. for measuring convergence).
loglike[rep] <- sum(dnbinom(frequency.matrix,size=outer(r.vec/lambda.vec,rep(1,species)),prob=1-outer(lambda.vec,rep(1,species)),log=TRUE))+sum(dgamma(r.vec,shape=theta[rep],rate=eps[rep],log=TRUE))+sum(dbeta(lambda.vec,alpha[rep],beta[rep],log=TRUE))

} # close MCMC loop

#### Print out what proportion of proposed moves are accepted. Should be aiming for about 25% of moves accepted.
theta_moves = round(1-mean(theta[-reps]==theta[-1]),digits=3)
eps_moves = round(1-mean(eps[-reps]==eps[-1]),digits=3)
alpha_moves = round(1-mean(alpha[-reps]==alpha[-1]),digits=3)
beta_moves = round(1-mean(beta[-reps]==beta[-1]),digits=3)

cat(paste("theta moves = ", theta_moves, "\neps moves = ", eps_moves, "\nalpha moves = ", alpha_moves, "\nbeta moves = ", beta_moves,"\n", sep=""))

#### Attach new estimated lambda.vec and r.vec rates into the data frame
cluster.frequency.wide$lambda <- lambda.vec
cluster.frequency.wide$r <- r.vec
rm(lambda.vec,r.vec,species)

#### Integrate over entire posterior distribution of unknown parameters. This approach seems to give more accurate results than the simple maximum-likelihood approach above.

# Find integrated cumulative posterior r distribution;
xVec <- seq(0,max(cluster.frequency.wide$r),l=n)
xMat <- outer(rep(1,reps),xVec)
thetaMat <- outer(theta,rep(1,n))
epsMat <- outer(eps,rep(1,n))
sampleVec <- cluster.frequency.wide$r
sampleMat <- outer( rep(1,reps), sampleVec )


post.r.dist <- colMeans(dgamma(xMat,shape=thetaMat,rate=epsMat))
cluster.frequency.wide$F.r <- colMeans(pgamma(sampleMat,shape=thetaMat,rate=epsMat))

pdf(paste(getwd(),"/plots/creation_rate_dist.pdf",sep=""))
r.plot <- hist( cluster.frequency.wide$r , xlab="r_i", ylab="probability", plot = TRUE, freq = FALSE)
		   lines( xVec, post.r.dist, type = "l" )
rm(r.plot)
bequiet <- dev.off()

# Find integrated posterior lambda distribution;
xVec <- seq(0,max(cluster.frequency.wide$lambda),l=n)
xMat <- outer(rep(1,reps),xVec)
alphaMat <- outer(alpha,rep(1,n))
betaMat <- outer(beta,rep(1,n))
sampleVec <- cluster.frequency.wide$lambda
sampleMat <- outer( rep(1,reps), sampleVec )

post.lambda.dist <- colMeans(dbeta(xMat,alphaMat,betaMat))
cluster.frequency.wide$F.lambda <- colMeans(pbeta(sampleMat,alphaMat,betaMat))

pdf(paste(getwd(),"/plots/transfer_rate_dist.pdf",sep=""))
lambda.plot <- hist(  cluster.frequency.wide$lambda , xlab="lambda_i", ylab="probability" , plot = TRUE, freq = FALSE)
	    lines(xVec, post.lambda.dist , type = "l")
rm(lambda.plot)
bequiet <- dev.off()

# Exporting a data frame with the distributions for r and lambda
parameter.data <- data.frame(alpha, beta, theta, eps)
write.table(parameter.data, file = paste(getwd(),"\\tables\\parameter_data.txt",sep=""))
rm(parameter.data)

#### Format the data frame so that we can start evaluating likelihoods of gene lengths.

cluster.frequency <- reshape(cluster.frequency.wide,
	varying = names(cluster.frequency.wide)[2:(dim(cluster.frequency.wide)[2]-4)],
	v.names = "frequency",
	timevar = "species",
	times = names(cluster.frequency.wide)[2:(dim(cluster.frequency.wide-2)[2]-4)],
	direction = "long")
rownames(cluster.frequency) <- NULL

# Removing the unnecessary 'id' column
cluster.frequency$id <- NULL

#### Associating likelihood with the appearance of a cluster, P(x|lambda,r)
# Cluster likelihood function to see likelihood of individual clusters.

clus_lik <- function(data.frame){
	funclambda <- data.frame[,"lambda"]
	funcr <- data.frame[,"r"]
	ratio <- funcr/funclambda
	obs <- data.frame[,"frequency"]
  	lik <- dnbinom(obs,size=ratio,prob=1-funclambda,log=FALSE)
	return(lik)
	}

# Adding on the likelihood column to the data frame
cluster.frequency$likelihood <- clus_lik(cluster.frequency)
likelihood.table <- matrix(cluster.frequency$likelihood, ncol = nspecies, nrow = dim(cluster.frequency)[1], byrow = F)
names(likelihood.table) <- levels(factor(cluster.frequency$species))
write.table(likelihood.table, file  = paste(getwd(),"/tables/likelihood_table.txt", sep = ""))

#### Forming ranking column for tbbm -----------------------------------------------------------------------------
attach(cluster.frequency)
ranking <- cluster.frequency.wide$r + 1/(cluster.frequency.wide$lambda)
ranking.column <- data.frame(ranking)
write.table(ranking.column, file = paste(getwd(),"\\tables\\ranking_column.txt", sep = ""))
detach(cluster.frequency)

#-----------------------------------------------------------------------------------------------------------------

# Exporting the cluster.frequency data frame for use by other scripts

write.table(cluster.frequency, file = paste(getwd(),"\\tables\\cluster_frequency.txt",sep=""))

# clearing loaded objects
rm(list = ls())

