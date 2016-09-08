# bamc
# bayesian analysis metric calculator
# author: Edward White
# date: 12/08/2016

# purpose:
# takes the distributions for lambda, r, and cluster length and uses them to calculate f(mu|x), which can then be used to assign # probability values to any given tree topology (with some modification).
# x = cluster length, mu = (lambda,r), i.e. our parameters taken from mcmc analysis


#### ----------------------------------------------------------------------------------------

# loading in dataframes needed for analysis
parameter.data <- read.table(paste(getwd(),"/tables/parameter_data.txt",sep=""), header = TRUE)
cluster.frequency.all <- read.table(paste(getwd(),"/tables/cluster_frequency.txt",sep=""), header = TRUE)
cluster.frequency.matrix <- matrix(cluster.frequency.all$frequency, nrow = length(levels(factor(cluster.frequency.all$cluster))), ncol = length(levels(factor(cluster.frequency.all$species))))
clusters <- length(levels(factor(cluster.frequency.all$cluster)))
cluster.frequency <- as.data.frame(cbind(cluster.frequency.matrix, cluster.frequency.all$lambda[1:clusters], cluster.frequency.all$r[1:clusters]))
names(cluster.frequency) <- c(levels(factor(cluster.frequency.all$species)), "lambda", "r")
species <- length(levels(factor(cluster.frequency.all$species)))
rm(cluster.frequency.all)

# Loading up the correct number of clusters
i <- read.table(paste(getwd(),"/tables/cluster_number.txt", sep = ""))[1,1]
cluster.frequency <- cluster.frequency[1:i,]

# extracting data into form for calculations
attach(parameter.data)
n <- dim(cluster.frequency)[1]
m <- dim(parameter.data)[1]
alphaMat <- outer(alpha,rep(1,n))
betaMat <- outer(beta,rep(1,n))
sampleVec.lambda <- cluster.frequency$lambda
sampleMat.lambda <- outer( rep(1,m), sampleVec.lambda )
thetaMat <- outer(theta,rep(1,n))
epsMat <- outer(eps,rep(1,n))
sampleVec.r <- cluster.frequency$r
sampleMat.r <- outer( rep(1,m), sampleVec.r )
parameter.vector <- colMeans(parameter.data)
detach(parameter.data)

# bayes metric calculator

# calculating P(x|mu)

num <- dnbinom( as.matrix(cluster.frequency[,1:species]),size = as.vector(as.matrix((cluster.frequency[,species+2]/cluster.frequency[,species+1]))), prob = as.vector(as.matrix(1 - cluster.frequency[,species+1])), log = FALSE)

# calculating P(mu)

param.prob <- (colMeans(pbeta(apply(sampleMat.lambda  + 0.1, c(1,2), function(x) min(1,x)),alphaMat,betaMat) - pbeta(apply(sampleMat.lambda  - 0.1, c(1,2), function(x) max(0,x)),alphaMat,betaMat)))*(colMeans(pgamma(sampleMat.r + 0.5,shape=thetaMat,rate=epsMat) - pgamma(apply(sampleMat.r - 0.5, c(1,2),function(x) min(0,x)),shape=thetaMat,rate=epsMat)))
	
# calculating P(x) using ecdf
data.vec <- as.numeric(data.matrix(cluster.frequency[,1:species]))
f <- ecdf(data.vec)
dnom <- apply(cluster.frequency[,1:species] + 2, c(1,2), f) - apply(cluster.frequency[,1:species] - 2, c(1,2), f)

# putting results together

results.matrix <- matrix(mapply("/", num, dnom), nrow = n, ncol = species)
results.matrix <- cbind(results.matrix, param.prob)

# exporting the data
	
bayes.data <- data.frame(results.matrix)
names(bayes.data) <- c(names(cluster.frequency)[1:species], "parameter_prob")
write.table( bayes.data, file = paste(getwd(),"/tables/bayes_data.txt",sep="") )

print("bamc complete")

#clearing loaded objects
rm(list = ls())