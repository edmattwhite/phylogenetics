 
# tbbm
# tree builder by bayesian metric
# author: Edward White
# date: 08/08/2016

# Purpose:
# Take as input the posterior distributions from MCBA and then use these to perform Bayesian analysis on phyogenetic trees.
# FIND A WAY OF OVERLAYING THE CLUSTERS TO WORK IN HYBRIDISATION

#### ------------------------------------------------------------------

library(abind)
library(ape)

# Loading up dataframes for the posterior pdf.
sieve.data <- read.table(paste(getwd(),"/comparativeAnalysis_table_clusters_counts.txt",sep=""), header = T)
cluster.col <- sieve.data[1,]

ranking.column <- read.table(file = paste(getwd(), "/tables/ranking_column.txt", sep =""), header = TRUE) 
names(ranking.column) <- "ranking"
prob.data <- read.table(paste(getwd(),"/tables/bayes_data.txt", sep = ""))

# Loading up the correct number of clusters
i <- read.table(paste(getwd(),"/tables/cluster_number.txt", sep = ""))[1,1]
sieve.data <- sieve.data[1:i,]
prob.data <- prob.data[1:i,]
ranking.column <- ranking.column[1:dim(sieve.data)[1],] 
print(paste("Clusters used:", dim(sieve.data)[1], sep = " "))

#### Setup for tree building

species <- dim(sieve.data)[2] - 1
# number of clusters to be used for the analysis
clusters <- dim(sieve.data)[1]
# Defining numper of time epochs 
gens <- species
# Setting up column to find f interval of ALL species by row
sieve.data$delta.f <- apply(sieve.data[,2:(species + 1)], 1 , function(x) (max(x)-min(x)) )
# Array to store blueprints of agreeing clusters storing epoch of time at which a cluster is sieved out.
tree.blueprint.array <- array( data = 0 , dim = c(species,species,gens,clusters) ) 
top.delta.table <- sort.list(sieve.data$delta.f[1:50], decreasing=TRUE)
write.table(top.delta.table, file = paste(getwd(),"/tables/identifying_clusters.txt",sep="")) 
sieve.data$ranking <- ranking.column 
attach(sieve.data)
sieve.data <- sieve.data[order(ranking),]
sieve.data <- sieve.data[1:clusters,]
cluster.col <- sieve.data[,1]
sieve.data <- sieve.data[,-1]
detach(sieve.data)
param.prob <- prob.data[,dim(prob.data)[2]]

# Starting resolution of filter
f.start <- 0.5*apply(sieve.data[,1:species], 1, function(x) max(apply( embed(x,2), 1, function(y) max(y) - min(y))))
# Ending resolution of filter 
f.end <- 0.5*apply(sieve.data[,1:species], 1, function(x) min(apply( embed(x,2), 1, function(y) max(y) - min(y))))

#### Function for grouping clusters according to delta.f
# Sorts data into sets with radius of each set \leq size.
# Will attempt to sort data into groups with smallest diameter and largest number of entries.
# Inputs: data = vector of data to sort, size = radius of largest sets allowed. 
# Outputs: sort.matrix

data.sieve <- function( data, size ){

	# Need to order first 
	old.ordering <- names(data)
	correct.order.positions <- vector(mode = "numeric", length = length(data))
	data <- sort(data, decreasing = FALSE)	
	new.ordering <- names(data)	
	
	# Initialising a matrix to specify groupings.
	sort.matrix <- matrix( 0 , ncol = length(data), nrow = length(data) )

	# If a point is a singleton (completely isolated by the sieve) we can discard it immediately.
	locality.count <- lapply(data, function(x) length(data[x-size<data & data<x+size]))
	singleton.positions <- lapply(locality.count, function(x) x==1)
	singleton.positions <- singleton.positions[old.ordering]
	
	# Now looking at only those elements with neighbours:
	grouped.positions <- mapply("-", locality.count, singleton.positions)
	
	# Finding peaks of the vector - corresponding to those points which have the largest number of sufficiently close neighbours.
	comparison.table <- embed(as.vector(grouped.positions), 3)  
	index.first.max <- max.col(comparison.table, ties.method = "first" ) == 2 
	index.second.max <- max.col(comparison.table, ties.method = "last" ) == 2
	index.max <- index.first.max + index.second.max
	index.max <- replace( index.max, index.max == 2 , TRUE) 
	if(grouped.positions[[1]]!=0){
		index.max <- c(TRUE, index.max)
	}else{ 
		index.max <- c(FALSE, index.max)
	}
	if(grouped.positions[[length(grouped.positions)]]!=0){
		index.max <- c(index.max, TRUE)
	}else{ 
		index.max <- c(index.max, FALSE)
	}
	index.max <- as.logical(index.max)

	# Using the peaks located above to sort the elements into their groups
	neighbour.data <- lapply( data, function(x) data[x-size<data & data<x+size] )
	group.data <- neighbour.data[index.max]
	position.vector <- vector(mode = "numeric", length = length(data))
	tics <- length(group.data)
	if( tics != 0){
		for( tic in 1:tics ){
			position.vector <- match(data, group.data[[tic]], nomatch = 0)
			position.vector <- replace(position.vector, position.vector>1, 1)
			names(position.vector) <- new.ordering
			position.vector <- position.vector[old.ordering]
			sort.matrix <- sort.matrix + outer(position.vector,position.vector)
		} 
	}
	sort.matrix <- sort.matrix + diag(species)
	sort.matrix <- matrix(1, ncol = species, nrow = species)*(sort.matrix>=1)
	return(sort.matrix)
}

#### TREE BUILDER ------------------------------------------------------
# Algorithm designed to build trees forwards in time based on the posterior distributions obtained using bamc. We start the tree with a single common ancestor.
# Then sift through every cluster type, trying to obtain an interval containing as many species as possible whilst having p of a measurement occuring in that interval. 
# We "peel off" each cluster as it diverges from the other species, recording when this peeling occurs. 
f <- f.start

for(gen in 1:gens){

# Spit out current iteration to the console periodically, just so we can see it is working.
  print(paste("tbbm completed iterations:",gen, "of", gens, sep = " "))
  flush.console()


# Getting the sorting matrix for this f size, whilst checking the sieve hasn't already separated out all these points - if it has then we can just paste over the result from the previous generation
for( i in 1: clusters){
	if( gen == 1 ){
	tree.blueprint.array[,,gen,i] <- i*data.sieve(data = sieve.data[i,1:species], size = f[i])
	}else{ 
		if(sum(tree.blueprint.array[,,gen-1,i] - diag(species))==0){
		tree.blueprint.array[,,gen,i] <- i*diag(species)
		}else{ 
		tree.blueprint.array[,,gen,i] <- i*data.sieve(data = sieve.data[i,1:species], size = f)
	}	}
}

# Decreasing filter size of sieve
f <- f.end + (f.start - f.end)*(gens - 1 - gen)/(gens-1)

} # Close tree building loop

# Assembling tree by superimposing tree for each cluster to give edge weights (collapse tree.blueprint.array to 2 dimensions), then 
tree.superimposed.array <- apply(tree.blueprint.array, c(1,2), sum)
weighted.adjacency <- matrix(1, ncol = species, nrow = species)*max(tree.superimposed.array) - tree.superimposed.array
d <- weighted.adjacency/(0.5*clusters*(clusters+1)*gens)
rownames(d) <- names(sieve.data)[1:species]
d <- as.dist(d)
names.list <- names(sieve.data)[1:species]

second.sort.hc <- hclust(d, method = "average", members = NULL)
second.sort.hcd <- as.dendrogram(second.sort.hc)

# trying to spot hybridisation in the network
sum.array <- matrix(data = 0, ncol = species, nrow = species)
hybrid.list <- list()

for(i in 1:species){
	hybrid.indicator.slice <- apply(tree.blueprint.array[,,2,], 2, function(x) (x*tree.blueprint.array[,i,2,]) > 0)
	hybrid.indicator.slice[,i] <- 0 
	for(j in 1:species){
		for(k in 1:j){
			sum.array[j,k] <- sum((hybrid.indicator.slice[,j] + hybrid.indicator.slice[,k])>0)
			sum.array[j,j] <- 0 
		}
	}
	index.position <- apply(sum.array, c(1,2), function(y) y==max(sum.array))
	index.1 <- sum(apply(index.position, 2, function(y) y*c(1:species)))
	index.2 <- sum(apply(index.position, 1, function(y) y*c(1:species)))
	hybrid.prop <- max(sum.array)/(clusters*species)
	hybrid.list[[i]] <- paste("Proposed hybrid parents for", names.list[i], "are", names.list[index.1], "and", names.list[index.2], "proportion of clusters accounted for", hybrid.prop, sep = " ")
}
l <- length(hybrid.list)
text <- ""
hybrid_file <- file(paste(getwd(), "/tables/hybrid_list.txt", sep = ""))
for(j in 1:l){text <- paste(text, toString(hybrid.list[[j]]), sep = "\r\n")}
writeLines(text, hybrid_file)
close(hybrid_file)
	
# "partnering" the weighted.adjacency matrix to produce a phylogenetic tree in string format

tree.list <- list()

min.vec <- 0
weighted.adjacency <- tree.superimposed.array*lower.tri(tree.superimposed.array, diag = F)
i <- 0
restore <- tree.blueprint.array
while(length(names.list) > 1){

	# This bit extracts the pairs

	i <- i+1
	tic <- i
	logical.position <- apply(weighted.adjacency, c(1,2), function(y) y==max(weighted.adjacency[weighted.adjacency!=0]))
	index.1 <- sum(apply(logical.position, 2, function(y) y*c(1:length(names.list))))
	index.2 <- sum(apply(logical.position, 1, function(y) y*c(1:length(names.list))))
	
	# This bit assigns a probability to the pairing - approximate probability of having matching parameters for all clusters and prints which clusters disagree with the pairing 
	indicator.slice <- tree.blueprint.array[index.1,index.2,2,] > 0
	difference.slice <- ((rep(1, clusters) - tree.blueprint.array[index.1,index.2,floor(gens/2),])*cluster.col)[(rep(1, clusters) - tree.blueprint.array[index.1,index.2,floor(gens/2),]) > 0]
	f.identical.species.slice <- (prob.data[,index.1]*param.prob)[indicator.slice]
	f.identical.species <- mean(f.identical.species.slice)
	tree.list[[i]] <- c(names.list[c(index.1,index.2)], f.identical.species, difference.slice)

	# This bit forms a new 'placeholder' species
	names.list <- append(names.list, paste("(",names.list[index.1],"+",names.list[index.2],")",sep=""))
	names.list <- names.list[-c(index.1,index.2)]
	if(length(names.list)==2) break
	new.array.slice <- tree.blueprint.array[index.1,,,]*tree.blueprint.array[index.2,,,]
	new.array.slice <- new.array.slice[-c(index.1,index.2),,]
	new.array.slice <- abind(new.array.slice, array(data = 1, dim = dim(new.array.slice)[-1]), along = 1)
	tree.blueprint.array <- abind(tree.blueprint.array, array(data = 0, dim = dim(tree.blueprint.array)[-1]), along = 1)
	tree.blueprint.array <- abind(tree.blueprint.array, array(data = 0, dim = dim(tree.blueprint.array)[-2]), along = 2)
	tree.blueprint.array <- tree.blueprint.array[-c(index.1,index.2),,,]
	tree.blueprint.array <- tree.blueprint.array[,-c(index.1,index.2),,]
	tree.blueprint.array[length(names.list),,,] <- new.array.slice
	tree.blueprint.array[,length(names.list),,] <- new.array.slice
	matrix <- apply(tree.blueprint.array, c(1,2), sum)
	weighted.adjacency <- matrix*lower.tri(matrix) 
	weighted.adjacency <- weighted.adjacency*lower.tri(weighted.adjacency, diag = F)
	prob.data <- cbind(prob.data[,-c(index.1,index.2,dim(prob.data)[2])],(prob.data[,index.1]+prob.data[,index.2])/2,param.prob)
	names(prob.data) <- c(names.list, "param_prob")
}	

indicator.slice <- tree.blueprint.array[1,2,2,] > 0
f.identical.species <- mean((prob.data[,names.list[1]]*prob.data[,tree.list[[i]][2]]*param.prob)[indicator.slice])
tree.list[[tic+1]] <- c(names.list[c(1,2)],f.identical.species)
text <- toString(tree.list[[1]])
tree_file <- file(paste(getwd(),"/tables/tree_list.txt", sep = ""))
for(j in 1:tic+1){text <- paste(text, toString(tree.list[[j]]), sep = "\n")}
writeLines(text, tree_file)
close(tree_file)

# plotting
pdf(paste(getwd(),"/plots/bayesian_tree_rooted.pdf",sep=""))
par(mar = c(5,4,4,2) + 5)
plot(second.sort.hcd, main = "Phylogenetic Tree by Bayesian Analysis", xlab = "Relative distance from perfect match", ylab = "Species", type = "rectangle", horiz = T)
bequiet <- dev.off()
pdf(paste(getwd(),"/plots/bayesian_tree_unrooted.pdf",sep=""))
plot(as.phylo(second.sort.hc), cex = 0.9, label.offset = 1, type = "unrooted")
bequiet <- dev.off()

# saving distance matrix
write(d, file = paste(getwd(),"/tables/distance_matrix.txt",sep = ""), ncolumns = species)

# clearing loaded objects 
rm(list = ls())