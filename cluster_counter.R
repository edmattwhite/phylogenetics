
# cluster counter
# author: Edward White
# date: 31/08/2016

# Purpose: find and save the number of clusters to use for tbbm and bamc 

cluster.table <- read.table(paste(getwd(),"/comparativeAnalysis_table_clusters_counts.txt",sep=""), header = T)
cluster.table <- cluster.table[,-1]

for(i in 1:dim(cluster.table)[1]){
	if(sum(cluster.table[i,]) < 0.001*sum(cluster.table[1,])){break}
}

if(i > 1000){
i <- 1000
}

num <- i 
write.table(num, file = paste(getwd(),"/tables/cluster_number.txt", sep = ""))
