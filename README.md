# phylogenetics
Project work for producing phylogenetic trees based on cluster frequency analysis

**Bayesian Tree Model Notes**

*mcmc – Markov Chain Monte Carlo*

This script uses Markov Chain Monte Carlo with a Metropolis (likelihood ratio based) update step to converge on distribution parameters for the BDI model – λ, r, θ, ε, α, β, taking the clusters with 1000 greatest frequencies as the data. Produces estimates of lambda and r for each cluster, as well as estimates of the parameters for the distributions that these are drawn from.  As yet there is no formal measure of convergence, nor is the distribution heated in any way to encourage movement. This hasn’t appeared to be a problem however – to review this check plots.

*bamc – Bayesian Analysis Metric Calculator*

Building on the estimates from mcmc, this script produces a table which is later used to assign probabilities to cluster matches in tbbm. It is here that we cut down to the number of clusters that will be used for the construction of the tree – as you can see this script is fairly simple. One thing to note is that when the unconditional probability of a specific data point is calculated we use an empirical approximation. The exact integral can be performed, however it is VERY computationally expensive – as a 2D integral over an unbounded region with very small function values it is also prone to large computational errors. It is because of this approximation that probabilities have to be tweaked in the next section, by smoothing out and giving each cluster frequency band it can lie within when assigning probability through the ecdf.

*tbbm - Tree Builder by Bayesian Metric*

This script is where the bulk of the computation happens, taking the data that we had previously and turning it into a phylogenetic tree.
  
  **Part 1)	Sorting the cluster frequency data**
  Here we take the data and sort it into “generation groups” for each cluster. This is done by taking the frequency data, seeing which lie closest to each other and then “grouping them”. In this way we build up an array holding the information of how closely grouped the frequency data is, and which species are “closest” to each other. It is worth noting that not all clusters are counted equally in this step – we take those with the highest birth rates r and those with the lowest transfer rates lambda and give these the highest weighting working down – weightings are assigned uniformly across generations per cluster. Take number of generations = number of species as this is of the correct order of “attachments” that must occur. This sorting is neither a purely quantitative or stochastic sort, rather it takes elements from both in an attempt to produce accurate trees with minimal further calculation. 
  
  **Part 2) Analysing the sorted array**
   Fed into dendrogram plotting package after summing over clusters and generations, and normalising to give a comparable distance measure – this is what produces the plots. The distance scale on the plot axis represents how closely we have to “look” before two species are distinguished: a score of 0 would imply that two species are the same, whilst a score of 1 would imply that two species are so distinct that they could never be ‘paired’ using this method. Slightly more sophisticated procedure take the species and cluster them together one by one based on how many clusters they each share, assigning a probability to each pairing, and forming from each pair a placeholder species which is used as the representative when further associating a pairing with another species or even another pair. Hybrid spotting is then performed by taking each species and pairing all the others to see whether they can account for the cluster profile of a daughter species between them. Most likely parents for each cluster are then recorded, as well as the proportion of the child’s cluster profile that they can account for.

*Probability Assignment to Pairing:*

  Calculate probability assigned by defining a sister species event for a specific cluster as:
  
  > S = { (λ’, r’): |λ’ – λ0 |< 0.1 & |r’ – r0| < 0.5}
  
  Where λ0 is the ‘true’ value of λ estimated by mcmc, (likewise for r0). Probabilities then given are probabilities that λ and r lie in the set S based solely on the data from the two samples (x1, x2) to be paired. We then take the mean of these probabilities over all agreeing clusters to get a measure of how “close” two species are.  

*Using the scripts:*

*	Best done in Windows – use batch file run_analysis and it will do everything for you automatically (requires you to find Rscript.exe, but this is fairly easy)

*	Otherwise have to use setwd(“path”) for the location of the cluster_analysis folder

*	Input file has to be in location cluster_analysis/comparativeAnalysis_table_clusters_counts.txt, and must be a tab delimited text file. 

*	Make sure you have moved or closed pdf plots before running the scripts

*	Outputs are plots and several tables of note:
*	
    ..*	hybrid_list tries to find hybrids in the data – and gives candidate parents as well as the proportion of its cluster profile inherited from the parents.

    ..*	tree_list gives the text based generation by generation construction of the tree, with the matching probabilities for each step, as well as a list of the clusters that don’t agree with the pairing.

    ..*	distance_matrix gives the distance from perfect match for each pairing possible

    ..*	cluster_number contains the number of clusters used for the construction of the tree.
