## Load all packages
install.packages("gengraph")
lapply(c("phybreak", "coda", "gplots", "phytools", "igraph", "seqinr", "gengraph", "adegenet"), require, character.only = TRUE)
genetic_data <- read.dna(file = "../data/syn_geneticData3.fasta",
                         format = "fasta")
distance <- dist.dna(genetic_data, model = "N", as.matrix = TRUE)
heatmap.2(distance, trace = "none", key.xlab = "Number of SNPs")

plot(density(distance))

## Make a graph object from the distance matrix
mygraph <- graph.adjacency(distance, weighted = TRUE)

## Compute and plot the minimum spanning tree
mstgraph <- as.undirected(minimum.spanning.tree(mygraph))

## Create a layout where nodes with more distance are further away from each other
layout <- layout_with_fr(mstgraph, weights = 1 / E(mstgraph)$weight)

plot(mstgraph,
  layout = layout,
  edge.label = edge_attr(mstgraph, "weight"), # Add distance as edge label
  edge.label.color = "black", # Color of edge labels
  vertex.size = 2, # Size of the nodes
  vertex.label.dist = 1,
  vertex.label.cex = 1.5
) # Add a little space between the node and its label
title(main = "Minimal spanning tree")


###### For subcluster  ######

make_cluster <- function(data_dist, c_thres) {
  if(c_thres<0) {
    stop("The chosen threshold can not be negative ")
  } else
  {
    hist(data_dist[upper.tri(data_dist, diag = FALSE)])
    clusters <- gengraph(data_dist, cutoff = c.thres)
    clusters
    clusters$graph
    plot(clusters$graph)
  }
  return(clusters)
} 

make_subcluster <- function(cl_data, mem_number, seq_data) {
  #' @title Make subclusters 
  #' @description This function will make define the subcluster with the 
  #' membership vale from the input. It returns the distance matrix and 
  #' dataframe with genetic data.
  #' @param cl_data The object with clustering information 
  #' @param mem_number Membership number which has to be in the clustering information
  #' @param seq_data dataframe with the original genetic data 
  #' @return The genetic subcluster data.
  
  if(mem_number>max(cl_data$clust$membership)) {
    stop("The chosen cluster does not exist")
  } else
  {
    
    id1 <- c(labels(cl_data$clust$membership[cl_data$clust$membership == mem_number]))
    id2 <- match(id1, labels(seq_data))
    sub_cluster <- seq_data[c(id2), ]
  }
  
  return(sub_cluster)
}

cl_seq <- make_cluster(distance, 20)
cl_seq
cl_group <- make_subcluster(cl_seq, mem_number = 15, genetic_data)

# Transform the subcluster of genetic data into a distance matrix 
distance_subcluster <- dist.dna(cl_group, model = "N", as.matrix = TRUE)

#Visualize the subcluster of genetic data
mygraph_subcluster <- graph.adjacency(distance_subcluster, weighted = TRUE)

## Compute and plot the minimum spanning tree
mstgraph_subcluster <- as.undirected(minimum.spanning.tree(mygraph_subcluster))
layout_sub <- layout_with_fr(mstgraph_subcluster, weights = 1 / E(mstgraph_subcluster)$weight)

plot(mstgraph_subcluster,
  layout = layout_sub,
  edge.label = edge_attr(mstgraph_subcluster, "weight"), # Add distance as edge label
  edge.label.color = "black", # Color of edge labels
  vertex.size = 2, # Size of the nodes
  vertex.label.dist = 1,
  vertex.label.cex = 1.5
) # Add a little space between the node and its label
title(main = "Minimal spanning tree")

## Compute and plot the Neighbour joining tree
njtree <- nj(distance)

# Plot without using a starting edge and no edge length for better overview
plot(njtree, type = "unrooted", use.edge.length = FALSE) 
title(main = "Neighbour joining tree")
edgelabels(round(njtree$edge.length, 1), frame = "none") # Add distance as edge labels


##############
meta_data <- data.frame(matrix(0, 100, 5))
meta_data$date <- as.Date(seq(as.Date("2024-01-03"), as.Date("2024-12-05"), by = 3)[1:100])
place_vector <- c("Utrecht", "Amsterdam", "Groningen", "Apeldoorn", "Texel", "Middelburg",
                  "Maastricht", "Oss", "Rheden", "Hoorn")
meta_data$place <- c(rep(as.character(place_vector), 10))

data_object <- phybreakdata(
  sequences = genetic_data,
  sample.times = as.Date(meta_data$date),
  host.names = meta_data$place
)

MCMC_init <- phybreak(data_object)

MCMC_init2 <- phybreak(data_object, prior.sample.mean.mean = 10, prior.sample.mean.sd = 1)
# (not runnable as code in this notebook)


## Perform a burn-in of the MCMC
MCMC_state <- burnin_phybreak(MCMC_init, ncycles = 100)

## To run the MCMC chain with 10,000 cycles
## The `nchains = 3` argument tells the model to use 3 parallel chains, in order
## to find the most optimal state

MCMC_state <- sample_phybreak(MCMC_state, nsample = 100, nchains = 3)

## To load the prerun MCMC chain provided as .rds file
print(MCMC_state)

mcmc <- get_mcmc(MCMC_state)
effectiveSize(mcmc)
plot(mcmc[, c("introductions")])

summary(mcmc[, c("mS")])

transtree(MCMC_state, method = "mpc")

## Show all infectors of Farm Utrecht, and their support
infectorsets(MCMC_state, which.hosts = "Utrecht")

plotTrans(MCMC_state, plot.which = "mpc", samplenr = 6)

## Plot the phylogenetic tree with maximum parent credibility
plotPhylo(MCMC_state, plot.which = "mpc", samplenr = 0)

plotPhyloTrans(MCMC_state, plot.which = "mpc")
