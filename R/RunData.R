## Sometimes the package gengraph has not loaded therefore the code for installing it 
# install.packages("gengraph")

## Load all packages
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
    clusters <- gengraph(data_dist, cutoff = c_thres)
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

cl_seq <- make_cluster(distance, c_thres=20)
cl_seq
cl_group <- make_subcluster(cl_seq, mem_number = 2, genetic_data)

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


