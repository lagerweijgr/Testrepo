## Load all packages
lapply(c("phybreak", "coda", "gplots", "phytools", "igraph","seqinr"), require, character.only = TRUE)


# setwd("C:/Users/lager008/OneDrive - Universiteit Utrecht/Phylobreak")
setwd("C:/Users/lager008//Testrepo/")

genetic.data <- read.dna(file = "Syn_geneticData3.fasta" , format = "fasta")
distance<- dist.dna(genetic.data, model = "N", as.matrix = TRUE) 
heatmap.2(distance, trace = "none", key.xlab = "Number of SNPs")

plot(density(distance))


library(gengraph)
library(adegenet)

cluster <- function(data_dist,c.thres)
{
  
  hist(data_dist[upper.tri(data_dist, diag = F)])
  # install.packages("adegenet")
  # install.packages("gengraph")
  # gengraph(distance)
  clusters <- gengraph(data_dist, cutoff = c.thres)
  # plot(clusters)
  clusters
  clusters$graph
  plot(clusters$graph)
  names(clusters$clust$membership[clusters$clus$membership == 1])
  names(clusters$clust$membership[clusters$clus$membership == 4])
  
  
  ## bigger clusters
  clustering <- clusters$clust$csize
  big.cluster <- which(clustering >=10 )
  
  return(clusters)
}
clusterin.data <- function(cl.data,mem.number,seq.data)
{
  id1 <- c(labels(cl.data$clust$membership[cl.data$clust$membership==mem.number]))
  id2 <- match(id1,labels(seq.data))
  sub.cluster <- seq.data[c(id2),]
  return(sub.cluster)
}

cl.seq <- cluster(distance,20)
cl.seq
cl.group <- clusterin.data(cl.seq,1,genetic.data)


distance.subcluster <- dist.dna(cl.group, model = "N", as.matrix = TRUE) 

## Make a graph object from the distance matrix
mygraph <- graph.adjacency(distance, weighted = TRUE)

## Compute and plot the minimum spanning tree
mstgraph <- as.undirected(minimum.spanning.tree(mygraph))

## Create a layout where nodes with more distance are further away from each other
layout <- layout_with_fr(mstgraph, weights=1/E(mstgraph)$weight)

plot(mstgraph,
     layout = layout,
     edge.label = edge_attr(mstgraph,"weight"),   # Add distance as edge label
     edge.label.color = "black",                  # Color of edge labels
     vertex.size = 2,                             # Size of the nodes
     vertex.label.dist= 1,
     vertex.label.cex = 1.5)                        # Add a little space between the node and its label
title(main = "Minimal spanning tree")



####### For subcluster  ###### 
mygraph.subcluster <- graph.adjacency(distance.subcluster, weighted = TRUE)

## Compute and plot the minimum spanning tree
mstgraph.subcluster <- as.undirected(minimum.spanning.tree(mygraph.subcluster))
layout.sub<- layout_with_fr(mstgraph.subcluster, weights=1/E(mstgraph.subcluster)$weight)

plot(mstgraph.subcluster,
     layout = layout.sub,
     edge.label = edge_attr(mstgraph.subcluster,"weight"),   # Add distance as edge label
     edge.label.color = "black",                  # Color of edge labels
     vertex.size = 2,                             # Size of the nodes
     vertex.label.dist= 1,
     vertex.label.cex = 1.5)                        # Add a little space between the node and its label
title(main = "Minimal spanning tree")

## Compute and plot the Neighbour joining tree
njtree <- nj(distance)

plot(njtree, type = "unrooted", use.edge.length = F)  # Plot without using a starting edge and no edge length for better overview
title(main = "Neighbour joining tree")
edgelabels(round(njtree$edge.length, 1), frame = "none") # Add distance as edge labels


##############
meta.data <- data.frame(matrix(0,100,5))
meta.data$date <- as.Date(seq(as.Date("2024-01-03"),as.Date("2024-12-05"),by=3)[1:100])
meta.data$place <- c(rep(as.character(c("Utrecht","Amsterdam","Groningen","Apeldoorn","Texel","Middelburg","Maastricht","Oss","Rheden","Hoorn")),10))

data.object <- phybreakdata(sequences = genetic.data,
                            sample.times = as.Date(meta.data$date),
                            host.names = meta.data$place)

MCMC_init <- phybreak(data.object)

MCMC_init2 <- phybreak(data.object, prior.sample.mean.mean = 10, prior.sample.mean.sd = 1)
# (not runnable as code in this notebook)


## Perform a burn-in of the MCMC
MCMC_state <- burnin_phybreak(MCMC_init, ncycles = 100)

## To run the MCMC chain with 10,000 cycles
## The `nchains = 3` argument tells the model to use 3 parallel chains, in order
## to find the most optimal state

MCMC_state <- sample_phybreak(MCMC_state, nsample = 100, nchains = 3)

## To load the prerun MCMC chain provided as .rds file
# MCMC_state <- readRDS("mcmc_samples.rds")
print(MCMC_state)

mcmc <- get_mcmc(MCMC_state)
effectiveSize(mcmc)
plot(mcmc[,c("introductions")])

summary(mcmc[,c("mS")])

transtree(MCMC_state, method = "mpc")

## Show all infectors of Farm 8a, and their support
infectorsets(MCMC_state, which.hosts = "F8a")

plotTrans(MCMC_state, plot.which = "mpc", samplenr = 6)

## Plot the phylogenetic tree with maximum parent credibility
plotPhylo(MCMC_state, plot.which = "mpc", samplenr = 0)


plotPhyloTrans(MCMC_state, plot.which = "mpc")

## Select three colorblind friendly colors
colors <- c("#000000", "#E69F00", "#56B4E9")

## Color the farms (hosts) according to the genetic clusters.
host.color <- sapply(unique(data.object$sample.hosts), function(n){
  cluster <- metadata$gen.cluster[metadata$farm == n][1]
  if (cluster == "A") color <- colors[1]
  else if (cluster == "B") color <- colors[2]
  else if (cluster == "D") color <- colors[3]
  return(color)
})

## Plot the transmission tree (the MPC tree) for which the labels are colored.
plotTrans(MCMC_state, plot.which = "mpc", label.col = host.color)
