#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --ntasks=31
#SBATCH --mem=100GB
#SBATCH --job-name=abm-net-02-dataprep
#SBATCH --mail-type=all
#SBATCH --mail-user=george.vegayon@utah.edu

library(igraph)
library(network)
# library(netplot)
library(intergraph)
library(ergm)
library(data.table)

ncores <- 30

# Read the Simulated data rds file from data/
networks <- readRDS("data/Simulated_1000_networks.rds")
networks <- unclass(networks)

# # Taking a sample of 100 networks
# set.seed(123)
# networks <- networks[sample(1:1000, 100)]

# Computing statistics using ERGM
S_ergm <- parallel::mclapply(networks, \(n) {
  summary_formula(
    n ~ edges + nodematch("gender") + nodematch("grade") + triangles + balance +
      twopath + gwdegree(decay = .25, fixed = TRUE) + isolates
  ) |> as.list() |> as.data.table()
}, mc.cores = ncores) |> rbindlist()

head(S_ergm)

# Computing statistics based on igraph
S_igraph <- parallel::mclapply(networks, \(n) { 
  n <- asIgraph(n)
  
  data.table(
    modularity      = modularity(cluster_fast_greedy(n)),
    transitivity    = transitivity(n),
    density         = igraph::edge_density(n),
    diameter        = diameter(n),
    avg_path_length = igraph::mean_distance(n),
    avg_degree      = mean(degree(n)),
    avg_betweenness = mean(betweenness(n)),
    avg_closeness   = mean(closeness(n)),
    avg_eigenvector = mean(eigen_centrality(n)$vector),
    components      = components(n)$no
  )
}, mc.cores = ncores) |> rbindlist()

setnames(S_ergm, new = paste0("ergm_", names(S_ergm)))
setnames(S_igraph, new = paste0("igraph_", names(S_igraph)))

# Combining the datasets
S <- cbind(S_igraph, S_ergm)
S[, netid := seq_len(.N)]

# Reading simulation results ---------------------------------------------------
simres <- readRDS("data/01-abm-simulation.rds")

simres <- lapply(simres, \(x) {

  if (inherits(x, "error"))
    return(NULL)

  # Computing the peak prevalence
  peak_preval <- with(x$history, counts[state == "Infected"])
  peak_time   <- which.max(peak_preval)
  peak_preval <- peak_preval[peak_time]
  rt          <- with(x$repnum, sum(avg * n, na.rm = TRUE)/
    sum(n, na.rm = TRUE))

  dispersion  <- with(x$repnum, sum(sd * n, na.rm = TRUE)/
    sum(n, na.rm = TRUE)) # Need to work on this

  dispersion  <- 1/dispersion^2

  gentime     <- with(x$gentime, sum(avg * n, na.rm = TRUE)/
    sum(n, na.rm = TRUE))

  # Final prevalence
  final_preval <- with(
    x$history, tail(counts[state == "Removed"], 1)
    )

  # Return a data.table
  data.table(
    netid        = x$netid,
    peak_time    = peak_time,
    peak_preval  = peak_preval,
    rt           = rt,
    dispersion   = dispersion,
    gentime      = gentime,
    final_preval = final_preval
  )

}) |> rbindlist()

# Merge the datasets
S <- merge(S, simres, by = "netid", all.x = TRUE)

fwrite(
  S,
  file = "data/02-dataprep-network-stats.csv.gz"
)


