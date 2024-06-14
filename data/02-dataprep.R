#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --mem=128GB
#SBATCH --job-name=02-dataprep
#SBATCH --time=24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=george.vegayon@utah.edu

library(network, lib.loc = "/uufs/chpc.utah.edu/common/home/vegayon-group1/george/R/x86_64-pc-linux-gnu-library/4.2/")
library(ergm, lib.loc = "/uufs/chpc.utah.edu/common/home/vegayon-group1/george/R/x86_64-pc-linux-gnu-library/4.2/")
library(igraph)
library(slurmR)
library(intergraph)
library(data.table)

ncores <- 30

# Set the slurmR options
SB_OPTS <- list(
  account    = "vegayon-np",
  partition  = "vegayon-shared-np",
  mem        = "8G"
  )

# Read the Simulated data rds file from data/
message("Loading the networks files")
networks <- list.files(
  "data/graphs", pattern = "[0-9]+-(ergm|sf|swp[0-9]{2}|degseq|er)\\.rds$",full.names = TRUE) 

# File name for each type
fn_ergm   <- "data/02-dataprep-ergm.csv.gz"
fn_igraph <- "data/02-dataprep-igraph.csv.gz"
fn_sim    <- "data/02-dataprep-sim.csv.gz"

# Computing statistics using ERGM
if (!file.exists(fn_ergm)) {
  message("Computing statistics using ERGM")
  S_ergm <- Slurm_lapply(seq_along(networks), \(i) {

    n <- networks[[i]]

    ng <- gsub("graphs/", "graphs_stats/", n)
    ng <- gsub(".rds$", "-ergm.csv.gz", ng)

    # If the file exists, load it and return it
    if (file.exists(ng)) {
      return(fread(ng))
    }

    # Reading the network
    n <- readRDS(n)

    res <- summary_formula(
      n ~ edges + nodematch("gender") + nodematch("grade") + triangles + balance +
        twopath + gwdegree(decay = .25, fixed = TRUE) + isolates
    ) |> as.list() |> as.data.table()

    res[, netfile := networks[[i]]]

    # Saving
    fwrite(res, ng, compress = "auto")

    res
  }, njobs = ncores, sbatch_opt =  SB_OPTS, export = "networks", job_name = "02-dataprep-ergm") |> rbindlist()

  # Pre-appending `ergm_` to all column names, except netfile
  setnames(S_ergm, new = paste0("ergm_", names(S_ergm)))
  setnames(S_ergm, old = "ergm_netfile", new = "netfile")

  fwrite(S_ergm, fn_ergm, compress = "auto")

  message("Done computing statistics using ERGM")
} else {
  message("Loading statistics using ERGM")
  S_ergm <- fread(fn_ergm)
}

head(S_ergm)

# Computing statistics based on igraph ----------------------------------------

if (!file.exists(fn_igraph)) {

  message("Computing statistics based on igraph")
  S_igraph <- Slurm_lapply(networks, \(i) { 

    ng <- gsub("graphs/", "graphs_stats/", i)
    ng <- gsub(".rds$", "-igraph.csv.gz", ng)

    # If the file exists, load it and return it
    if (file.exists(ng)) {
      return(fread(ng))
    }

    # Reading the network and converting it into igraph
    n <- readRDS(i)
    n <- intergraph::asIgraph(n)
    
    n <- data.table(
      netfile         = i,
      modularity      = modularity(cluster_fast_greedy(n)),
      transitivity    = transitivity(n),
      density         = igraph::edge_density(n),
      diameter        = diameter(n),
      avg_path_length = igraph::mean_distance(n),
      avg_degree      = mean(degree(n), na.rm = TRUE),
      avg_betweenness = mean(betweenness(n), na.rm = TRUE),
      avg_closeness   = mean(closeness(n), na.rm = TRUE),
      avg_eigenvector = mean(eigen_centrality(n)$vector, na.rm = TRUE),
      components      = components(n)$no
    )

    # Saving
    fwrite(n, ng, compress = "auto")

    n
    
  }, njobs = ncores, sbatch_opt =  SB_OPTS, export = "networks", job_name = "02-dataprep-igraph") |> rbindlist()

  head(S_igraph)

  setnames(S_igraph, new = paste0("igraph_", names(S_igraph)))
  setnames(S_igraph, old = "igraph_netfile", new = "netfile")

  fwrite(S_igraph, fn_igraph, compress = "auto")

  message("Done computing statistics based on igraph")
} else {
  message("Loading statistics based on igraph")
  S_igraph <- fread(fn_igraph)
}

# Combining the datasets
S <- merge(S_igraph, S_ergm, by = "netfile", all = TRUE)

S[, nettype := gsub(".+-([a-z0-9]+)\\.rds", "\\1", netfile)]
S[, netid   := gsub(".+/([0-9]+-[a-z0-9]+)\\.rds", "\\1", netfile)]


# Reading simulation results ---------------------------------------------------

if (!file.exists(fn_sim)) {

  message("Processing the simulation results")
  simfiles <- list.files("data/sims", pattern = "-sim-[0-9]+\\.rds$", full.names = TRUE)

  simres <- Slurm_lapply(simfiles, \(fn) {

    x <- readRDS(fn)
    
    res <- tryCatch({
      history <- data.table(x$history)

      # Index of the max prevalence
      peak_idx    <- which.max(x$incidence$Exposed)
      peak_preval <- x$incidence$Exposed[peak_idx]
      peak_time   <- as.integer(rownames(x$incidence)[peak_idx])
      
      # Rt in the peak
      rt_idx     <- with(x$repnum, which.min(abs(peak_time - date)))
      rt         <- x$repnum$avg[rt_idx]
      dispersion <- 1/x$repnum$sd[rt_idx]^2

      # Mean Rt
      r_mean     <- with(x$repnum, sum(avg * n, na.rm = TRUE)/
        sum(n, na.rm = TRUE))

      #  Rt per day
      rt_0 <- with(x$repnum[x$repnum$date == 0,,drop=FALSE], avg)
      rt_1 <- with(x$repnum[x$repnum$date == 1,,drop=FALSE], avg)
      rt_2 <- with(x$repnum[x$repnum$date == 2,,drop=FALSE], avg)
      rt_3 <- with(x$repnum[x$repnum$date == 3,,drop=FALSE], avg)
      rt_4 <- with(x$repnum[x$repnum$date == 4,,drop=FALSE], avg)
      rt_5 <- with(x$repnum[x$repnum$date == 5,,drop=FALSE], avg)

      gentime <- with(x$gentime, sum(avg * n, na.rm = TRUE)/
        sum(n, na.rm = TRUE))

      # Final prevalence
      final_preval <- with(
        x$history, tail(counts[state == "Removed"], 1)
        )

      # Return a data.table
      data.table(
        simfile           = if (length(fn)) fn else NA,
        simid             = if (length(x$simid)) x$simid else NA,
        peak_time         = if (length(peak_time)) peak_time else NA,
        peak_preval       = if (length(peak_preval)) peak_preval else NA,
        rt                = if (length(rt)) rt else NA,
        rt_mean           = if (length(r_mean)) r_mean else NA,
        rt_0              = if (length(rt_0)) rt_0 else NA,
        rt_1              = if (length(rt_1)) rt_1 else NA,
        rt_2              = if (length(rt_2)) rt_2 else NA,
        rt_3              = if (length(rt_3)) rt_3 else NA,
        rt_4              = if (length(rt_4)) rt_4 else NA,
        rt_5              = if (length(rt_5)) rt_5 else NA,
        dispersion        = if (length(dispersion)) dispersion else NA,
        gentime           = if (length(gentime)) gentime else NA,
        final_preval      = if (length(final_preval)) final_preval else NA,
        transmission_rate = if (length(x$param$transmission_rate)) x$param$transmission_rate else NA,
        inc_days          = if (length(x$param$inc_days)) x$param$inc_days else NA,
        recovery_rate     = if (length(x$param$recovery_rate)) x$param$recovery_rate else NA
      )
    }, error=function(e) e)

    if (inherits(res, "error"))
      return(NULL)

    res

  }, njobs = ncores, sbatch_opt =  SB_OPTS, export = "networks", job_name = "02-dataprep-sims") |> rbindlist()

  fwrite(simres, fn_sim, compress = "auto")

  message("Done processing the simulation results")

} else {

  message("Loading the simulation results")
  simres <- fread(fn_sim)
  
}

simres[, netid := gsub(".+/([0-9]+-[a-z0-9]+)-sim.+", "\\1", simfile)]

print(head(S))

print(head(simres))

# Merge the datasets
S <- merge(simres, S, by = "netid", all.x = TRUE, all.y = FALSE)

fwrite(
  S,
  file = "data/02-dataprep-network-stats.csv.gz"
)

message("Done!")
