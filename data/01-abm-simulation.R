#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=abm-simulation-main
#SBATCH --output=abm-simulation-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --mem=32G

library(epiworldR)
library(slurmR)
library(network)

# Set the seed for reproducibility
set.seed(1231)

# Listing all the files under data/graphs and keeping full path
# to the files
graph_files <- list.files("data/graphs", full.names = TRUE)

# Set the parameters
nsims <- 20000
Njobs <- 40L

# Set the slurmR options
SB_OPTS <- list(
    account    = "vegayon-np",
    partition  = "vegayon-shared-np",
    mem        = "8G"
    )


# Model parameters
incubation_days    <- 7
recovery_rates     <- 1/7
contact_rates      <- 14 # Average observed
reproductive       <- 2
transmission_rates <- recovery_rates/(
  contact_rates/reproductive + recovery_rates - 1
  ) # Analytical sol


# Group infectiousness and recovery_rate into a list of length
# nsims featuring one element per simulation
netfiles <- sample(graph_files, nsims, replace = TRUE)
seeds    <- sample.int(.Machine$integer.max, nsims, replace = TRUE)

params <- Map(\(t_rate, r_rate, i, netfile, simid) {
  list(
    netfile           = netfile,
    transmission_rate = t_rate,
    recovery_rate     = r_rate,
    inc_days          = i,
    seed              = seeds[simid],
    simid             = simid
    )
  },
  t_rate  = transmission_rates,
  r_rate  = recovery_rates,
  i       = incubation_days,
  netfile = netfiles,
  simid   = 1:nsims
  )

# Sampling

# opts_slurmR$debug_on()

res <- Slurm_lapply(params, FUN = \(param) {

  message(
    "Initiating model for simid ", param$simid, " and netfile ", param$netfile,
    "...", appendLF = FALSE
    )

  # Generating filename using param$netfile
  fn <- gsub(
    pattern = "\\.rds$",
    replacement = sprintf("-sim-%04i.rds", param$simid),
    x = param$netfile
    )

  fn <- gsub(
    pattern = "data/graphs/",
    replacement = "data/sims/",
    x = fn
  )

  if (file.exists(fn)) {
    message("File ", fn, " already exists, skipping...")
    return(NULL)
  }

  ans <- tryCatch({
    # Creating the SEIR model
    model <- ModelSEIR(
      name              = "SEIR",
      prevalence        = 0.002, # One seed
      transmission_rate = param$transmission_rate,
      incubation_days   = param$inc_days,
      recovery_rate     = param$recovery_rate
    )

    verbose_off(model)

    # Adding the network to model
    net <- readRDS(param$netfile)

    nnodes <- network::network.size(net)
    nedges <- network::network.edgecount(net)
    net <- network::as.edgelist(net)

    if (any(is.na(net))) {
      message("NA in net for simid ", param$simid, " and netfile ", param$netfile)
      print(net)
      return(NULL)
    }

    if (nrow(net) != nedges) {
      message("Missmatch of edgelist with num of edges for simid ", param$simid, " and netfile ", param$netfile)
      print(net)
      return(NULL)
    }

    if (inherits(net, "character")) {
      message("Error for simid ", param$simid, " and netfile ", param$netfile, ". Seems to be a character.")
      print(net)
      return(NULL)
    }

    agents_from_edgelist(
      model,
      size   = nnodes,
      source = net[, 1] - 1L,
      target = net[, 2] - 1L,
      directed = FALSE
      )

    # Running the simulation
    run(model, ndays = 100, seed = param$seed)

    # Get the results
    saveRDS(list(
      simid     = param$simid,
      netfile   = param$netfile,
      history   = get_hist_total(model),
      repnum    = plot_reproductive_number(model, plot = FALSE),
      incidence = plot_incidence(model, plot = FALSE),
      gentime   = plot_generation_time(model, plot = FALSE),
      params    = param
    ), fn, compress = FALSE)

    fn

  }, error = function(e) e)

  if (inherits(ans, "error")) {
    message(
      "Error for netid ", param$netid, " and simid ", param$simid,
      appendLF = FALSE
      )
    return(ans)
  } 

  message("done")

  ans

}, njobs = Njobs, sbatch_opt = SB_OPTS, job_name = "abm-simulation-lapply",
plan = "collect")

message("Saving the results under data/")

# Saving the results under data/
saveRDS(
  res, file = "data/01-abm-simulation-fn.rds", 
  compress = FALSE
  )

