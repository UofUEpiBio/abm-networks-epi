library(data.table)
library(slurmR)
library(epiworldR)

# Listing the files under data/sims
simfiles <- list.files(
  "data/sims", pattern = "-sim-[0-9]+\\.rds$",
  full.names = TRUE
  )

ncores <- 40

# Set the slurmR options
SB_OPTS <- list(
  account    = "vegayon-np",
  partition  = "vegayon-shared-np",
  mem        = "8G"
  )

# Extracting the data
res <- Slurm_lapply(simfiles, \(fn) {

  dat <- readRDS(fn)

  # Extracting the reproductive number
  dat <- data.table(dat$repnum)
  dat[, simfile := fn]

}, njobs = ncores, sbatch_opt =  SB_OPTS, job_name = "03-rt-data") |>
  rbindlist()

# Getting the nettype
res[, nettype := gsub(".+/[0-9]+-(er|sf|ergm|degseq|swp[0-9]{2})-sim-.+", "\\1", simfile)]

res[, table(nettype)]

# Getting the simulation id
res[, simid := gsub(".+/[0-9]+-(er|sf|ergm|degseq|swp[0-9]{2})-sim-([0-9]+)\\.rds", "\\2", simfile)]

# Saving the data
fwrite(res, "data/03-rt-data.csv.gz", compress = "auto")
