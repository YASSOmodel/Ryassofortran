# Amount of rows to run model over
n_runs <- 4

# Time for each run
time <- c(1.01, 2.03, 3.01, 4.00)

# Weather for each run
prec <- rep(814.585, 4)
temp <- c(-1.1, -1.2, 1, 5.6, 10.8, 14.7, 16.8, 16.2, 13, 8.5, 4.5, 1.4)
temp <- matrix(temp, nrow = n_runs, ncol = length(temp), byrow = TRUE)

# Litter input for each run
litter <- matrix(rep(0, len = n_runs * 5), nrow = n_runs)

# Woody size for each run
wsize <- rep(0, 4)

# Leaching for each run
leac <- rep(0, 4)

# Initial carbon for each run
init <- c(40.5, 30.5, 20.3, 10.3, 0)

# Create a list, typeset for fortran
sample_data_run <- list(
  "n_runs" = as.integer(n_runs),
  "time" = as.double(time),
  "temp" = as.matrix(temp),
  "prec" = as.double(prec),
  "init" = as.double(init),
  "litter" = as.matrix(litter),
  "wsize" = as.double(wsize),
  "leac" = as.double(leac)
)

usethis::use_data(sample_data_run, overwrite = TRUE)
