# Parameters
sample_parameters <- c(0.51, 5.19, 0.13, 0.1,
                       0.5, 0, 1, 1, 0.99, 0,
                       0, 0, 0, 0, 0.163, 0,
                       0, 0, 0, 0, 0,
                       0.158, -0.002, 0.17, -0.005, 0.067, 0,
                       -1.44, -2, -6.9,
                       0.0042, 0.0015,
                       -2.55, 1.24, 0.25)

parameter_names <- c('aA', 'aW', 'aE', 'aN',
                     'pWA', 'pEA', 'pNA', 'pAW', 'pEW', 'pNW',
                     'pAE', 'pWE', 'pNE', 'pAN', 'pWN', 'pEN',
                     'w1', 'w2', 'w3', 'w4', 'w5',
                     'b1', 'b2', 'bN1', 'bN2', 'bH1', 'bH2',
                     'g', 'gN', 'gH',
                     'pH', 'aH',
                     'th1', 'th2', 'r')

names(sample_parameters) <- parameter_names

usethis::use_data(sample_parameters)
