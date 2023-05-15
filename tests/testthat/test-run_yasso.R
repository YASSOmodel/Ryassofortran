test_that("Simulated carbon is correct", {

  # Expected output
  real_out <- as.matrix(structure(c(
    40.5, 35.2177536379588, 17.5375289220381, 8.17949877248404,
    4.6057473374274, 30.5, 3.7080866621023, 1.82968250966202, 0.840045164509876,
    0.466335508642046, 20.3, 13.393577282265, 5.92172186727114, 1.94116296851444,
    0.619738877088547, 10.3, 20.0008636090304, 25.6552313829448,
    25.0740864612265, 20.4820925829699, 0, 0.493026914441875, 0.954475314649129,
    1.28708059992064, 1.51383849850536), .Dim = c(5L, 5L)))

  # Simulated output
  simulated_out <- run_yasso(
    par = sample_parameters,
    n_runs = sample_data_run$n_runs,
    time = sample_data_run$time,
    temp = sample_data_run$temp,
    prec = sample_data_run$prec,
    init = sample_data_run$init,
    litter = sample_data_run$litter,
    wsize = sample_data_run$wsize,
    leac = sample_data_run$leac
  )

  # Compare
  expect_equal(simulated_out, real_out, tolerance = 1e-4)
})

# Commenting out this steady state test as it is not constructed as intended.
# Should only be run for a single time step, and initial state should not impact
# the steady state output.
#
# test_that("Steady state predictions are correct", {
#
#   # Expected output
#   real_out <- as.matrix(structure(c(
#     40.5, 62.9975581347448, 71.7491909991854, 80.500823863626,
#     89.2524567280666, 30.5, 6.95824142758275, 7.95622030996732, 8.95419919235189,
#     9.95217807473646, 20.3, 25.6575789594938, 28.4730058325039, 31.2884327055141,
#     34.1038595785242, 10.3, 299.547585467295, 329.938996687453, 360.330407907612,
#     390.72181912777, 0, 406.312700833162, 459.969829898441, 513.62695896372,
#     567.284088028999), .Dim = c(5L, 5L)))
#
#   # Simulated output
#   simulated_out <- run_yasso(
#     par = sample_parameters,
#     n_runs = sample_data_run$n_runs,
#     time = sample_data_run$time,
#     temp = sample_data_run$temp,
#     prec = sample_data_run$prec,
#     init = sample_data_run$init,
#     litter = sample_data_run$litter + matrix(1:20, nrow = 4),
#     wsize = sample_data_run$wsize,
#     leac = sample_data_run$leac,
#     sspred = 1L
#   )
#
#   # Compare
#   expect_equal(simulated_out, real_out, tolerance = 1e-4)
# })
