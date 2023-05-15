test_that("Simulated carbon is correct for calibration", {

  # Expected output
  real_out <- as.matrix(structure(c(
    35.2177536379588, 24.4092013024, 17.7015601492804,
    13.2497972896469, 3.7080866621023, 2.55822853619419, 1.84706480081234,
    1.37558370901776, 13.393577282265, 8.85051681596413, 5.99194784098137,
    4.08418742191974, 20.0008636090304, 23.8914783565135, 25.6244546620003,
    26.1308232027941, 0.493026914441875, 0.763746117881212, 0.949634863889708,
    1.08913120892361), .Dim = 4:5))

  # Simulated output
  simulated_out <- calibrate_yasso(
    par = sample_parameters,
    n_runs = sample_data_cal$n_runs,
    time = sample_data_cal$time,
    temp = sample_data_cal$temp,
    prec = sample_data_cal$prec,
    init = sample_data_cal$init,
    litter = sample_data_cal$litter,
    wsize = sample_data_cal$wsize,
    leac = sample_data_cal$leac
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
#     62.9975581347448, 71.7491909991854, 80.500823863626,
#     89.2524567280666, 6.95824142758275, 7.95622030996732, 8.95419919235189,
#     9.95217807473646, 25.6575789594938, 28.4730058325039, 31.2884327055141,
#     34.1038595785242, 299.547585467295, 329.938996687453, 360.330407907612,
#     390.72181912777, 406.312700833162, 459.969829898441, 513.62695896372,
#     567.284088028999), .Dim = 4:5))
#
#   # Simulated output
#   simulated_out <- calibrate_yasso(
#     par = sample_parameters,
#     n_runs = sample_data_cal$n_runs,
#     time = sample_data_cal$time,
#     temp = sample_data_cal$temp,
#     prec = sample_data_cal$prec,
#     init = sample_data_cal$init,
#     litter = sample_data_cal$litter + matrix(1:20, nrow = 4),
#     wsize = sample_data_cal$wsize,
#     leac = sample_data_cal$leac,
#     sspred = 1L
#   )
#
#   # Compare
#   expect_equal(simulated_out, real_out, tolerance = 1e-4)
# })
