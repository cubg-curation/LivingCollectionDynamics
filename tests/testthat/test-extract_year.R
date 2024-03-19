test_that("extract_year works", {
  expect_equal(extract_year(c('2023-04-28', '28/04/2023','26th April 2023','Apr. 2023', '2023-04', '2023')), rep(2023,6))
  expect_equal(extract_year(c('1000', '0000')), c(1000,0))

})
