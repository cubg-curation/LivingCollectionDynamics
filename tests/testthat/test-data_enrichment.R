test_that("is_autonym works", {
  expect_equal(is_autonym('Codiaeum variegatum var. variegatum'), TRUE)
  expect_equal(is_autonym('Codiaeum variegatum aff. variegatum'), FALSE)
  expect_equal(is_autonym('Crinum pedunculatum f. purple'), FALSE)
  expect_equal(is_autonym("Haworthia retusa var. retusa 'Jolly Green Giant'"), FALSE)
  expect_equal(is_autonym("Indet. Poaceae"), FALSE)
  expect_equal(is_autonym("Crataegus × media var. media"), FALSE)
  expect_equal(is_autonym("Haworthia coarctata var. coarctata forma chalwinii"), FALSE)
  expect_equal(is_autonym("Acer saccharum ssp. saccharum"), TRUE)
  expect_equal(is_autonym("Rhododendron augustinii subsp. augustinii [Tower Court form]"), FALSE)
  expect_equal(is_autonym("Camellia sinensis var. sinensis f. rosea"), FALSE)
  expect_equal(is_autonym("Grevillea exul subsp. exul var. exul"), TRUE)
  expect_equal(is_autonym('Codiaeum variegatum nothosubsp. variegatum'), TRUE)
})

test_that("get_accepted_plant() works", {
  expect_equal(get_accepted_plant('582307-1'), c("Campomanesia thea", "44106-2"))
  expect_equal(get_accepted_plant('546134-4'), c(NA,NA))

})
