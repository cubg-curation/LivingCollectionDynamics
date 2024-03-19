test_that("sanitise_name works", {
  expect_equal(sanitise_name('hypericum adenotrichum var. myriotrichum'), 'Hypericum adenotrichum var. myriotrichum')
  expect_equal(sanitise_name('TRIGONELLA smyrnaea'), 'Trigonella smyrnaea')
  expect_equal(sanitise_name('× orchiserapias bevilacquae'), '× Orchiserapias bevilacquae')
  expect_equal(sanitise_name('× ORCHiserapias bevilacquae'), '× Orchiserapias bevilacquae')
  expect_equal(sanitise_name('+ ORCHiserapias bevilacquae'), '+ Orchiserapias bevilacquae')
  expect_equal(sanitise_name('x ORCHiserapias bevilacquae'), '× Orchiserapias bevilacquae')
  expect_equal(sanitise_name('Aruncus dioicus var acuminatus'), 'Aruncus dioicus var. acuminatus')
  expect_equal(sanitise_name('Osmanthus fragrans f aurantiacum'), 'Osmanthus fragrans f. aurantiacum')
  expect_equal(sanitise_name('Picea abies v columnaris'), 'Picea abies var. columnaris')
  expect_equal(sanitise_name('Zaleya galericulata subsp galericulata'), 'Zaleya galericulata subsp. galericulata')
  expect_equal(sanitise_name('Dryopteris x furadensis'), 'Dryopteris × furadensis')
  expect_equal(sanitise_name('Halimium X pauanum'), 'Halimium × pauanum')



})
