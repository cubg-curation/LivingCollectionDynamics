test_that("author_check works", {
  expect_equal(author_check('(G.Forst.) Steud.','(G.Forst.) Endl.'), 'Partial')
  expect_equal(author_check('Wall. C.B.Clarke','C.B.Clarke'),'Partial')

})
