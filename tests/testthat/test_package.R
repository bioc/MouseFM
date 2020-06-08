# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


test_that("Package works properly", {
  
  
  expect_equal(length(prio("chr1",
    start = 5000000, end = 5001000, strain1 = "C57BL_6J",
    strain2 = "AKR_J"
  )$genotypes), 0)

  
  expect_equal(length(prio("chr1",
    start = 5000000, end = 6000000, strain1 = "C57BL_6J",
    strain2 = "AKR_J", return_obj = "granges"
  )$genotypes), 683)
  
  
  expect_equal(length(vis_reduction_factors(prio("chr1",
                                                 start = 5000000, end = 6000000, strain1 = "C57BL_6J",
                                                 strain2 = "AKR_J"
  )$genotypes, 
  prio("chr1",
       start = 5000000, end = 6000000, strain1 = "C57BL_6J",
       strain2 = "AKR_J"
  )$reduction, 3)), 3)
  
  
  expect_equal(length(vis_reduction_factors(prio("chr1",
                           start = 5000000, end = 6000000, strain1 = "C57BL_6J",
                           strain2 = "AKR_J", return_obj = "granges"
  )$genotypes, 
  prio("chr1",
       start = 5000000, end = 6000000, strain1 = "C57BL_6J",
       strain2 = "AKR_J", return_obj = "granges"
  )$reduction, 3)), 3)

  
  expect_equal(nrow(fetch("chr1", start = 5000000, end = 5001000)), 33)

  
  expect_equal(nrow(fetch("chr1", start = 5000000, end = 5000000)), 0)

  
  expect_equal(length(fetch("chr1", start = 5000000, end = 5000000, return_obj = "granges")), 0)

  
  expect_equal(nrow(finemap("chr1",
    start = 5000000, end = 6000000,
    strain1 = c("C57BL_6J"), strain2 = c("129S1_SvImJ", "129S5SvEvBrd", "AKR_J")
  )), 830)

  
  expect_equal(length(finemap("chr1",
    start = 5000000, end = 6000000,
    strain1 = c("C57BL_6J"), strain2 = c(
      "129S1_SvImJ", "129S5SvEvBrd",
      "AKR_J"
    ), return_obj = "granges"
  )), 830)
})
