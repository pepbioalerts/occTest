require (occTest)

# setup test data
occTest_output <- readRDS (system.file('ext/output_occTest.rds',package = 'occTest'))
filtered_dataset <- occFilter (occTest_output)
length (filtered_dataset)
names (filtered_dataset)
class(filtered_dataset)

# run test to check outputs of occFilter are correct
test_that("occFilter outputs", {
  expect_equal(length (filtered_dataset), 3)
  expect_equal( c('occFilter','list'),class(filtered_dataset))
  expect_in(names(filtered_dataset),c('filteredDataset','summaryStats','rule'))
  expect_in('data.frame',class(filtered_dataset$filteredDataset))
  expect_equal('data.frame',class(filtered_dataset$summaryStats))
  expect_equal('data.frame',class(filtered_dataset$rule))
})