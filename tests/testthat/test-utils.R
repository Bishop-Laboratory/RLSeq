test_that("getRLFSAnno works", {
  expect_s4_class(getRLFSAnno("hg38"), "GRanges")
})

test_that("urlExists works", {
  expect_true(urlExists("www.google.com"))
})

test_that("getChromSizes works", {
  expect_s3_class(getChromSizes("hg38"), class="tbl")
})

test_that("checkRLFSAnno works", {
  expect_true(checkRLFSAnno("hg38"))
})

test_that("checkGenes works", {
  expect_true(checkGenes("hg38"))
})

test_that("checkHomer works", {
  expect_true(checkHomer("hg38"))
})
