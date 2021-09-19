
test_that(desc = "Test report", {
  expect_null(
    RLSeq::report(RLSeq::SRX1025890, outputFile = "report.html")
  )
  file.remove("report.html")
})
