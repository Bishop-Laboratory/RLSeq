
test_that(desc = "Test makeReport", {
  expect_null(
    RSeqR::makeReport(RSeqR::SRX1025890, outputFile = "report.html")
  )
  file.remove("report.html")
})

