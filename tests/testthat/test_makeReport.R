
test_that(desc = "Test makeReport", {
  expect_null(
    RSeqR::makeReport(RSeqR::SRX1025890, output_file = "report.html")
  )
  file.remove("report.html")
})

