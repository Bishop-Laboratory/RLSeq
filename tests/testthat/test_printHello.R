
test_that(desc = "Test someFunction", {
  expect_null(
    RSeqR::printHello("world")
  )
  expect_output(
    RSeqR::printHello("world"), "Hello world"
  )
})


