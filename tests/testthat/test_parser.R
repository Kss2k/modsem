devtools::load_all()

testthat::expect_error({
  modsemify("y ~ x + z + (x + z | cluster)", parentheses.as.string = FALSE)
})

testthat::expect_no_condition({
  modsemify("y ~ x + z + (x + z | cluster)", parentheses.as.string = TRUE)
})

testthat::expect_true(!PARSER_SETTINGS$parentheses.as.string)
