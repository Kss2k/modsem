devtools::load_all()

testthat::expect_error({
  modsemify("y ~ x + z + (x + z | cluster)", parentheses.as.string = FALSE)
})

testthat::expect_no_condition({
  modsemify("y ~ x + z + (x + z | cluster)", parentheses.as.string = TRUE)
})

testthat::expect_no_condition({
  modsemify("y ~ x + z + # comment with a Bunch of stuff \\n here...*ø
            u + w # line continues here")
})

testthat::expect_no_condition({
  modsemify("
    level: 1
    y ~ x + z
    level: 2
    u ~ w
  ")
})

testthat::expect_no_condition({
  modsemify("y ~ x # trailing comment without a newline at the end")
})

testthat::expect_no_condition({
  modsemify("level: 1 \r\ny ~ x\nlevel: 2\t\r\nu ~ w")
})

testthat::expect_no_condition({
  modsemify("
    level: within
    y ~ x
    level: between
    u ~ w
  ")
})

testthat::expect_true(grepl("a#b", modsemify("lab := 'a#b'\n")$rhs))

testthat::expect_true(!PARSER_SETTINGS$parentheses.as.string)
