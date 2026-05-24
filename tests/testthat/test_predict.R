devtools::load_all()

tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
  BEH ~ INT:PBC
'


testthat::expect_no_condition({
  TPB2 <- TPB
  TPB2[23:150, "att1"] <- NA

  fit_fiml <- modsem(tpb, TPB2, method = "lms", nodes = 32, missing = "fiml")

  ebm     <- modsem_predict(fit_fiml, method = "EBM", se = TRUE)
  ml      <- modsem_predict(fit_fiml, method = "ML", type = "lv")
  ml_all  <- modsem_predict(fit_fiml, method = "ML", type = "all")

  # newdata: use first 20 rows (may contain missing)
  ebm_new <- modsem_predict(fit_fiml, newdata = TPB2[1:20, ], method = "EBM")
  ml_new  <- modsem_predict(fit_fiml, newdata = TPB2[1:20, ], method = "ML")
})


testthat::expect_no_condition({
  fit_qml <- modsem(tpb, TPB, method = "qml")

  ebm     <- modsem_predict(fit_qml, method = "EBM", se = TRUE)
  ml      <- modsem_predict(fit_qml, method = "ML")
  ml_all  <- modsem_predict(fit_qml, method = "ML", type = "all")

  # newdata: use last 20 rows
  ebm_new <- modsem_predict(fit_qml, newdata = TPB[181:200, ], method = "EBM", type = "lv")
  ml_new  <- modsem_predict(fit_qml, newdata = TPB[181:200, ], method = "ML")
})


tpb.comp <- '
  ATT <~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH <~ b1 + b2

  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
'

testthat::expect_no_condition({
  fit_comp <- modsem(tpb.comp, TPB, method = "lms")

  ebm     <- modsem_predict(fit_comp, method = "EBM")
  ml      <- modsem_predict(fit_comp, method = "ML", type = "lv")

  ebm_new <- modsem_predict(fit_comp, newdata = TPB[1:20, ], method = "EBM", type = "all")
  ml_new  <- modsem_predict(fit_comp, newdata = TPB[1:20, ], method = "ML")
})


testthat::test_that("modsem_predict.modsem_da formats SE attributes like lavaan for one pattern", {
  fit_qml <- modsem(tpb, TPB[1:80, ], method = "qml", calc.se = FALSE)

  lv  <- modsem_predict(fit_qml, method = "EBM", type = "lv",  se = TRUE)
  ov  <- modsem_predict(fit_qml, method = "EBM", type = "ov",  se = TRUE)
  all <- modsem_predict(fit_qml, method = "EBM", type = "all", se = TRUE)

  testthat::expect_equal(dim(attr(lv, "se")),  c(1L, ncol(lv)))
  testthat::expect_equal(dim(attr(ov, "se")),  c(1L, ncol(ov)))
  testthat::expect_equal(dim(attr(all, "se")), c(1L, ncol(all)))

  testthat::expect_true(is.matrix(attr(lv, "acov")))
  testthat::expect_true(is.matrix(attr(ov, "acov")))
  testthat::expect_true(is.matrix(attr(all, "acov")))

  testthat::expect_null(attr(lv, "vcov"))
  testthat::expect_null(attr(lv, "pattern"))
  testthat::expect_null(attr(attr(lv, "acov"), "patterns"))
})


testthat::test_that("modsem_predict.modsem_da formats SE attributes by missingness pattern", {
  TPB2 <- TPB[1:80, ]
  TPB2[1:10, "att1"] <- NA

  fit_fiml <- modsem(
    tpb, TPB2, method = "lms", nodes = 32, missing = "fiml",
    calc.se = FALSE
  )

  lv <- modsem_predict(fit_fiml, method = "EBM", type = "lv", se = TRUE)

  testthat::expect_equal(dim(attr(lv, "se")), dim(lv))
  testthat::expect_null(attr(lv, "vcov"))
  testthat::expect_null(attr(lv, "pattern"))

  acov <- attr(lv, "acov")
  testthat::expect_true(is.list(acov))
  testthat::expect_equal(names(acov), c("pattern1", "pattern2"))
  testthat::expect_true(all(vapply(acov, is.matrix, logical(1))))

  patterns <- attr(acov, "patterns")
  testthat::expect_equal(dim(patterns)[1], 2L)
  testthat::expect_equal(rownames(patterns), c("pattern1", "pattern2"))
  testthat::expect_true("att1" %in% colnames(patterns))
  testthat::expect_true(any(!patterns[, "att1"]))

  lv_grouped <- modsem_predict(
    fit_fiml, method = "EBM", type = "lv", se = TRUE,
    drop.list.single.group = FALSE
  )

  acov_grouped <- attr(lv_grouped, "acov")
  testthat::expect_equal(names(acov_grouped), "group1")
  testthat::expect_equal(names(acov_grouped$group1), c("pattern1", "pattern2"))
  testthat::expect_equal(
    attr(acov_grouped$group1, "patterns"),
    patterns
  )
})
