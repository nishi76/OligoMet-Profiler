# test-custom-chemistry-validation.R -- the Custom Chemistry table's guard
# rail. An editable DT table accepts any string; is_valid_formula_string()
# is what stands between a typo'd formula and it silently reaching
# build_dictionary(overrides = ...), which would otherwise propagate a
# wrong mass into every downstream metabolite without any error at all.

test_that("valid element formulas pass", {
  expect_true(is_valid_formula_string("C7H12O4"))
  expect_true(is_valid_formula_string("H2O"))
  expect_true(is_valid_formula_string("Na"))
  expect_true(is_valid_formula_string(" C7H12O4 "))  # surrounding whitespace is trimmed
})

test_that("malformed or nonsense formulas are rejected", {
  expect_false(is_valid_formula_string("not a formula"))
  expect_false(is_valid_formula_string(""))
  expect_false(is_valid_formula_string("  "))
  expect_false(is_valid_formula_string("c7h12o4"))     # wrong case
  expect_false(is_valid_formula_string("Q5"))           # not a real element
  expect_false(is_valid_formula_string("C7H12O4!"))     # stray punctuation
  expect_false(is_valid_formula_string("C7 H12"))       # internal space
  expect_false(is_valid_formula_string(NA))
})

test_that("build_dictionary() overrides only accept formulas that already validate", {
  # This mirrors .overrides_from_table() in app.R: a row is only turned into
  # an override once its Formula cell passes is_valid_formula_string(), so
  # build_dictionary() itself is never asked to interpret a malformed string.
  ok_dict <- build_dictionary(overrides = list(myBase = list(formula = "C10H13N5O4", kind = "base")))
  expect_equal(unname(ok_dict[["myBase"]]$formula[["C"]]), 10)
  expect_equal(ok_dict[["myBase"]]$kind, "base")
})
