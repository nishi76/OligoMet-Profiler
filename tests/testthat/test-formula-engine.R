# test-formula-engine.R -- the formula/mass engine against published truth.
#
# These are the highest-stakes tests in the suite: chemistry_dict.R and
# mass_isotope.R sit underneath every mass this app reports, so a wrong
# residue formula or a broken condensation step propagates silently into
# every downstream m/z (workbook, report, PRM list, spectral libraries).

test_that("the formula engine reproduces published drug label formulas", {
  ref <- validate_reference(verbose = FALSE)
  expect_true(isTRUE(ref$ok),
              info = "validate_reference() no longer matches the reference oligos' published molecular formulas")
})

test_that("parse_formula() / format_formula() round-trip", {
  for (fs in c("C164H221N44O97P15S15", "C8H10N4O2", "C6H12O6", "S8", "H2O")) {
    f <- parse_formula(fs)
    expect_equal(format_formula(f), fs)
  }
})

test_that("parse_formula() ignores unknown tokens gracefully rather than erroring", {
  expect_equal(unname(parse_formula("")[["C"]]), 0)
})

test_that("charge envelope m/z back-calculates to the same neutral mass at every charge", {
  sp <- parse_input(INOTERSEN_TRIPLET)
  mets <- generate_metabolites(sp, list(oligo_name = "inotersen", endo = FALSE))
  info <- metabolite_mass_info(mets[[1]])

  zs <- 4:12
  neutral <- vapply(zs, function(z) {
    mz <- charge_envelope(info$mono_mass, z)$mz[1]
    mz * z + z * .PROTON
  }, numeric(1))
  expect_lt(max(abs(neutral - info$mono_mass)), 1e-6)
})

test_that("isotope pattern never places a peak below the monoisotopic mass, for either engine", {
  sp <- parse_input(INOTERSEN_TRIPLET)
  mets <- generate_metabolites(sp, list(oligo_name = "inotersen", endo = FALSE))
  info <- metabolite_mass_info(mets[[1]])

  check_one <- function(formula_vec, mono) {
    for (engine in c(TRUE, FALSE)) {
      p <- isotope_pattern(formula_vec, n_top = 60, use_envipat = engine)
      if (is.null(p) || nrow(p) == 0) next
      expect_gte(min(p$mass), mono - 0.05,
                 label = sprintf("min isotope mass (use_envipat=%s)", engine))
      closest <- min(abs(p$mass - mono))
      expect_lt(closest, 0.05,
                label = sprintf("closest peak to monoisotopic mass (use_envipat=%s)", engine))
    }
  }
  check_one(info$formula_vec, info$mono_mass)
  for (fs in c("C8H10N4O2", "C6H12O6", "S8")) {
    check_one(parse_formula(fs), formula_mass(parse_formula(fs)))
  }
})

test_that("each PS -> PO oxidation event shifts mass by exactly one S->O swap", {
  sp <- parse_input(INOTERSEN_TRIPLET)
  mets <- generate_metabolites(sp, list(oligo_name = "inotersen", endo = FALSE))
  ox <- ps_oxidation_series(mets[[1]], max_oxid = 6)
  deltas <- diff(ox$mono_mass)
  expect_equal(unique(round(deltas, 6)), round(-.PS_TO_PO_SHIFT, 6))
})

test_that("build_dictionary() overrides replace a formula without touching sibling entries", {
  base_dict <- build_dictionary()
  ov_dict <- build_dictionary(overrides = list(cEt = list(formula = "C7H12O4", name = "test override")))
  expect_equal(ov_dict[["cEt"]]$formula[["C"]], 7)
  expect_equal(ov_dict[["cEt"]]$formula[["H"]], 12)
  # An unrelated entry must be untouched by the override.
  expect_equal(ov_dict[["d"]]$formula, base_dict[["d"]]$formula)
})

test_that("a malformed custom-chemistry formula string does not silently assemble to an empty formula", {
  # A typo'd formula (lowercase element, or a bare non-element token) must
  # not parse to something chemically meaningless without at least being
  # detectable by is_valid_formula_string() -- see also
  # test-custom-chemistry-validation.R for the validator itself.
  expect_true(is_valid_formula_string("C7H12O4"))
  expect_false(is_valid_formula_string("not a formula"))
  expect_false(is_valid_formula_string(""))
})
