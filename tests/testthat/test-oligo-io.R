# test-oligo-io.R -- notation parsing and round-trips.

test_that("triplet notation round-trips through parse_input()/format_triplet()", {
  triplet <- "Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm"
  sp <- parse_input(triplet)
  expect_equal(format_triplet(sp), triplet)
  expect_equal(sp$n, 16L)
})

test_that("three-line entry (bases/sugars/linkages) matches the SEQUENCE_GUIDE worked example", {
  sp <- parse_three_line("TSASTTTSATAATGSTGG", "eeeeeeeeeeeeeeeeee", "sssssssssssssssss")
  expect_equal(sp$n, 18L)
  expect_true(all(sp$sugars == "e"))
  expect_true(all(sp$linkages[1:17] == "s"))
})

test_that("format_biopharma_fasta() output round-trips back through parse_input() (auto-detected)", {
  sp <- parse_input("Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm")
  fasta <- format_biopharma_fasta(sp, "test_oligo")
  expect_true(grepl("^>test_oligo", fasta))

  # parse_fasta() (R/oligo_io.R) must recover the same spec that was
  # exported -- the whole point of adding an import path alongside the
  # existing export is that the two are symmetric. The FASTA record
  # (header + sequence line together) is passed straight to parse_input(),
  # which must auto-detect it as FASTA from the leading '>'.
  reparsed <- parse_input(fasta)
  expect_equal(reparsed$notation, "fasta")
  expect_equal(reparsed$bases, sp$bases)
  expect_equal(reparsed$sugars, sp$sugars)
  expect_equal(reparsed$linkages, sp$linkages)
})

test_that("a FASTA header's terminal-conjugate comment round-trips through parse_fasta()", {
  sp <- parse_input("Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm")
  sp$conj5 <- "GalNAc3"
  sp$conj3 <- "cholesterol"
  fasta <- format_biopharma_fasta(sp, "test_oligo")
  reparsed <- parse_input(fasta)
  expect_equal(reparsed$conj5, "GalNAc3")
  expect_equal(reparsed$conj3, "cholesterol")
})

test_that("validate_spec() rejects an unknown code", {
  bad_spec <- list(bases = c("A", "Z"), sugars = c("d", "d"),
                    linkages = c("s", NA), conj5 = "none", conj3 = "none")
  expect_error(parse_input(bad_spec), "Unknown base code")
})
