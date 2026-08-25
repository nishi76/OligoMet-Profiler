# =============================================================================
# about.R
# Authorship and the research-use-only disclaimer, in one place.
#
# Every surface that carries attribution -- the Shiny dashboard, the HTML/PDF
# report, the Excel workbook, the MGF/MSP spectral libraries, the README and
# the guides -- draws its text from the constants below, so the wording stays
# identical everywhere and only has to be changed once.
#
# Deliberately NOT added to the Thermo acquisition CSVs (MS1 inclusion list,
# MS2 PRM target list, fragment reference): those are strict column-layout
# files that the Method Editor and Skyline parse positionally, and a comment
# or banner line would break the import.
# =============================================================================

## ---- Authorship -------------------------------------------------------------
OLIGOMET_AUTHOR <- "Nishikant Wase, PhD"
OLIGOMET_AUTHOR_EMAIL <- "nishikant.wase@gmail.com"
OLIGOMET_AUTHOR_ROLE <- "Author and developer"
OLIGOMET_URL <- "https://github.com/nishi76/OligoMet-Profiler"

# "Nishikant Wase, PhD <nishikant.wase@gmail.com>"
oligomet_author_line <- function(email = TRUE) {
  if (email) sprintf("%s <%s>", OLIGOMET_AUTHOR, OLIGOMET_AUTHOR_EMAIL)
  else OLIGOMET_AUTHOR
}

## ---- Disclaimer -------------------------------------------------------------
# One-line form, for places with no room for more (file headers, footers).
OLIGOMET_DISCLAIMER_SHORT <- paste0(
  "FOR RESEARCH USE ONLY. Not for diagnostic, clinical, or regulatory ",
  "submission use. All values are computed predictions, not measurements. ",
  "Provided without warranty; the author accepts no liability for their use.")

# Full statement. Plain text, one element per paragraph, so each surface can
# wrap or format it however it likes.
OLIGOMET_DISCLAIMER <- c(
  paste0(
    "Research use only. OligoMetProfiler is a research tool. It is not a ",
    "medical device, and it is not intended or validated for diagnostic ",
    "use, for clinical decision-making, for patient care, for quality ",
    "control release testing, or for inclusion in a regulatory submission. ",
    "Any such use is outside the scope of this software and is done ",
    "entirely at the user's own risk."),
  paste0(
    "Everything this software reports is a prediction. Metabolite ",
    "libraries, molecular formulas, monoisotopic and average masses, ",
    "charge envelopes, isotope patterns, MS/MS fragment ions, oxidation ",
    "series, acquisition target lists and spectral libraries are computed ",
    "from a chemistry dictionary and a set of assumptions about how ",
    "oligonucleotides are metabolised and fragment. They are not ",
    "measurements, and they are not evidence that any species exists in a ",
    "sample. Every assignment must be confirmed experimentally before it ",
    "is acted on."),
  paste0(
    "Some formulas are best estimates. Dictionary entries flagged ",
    "verify = TRUE have not been checked against an independently measured ",
    "mass. They are internally consistent, but the user is responsible for ",
    "confirming any of them against a known standard before relying on a ",
    "result that depends on one."),
  paste0(
    "No warranty and no liability. This software is provided \"as is\", ",
    "without warranty of any kind, express or implied, including but not ",
    "limited to the warranties of merchantability, fitness for a ",
    "particular purpose and non-infringement, as set out in the MIT ",
    "licence under which it is distributed. In no event shall the author ",
    "be liable for any claim, damages or other liability -- including any ",
    "direct, indirect, incidental, consequential or special loss, any loss ",
    "of data, any wasted instrument time or materials, or any erroneous ",
    "scientific conclusion -- arising from or in connection with this ",
    "software or its use. The user is solely responsible for verifying that ",
    "the software is fit for their intended purpose and for validating ",
    "every result it produces."),
  paste0(
    "No affiliation or endorsement. The author is not affiliated with, and ",
    "this software is not endorsed by, any instrument vendor, software ",
    "vendor or pharmaceutical company named anywhere in this package. ",
    "Product, instrument and drug names are used for identification and ",
    "interoperability only and remain the trademarks of their respective ",
    "owners. The bundled reference drugs are worked examples drawn from ",
    "published literature and public regulatory filings; they are not ",
    "supplied, reviewed or approved by their manufacturers."),
  paste0(
    "Opinions and outputs are the author's own and do not represent any ",
    "employer or institution.")
)

## ---- Convenience -------------------------------------------------------------
# Print authorship, version and the full disclaimer to the console.
# Wrapped at `width` so it stays readable in a terminal.
oligomet_about <- function(width = 76) {
  ver <- tryCatch(as.character(utils::packageVersion("OligoMetProfiler")),
                  error = function(e) "(source)")
  cat("OligoMetProfiler ", ver, "\n", sep = "")
  cat(OLIGOMET_AUTHOR_ROLE, ": ", oligomet_author_line(), "\n", sep = "")
  cat(OLIGOMET_URL, "\n\n", sep = "")
  cat("DISCLAIMER\n")
  for (p in OLIGOMET_DISCLAIMER) {
    cat(paste(strwrap(p, width = width), collapse = "\n"), "\n\n", sep = "")
  }
  invisible(NULL)
}

# Comment-prefixed disclaimer block for text file headers (MGF, MSP).
.disclaimer_comment_block <- function(prefix = "# ", width = 76) {
  lines <- c(
    paste0("Generated by OligoMetProfiler -- ", OLIGOMET_URL),
    paste0(OLIGOMET_AUTHOR_ROLE, ": ", oligomet_author_line()),
    "",
    strwrap(OLIGOMET_DISCLAIMER_SHORT, width = width))
  paste0(prefix, lines)
}
