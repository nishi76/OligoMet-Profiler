# inspect_workbook.R -- understand the ION337 workbook charge-envelope layout
library(openxlsx)
.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
REFERENCE_WORKBOOK <- Sys.getenv("ION337_WORKBOOK_PATH", file.path(.pkg_root, "tests", "ION337_Workbook.xlsx"))
if (!file.exists(REFERENCE_WORKBOOK)) {
  stop("Reference workbook not found at ", REFERENCE_WORKBOOK, ". ",
       "Place ION337_Workbook.xlsx in tests/, or set the ",
       "ION337_WORKBOOK_PATH environment variable to its location. ",
       "This file is not distributed with the package.")
}
wb <- readWorkbook(REFERENCE_WORKBOOK, sheet = 1)
cat("Columns:", paste(names(wb), collapse = " | "), "\n\n")

# Print rows 1-35 (parent block + charge states) with all columns
cat("=== Rows 1-35 (parent block) ===\n")
for (i in 1:35) {
  vals <- sapply(wb[i, ], function(x) if (is.na(x)) "" else as.character(x))
  nonzero <- vals[vals != ""]
  if (length(nonzero) > 0) {
    cat(sprintf("row %2d: %s\n", i,
        paste(paste0("[", names(nonzero), "]=", nonzero), collapse = "  ")))
  }
}
