# =============================================================================
# Tests: Input validation (validate_excel_input from mod_preprocessing.R)
# =============================================================================

# The validate_excel_input function is defined inside preprocessingServer.
# We re-implement the same logic here for testing since it's a local function.
# This tests the validation rules themselves.

validate_excel_input_standalone <- function(file_path) {
  errors <- character(0)
  if (!file.exists(file_path)) return("File not found.")
  sheets <- readxl::excel_sheets(file_path)
  if (length(sheets) < 3) {
    return(paste0("Excel file must have at least 3 sheets (found ", length(sheets), ")."))
  }
  mat <- tryCatch(readxl::read_excel(file_path, sheet = 1), error = function(e) NULL)
  meta <- tryCatch(readxl::read_excel(file_path, sheet = 2), error = function(e) NULL)
  annot <- tryCatch(readxl::read_excel(file_path, sheet = 3), error = function(e) NULL)
  if (is.null(mat)) errors <- c(errors, "Cannot read Sheet 1 (Matrix).")
  if (is.null(meta)) errors <- c(errors, "Cannot read Sheet 2 (Metadata).")
  if (is.null(annot)) errors <- c(errors, "Cannot read Sheet 3 (Annotation).")
  if (length(errors) > 0) return(paste(errors, collapse = "\n"))
  if (ncol(mat) < 2) errors <- c(errors, "Matrix sheet must have at least 2 columns.")
  if (!"label" %in% colnames(meta)) errors <- c(errors, "Metadata must contain 'label' column.")
  if (ncol(annot) < 2) errors <- c(errors, "Annotation must have at least 2 columns.")
  if (length(errors) > 0) return(paste(errors, collapse = "\n"))
  return(NULL)
}

test_that("Valid Excel file returns NULL (no errors)", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  path <- create_test_excel()
  result <- validate_excel_input_standalone(path)
  expect_null(result)
})

test_that("Missing file returns error", {
  skip_if_not_installed("readxl")

  result <- validate_excel_input_standalone("/nonexistent/path/file.xlsx")
  expect_type(result, "character")
  expect_true(grepl("not found", result))
})

test_that("Too few sheets errors", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  # Create Excel with only 2 sheets
  path <- file.path(tempdir(), "two_sheets.xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Sheet1")
  openxlsx::writeData(wb, "Sheet1", data.frame(a = 1:5))
  openxlsx::addWorksheet(wb, "Sheet2")
  openxlsx::writeData(wb, "Sheet2", data.frame(b = 1:5))
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  result <- validate_excel_input_standalone(path)
  expect_type(result, "character")
  expect_true(grepl("3 sheets", result))
})

test_that("Missing 'label' column errors", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  mat   <- create_test_matrix()
  annot <- create_test_annotation()
  # Metadata without 'label' column
  meta_bad <- data.frame(
    sample_name = paste0("S", 1:8),
    condition   = rep(c("A", "B"), each = 4),
    stringsAsFactors = FALSE
  )

  path <- file.path(tempdir(), "no_label.xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Matrix")
  openxlsx::writeData(wb, "Matrix", mat)
  openxlsx::addWorksheet(wb, "Metadata")
  openxlsx::writeData(wb, "Metadata", meta_bad)
  openxlsx::addWorksheet(wb, "Annotation")
  openxlsx::writeData(wb, "Annotation", annot)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  result <- validate_excel_input_standalone(path)
  expect_type(result, "character")
  expect_true(grepl("label", result))
})

test_that("Matrix with single column errors", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  meta  <- create_test_metadata()
  annot <- create_test_annotation()

  # Matrix with only 1 column (no sample data)
  mat_bad <- data.frame(ProteinID = paste0("P", 1:5))

  path <- file.path(tempdir(), "single_col_matrix.xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Matrix")
  openxlsx::writeData(wb, "Matrix", mat_bad)
  openxlsx::addWorksheet(wb, "Metadata")
  openxlsx::writeData(wb, "Metadata", meta)
  openxlsx::addWorksheet(wb, "Annotation")
  openxlsx::writeData(wb, "Annotation", annot)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  result <- validate_excel_input_standalone(path)
  expect_type(result, "character")
  expect_true(grepl("at least 2 columns", result))
})

# =============================================================================
# Tests: Path sanitization for special characters
# =============================================================================

test_that("File path with special characters can be sanitized for cache_root", {
  # Simulate the sanitization logic from mod_preprocessing.R
  file_paths <- c(
    "/data/experiment & results/patient(1)/data.xlsx",
    "/data/my+project [2024]/test.xlsx",
    "/data/résumé data/café.xlsx",
    "/data/normal_path/file.xlsx"
  )

  for (fp in file_paths) {
    file_id <- tools::file_path_sans_ext(basename(fp))
    safe_id <- gsub("[^a-zA-Z0-9._-]", "_", file_id)

    # Should not contain problematic characters
    expect_false(grepl("[&+()\\[\\]]", safe_id),
                 info = paste("Failed for:", fp))
    # Should not be empty
    expect_true(nchar(safe_id) > 0, info = paste("Empty for:", fp))
  }
})

test_that("Group names with special characters are sanitized for directory names", {
  # Simulate get_base_name() sanitization from mod_analysis.R
  test_cases <- list(
    list(case = "Treatment (high dose)", control = "Control & Placebo",
         expect_no = "[&()]"),
    list(case = "Group+A", control = "Group-B",
         expect_no = "[+]"),
    list(case = "Type [1]", control = "Type [2]",
         expect_no = "[\\[\\]]"),
    list(case = "50% Response", control = "Normal",
         expect_no = "[%]")
  )

  for (tc in test_cases) {
    raw_name <- paste0(tc$case, "_vs_", tc$control)
    safe_name <- gsub("[^a-zA-Z0-9._-]", "_", raw_name)

    expect_false(grepl(tc$expect_no, safe_name),
                 info = paste("Failed for:", raw_name))
    expect_true(nchar(safe_name) > 0)
  }
})

test_that("Excel file in directory with special characters can be read", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  # Create a directory with special characters
  special_dir <- file.path(tempdir(), "test&data (project) [v2]")
  dir.create(special_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(special_dir)) {
    skip("OS does not support special characters in directory names")
  }

  path <- create_test_excel(dir = special_dir)
  result <- validate_excel_input_standalone(path)
  expect_null(result)

  # Cleanup
  unlink(special_dir, recursive = TRUE)
})

test_that("Annotation with single column errors", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  mat  <- create_test_matrix()
  meta <- create_test_metadata()

  # Annotation with only 1 column
  annot_bad <- data.frame(ProteinID = paste0("P", 1:5))

  path <- file.path(tempdir(), "single_col_annot.xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Matrix")
  openxlsx::writeData(wb, "Matrix", mat)
  openxlsx::addWorksheet(wb, "Metadata")
  openxlsx::writeData(wb, "Metadata", meta)
  openxlsx::addWorksheet(wb, "Annotation")
  openxlsx::writeData(wb, "Annotation", annot_bad)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  result <- validate_excel_input_standalone(path)
  expect_type(result, "character")
  expect_true(grepl("at least 2 columns", result))
})
