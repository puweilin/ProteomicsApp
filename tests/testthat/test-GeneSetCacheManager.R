# =============================================================================
# Tests: GeneSetCacheManager
# =============================================================================

test_that("GeneSetCacheManager initialization creates cache dir", {
  cache_dir <- file.path(tempdir(), paste0("test_gscm_", as.integer(Sys.time())))
  mgr <- GeneSetCacheManager$new(cache_dir = cache_dir)

  expect_true(dir.exists(mgr$cache_dir))
  expect_equal(mgr$cache_dir, cache_dir)

  # Cleanup
  unlink(cache_dir, recursive = TRUE)
})

test_that("is_cache_exists() returns FALSE for empty cache", {
  cache_dir <- file.path(tempdir(), paste0("test_gscm_empty_", as.integer(Sys.time())))
  mgr <- GeneSetCacheManager$new(cache_dir = cache_dir)

  expect_false(mgr$is_cache_exists())

  # Cleanup
  unlink(cache_dir, recursive = TRUE)
})

test_that("is_cache_exists() returns TRUE for complete mock cache", {
  mgr <- create_test_cache_manager()
  expect_true(mgr$is_cache_exists())
})

test_that("get_cache_info() returns informative string when cache exists", {
  mgr <- create_test_cache_manager()
  info <- mgr$get_cache_info()

  expect_type(info, "character")
  expect_true(nchar(info) > 0)
  expect_true(grepl("Cache built", info))
})

test_that("get_cache_info() returns message when no cache metadata", {
  cache_dir <- file.path(tempdir(), paste0("test_gscm_noinfo_", as.integer(Sys.time())))
  mgr <- GeneSetCacheManager$new(cache_dir = cache_dir)

  info <- mgr$get_cache_info()
  expect_type(info, "character")
  expect_true(grepl("not built", info))

  # Cleanup
  unlink(cache_dir, recursive = TRUE)
})

test_that("clean_msigdbr_name() cleans names correctly", {
  mgr <- create_test_cache_manager()
  clean_fn <- mgr$.__enclos_env__$private$clean_msigdbr_name

  result <- clean_fn("WP_SELENIUM_MICRONUTRIENT_NETWORK", "WP")
  expect_equal(result, "selenium micronutrient network")

  result2 <- clean_fn("GOBP_CELL_CYCLE_PROCESS", "GOBP")
  expect_equal(result2, "cell cycle process")
})

test_that("Custom cache dir is respected", {
  custom_dir <- file.path(tempdir(), "my_custom_cache")
  mgr <- GeneSetCacheManager$new(cache_dir = custom_dir)

  expect_equal(mgr$cache_dir, custom_dir)
  expect_true(dir.exists(custom_dir))

  # Cleanup
  unlink(custom_dir, recursive = TRUE)
})

test_that("CACHE_TTL_DAYS default is 30", {
  cache_dir <- file.path(tempdir(), paste0("test_ttl_", as.integer(Sys.time())))
  mgr <- GeneSetCacheManager$new(cache_dir = cache_dir)

  expect_equal(mgr$CACHE_TTL_DAYS, 30)

  # Cleanup
  unlink(cache_dir, recursive = TRUE)
})

test_that("get_term2gene() returns cached data", {
  mgr <- create_test_cache_manager()

  result <- mgr$get_term2gene("GO_BP")
  expect_type(result, "list")
  expect_true("TERM2GENE" %in% names(result))
  expect_true("TERM2NAME" %in% names(result))
  expect_s3_class(result$TERM2GENE, "data.frame")
})

test_that("get_term2gene() errors for missing database", {
  cache_dir <- file.path(tempdir(), paste0("test_gscm_missing_", as.integer(Sys.time())))
  mgr <- GeneSetCacheManager$new(cache_dir = cache_dir)

  expect_error(mgr$get_term2gene("GO_BP"), "Cache file not found")

  # Cleanup
  unlink(cache_dir, recursive = TRUE)
})

test_that("is_cache_expired() returns TRUE when no metadata", {
  cache_dir <- file.path(tempdir(), paste0("test_expire_", as.integer(Sys.time())))
  mgr <- GeneSetCacheManager$new(cache_dir = cache_dir)

  expect_true(mgr$is_cache_expired())

  # Cleanup
  unlink(cache_dir, recursive = TRUE)
})

test_that("check_for_updates() returns structured list", {
  mgr <- create_test_cache_manager()

  result <- mgr$check_for_updates()
  expect_type(result, "list")
  expect_true(all(c("has_update", "is_expired", "cache_age_days", "needs_refresh") %in%
                  names(result)))
})
