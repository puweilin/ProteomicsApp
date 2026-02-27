# =============================================================================
# GeneSetCacheManager: Local gene set cache for offline enrichment analysis
# =============================================================================
# Builds and caches TERM2GENE / TERM2NAME data from:
#   - GO (BP/MF/CC): org.Hs.eg.db + GO.db via GOALL (ancestor-propagated)
#     -> identical coverage to enrichGO()
#   - Reactome: reactome.db (PATHID2EXTID + PATHID2NAME)
#     -> identical coverage to enrichPathway()
#   - WikiPathways: msigdbr (C2:CP:WIKIPATHWAYS)
#     -> enrichWP online API is broken (404), msigdbr is the best available source
#   - KEGG: KEGGREST API
#     -> identical coverage to enrichKEGG()
#
# Enrichment then uses enricher() / GSEA() with cached data -- fully offline
# after the initial build.
# =============================================================================

library(R6)

GeneSetCacheManager <- R6Class("GeneSetCacheManager",
  public = list(
    cache_dir = NULL,
    cache_metadata = NULL,
    CACHE_TTL_DAYS = 30,

    # ------------------------------------------------------------------
    # initialize
    # ------------------------------------------------------------------
    initialize = function(cache_dir = NULL) {
      if (is.null(cache_dir)) {
        cache_dir <- file.path(tools::R_user_dir("ProteomicsApp", "cache"),
                               "geneset_cache")
      }
      self$cache_dir <- cache_dir
      if (!dir.exists(self$cache_dir)) {
        dir.create(self$cache_dir, recursive = TRUE)
      }
      # Load existing metadata if present
      meta_path <- file.path(self$cache_dir, "cache_meta.rds")
      if (file.exists(meta_path)) {
        self$cache_metadata <- readRDS(meta_path)
      }
      message("GeneSetCacheManager initialized. Cache dir: ", self$cache_dir)
    },

    # ------------------------------------------------------------------
    # build_cache
    # ------------------------------------------------------------------
    build_cache = function(dbs = c("GO_BP", "GO_MF", "GO_CC",
                                   "KEGG", "Reactome", "Wiki"),
                           progress_callback = NULL) {

      # Check required packages per database
      go_dbs <- intersect(dbs, c("GO_BP", "GO_MF", "GO_CC"))
      if (length(go_dbs) > 0) {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
          stop("Package 'org.Hs.eg.db' is required for GO cache.")
        if (!requireNamespace("GO.db", quietly = TRUE))
          stop("Package 'GO.db' is required for GO cache.")
      }
      if ("Reactome" %in% dbs) {
        if (!requireNamespace("reactome.db", quietly = TRUE))
          stop("Package 'reactome.db' is required for Reactome cache.")
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
          stop("Package 'org.Hs.eg.db' is required for Reactome Entrez->Symbol conversion.")
      }
      if ("Wiki" %in% dbs) {
        if (!requireNamespace("msigdbr", quietly = TRUE))
          stop("Package 'msigdbr' is required for WikiPathways cache.")
      }
      if ("KEGG" %in% dbs) {
        if (!requireNamespace("KEGGREST", quietly = TRUE))
          stop("Package 'KEGGREST' is required for KEGG cache.")
      }
      if (!requireNamespace("clusterProfiler", quietly = TRUE))
        stop("Package 'clusterProfiler' is required.")

      db_counts <- list()

      for (i in seq_along(dbs)) {
        db <- dbs[i]
        message("--- Building cache for: ", db, " (", i, "/", length(dbs), ") ---")
        if (!is.null(progress_callback)) {
          progress_callback(db, i, length(dbs))
        }

        result <- tryCatch(
          private$build_single_db(db),
          error = function(e) {
            warning("Failed to build cache for ", db, ": ", e$message)
            NULL
          }
        )

        if (!is.null(result)) {
          saveRDS(result, file.path(self$cache_dir, paste0(db, ".rds")))
          db_counts[[db]] <- nrow(result$TERM2GENE)
          message("  Saved ", db, ".rds  (", nrow(result$TERM2GENE),
                  " gene-term pairs, ",
                  length(unique(result$TERM2GENE$term)), " terms)")
        }
      }

      # Save metadata with version info for all data sources
      self$cache_metadata <- list(
        date_created    = Sys.time(),
        org_hs_version  = as.character(packageVersion("org.Hs.eg.db")),
        go_db_version   = as.character(packageVersion("GO.db")),
        reactome_db_version = tryCatch(
          as.character(packageVersion("reactome.db")),
          error = function(e) NA_character_
        ),
        msigdbr_version = tryCatch(
          as.character(packageVersion("msigdbr")),
          error = function(e) NA_character_
        ),
        kegg_date       = Sys.Date(),
        db_counts       = db_counts
      )
      saveRDS(self$cache_metadata,
              file.path(self$cache_dir, "cache_meta.rds"))
      message("--- Cache build complete ---")
      invisible(self$cache_metadata)
    },

    # ------------------------------------------------------------------
    # get_term2gene
    # ------------------------------------------------------------------
    get_term2gene = function(db) {
      rds_path <- file.path(self$cache_dir, paste0(db, ".rds"))
      if (!file.exists(rds_path)) {
        stop("Cache file not found for '", db,
             "'. Run build_cache() first.")
      }
      readRDS(rds_path)
    },

    # ------------------------------------------------------------------
    # is_cache_exists
    # ------------------------------------------------------------------
    is_cache_exists = function() {
      required <- c("GO_BP", "GO_MF", "GO_CC", "KEGG",
                     "Reactome", "Wiki", "cache_meta")
      all(file.exists(
        file.path(self$cache_dir, paste0(required, ".rds"))
      ))
    },

    # ------------------------------------------------------------------
    # is_cache_expired
    # ------------------------------------------------------------------
    is_cache_expired = function() {
      if (is.null(self$cache_metadata)) return(TRUE)
      age <- as.numeric(difftime(Sys.time(),
                                 self$cache_metadata$date_created,
                                 units = "days"))
      age > self$CACHE_TTL_DAYS
    },

    # ------------------------------------------------------------------
    # check_for_updates
    # ------------------------------------------------------------------
    check_for_updates = function() {
      if (is.null(self$cache_metadata)) {
        return(list(
          has_update     = FALSE,
          is_expired     = TRUE,
          cache_age_days = NA,
          needs_refresh  = TRUE
        ))
      }

      age <- as.numeric(difftime(Sys.time(),
                                 self$cache_metadata$date_created,
                                 units = "days"))

      # Check if any annotation DB version has changed
      has_update <- FALSE
      version_changes <- list()

      check_pkg <- function(pkg_name, cached_field) {
        current <- tryCatch(as.character(packageVersion(pkg_name)),
                            error = function(e) NA_character_)
        cached  <- self$cache_metadata[[cached_field]]
        if (is.null(cached)) cached <- NA_character_
        if (!is.na(current) && !is.na(cached) && current != cached) {
          version_changes[[pkg_name]] <<- list(cached = cached, current = current)
          has_update <<- TRUE
        }
      }

      check_pkg("org.Hs.eg.db",  "org_hs_version")
      check_pkg("GO.db",          "go_db_version")
      check_pkg("reactome.db",    "reactome_db_version")
      check_pkg("msigdbr",        "msigdbr_version")

      is_expired <- age > self$CACHE_TTL_DAYS

      list(
        has_update      = has_update,
        is_expired      = is_expired,
        cache_age_days  = round(age, 1),
        version_changes = version_changes,
        needs_refresh   = has_update || is_expired
      )
    },

    # ------------------------------------------------------------------
    # get_cache_info
    # ------------------------------------------------------------------
    get_cache_info = function() {
      if (is.null(self$cache_metadata)) {
        return("Gene set cache not built yet.")
      }

      age <- round(as.numeric(difftime(Sys.time(),
                                       self$cache_metadata$date_created,
                                       units = "days")), 1)
      update <- self$check_for_updates()

      status <- if (update$needs_refresh) {
        if (update$has_update) "Update available"
        else "Expired"
      } else {
        "Up to date"
      }

      counts_str <- paste(
        vapply(names(self$cache_metadata$db_counts), function(db) {
          paste0(db, ": ", self$cache_metadata$db_counts[[db]])
        }, character(1)),
        collapse = ", "
      )

      # Helper to safely retrieve metadata field
      get_meta <- function(field) {
        val <- self$cache_metadata[[field]]
        if (is.null(val) || is.na(val)) "N/A" else val
      }

      paste0(
        "Cache built: ",
        format(self$cache_metadata$date_created, "%Y-%m-%d %H:%M"), "<br>",
        "Age: ", age, " days<br>",
        "Sources: org.Hs.eg.db ", get_meta("org_hs_version"),
        ", GO.db ", get_meta("go_db_version"),
        ", reactome.db ", get_meta("reactome_db_version"),
        ", msigdbr ", get_meta("msigdbr_version"), "<br>",
        "Gene-term pairs: ", counts_str, "<br>",
        "Status: <b>", status, "</b>"
      )
    }
  ),

  # ====================================================================
  # Private helpers
  # ====================================================================
  private = list(

    # Convert msigdbr gs_name to lowercase human-readable description.
    # "WP_SELENIUM_MICRONUTRIENT_NETWORK"
    #   -> "selenium micronutrient network"
    clean_msigdbr_name = function(name, prefix) {
      cleaned <- sub(paste0("^", prefix, "_"), "", name)
      cleaned <- tolower(gsub("_", " ", cleaned))
      cleaned
    },

    # Build TERM2GENE + TERM2NAME for a single database
    build_single_db = function(db) {
      TERM2GENE <- NULL
      TERM2NAME <- NULL

      if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
        # =============================================================
        # GO via org.Hs.eg.db + GO.db (GOALL = ancestor-propagated)
        # This replicates exactly what enrichGO() does internally via
        # clusterProfiler:::get_GO_data()
        # =============================================================
        ont_map <- c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC")
        ont <- ont_map[[db]]

        message("  Building GO (", ont, ") from org.Hs.eg.db via GOALL...")

        # 1. Get all GO terms for this ontology from GO.db
        goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
        goterms <- goterms[goterms == ont]
        message("  Total ", ont, " GO terms in GO.db: ", length(goterms))

        # 2. Map each GO term -> gene symbols via GOALL (ancestor-propagated)
        #    This is the same approach as clusterProfiler:::get_GO_data()
        go2gene <- suppressMessages(
          AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys     = names(goterms),
            column   = "SYMBOL",
            keytype  = "GOALL",
            multiVals = "list"
          )
        )

        # 3. Convert to data.frame (stack the list)
        goAnno <- stack(go2gene)
        colnames(goAnno) <- c("gene", "go_id")
        goAnno <- unique(goAnno[!is.na(goAnno$gene), ])

        message("  Mapped ", nrow(goAnno), " gene-GO pairs across ",
                length(unique(goAnno$go_id)), " terms")

        # 4. Get GO term descriptions from GO.db
        go_terms_info <- AnnotationDbi::Term(GO.db::GOTERM)
        goAnno$term_name <- go_terms_info[goAnno$go_id]

        # Use GO ID as term key (matches enrichGO output format)
        TERM2GENE <- data.frame(
          term = goAnno$go_id,
          gene = goAnno$gene,
          stringsAsFactors = FALSE
        )
        TERM2NAME <- unique(data.frame(
          term = goAnno$go_id,
          name = goAnno$term_name,
          stringsAsFactors = FALSE
        ))
        TERM2NAME <- TERM2NAME[!is.na(TERM2NAME$name), ]

      } else if (db == "Reactome") {
        # =============================================================
        # Reactome via reactome.db
        # This replicates exactly what ReactomePA::enrichPathway() does
        # internally via get_Reactome_DATA()
        # =============================================================
        message("  Building Reactome from reactome.db...")

        # 1. Get human gene universe (same as ReactomePA:::getALLEG("human"))
        getALLEG <- getFromNamespace("getALLEG", "ReactomePA")
        alleg <- getALLEG("human")

        # 2. Get pathway -> Entrez ID mappings
        PATHID2EXTID <- as.list(reactome.db::reactomePATHID2EXTID)

        # 3. Get pathway names
        PATHID2NAME_raw <- as.list(reactome.db::reactomePATHID2NAME)
        # reactome.db stores multiple names per pathway; take first
        PI <- names(PATHID2NAME_raw)
        PATHID2NAME_raw <- lapply(PATHID2NAME_raw, function(x) x[1])
        names(PATHID2NAME_raw) <- PI

        # 4. Filter to human pathways (R-HSA-*) that have mapped genes
        EXTID2PATHID <- as.list(reactome.db::reactomeEXTID2PATHID)
        EXTID2PATHID <- EXTID2PATHID[names(EXTID2PATHID) %in% alleg]
        valid_pathways <- unique(unlist(EXTID2PATHID))

        PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% names(PATHID2NAME_raw)]
        PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% valid_pathways]
        PATHID2EXTID <- lapply(PATHID2EXTID, function(x) intersect(x, alleg))

        # 5. Clean pathway names: strip "Homo sapiens: " prefix
        path_names <- unlist(PATHID2NAME_raw[names(PATHID2EXTID)])
        path_names <- sub("^[^:]+:\\s*", "", path_names)

        message("  Human Reactome pathways: ", length(PATHID2EXTID))

        # 6. Build long-format data.frame (pathway_id -> entrez)
        df_list <- lapply(names(PATHID2EXTID), function(pid) {
          data.frame(
            PathwayID = pid,
            ENTREZID  = PATHID2EXTID[[pid]],
            stringsAsFactors = FALSE
          )
        })
        df_all <- do.call(rbind, df_list)

        # 7. Convert Entrez -> Symbol
        gene_map <- clusterProfiler::bitr(
          unique(df_all$ENTREZID),
          fromType  = "ENTREZID",
          toType    = "SYMBOL",
          OrgDb     = "org.Hs.eg.db"
        )
        df_all <- merge(df_all, gene_map, by = "ENTREZID")

        message("  Mapped ", nrow(df_all), " gene-pathway pairs")

        TERM2GENE <- data.frame(
          term = df_all$PathwayID,
          gene = df_all$SYMBOL,
          stringsAsFactors = FALSE
        )
        TERM2NAME <- unique(data.frame(
          term = names(PATHID2EXTID),
          name = path_names[names(PATHID2EXTID)],
          stringsAsFactors = FALSE
        ))

      } else if (db == "Wiki") {
        # =============================================================
        # WikiPathways via msigdbr (C2:CP:WIKIPATHWAYS)
        # enrichWP online API is broken (WikiPathways GMT returns 404),
        # so msigdbr is the best available offline source.
        # Falls back to bundled RDS if msigdbr download (Zenodo) fails.
        # =============================================================
        message("  Building WikiPathways from msigdbr...")
        df <- tryCatch({
          raw <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
          raw[raw$gs_subcat == "CP:WIKIPATHWAYS", c("gs_name", "gene_symbol")]
        }, error = function(e) {
          message("  msigdbr download failed (", e$message,
                  "), using bundled WikiPathways data...")
          # Fallback: bundled RDS shipped with the package
          fallback <- system.file("extdata", "wikipathways_msigdb.rds",
                                  package = "ProteomicsApp")
          if (fallback == "") {
            # Not installed as package — try relative path
            fallback <- file.path("inst", "extdata", "wikipathways_msigdb.rds")
          }
          if (file.exists(fallback)) {
            readRDS(fallback)
          } else {
            stop("Neither msigdbr nor bundled WikiPathways data available.")
          }
        })

        readable <- private$clean_msigdbr_name(df$gs_name, "WP")

        TERM2GENE <- data.frame(
          term = readable,
          gene = df$gene_symbol,
          stringsAsFactors = FALSE
        )
        TERM2NAME <- unique(data.frame(
          term = readable,
          name = readable,
          stringsAsFactors = FALSE
        ))

      } else if (db == "KEGG") {
        # =============================================================
        # KEGG via KEGGREST API
        # =============================================================
        message("  Building KEGG from KEGGREST API...")
        link_data <- KEGGREST::keggLink("pathway", "hsa")
        name_data <- KEGGREST::keggList("pathway", "hsa")

        if (length(link_data) == 0 || length(name_data) == 0) {
          stop("Failed to fetch data from KEGG API.")
        }

        gene_ids <- gsub("hsa:", "", names(link_data))
        path_ids <- gsub("path:", "", unname(link_data))

        # Entrez -> Symbol
        gene_map <- clusterProfiler::bitr(
          unique(gene_ids),
          fromType = "ENTREZID",
          toType   = "SYMBOL",
          OrgDb    = "org.Hs.eg.db"
        )

        df_map <- data.frame(ENTREZID = gene_ids, PathwayID = path_ids,
                             stringsAsFactors = FALSE)
        df_map <- merge(df_map, gene_map, by = "ENTREZID")

        # Clean pathway names
        path_names <- name_data[df_map$PathwayID]
        path_names <- gsub(" - Homo sapiens \\(human\\)", "", path_names)
        path_names <- trimws(path_names)

        TERM2GENE <- data.frame(
          term = df_map$PathwayID,
          gene = df_map$SYMBOL,
          stringsAsFactors = FALSE
        )
        TERM2NAME <- unique(data.frame(
          term = df_map$PathwayID,
          name = path_names,
          stringsAsFactors = FALSE
        ))

      } else {
        stop("Unknown database: ", db)
      }

      # Remove duplicates
      TERM2GENE <- unique(TERM2GENE)
      TERM2NAME <- unique(TERM2NAME)

      list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)
    }
  )
)
