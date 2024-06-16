#' Pyscenic workflow
#'
#' @param counts Single cell gene expression counts matrix, (rows=genes x
#' columns=cells). It's okay to provide a csv or loom file. If provided as a
#' file, the matrix must be (rows=cells x columns=genes), otherwise, you should
#' specify the transpose argument.
#' @param tf_list Transcription factors file (TXT; one TF per line). See
#' <https://resources.aertslab.org/cistarget/tf_lists/>.
#' @param motif2tf Motif annotations file. See
#' <https://resources.aertslab.org/cistarget/motif2tf/>.
#' @param motif_ranks The regulatory feature databases file. Two file
#' formats can be supported: feather or db (legacy). See
#' <https://resources.aertslab.org/cistarget/databases/>.
#' @param method The algorithm for gene regulatory network reconstruction, one
#' of "genie3" or "grnboost2". Default: `grnboost2`.
#' @param mode The mode to be used for computing. One of
#' "custom_multiprocessing", "dask_multiprocessing", "dask_cluster". Default:
#' `custom_multiprocessing`.
#' @param pruning A boolean value indicates wether perform pruning when finding
#' enriched motifs. Default: `TRUE`.
#' @param all_modules A boolean value indicates whether including both positive
#' and negative regulons in the analysis.
#' @param chunk_size The size of the module chunks assigned to a node in the
#' dask graph (default: `100`).
#' @param min_orthologous_identity Minimum orthologous identity to use when
#' annotating enriched motifs (default: `0.0`).
#' @param max_similarity_fdr Maximum FDR in motif similarity to use when
#' annotating enriched motifs (default: `0.001`).
#' @param thresholds The first method to create the TF-modules based on the best
#' targets for each transcription factor (default: `c(0.75, 0.90)`).
#' @param top_n_targets The second method is to select the top targets for a
#' given TF. (default: `50`)
#' @param top_n_regulators  The alternative way to create the TF-modules is to
#' select the best regulators for each gene. (default: `c(5, 10, 50)`).
#' @param min_genes The minimum number of genes in a module (default: `20`).
#' @param mask_dropouts A boolean value indicates whether cell dropouts (cells
#' in which expression of either TF or target gene is 0) are masked when
#' calculating the correlation between a TF-target pair. This affects which
#' target genes are included in the initial modules, and the final pruned
#' regulon (by default only positive regulons are kept (see --all_modules
#' option)). The default value in pySCENIC 0.9.16 and previous versions was to
#' mask dropouts when calculating the correlation; however, all cells are now
#' kept by default, to match the R version.
#' @param rank_threshold The rank threshold used for deriving the target genes
#' of an enriched motif (default: `5000`).
#' @param auc_threshold The threshold used for calculating the AUC of a feature
#' as fraction of ranked genes (default: `0.05`).
#' @param nes_threshold The Normalized Enrichment Score (NES) threshold for
#' finding enriched features (default: `3.0`).
#' @param weights Use weights associated with genes in recovery analysis. Is
#' only relevant when `ctx_ofile` is supplied as json format.
#' @param counts_ofile Output file (must end With `.loom`) of the counts matrix.
#' If `NULL`, a temporary file will be used and removed when function exit. If
#' you want to save this file, just specify this argument.
#' @param grn_ofile Output file of the TF-target genes (CSV).
#' @param regulon_ofile Output file of the enriched motifs and target genes
#' (csv, tsv).
#' @param aucell_ofile Output file/stream, a matrix of AUC values (must end with
#' .loom). the loom file while contain the original expression matrix and the
#' calculated AUC values as extra column attributes.
#' @param transpose Transpose the expression matrix if counts is supplied as a
#' csv (rows=genes x columns=cells).
#' @param gene_atrr The name of the row attribute that specifies the gene
#' symbols in the loom file.
#' @param cell_atrr The name of the column attribute that specifies the
#' identifiers of the cells in the loom file.
#' @param threads The number of workers to use. Only valid if using
#' dask_multiprocessing, custom_multiprocessing or local as mode. (default:
#' `1`).
#' @param seed Seed for the expression matrix ranking step. The default is to
#' use a random seed.
#' @param override A boolean value indicates whether overriding the
#' `counts_ofile` or `grn_ofile` if they exist. Since both process are
#' time-consuming.
#' @param envpath A character to define the `PATH` environment variables.
#' @seealso [pyscenic()][biosys::pyscenic]
#' @references <https://github.com/aertslab/pySCENIC>
#' @export
run <- function(counts, tf_list, motif2tf, motif_ranks,
                # pyscenic grn ------------------------
                method = NULL,
                # pyscenic ctx ------------------------
                mode = NULL,
                pruning = TRUE,
                all_modules = FALSE,
                chunk_size = 100L,
                min_orthologous_identity = 0,
                max_similarity_fdr = 0.001,
                thresholds = c(0.75, 0.90),
                top_n_targets = 50L,
                top_n_regulators = c(5, 10, 50),
                min_genes = 20L,
                mask_dropouts = FALSE,
                # pyscenic aucell ---------------------
                weights = FALSE,
                # motif enrichment arguments ----------
                # For both pyscenic `ctx` and `aucell`
                rank_threshold = 5000L,
                auc_threshold = 0.05,
                nes_threshold = 3.0,
                # output arguments --------------------
                counts_ofile = NULL,
                grn_ofile = "grn_adj.csv",
                regulon_ofile = "regulons.csv",
                aucell_ofile = "aucell.loom",
                # common arguments for loom counts matrix ----
                transpose = FALSE,
                gene_atrr = "GeneID",
                cell_atrr = "CellID",
                # common arguments --------------------
                threads = 1L, seed = NULL,
                override = FALSE, envpath = NULL) {
    method <- match.arg(method, c("grnboost2", "genie3"))
    mode <- match.arg(
        mode,
        c("custom_multiprocessing", "dask_multiprocessing", "dask_cluster")
    )
    assert_bool(pruning)
    assert_bool(all_modules)
    assert_(chunk_size, function(x) is_number(x) && x >= 0, "a number (>= 0)")
    assert_(
        min_orthologous_identity,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    assert_(
        max_similarity_fdr,
        function(x) is_number(x) && x >= 0 && x <= 1, "a number ([0, 1])"
    )
    assert_(
        top_n_targets,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    assert_(
        min_genes,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    assert_bool(mask_dropouts)
    assert_bool(weights)
    assert_(
        rank_threshold,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    assert_(
        auc_threshold,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    assert_(
        nes_threshold,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    assert_string(grn_ofile, empty_ok = FALSE)
    assert_(regulon_ofile, function(x) {
        rlang::is_string(x) && (endsWith(x, ".csv") || endsWith(x, ".tsv"))
    }, "a string ends with `.csv` or `.tsv`", empty_ok = FALSE)
    assert_(aucell_ofile, function(x) {
        rlang::is_string(x) && endsWith(x, ".loom")
    }, "a string ends with `.loom`", empty_ok = FALSE)
    assert_bool(transpose)
    assert_string(gene_atrr, empty_ok = FALSE)
    assert_string(cell_atrr, empty_ok = FALSE)
    assert_(
        threads,
        function(x) is_number(x) && x >= 0, "a number (>= 0)"
    )
    threads <- as.integer(threads)
    if (inherits(counts, what = c("matrix", "Matrix"))) {
        assert_(counts_ofile, function(x) {
            rlang::is_string(x) && endsWith(x, ".loom")
        }, "a string ends with `.loom`", null_ok = TRUE)
        if (is.null(counts_ofile)) {
            counts_mat_file <- tempfile("pyscenic", fileext = ".loom")
            on.exit(file.remove(counts_mat_file))
        } else {
            counts_mat_file <- counts_ofile
        }
        # re-using counts matrix file or not
        if (is.null(counts_ofile) || !file.exists(counts_mat_file) ||
            override) {
            assert_pkg("loomR")
            con <- loomR::create(
                # the output file name, which is also the input of pyscenic
                counts_mat_file,
                # row are genes, column are cells
                # the internal will transpose the matrix
                data = counts,
                gene.attrs = structure(
                    list(rownames(counts)),
                    names = gene_atrr
                ), # -- gene_attribute
                cell.attrs = structure(
                    list(colnames(counts)),
                    names = cell_atrr
                ) # -- cell_id_attribute
            )
            con$close_all()
        } else {
            cli::cli_alert_info(
                "Re-using counts matrix file from: {.path {counts_mat_file}}"
            )
        }
    } else if (rlang::is_string(counts) && counts != "") {
        counts_mat_file <- counts
    } else {
        cli::cli_abort(
            "{.arg counts} must be a string of file path or a matrix"
        )
    }
    # prepare seed -------------------------------------
    seed <- as.integer(seed)
    if (is_scalar(seed)) {
        set_seed(seed, TRUE)
        seed <- random_seed(2L)
    } else if (length(seed) >= 2L) {
        seed <- seed[seq_len(2L)]
    } else if (anyNA(seed)) {
        cli::cli_abort("{.arg seed} cannot be `NA`")
    } else if (length(seed) == 0L) {
        set_seed(add = TRUE)
        seed <- random_seed(2L)
    }

    # GRN inference using the GRNBoost2 algorithm -----------
    # https://github.com/aertslab/pySCENIC/issues/525#issuecomment-2041298258
    # pip install dask-expr==0.5.3
    if (!file.exists(grn_ofile) || override) {
        biosys::pyscenic(
            "grn",
            "--seed", seed[1L],
            "--num_workers", threads,
            if (transpose) "--transpose",
            "-m", method,
            # output file ----------------------------------------
            "--output", grn_ofile,
            # expression matrix file -----------------------------
            "--gene_attribute", gene_atrr,
            "--cell_id_attribute", cell_atrr,
            "--sparse",
            counts_mat_file,
            # the list of transcription factors ------------------
            tf_list
        )$setup_envpath(envpath)$run()
    } else {
        cli::cli_alert_info(
            "Re-using pyscenic grn output from: {.path {grn_ofile}}"
        )
    }

    # Regulon prediction aka cisTarget from CLI ------------
    biosys::pyscenic(
        "ctx",
        # output file ----------------------------------------
        "--output", regulon_ofile,
        "--num_workers", threads,
        if (!pruning) "--no_pruning",
        "--chunk_size", chunk_size,
        "--mode", mode,
        if (all_modules) "--all_modules",
        if (transpose) "--transpose",
        # motif enrichment arguments -----------------------
        "--rank_threshold", rank_threshold,
        "--auc_threshold", auc_threshold,
        "--nes_threshold", nes_threshold,
        # motif annotation arguments -------------------------
        "--min_orthologous_identity", min_orthologous_identity,
        "--max_similarity_fdr", max_similarity_fdr,
        # motif2tf file --------------------------------------
        "--annotations_fname", motif2tf,
        # module generation arguments ------------------------
        "--thresholds", thresholds,
        "--top_n_targets", top_n_targets,
        "--top_n_regulators", top_n_regulators,
        "--min_genes", min_genes,
        if (mask_dropouts) "--mask_dropouts",
        # expression matrix file -----------------------------
        "--gene_attribute", gene_atrr,
        "--cell_id_attribute", cell_atrr,
        "--sparse",
        "--expression_mtx_fname", counts_mat_file,
        # output from pyscenic grn ---------------------------
        grn_ofile,
        # followed by motif ranking databases ----------------
        # it could be multiple files
        motif_ranks
    )$setup_envpath(envpath)$run()

    # Cellular enrichment ------------------------------------
    biosys::pyscenic(
        "aucell",
        "--seed", seed[2L],
        "--num_workers", threads,
        if (transpose) "--transpose",
        # Use weights associated with genes in recovery analysis. Is only
        # relevant when gene signatures are supplied as json format.
        if (weights) "--weights",
        # output file ----------------------------------------
        "--output", aucell_ofile,

        # motif enrichment arguments -----------------------
        "--rank_threshold", rank_threshold,
        "--auc_threshold", auc_threshold,
        "--nes_threshold", nes_threshold,

        # expression matrix file -----------------------------
        "--gene_attribute", gene_atrr,
        "--cell_id_attribute", cell_atrr,
        "--sparse",
        counts_mat_file,
        # output from pyscenic ctx ---------------------------
        regulon_ofile
    )$setup_envpath(envpath)$run()
}
