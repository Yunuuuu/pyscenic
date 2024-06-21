#' Import pyscenic auc into R
#'
#' @param auc_scenic A loom file return by pyscenic aucell.
#' @return
#' A list of following elements:
#'  - `auc`: AUC value matrix.
#'  - `regulons`: A list of regulon signatures attached with attributes binary
#'              AUC value matrix ("binary_incidence") and thretholds
#'              ("threshold").
#' @export
import_auc <- function(auc_scenic) {
    assert_(auc_scenic, function(x) {
        rlang::is_string(x) && endsWith(x, ".loom")
    }, "a string ends with `.loom`", empty_ok = FALSE)

    # https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/importing_pySCENIC.html
    loom <- loomR::connect(auc_scenic)
    gene_ids <- loom$row.attrs$gene_names[]
    cell_ids <- loom$col.attrs$cell_names[]

    # get_regulons_AUC --------------------
    regulons_auc <- loom$col.attrs$RegulonsAUC[]
    rownames(regulons_auc) <- cell_ids
    regulons_auc <- t(regulons_auc)

    # get_regulons ------------------------
    regulons_incidence <- loom$row.attrs$Regulons[] # incid mat
    rownames(regulons_incidence) <- gene_ids

    # regulonsToGeneLists -----------------
    regulons <- apply(regulons_incidence, 2L,
        function(x) names(x)[x > 0L],
        simplify = FALSE
    )

    # get_dgem ------------------------
    # counts <- loom$matrix[, ]
    # dimnames(counts) <- list(cell_ids, gene_ids)

    # get_regulonThresholds ------------------------
    md <- load_meta_data(hdf5r::h5attr(loom, "MetaData"))
    regulonsAucThresholds <- vapply(
        md$regulonThresholds, .subset2, numeric(1L), "defaultThresholdValue"
    )
    names(regulonsAucThresholds) <- vapply(
        md$regulonThresholds, .subset2, character(1L), "regulon"
    )
    attr(regulons, "binary_incidence") <- regulons_incidence
    attr(regulons, "threshold") <- regulonsAucThresholds
    list(auc = regulons_auc, regulons = regulons)
}

load_meta_data <- function(meta.data) {
    if (is_base64_encoded(value = meta.data)) {
        rjson::fromJSON(decompress_gzb64(gzb64c = meta.data))
    } else {
        tryCatch(
            rjson::fromJSON(json_str = meta.data),
            function(cnd) {
                cli::cli_abort(
                    "The global MetaData attribute in the given loom is corrupted."
                )
            }
        )
    }
}

is_base64_encoded <- function(value) {
    grepl("^([A-Za-z0-9+/]{4})*([A-Za-z0-9+/]{3}=|[A-Za-z0-9+/]{2}==)?$",
        x = value, perl = TRUE
    )
}

decompress_gzb64 <- function(gzb64c) {
    rawToChar(
        memDecompress(
            from = base64enc::base64decode(what = gzb64c),
            type = "gzip", asChar = FALSE
        ),
        multiple = FALSE
    )
}
