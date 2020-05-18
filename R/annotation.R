# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck(/path/to/project)


source("R/granges_conversion.R")


#'Annotate consequences
#'@description Request variant consequences from Variant Effect Predictor
#'(VEP) via Ensembl Rest Service.
#'@param geno Data frame or GenomicRanges::GRanges object including columns rsid, ref, alt.
#'@param species Species name, e.g. mouse (GRCm38) or human (GRCh38).
#'@return Data frame.
#'@examples geno = finemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("AKR_J", "A_J", "BALB_cJ"))
#'
#'df = annotate_consequences(geno[seq_len(10),], "mouse")
#'
#'geno.granges = finemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("AKR_J", "A_J", "BALB_cJ"), return_obj="granges")
#'
#'df2 = annotate_consequences(geno.granges[seq_len(10),], "mouse")
#'@export
annotate_consequences = function(geno, species) {
    if ("GRanges" %in% methods::is(geno)) {
        geno = GRanges2df(geno)
    }

    df.split = df_split(geno[!is.na(geno$rsid),], 200)

    res = lapply(df.split, function(x) {
        message(paste0("Query ", nrow(x), " variants..."))
        ensembl_rest_vep(x, "mouse")
    })

    return(as.data.frame(data.table::rbindlist(res)))
}


#'Request variant consequences from Variant Effect Predictor (VEP)
#'via Ensembl Rest Service
#'@param geno Data frame including columns rsid, ref, alt.
#'@param species Species name, e.g. mouse or human.
#'@return Data frame.
#'@keywords internal
ensembl_rest_vep = function(geno, species) {
    # Check if there is an internet connection
    if (!curl::has_internet())
        stop("No internet connection detected...")


    # Request
    server = "https://rest.ensembl.org"
    ext = paste0("/vep/", species, "/id")
    r = httr::POST(
        paste(server, ext, sep = ""),
        httr::content_type("application/json"),
        httr::accept("application/json"),
        body = paste0('{ "ids" : ["', paste0(geno$rsid, collapse =
                                                 '","'), '" ] }')
    )


    httr::stop_for_status(r)


    res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))


    # Extract consequences
    cons.list = lapply(seq_len(length(res$input)),
                       function(i) {
                           if (!is.null(res$transcript_consequences[[i]])) {
                               res$transcript_consequences[[i]]$snp = res$input[[i]]
                               res$transcript_consequences[[i]]$consequence_terms =
                                   lapply(seq_len(
                                       length(res$transcript_consequences[[i]]$consequence_terms)
                                   ),
                                   function(j) {
                                       paste0(res$transcript_consequences[[i]]$consequence_terms[[j]],
                                              collapse = ",")
                                   })
                               return(res$transcript_consequences[[i]])
                           }
                           else if (!is.null(res$intergenic_consequences[[i]])) {
                               res$intergenic_consequences[[i]]$snp = res$input[[i]]
                               res$intergenic_consequences[[i]]$consequence_terms =
                                   lapply(seq_len(
                                       length(res$intergenic_consequences[[i]]$consequence_terms)
                                   ),
                                   function(j) {
                                       paste0(res$intergenic_consequences[[i]]$consequence_terms[[j]],
                                              collapse = ",")
                                   })
                               return(res$intergenic_consequences[[i]])
                           }
                       })


    final = data.table::rbindlist(cons.list, fill = TRUE)


    # Reformat
    final[final == "NULL"] = NA
    cols = c(
        "snp",
        "variant_allele",
        "consequence_terms",
        "impact",
        "transcript_id",
        "gene_id",
        "gene_symbol",
        "gene_symbol_source",
        "strand",
        "distance",
        "cdna_end",
        "cdna_start",
        "codons",
        "cds_end",
        "amino_acids",
        "protein_start",
        "cds_start",
        "protein_end",
        "sift_score",
        "sift_prediction"
    )
    cols.missing = cols[!cols %in% colnames(final)]
    if (length(cols.missing) > 0)
        final[, cols.missing] = NA
    final = subset(final, select = cols)

    final = as.data.frame(matrix(unlist(final), ncol = ncol(final)))
    colnames(final) = cols


    # Filter for alleles
    final = final[final$snp %in% geno$rsid &
                      (final$variant_allele %in% geno$ref |
                           final$variant_allele %in% geno$alt),]


    return(final)
}


#'Splits data frame df into subsets with maximum n rows
#'@param df Data frame.
#'@param n Max number of rows per subset.
#'@return List of data frames.
#'@keywords internal
df_split = function(df, n) {
    l = nrow(df)
    s = ceiling(l / n)

    return(lapply(seq_len(s), function(i) {
        to = min(l, i * 200)
        from = i * 200 - n + 1
        return(df[from:to,])
    }))
}
