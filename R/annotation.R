# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


# source("R/granges_conversion.R")


#' Annotate with consequences
#' @description Request variant consequences from Variant Effect Predictor (VEP)
#' via Ensembl Rest Service. Not recommended for large queries.
#' @param geno Data frame or GenomicRanges::GRanges object including columns
#' rsid, ref, alt.
#' @param species Species name, e.g. mouse (GRCm38) or human (GRCh38).
#' @return Data frame.
#' @examples
#' geno = finemap("chr1",
#'   start = 5000000, end = 6000000,
#'   strain1 = c("C57BL_6J"), strain2 = c("AKR_J", "A_J", "BALB_cJ")
#' )
#'
#' df = annotate_consequences(geno[seq_len(10), ], "mouse")
#'
#' geno.granges = finemap("chr1",
#'     start = 5000000, end = 6000000,
#'     strain1 = c("C57BL_6J"), strain2 = c("AKR_J", "A_J", "BALB_cJ"),
#'     return_obj = "granges"
#' )
#'
#' df2 = annotate_consequences(geno.granges[seq_len(10), ], "mouse")
#' @export
#' @importFrom methods is
#' @importFrom data.table rbindlist
annotate_consequences = function(geno, species) {
  if ("GRanges" %in% is(geno)) {
    geno = GRanges2df(geno)
  }

  if (length(is.na(geno$rsid)) > 0) {
    message("Please consider: only SNVs with existing rsID can be annotated")
  }
  df.split = df_split(geno[!is.na(geno$rsid), ], 200)

  res = lapply(df.split, function(x) {
    message(paste0("Query ", nrow(x), " variants..."))
    ensembl_rest_vep(x, "mouse")
  })

  return(as.data.frame(rbindlist(res)))
}


#' Request variant consequences from Variant Effect Predictor (VEP)
#' via Ensembl Rest Service
#' @param geno Data frame including columns rsid, ref, alt.
#' @param species Species name, e.g. mouse or human.
#' @return Data frame.
#' @keywords internal
#' @importFrom curl has_internet
#' @importFrom httr POST content_type accept stop_for_status
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom data.table rbindlist
ensembl_rest_vep = function(geno, species) {
  # Check if there is an internet connection
  if (!has_internet()) {
    stop("No internet connection detected...")
  }


  # Request
  server = "https://rest.ensembl.org"
  ext = paste0("/vep/", species, "/id")
  r = POST(
    paste(server, ext, sep = ""),
    content_type("application/json"),
    accept("application/json"),
    body = paste0('{ "ids" : ["', paste0(geno$rsid,
      collapse =
        '","'
    ), '" ] }')
  )


  stop_for_status(r)


  res = fromJSON(toJSON(content(r)))


  # Extract consequences
  cons.list = lapply(
    seq_len(length(res$input)),
    function(i) {
      if (!is.null(res$transcript_consequences[[i]])) {
        res$transcript_consequences[[i]]$snp = res$input[[i]]
        res$transcript_consequences[[i]]$consequence_terms =
          lapply(
            seq_len(
              length(res$transcript_consequences[[i]]$consequence_terms)
            ),
            function(j) {
              paste0(res$transcript_consequences[[i]]$consequence_terms[[j]],
                collapse = ","
              )
            }
          )
        return(res$transcript_consequences[[i]])
      }
      else if (!is.null(res$intergenic_consequences[[i]])) {
        res$intergenic_consequences[[i]]$snp = res$input[[i]]
        res$intergenic_consequences[[i]]$consequence_terms =
          lapply(
            seq_len(
              length(res$intergenic_consequences[[i]]$consequence_terms)
            ),
            function(j) {
              paste0(res$intergenic_consequences[[i]]$consequence_terms[[j]],
                collapse = ","
              )
            }
          )
        return(res$intergenic_consequences[[i]])
      }
    }
  )


  final = rbindlist(cons.list, fill = TRUE)


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
  if (length(cols.missing) > 0) {
    final[, cols.missing] = NA
  }
  final = subset(final, select = cols)

  final = as.data.frame(matrix(unlist(final), ncol = ncol(final)))
  colnames(final) = cols


  # Filter for alleles
  final = final[final$snp %in% geno$rsid &
    (final$variant_allele %in% geno$ref |
      final$variant_allele %in% geno$alt), ]


  return(final)
}


#' Splits data frame df into subsets with maximum n rows
#' @param df Data frame.
#' @param n Max number of rows per subset.
#' @return List of data frames.
#' @keywords internal
df_split = function(df, n) {
  l = nrow(df)
  s = ceiling(l / n)

  return(lapply(seq_len(s), function(i) {
    to = min(l, i * 200)
    from = i * 200 - n + 1
    return(df[from:to, ])
  }))
}


#' Annotate with genes
#' @description Request mouse genes from Ensembl Biomart.
#' @param geno Data frame or GenomicRanges::GRanges object including columns
#' chr, pos.
#' @param flanking Size of flanking sequence to be included.
#' @return Data frame.
#' @examples
#' geno = finemap("chr1",
#'   start = 5000000, end = 6000000,
#'   strain1 = c("C57BL_6J"), strain2 = c("AKR_J", "A_J", "BALB_cJ")
#' )
#'
#' genes = annotate_mouse_genes(geno, 50000)
#' @export
#' @importFrom curl has_internet
#' @importFrom methods is
#' @importFrom GenomicRanges start end start<- end<- intersect
#' @importFrom biomaRt useMart listDatasets useDataset getBM
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges subsetByOverlaps
annotate_mouse_genes = function(geno, flanking = NULL) {

    # Check if there is an internet connection
    if (!has_internet()) {
        stop("No internet connection detected...")
    }

    if (!("GRanges" %in% is(geno))) {
        geno = df2GRanges(geno)
    }

    if (is.numeric(flanking)) {
        start(geno) = start(geno) - flanking
        end(geno) = end(geno) + flanking
    }

    # biomaRt::listMarts()
    # Archives: https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#using-archived-versions-of-ensembl
    m = useMart("ENSEMBL_MART_ENSEMBL", host = "https://nov2020.archive.ensembl.org")
    datasets = listDatasets(m)
    # head(datasets[grep ("musculus", datasets$dataset),])

    if (!startsWith(datasets$version[datasets$dataset ==
        "mmusculus_gene_ensembl"], ref_genome())) {
        stop(
            "Reference genome version of BiomaRt is different from the one used 
            in the R package. Contact maintainer."
        )
    }

    ds = useDataset("mmusculus_gene_ensembl", mart = m)

    # filters = listFilters(ds)
    # attributes = listAttributes(ds)


    # Request Biomart
    res = getBM(
        attributes = c(
            "external_gene_name",
            "external_gene_source",
            "ensembl_gene_id",
            "description",
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
            "gene_biotype"
        ),
        filters = "chromosome_name",
        values = seqlevels(geno),
        mart = ds
    )


    # Convert Biomart output to GRanges object
    res$strand = vapply(res$strand, function(x) {
        if (x == 1) {
            "+"
        } else {
            "-"
        }
    }, character(1))
    res.granges = df2GRanges(
        res,
        chr_name = "chromosome_name",
        start_name = "start_position",
        end_name = "end_position"
    )


    # Overlap
    inters = intersect(res.granges,
        geno,
        ignore.strand = TRUE
    )

    geno.subset = subsetByOverlaps(geno, inters)
    res.subset = subsetByOverlaps(res.granges, inters)


    # Reformat
    res.subset = as.data.frame(res.subset)
    colnames(res.subset) = c(
        "chr",
        "start",
        "end",
        "width",
        "strand",
        "symbol",
        "symbol_source",
        "ensgid",
        "description",
        "biotype"
    )
    res.subset = res.subset[, c(
        "chr",
        "start",
        "end",
        "symbol",
        "symbol_source",
        "ensgid",
        "description",
        "biotype",
        "strand"
    )]

    return(res.subset)
}
