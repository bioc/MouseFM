# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


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
#' @description Request mouse gene annotations from Ensembl BioMart.
#' @param geno Data frame or \code{GenomicRanges::GRanges} object including
#' chromosome and genomic position information.
#' @param flanking Numeric. Size of the flanking region to include on both
#' sides of each interval.
#' @param host Character string. Ensembl host passed to
#' \code{biomaRt::useMart()}. If the selected host does not match the
#' reference genome used by MouseFM, another Ensembl archive host should be
#' selected, for example from \code{biomaRt::listEnsemblArchives()}.
#' @param biomart Character string. Name of the BioMart database passed to
#' \code{biomaRt::useMart()}.
#' @param dataset_name Character string. Name of the BioMart dataset to use.
#' Default is \code{"mmusculus_gene_ensembl"}.
#' @return Data frame with overlapping gene annotations.
#' @examplesIf interactive() && curl::has_internet()
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
#' @importFrom Seqinfo seqlevels
#' @importFrom IRanges subsetByOverlaps
annotate_mouse_genes = function(geno,
                                flanking = NULL,
                                host = "https://nov2020.archive.ensembl.org",
                                biomart = "ENSEMBL_MART_ENSEMBL",
                                dataset_name = "mmusculus_gene_ensembl") {
    
    
    # Check if there is an internet connection
    if (!has_internet()) {
        stop("No internet connection detected...", call. = FALSE)
    }
    
    if (!("GRanges" %in% is(geno))) {
        geno = df2GRanges(geno)
    }
    
    if (is.numeric(flanking)) {
        start(geno) = start(geno) - flanking
        end(geno) = end(geno) + flanking
    }
    
    m = tryCatch(
        useMart(biomart, host = host),
        error = function(e) {
            stop(
                sprintf(
                    paste(
                        "Could not connect to BioMart '%s' at host '%s'.",
                        "The selected Ensembl archive may be unavailable or slow.",
                        "Please try another host listed by biomaRt::listEnsemblArchives().",
                        "Original error: %s"
                    ),
                    biomart, host, conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )
    
    datasets = tryCatch(
        listDatasets(m),
        error = function(e) {
            stop(
                sprintf(
                    paste(
                        "Connected to BioMart host '%s', but could not retrieve datasets.",
                        "Please try another host listed by biomaRt::listEnsemblArchives().",
                        "Original error: %s"
                    ),
                    host, conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )
    
    dataset_version = datasets$version[datasets$dataset == dataset_name][1]
    expected_version = ref_genome()
    
    if (is.na(dataset_version) || !nzchar(dataset_version)) {
        stop(
            sprintf(
                paste(
                    "Could not determine the Ensembl dataset version for '%s' from host '%s'.",
                    "Please check the selected host or choose another one with",
                    "biomaRt::listEnsemblArchives()."
                ),
                dataset_name, host
            ),
            call. = FALSE
        )
    }
    
    if (!startsWith(dataset_version, expected_version)) {
        stop(
            sprintf(
                paste(
                    "Reference genome mismatch between MouseFM and the selected Ensembl host.",
                    "MouseFM expects a genome version starting with '%s',",
                    "but dataset '%s' on host '%s' reports version '%s'.",
                    "Please pass another 'host' to annotate_mouse_genes().",
                    "Available Ensembl archive hosts can be checked with",
                    "biomaRt::listEnsemblArchives()."
                ),
                expected_version, dataset_name, host, dataset_version
            ),
            call. = FALSE
        )
    }
    
    ds = useDataset(dataset_name, mart = m)
    
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
    
    inters = intersect(res.granges, geno, ignore.strand = TRUE)
    
    geno.subset = subsetByOverlaps(geno, inters)
    res.subset = subsetByOverlaps(res.granges, inters)
    
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
