# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


source("R/send_request.R")
source("R/make_query.R")
source("R/granges_conversion.R")


#' Finemapping of genetic regions
#' @description Finemapping of genetic regions in 37 inbred mice by taking
#' advantage of their very high homozygosity rate (>95%). For one ore more
#' chromosomal regions (GRCm38), this method extracts homozygous SNVs for which
#' the allele differs between two sets of strains (e.g. case vs controls) and
#' outputs respective causal SNV/gene candidates.
#' @param chr Vector of chromosome names.
#' @param start Optional vector of chromosomal start positions of target regions
#' (GRCm38).
#' @param end Optional vector of chromosomal end positions of target regions
#' (GRCm38).
#' @param strain1 First strain set with strains from avail_strains().
#' @param strain2 Second strain set with strains from avail_strains().
#' @param consequence Optional vector of consequence types.
#' @param impact Optional vector of impact types.
#' @param thr1 Number discordant strains in strain1. Between 0 and
#' length(strain1)-1. 0 by default.
#' @param thr2 Number discordant strains in strain2. Between 0 and
#' length(strain2)-1. 0 by default.
#' @param return_obj The user can choose to get the result to be returned as
#' data frame ("dataframe") or as a GenomicRanges::GRanges ("granges") object.
#' Default value is "dataframe".
#' @return Data frame or GenomicRanges::GRanges object containing result data.
#' @examples
#' geno = finemap("chr1",
#'     start = 5000000, end = 6000000,
#'     strain1 = c("C57BL_6J"), strain2 = c(
#'         "129S1_SvImJ", "129S5SvEvBrd",
#'         "AKR_J"
#'     )
#' )
#'
#' comment(geno)
#' @export
#' @importFrom scales comma
#' @importFrom data.table rbindlist
#' @importFrom stats setNames
finemap = function(chr,
                    start = NULL,
                    end = NULL,
                    strain1,
                    strain2,
                    consequence = NULL,
                    impact = NULL,
                    thr1 = 0,
                    thr2 = 0,
                    return_obj = "dataframe") {
    # Create URL and query data
    res = lapply(seq_len(length(chr)), function(i) {
        message(paste0(
            "Query ",
            chr[i],
            if (is.numeric(start[i]) &&
                is.numeric(end[i])) {
                paste0(":", comma(start[i]), "-", comma(end[i]))
            } else {
                ""
            }
        ))
        q = finemap_query(
            chr[i],
            start[i],
            end[i],
            strain1,
            strain2,
            consequence,
            impact,
            thr1,
            thr2
        )
        backend_request(q)
    })


    geno = as.data.frame(rbindlist(res))


    # Return if no results
    if (nrow(geno) == 0) {
        return(geno)
    }

    
    # Keep only input strains
    geno = geno[tolower(names(geno)) %in% c(
        "rsid",
        "chr",
        "pos",
        "ref",
        "alt",
        "consequences",
        tolower(unique(strain1)),
        tolower(unique(strain2))
    )]
    

    # Convert to respective data types
    geno[!names(geno) %in% c("rsid", "ref", "alt", "most_severe_consequence", "consequences")] =
        vapply(
            geno[!names(geno) %in% c("rsid", "ref", "alt", "most_severe_consequence", "consequences")],
            as.numeric, rep(numeric(1), nrow(geno))
        )


    # Add comments
    comment(geno) = comment(res[[1]])


    # Create GRanges container
    if (tolower(return_obj) == "granges") {
        geno$strand = "+"
        seq_lengths = setNames(
            as.list(avail_chromosomes()$length),
            avail_chromosomes()$chr
        )
        return(df2GRanges(geno, strand_name = "strand", 
                          seq_lengths = seq_lengths))
    } else {
        return(geno)
    }
}
