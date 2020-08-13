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


#' Fetch
#' @description Fetch homozygous genotypes for a specified chromosomal region in
#' 37 inbred mouse strains.
#' @param chr Vector of chromosome names.
#' @param start Optional vector of chromosomal start positions of target regions
#' (GRCm38).
#' @param end Optional vector of chromosomal end positions of target regions
#' (GRCm38).
#' @param consequence Optional vector of consequence types.
#' @param impact Optional vector of impact types.
#' @param return_obj The user can choose to get the result to be returned
#' as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges")
#' object. Default value is "dataframe".
#' @return Data frame or GenomicRanges::GRanges object containing result data.
#' @examples
#' geno = fetch("chr7", start = 5000000, end = 6000000)
#'
#' comment(geno)
#' @export
#' @importFrom scales comma
#' @importFrom data.table rbindlist
#' @importFrom stats setNames
fetch = function(chr,
                  start = NULL,
                  end = NULL,
                  consequence = NULL,
                  impact = NULL,
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
        q = finemap_query(chr[i], start[i], end[i],
            consequence = consequence,
            impact = impact
        )
        backend_request(q)
    })


    geno = as.data.frame(rbindlist(res))


    # Return if no results
    if (nrow(geno) == 0) {
        return(geno)
    }


    # Add comments
    comment(geno) = comment(res[[1]])


    # Convert to respective data types
    geno[!names(res) %in% c("rsid", "ref", "alt", "most_severe_consequence", "consequences")] =
        vapply(
            geno[!names(geno) %in% c("rsid", "ref", "alt", "most_severe_consequence", "consequences")],
            as.numeric, rep(numeric(1), nrow(geno))
        )


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
