# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck(/path/to/project)


source("R/send_request.R")
source("R/make_query.R")


#'MMUS Fetch
#'@description Detch genetic variant data for chromosomal region
#'@param chr Vector of chromosome names.
#'@param start Optional vector of chromosomal start positions of target regions (GRCm38).
#'@param end Optional vector of chromosomal end positions of target regions (GRCm38).
#'@param consequence Optional vector of consequence types.
#'@param impact Optional vector of impact types.
#'@param return_obj The user can choose to get the result to be returned
#'as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges") object. Default value is "dataframe".
#'@return Data frame or GenomicRanges::GRanges object containing result data.
#'@examples geno = mmusfetch("chr7", start=5000000, end=6000000)
#'
#'comment(geno)
#'@export
mmusfetch = function(chr,
                     start = NULL,
                     end = NULL,
                     consequence = NULL,
                     impact = NULL,
                     return_obj = "dataframe") {
    # Create URL and query data
    res = lapply(seq_len(length(chr)), function(i) {
        message(paste0("Query ", chr[i], if (is.numeric(start[i]) &&
                                             is.numeric(end[i]))
            paste0(":", start[i], "-", end[i])
            else
                ""))
        q = filter_query(chr[i], start[i], end[i], consequence, impact)
        genehopper_request(q)
    })


    geno = as.data.frame(data.table::rbindlist(res))


    # Add comments
    comment(geno) = comment(res[[1]])


    # Convert to respective data types
    geno[geno == "-" | geno == "."] = NA
    geno[!names(res) %in% c("rsid", "ref", "alt", "consequences")] =
        vapply(geno[!names(geno) %in% c("rsid", "ref", "alt", "consequences")],
               as.numeric, rep(numeric(1), nrow(geno)))


    if (tolower(return_obj) == "granges") {
        # Create GRanges container
        gres = GenomicRanges::makeGRangesFromDataFrame(
            geno,
            start.field = "pos",
            end.field = "pos",
            seqnames.field = "chr",
            keep.extra.columns = TRUE
        )

        GenomicRanges::strand(gres) = "+"
        comment(gres) = comment(geno)

        return(gres)

    } else {
        return(geno)
    }
}
