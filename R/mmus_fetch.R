# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()


source("R/send_request.R")
source("R/make_query.R")


#'MMUS Fetch
#'@description Detch genetic variant data for chromosomal region
#'@param chr Chromosome name.
#'@param start Optional chromosomal start position of target region (GRCm38). NA by default.
#'@param end Optional chromosomal end position of target region (GRCm38). NA by default.
#'@param consequence Vector containing consequence types. NA by default.
#'@param impact Vector containing impact types. NA by default.
#'@param return_obj The user can choose to get the result to be returned
#'as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges") object. Default value is "dataframe".
#'@return Data frame or GenomicRanges::GRanges object containing result data.
#'@examples mmusfetch("chr1", start=5000000, end=6000000)
#'@export
mmusfetch = function(chr, start = NA, end = NA, consequence = NA, impact = NA, return_obj = "dataframe"){

  # Check if there is an internet connection
  if (!curl::has_internet())
    stop("No internet connection detected...")


  # Create URL
  message("Build query...")
  q = filter_query(chr, start, end, consequence, impact)


  # Send query
  message("Retrieve data...")
  res = genehopper_request(q)


  # Convert to respective data types
  res[res == "-" | res == "."] = NA
  res[! names(res) %in% c("rsid", "ref", "alt", "consequences")] =
    sapply(res[! names(res) %in% c("rsid", "ref", "alt", "consequences")], as.numeric)

  if(tolower(return_obj) == "granges"){
    # Create GRanges container
    gres = GenomicRanges::makeGRangesFromDataFrame(res, start.field = "pos", end.field = "pos",
                                                   seqnames.field = "chr", keep.extra.columns = TRUE)

    GenomicRanges::strand(gres) = "+"
    comment(gres) = comment(res)

    return(gres)

  } else {
    return(res)
  }
}


