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
#'@param chr Chromosome name or vector of chromosome names.
#'@param start Optional Chromosomal start position of target region (GRCm38). Multiple positions can be passed as vector. NA by default.
#'@param end Optional Chromosomal end position of target region (GRCm38). Multiple positions can be passed as vector. NA by default.
#'@param consequence Vector containing consequence types. NA by default.
#'@param impact Vector containing impact types. NA by default.
#'@param return_obj The user can choose to get the result to be returned
#'as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges") object. Default value is "dataframe".
#'@return Data frame or GenomicRanges::GRanges object containing result data.
#'@examples mmusfetch("chr7", start=5000000, end=6000000)
#'@export
mmusfetch = function(chr, start = NA, end = NA, consequence = NA, impact = NA, return_obj = "dataframe"){

  # Check if there is an internet connection
  if (!curl::has_internet())
    stop("No internet connection detected...")


  # Create URL and query data
  res = lapply(1:length(chr), function(i){
    message(paste0("Query ", chr[i], if(is.numeric(start[i]) && is.numeric(end[i])) paste0(":", start[i], "-", end[i]) else ""))
    q = filter_query(chr[i], start[i], end[i], consequence, impact)
    genehopper_request(q)
  })


  geno = as.data.frame(data.table::rbindlist(res))


  # Add comments
  comment(geno) = comment(res[[1]])


  # Convert to respective data types
  geno[geno == "-" | geno == "."] = NA
  geno[! names(res) %in% c("rsid", "ref", "alt", "consequences")] =
    sapply(geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")], as.numeric)


  if(tolower(return_obj) == "granges"){
    # Create GRanges container
    gres = GenomicRanges::makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos",
                                                   seqnames.field = "chr", keep.extra.columns = TRUE)

    GenomicRanges::strand(gres) = "+"
    comment(gres) = comment(geno)

    return(gres)

  } else {
    return(geno)
  }
}


