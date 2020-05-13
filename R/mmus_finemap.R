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


#'MMUS Finemapping
#'@description Finemapping
#'@param chr Chromosome name.
#'@param start Optional chromosomal start position of target region (GRCm38). NA by default.
#'@param end Optional chromosomal end position of target region (GRCm38). NA by default.
#'@param strain1 First strain(s).
#'@param strain2 Second strain(s).
#'@param consequence Vector containing consequence types. NA by default.
#'@param impact Vector containing impact types. NA by default.
#'@param thr1 Number discordant strains in strain1. Between 0 and length(strain1)-1. 0 by default.
#'@param thr2 Number discordant strains in strain2. Between 0 and length(strain2)-1. 0 by default.
#'@param return_obj The user can choose to get the result to be returned
#'as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges") object. Default value is "dataframe".
#'@return Data frame or GenomicRanges::GRanges object containing result data.
#'@examples mmusfinemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("129S1_SvImJ", "129S5SvEvBrd", "AKR_J"))
#'
#'mmusfinemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("AKR_J", "A_J", "BALB_cJ"))
#'
#'mmusfinemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("NOD_ShiLtJ", "CBA_J", "KK_HiJ", "PWK_PhJ"), thr2=1)
#'@export
mmusfinemap = function(chr, start = NA, end = NA, strain1, strain2, consequence = NA, impact = NA, thr1 = 0, thr2 = 0, return_obj = "dataframe"){

  # Check if there is an internet connection
  if (!curl::has_internet())
    stop("No internet connection detected...")

  # Create URL
  message("Build query...")
  q = finemap_query(chr, start, end, strain1, strain2, consequence, impact, thr1, thr2)


  # Send query
  message("Retrieve data...")
  res = genehopper_request(q)


  # Convert to respective data types
  res[res == "-" | res == "."] = NA
  res[! names(res) %in% c("rsid", "ref", "alt", "consequences")] =
    sapply(res[! names(res) %in% c("rsid", "ref", "alt", "consequences")], as.numeric)


  # Keep only input strains
  res = res[tolower(names(res)) %in% c("rsid", "ref", "alt", "consequences",
                                       tolower(unique(strain1)), tolower(unique(strain2)))]

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

