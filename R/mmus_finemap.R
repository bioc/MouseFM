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
#'@param chr Vector of chromosome names.
#'@param start Optional vector of chromosomal start positions of target regions (GRCm38).
#'@param end Optional vector of chromosomal end positions of target regions (GRCm38).
#'@param strain1 First strain(s).
#'@param strain2 Second strain(s).
#'@param consequence Optional vector of consequence types.
#'@param impact Optional vector of impact types.
#'@param thr1 Number discordant strains in strain1. Between 0 and length(strain1)-1. 0 by default.
#'@param thr2 Number discordant strains in strain2. Between 0 and length(strain2)-1. 0 by default.
#'@param return_obj The user can choose to get the result to be returned
#'as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges") object. Default value is "dataframe".
#'@return Data frame or GenomicRanges::GRanges object containing result data.
#'@examples
#'mmusfinemap("chr1", start=5000000, end=6000000,
#'            strain1=c("C57BL_6J"), strain2=c("129S1_SvImJ", "129S5SvEvBrd", "AKR_J"))
#'
#'mmusfinemap(chr=c("chr7", "1"),strain1=c("C57BL_6J","C57L_J","CBA_J","NZB_B1NJ"),
#'            strain2=c("C3H_HEJ","MOLF_EiJ","NZW_LacJ","WSB_EiJ","SPRET_EiJ"))
#'
#'mmusfinemap("chr1", start=5000000, end=6000000,
#'            strain1=c("C57BL_6J"), strain2=c("NOD_ShiLtJ", "CBA_J", "KK_HiJ", "PWK_PhJ"), thr2=1)
#'@export
mmusfinemap = function(chr, start = NULL, end = NULL, strain1, strain2, consequence = NULL, impact = NULL, thr1 = 0, thr2 = 0, return_obj = "dataframe"){


  # Create URL and query data
  res = lapply(1:length(chr), function(i){
    message(paste0("Query ", chr[i], if(is.numeric(start[i]) && is.numeric(end[i])) paste0(":", start[i], "-", end[i]) else ""))
    q = finemap_query(chr[i], start[i], end[i], strain1, strain2, consequence, impact, thr1, thr2)
    #genehopper_request(q)
  })

  return()

  # geno = as.data.frame(data.table::rbindlist(res))
  #
  #
  # # Convert to respective data types
  # geno[geno == "-" | geno == "."] = NA
  # geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")] =
  #   sapply(geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")], as.numeric)
  #
  #
  # # Keep only input strains
  # geno = geno[tolower(names(geno)) %in% c("rsid", "ref", "alt", "consequences",
  #                                      tolower(unique(strain1)), tolower(unique(strain2)))]
  #
  # # Add comments
  # comment(geno) = comment(res[[1]])
  #
  #
  # if(tolower(return_obj) == "granges"){
  #   # Create GRanges container
  #   gres = GenomicRanges::makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos",
  #                                                  seqnames.field = "chr", keep.extra.columns = TRUE)
  #
  #   GenomicRanges::strand(gres) = "+"
  #   comment(gres) = comment(geno)
  #
  #   return(gres)
  #
  # } else {
  #   return(geno)
  # }
}

