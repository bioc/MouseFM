# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


#' Data frame to GenomicRanges::GRanges object
#' @description Wrapper for GenomicRanges::makeGRangesFromDataFrame().
#' @param geno Data frame.
#' @param chr_name Name of chromosome column. Default is 'chr'.
#' @param start_name Name of start position column. Default is 'pos.'
#' @param end_name Name of end position column. Default is 'pos'
#' @param strand_name Name of end position column. Default is NULL.
#' @param ref_version Reference genome version. Default is 'ref_genome()'.
#' @param seq_lengths List of sequence lengths with sequence name as key.
#' Default is NULL.
#' @param is_circular Whether genome is circular. Default is FALSE.
#' @return GenomicRanges::GRanges object.
#' @examples
#' geno = finemap("chr1",
#'   start = 5000000, end = 6000000,
#'   strain1 = c("C57BL_6J"), strain2 = c("AKR_J", "A_J", "BALB_cJ")
#' )
#'
#' geno$strand = "+"
#' seq_lengths = stats::setNames(
#'     as.list(avail_chromosomes()$length),
#'     avail_chromosomes()$chr
#' )
#' geno.granges = df2GRanges(geno,
#'     strand_name = "strand",
#'     seq_lengths = seq_lengths
#' )
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame strand<-
#' @importFrom GenomeInfoDb genome<- isCircular<- seqlevels seqlengths<-
df2GRanges = function(geno,
                       chr_name = "chr",
                       start_name = "pos",
                       end_name = "pos",
                       strand_name = NULL,
                       ref_version = ref_genome(),
                       seq_lengths = NULL,
                       is_circular = FALSE) {
  gres = makeGRangesFromDataFrame(
    geno,
    start.field = start_name,
    end.field = end_name,
    seqnames.field = chr_name,
    keep.extra.columns = TRUE
  )

  if (!is.null(strand_name)) {
    strand(gres) = geno[, strand_name]
  }


  # SeqInfo
  if (!is.null(ref_version)) {
    genome(gres) = ref_version
  }

  if (!is.null(is_circular)) {
    isCircular(gres) = rep(is_circular, length(seqlevels(gres)))
  }

  if (!is.null(seq_lengths)) {
    seqlengths(gres) =
      unlist(seq_lengths[seqlevels(gres)], use.names = FALSE)
  }


  comment(gres) = comment(geno)

  return(gres)
}


#' GenomicRanges::GRanges object to data frame
#' @description Wrapper for as.data.frame().
#' @param granges GenomicRanges::GRanges object
#' @return Data frame.
#' @examples
#' geno.granges = finemap("chr1",
#'     start = 5000000, end = 6000000,
#'     strain1 = c("C57BL_6J"), strain2 = c("AKR_J", "A_J", "BALB_cJ"),
#'     return_obj = "granges"
#' )
#'
#' geno = GRanges2df(geno.granges)
#' @export
GRanges2df = function(granges) {
  geno = as.data.frame(granges)
  
  names(geno)[1] = "chr"
  
  if(all(geno$start == geno$end)){
    geno = geno[ , -which(names(geno) %in% c("end"))]
    names(geno)[2] = "pos"
  }
  
  comment(geno) = comment(granges)

  return(geno)
}
