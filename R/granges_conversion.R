# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck(/path/to/project)


#'Data frame to GenomicRanges::GRanges object
#'@description Wrapper for GenomicRanges::makeGRangesFromDataFrame().
#'@param geno Data frame
#'@param chr_name Name of chromosome column. Default is 'chr'
#'@param start_name Name of start position column. Default is 'pos'
#'@param end_name Name of end position column. Default is 'pos'
#'@param strand Strand. Default is '+'.
#'@param ref_genome Reference genome version. Default is 'GRCm38'
#'@return GenomicRanges::GRanges object.
#'@examples geno = finemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("AKR_J", "A_J", "BALB_cJ"))
#'
#'geno.granges = df2GRanges(geno)
#'@export
df2GRanges = function(geno,
                      chr_name = "chr",
                      start_name = "pos",
                      end_name = "pos",
                      strand = "+",
                      ref_genome = "GRCm38") {
    gres = GenomicRanges::makeGRangesFromDataFrame(
        geno,
        start.field = start_name,
        end.field = end_name,
        seqnames.field = chr_name,
        keep.extra.columns = TRUE
    )

    GenomicRanges::strand(gres) = strand
    GenomeInfoDb::genome(gres) = ref_genome
    comment(gres) = comment(geno)

    return(gres)

}


#'GenomicRanges::GRanges object to data frame
#'@description Wrapper for as.data.frame().
#'@param granges GenomicRanges::GRanges object
#'@return Data frame.
#'@examples geno.granges = finemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("AKR_J", "A_J", "BALB_cJ"), return_obj="granges")
#'
#'geno = GRanges2df(geno.granges)
#'@export
GRanges2df = function(granges) {

    geno = as.data.frame(granges)
    comment(geno) = comment(granges)

    return(geno)

}
