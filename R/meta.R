# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


#' Reference genome version
#' @description Returns version of reference genome used in package MouseFM.
#' @return Vector.
#' @examples
#' ref_genome()
#' @export
ref_genome = function() {
  return("GRCm38")
}


#' Available strains
#' @description There are 37 strains available.
#' @return Data frame.
#' @examples
#' avail_strains()
#' @export
avail_strains = function() {
  df = data.frame(
    id = c(
      "129P2_OlaHsd",
      "129S1_SvImJ",
      "129S5SvEvBrd",
      "A_J",
      "AKR_J",
      "BALB_cJ",
      "BTBR",
      "BUB_BnJ",
      "C3H_HeH",
      "C3H_HeJ",
      "C57BL_10J",
      "C57BL_6J",
      "C57BL_6NJ",
      "C57BR_cdJ",
      "C57L_J",
      "C58_J",
      "CAST_EiJ",
      "CBA_J",
      "DBA_1J",
      "DBA_2J",
      "FVB_NJ",
      "I_LnJ",
      "KK_HiJ",
      "LEWES_EiJ",
      "LP_J",
      "MOLF_EiJ",
      "NOD_ShiLtJ",
      "NZB_B1NJ",
      "NZO_HlLtJ",
      "NZW_LacJ",
      "PWK_PhJ",
      "RF_J",
      "SEA_GnJ",
      "SPRET_EiJ",
      "ST_bJ",
      "WSB_EiJ",
      "ZALENDE_EiJ"
    ),
    strain = c(
      "129P2/OlaHsd",
      "129S1/SvImJ",
      "129S5/SvEvBrd",
      "A/J",
      "AKR/J",
      "BALB/cJ",
      "BTBR",
      "BUB/BnJ",
      "C3H/HeH",
      "C3H/HeJ",
      "C57BL/10J",
      "C57BL/6J",
      "C57BL/6NJ",
      "C57BR/cdJ",
      "C57L/J",
      "C58/J",
      "CAST/EiJ",
      "CBA/J",
      "DBA/1J",
      "DBA/2J",
      "FVB/NJ",
      "I/LnJ",
      "KK/HiJ",
      "LEWES/EiJ",
      "LP/J",
      "MOLF/EiJ",
      "NOD/ShiLtJ",
      "NZB/B1NJ",
      "NZO/HlLtJ",
      "NZW/LacJ",
      "PWK/PhJ",
      "RF/J",
      "SEA/GnJ",
      "SPRET/EiJ",
      "ST/bJ",
      "WSB/EiJ",
      "ZALENDE/EiJ"
    )
  )

  return(df)
}


#' Available consequences
#' @description Available consequence and impact types.
#' @return Data frame.
#' @examples
#' avail_consequences()$consequence
#'
#' unique(avail_consequences()$impact)
#' @export
avail_consequences = function() {
    df = data.frame(
        consequence = c(
            "splice_acceptor_variant",
            "splice_donor_variant",
            "stop_gained",
            "stop_lost",
            "splice_region_variant",
            "incomplete_terminal_codon_variant",
            "stop_retained_variant",
            "synonymous_variant",
            "missense_variant",
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant",
            "intron_variant",
            "NMD_transcript_variant",
            "non_coding_transcript_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "intergenic_variant",
            "start_lost"
        ),
        impact = c(
            "HIGH",
            "HIGH",
            "HIGH",
            "HIGH",
            "LOW",
            "LOW",
            "LOW",
            "LOW",
            "MODERATE",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "MODIFIER",
            "HIGH"
        ),
        so_description = c(
            "A splice variant that changes the 2 base region at the 3' end of an 
            intron",
            "A splice variant that changes the 2 base region at the 5' end of an 
            intron",
            "A sequence variant whereby at least one base of a codon is changed, 
            resulting in a premature stop codon, leading to a shortened 
            transcript",
            "A sequence variant where at least one base of the terminator codon 
            (stop) is changed, resulting in an elongated transcript",
            "A sequence variant in which a change has occurred within the region 
            of the splice site, either within 1-3 bases of the exon or 3-8 bases 
            of the intron",
            "A sequence variant where at least one base of the final codon of an 
            incompletely annotated transcript is changed",
            "A sequence variant where at least one base in the terminator codon is 
            changed, but the terminator remains",
            "A sequence variant where there is no resulting change to the encoded 
            amino acid",
            "A sequence variant, that changes one or more bases, resulting in a 
            different amino acid sequence but where the length is preserved",
            "A sequence variant that changes the coding sequence",
            "A transcript variant located with the sequence of the mature 
            miRNA",
            "A UTR variant of the 5' UTR",
            "A UTR variant of the 3' UTR",
            "A sequence variant that changes non-coding exon sequence in a 
            non-coding transcript",
            "A transcript variant occurring within an intron",
            "A variant in a transcript that is the target of NMD",
            "A transcript variant of a non coding RNA gene",
            "A sequence variant located 5' of a gene",
            "A sequence variant located 3' of a gene",
            "A sequence variant located in the intergenic region, between genes",
            "A codon variant that changes at least one base of the canonical start codon"
        ),
        so_term = c(
            "SO:0001574",
            "SO:0001575",
            "SO:0001587",
            "SO:0001578",
            "SO:0001630",
            "SO:0001626",
            "SO:0001567",
            "SO:0001819",
            "SO:0001583",
            "SO:0001580",
            "SO:0001620",
            "SO:0001623",
            "SO:0001624",
            "SO:0001792",
            "SO:0001627",
            "SO:0001621",
            "SO:0001619",
            "SO:0001631",
            "SO:0001632",
            "SO:0001628",
            "SO:0002012"
        ),
        display_term = c(
            "Splice acceptor variant",
            "Splice donor variant",
            "Stop gained",
            "Stop lost",
            "Splice region variant",
            "Incomplete terminal codon variant",
            "Stop retained variant",
            "Synonymous variant",
            "Missense variant",
            "Coding sequence variant",
            "Mature miRNA variant",
            "5 prime UTR variant",
            "3 prime UTR variant",
            "Non coding transcript exon variant",
            "Intron variant",
            "NMD transcript variant",
            "Non coding transcript variant",
            "Upstream gene variant",
            "Downstream gene variant",
            "Intergenic variant",
            "Start lost"
        )
    )
    return(df)
}


#' Available chromosomes
#' @description Available mouse chromosomes.
#' @return Data frame
#' @examples
#' avail_chromosomes()
#' @export
avail_chromosomes = function() {
  df = data.frame(
    chr = seq_len(19),
    length = c(
      195471971,
      182113224,
      160039680,
      156508116,
      151834684,
      149736546,
      145441459,
      129401213,
      124595110,
      130694993,
      122082543,
      120129022,
      120421639,
      124902244,
      104043685,
      98207768,
      94987271,
      90702639,
      61431566
    ),
    min_pos = c(
      3000185,
      3050115,
      3000104,
      3050235,
      3000846,
      3050051,
      3001115,
      3000222,
      3000063,
      3100281,
      3100189,
      3000048,
      3000177,
      3000306,
      3050014,
      3000350,
      3000030,
      3000033,
      3000287
    ),
    max_pos = c(
      195371784,
      182012539,
      159939640,
      156357848,
      151734582,
      149586443,
      145340897,
      129299441,
      124491297,
      130594765,
      121977829,
      120028129,
      120321311,
      124801447,
      103943534,
      98107534,
      94886990,
      90601947,
      61330097
    )
  )

  return(df)
}
