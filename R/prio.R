# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck("/path/to/project")


# source("R/send_request.R")
# source("R/make_query.R")
# source("R/granges_conversion.R")


#' Prioritization of inbred mouse strains for refining genetic regions
#' @description This method allows to select strain combinations which best
#' refine a specified genetic region (GRCm38). E.g. if a crossing experiment
#' with two inbred mouse strains 'strain1' and 'strain2' resulted in a QTL, the
#' outputted strain combinations can be used to refine the respective region in
#' further crossing experiments.
#' @param chr Vector of chromosome names.
#' @param start Optional vector of chromosomal start positions of target regions
#' (GRCm38).
#' @param end Optional vector of chromosomal end positions of target regions
#' (GRCm38).
#' @param strain1 First strain set with strains from avail_strains().
#' @param strain2 Second strain set with strains from avail_strains().
#' @param consequence Optional vector of consequence types.
#' @param impact Optional vector of impact types.
#' @param max_set_size Maximum set of strains.
#' @param min_strain_benef Minimum reduction factor (min) of a single strain.
#' @param return_obj The user can choose to get the result to be returned
#' as data frame ("dataframe") or as a GenomicRanges::GRanges ("granges")
#' object. Default value is "data frame".
#' @return Data frame
#' @examples
#' res = prio("chr1",
#'   start = 5000000, end = 6000000, strain1 = "C57BL_6J",
#'   strain2 = "AKR_J"
#' )
#'
#' comment(res$genotypes)
#' @export
#' @importFrom scales comma
#' @importFrom data.table rbindlist
#' @importFrom stats setNames complete.cases
prio = function(chr,
                 start = NULL,
                 end = NULL,
                 strain1 = NULL,
                 strain2 = NULL,
                 consequence = NULL,
                 impact = NULL,
                 min_strain_benef = 0.1,
                 max_set_size = 3,
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
      impact
    )
    backend_request(q)
  })


  geno = as.data.frame(rbindlist(res))


  # Return if no results
  if (nrow(geno) == 0) {
    return(list())
  }

  
  # Remove variants with missing genotypes
  geno.strains = geno[!(
    tolower(names(geno)) %in% c(
      "chr",
      "pos",
      "rsid",
      "ref",
      "alt",
      "most_severe_consequence",
      "consequences"
    )
  )]
  geno = geno[complete.cases(geno.strains),]
  

  # Convert to respective data types
  geno[!names(geno) %in% c("rsid", "ref", "alt", "most_severe_consequence", "consequences")] =
    vapply(
      geno[!names(geno) %in% c("rsid", "ref", "alt", "most_severe_consequence", "consequences")],
      as.numeric, rep(numeric(1), nrow(geno))
    )


  # Add comments
  comment(geno) = comment(res[[1]])

  
  # Extract additional strains
  geno.add_strains = geno[!(
    tolower(names(geno)) %in% c(
      "chr",
      "pos",
      "rsid",
      "ref",
      "alt",
      "most_severe_consequence",
      "consequences",
      tolower(strain1),
      tolower(strain2)
    )
  )]


  # Calculate reduction factors
  message("Calculate reduction factors...")
  red = comb(geno.add_strains, min_strain_benef, max_set_size)


  # Create final reduction factor data frame
  red = cbind(
    strain1 = rep(if (is.null(strain1)) {
      NA
    } else {
      paste(sort(strain1), collapse = ",")
    }, nrow(red)),
    strain2 = rep(if (is.null(strain2)) {
      NA
    } else {
      paste(sort(strain2), collapse = ",")
    }, nrow(red)),
    red[rev(order(red$min)), ]
  )
  row.names(red) = seq_len(nrow(red))


  # Create GRanges container
  if (tolower(return_obj) == "granges") {
    geno$strand = "+"
    l = setNames(
      as.list(avail_chromosomes()$length),
      avail_chromosomes()$chr
    )
    return(list(genotypes = df2GRanges(geno, strand_name = "strand", seq_lengths = l), reduction = red))
  }


  return(list(genotypes = geno, reduction = red))
}


#' Strain combination builder
#' @description Generate strain sets and calculate reduction factors
#' @param geno Data frame of genotypes for additional strains.
#' @param max_set_size Maximum set of strains. Default is 3.
#' @param min_strain_benef Minimum reduction factor (min) of a single strain.
#' Default is 0.1.
#' @return Data frame
#' @keywords internal
#' @importFrom gtools combinations
#' @importFrom scales comma
#' @importFrom tidyr unite
#' @importFrom data.table rbindlist setDF
comb = function(geno,
                 min_strain_benef = 0.1,
                 max_set_size = 3) {
  combination = NULL
  res.list = list()
  n_strains = ncol(geno)

  for (i in seq_len(max_set_size)) {
    combs = as.data.frame(combinations(ncol(geno), i, colnames(geno)))
    message(paste0(
      "Set size ", i, ": ",
      comma(nrow(combs)), " combinations"
    ))
    red = reduction(combs, geno)

    if (i == 1 && max_set_size > 1 && min_strain_benef > 0) {
      # Filter
      geno = geno[names(geno) %in%
        unlist(combs[red$min >= min_strain_benef, ])]

      message(paste0(
        "Set size ",
        i,
        ": continue with ",
        comma(ncol(geno)),
        " of ",
        comma(n_strains),
        " strains"
      ))
    }


    # Extract and check for good combos
    df = cbind(unite(combs, combination, sep = ","), red)
    df$n = i
    res.list[[i]] = df
  }

  final = rbindlist(res.list)
  setDF(final)

  return(final)
}


#' Reduction factor calculation
#' @description Generate strain sets and calculate reduction factors
#' @param combs Data frame of strain sets.
#' @param geno Data frame of genotypes for additional strains.
#' @return Data frame
#' @keywords internal
#' @importFrom dplyr as_tibble mutate group_by_all n
reduction = function(combs, geno) {
  res = t(apply(combs, 1, function(x) {
    geno.subset = as_tibble(geno[, x]) # For single column data frames
    red = mutate(group_by_all(geno.subset),
      reduction_factor = 1 - (n() / nrow(geno.subset))
    )

    return(c(
      mean(red$reduction_factor),
      min(red$reduction_factor),
      max(red$reduction_factor)
    ))
  }))

  res = as.data.frame(res)
  colnames(res) = c("mean", "min", "max")

  return(res)
}


#' Best strain combinations
#' @description Get best strain combinations
#' @param red Reduction factors data frame.
#' @param n_top Number of combinations to be returned.
#' @return Data frame
#' @examples
#' l = prio("chr1",
#'   start = 5000000, end = 6000000,
#'   strain1 = "C57BL_6J", strain2 = "AKR_J"
#' )
#'
#' get_top(l$reduction, 3)
#' @export
get_top = function(red, n_top) {
  red = red[rev(order(red$min)), ]
  top.n = red[seq_len(min(nrow(red), n_top)), ]

  return(top.n)
}


#' Visualize
#' @description Visualize reduction factors
#' @param geno Genotype data frame or GenomicRanges::GRanges object.
#' @param red Reduction factor data frame.
#' @param n_top Number if combinations to be returned.
#' @return Data frame
#' @examples
#' l = prio(c("chr1", "chr2"),
#'   start = c(5000000, 5000000),
#'   end = c(6000000, 6000000), strain1 = c("C3H_HeH"), strain2 = "AKR_J"
#' )
#'
#' plots = vis_reduction_factors(l$genotypes, l$reduction, 2)
#'
#' plots[[1]]
#' plots[[2]]
#' @export
#' @importFrom rlist list.cbind
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual guide_legend
#' theme_bw labs theme element_blank element_rect element_text scale_x_discrete
#' scale_y_discrete
#' @importFrom scales comma
vis_reduction_factors = function(geno, red, n_top) {
  if ("GRanges" %in% is(geno)) {
    geno = GRanges2df(geno)
  }
  
  
  pos = strain = allele = NULL

  
  strain1_vec = unlist(strsplit(red$strain1[1], ","))
  strain2_vec = unlist(strsplit(red$strain2[1], ","))
  strain1 = red$strain1[1]
  strain2 = red$strain2[1]
  top.n = get_top(red, n_top)

  
  geno = geno[order(geno$chr, geno$pos), ]
  geno$pos = seq_len(nrow(geno))


  # Add genotype of strain1 if it's a combo
  if (length(strain1_vec) > 1) {
    s1 = rowSums(geno[tolower(names(geno)) %in%
      tolower(strain1_vec)]) / length(strain1_vec)
    geno[, strain1] = s1
  }


  # Add genotype of strain2 if it's a combo
  if (length(strain2_vec) > 1) {
    s2 = rowSums(geno[tolower(names(geno)) %in%
      tolower(strain2_vec)]) / length(strain2_vec)
    geno[, strain2] = s2
  }


  # All strains except strain1 and strain2
  other_strains = names(geno)[!tolower(names(geno)) %in%
    tolower(c(
      "chr", "pos", "rsid", "ref", "alt",
      "consequences", strain1, strain2
    ))]


  # Replace 0/1 by strain1/strain2 name
  geno.mod = lapply(other_strains, function(s) {
    vapply(seq_len(nrow(geno)), function(i) {
      if (geno[i, strain1] == geno[i, s]) {
        strain1
      } else if (geno[i, strain2] == geno[i, s]) {
        strain2
      } else {
        "NA"
      }
    }, character(1))
  })


  # Reformat
  geno.mod = as.data.frame(list.cbind(geno.mod))
  colnames(geno.mod) = other_strains
  geno = cbind(
    geno[c("chr", "pos", "rsid", "ref", "alt", "consequences")],
    geno.mod
  )
  geno[, strain1] = strain1
  geno[, strain2] = strain2


  plots = list()


  for (i in seq_len(nrow(top.n))) {
    strain.comb = unlist(strsplit(top.n$combination[i], ","))

    geno.filtered = geno[(tolower(names(geno)) %in%
      c("chr", "pos", tolower(c(strain1, strain2, strain.comb))))]

    geno.melt = melt(
      geno.filtered,
      id.vars = c("chr", "pos"),
      variable.name = "strain",
      value.name = "allele"
    )

    geno.melt$strain = ordered(geno.melt$strain,
      levels = rev(c(strain1, strain2, strain.comb))
    )

    p = ggplot(data = geno.melt, aes(
      x = pos,
      y = strain,
      fill = factor(allele)
    )) +
      geom_tile() +
      scale_fill_manual(
        values = c("#3b5998", "#e9ebee", "#f6f7f9"),
        guide = guide_legend(title = "Allele")
      ) +
      theme_bw() +
      labs(
        y = "Strain",
        x = "SNP (Ordered by chromosomal position)",
        title = paste0(
          "Strain combination: ",
          paste0(strain.comb, collapse = ", ")
        ),
        subtitle = paste0(
          " Reduction factor: ",
          round(top.n$min[i], digits = 3),
          " | ",
          paste0(comma(nrow(geno)), " SNV(s)")
        )
      ) +
      theme(
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0))

    plots[[paste0("top", i)]] = p
  }

  return(plots)
}
