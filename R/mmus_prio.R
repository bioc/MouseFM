# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck(/path/to/project)


source("R/send_request.R")
source("R/make_query.R")


#'MMUS Prio
#'@description Prioritize mouse strains
#'@param chr Vector of chromosome names.
#'@param start Optional vector of chromosomal start positions of target regions (GRCm38).
#'@param end Optional vector of chromosomal end positions of target regions (GRCm38).
#'@param strain1 First strain.
#'@param strain2 Second strain.
#'@param consequence Optional vector of consequence types.
#'@param impact Optional vector of impact types.
#'@param max_set_size Maximum set of strains.
#'@param min_strain_benef Minimum reduction factor (min) of a single strain.
#'@return Dataframe
#'@examples res = mmusprio("chr1", start=5000000, end=6000000, strain1="C57BL_6J", strain2="AKR_J")
#'
#'comment(res$genotypes)
#'@export
mmusprio = function(chr, start = NULL, end = NULL, strain1, strain2, consequence = NULL, impact = NULL, min_strain_benef = 0.1, max_set_size = 3){


  stopifnot(is.vector(strain1), is.vector(strain2), length(strain1) == 1, length(strain2) == 1)


  # Create URL and query data
  res = lapply(seq_len(length(chr)), function(i){
    message(paste0("Query ", chr[i], if(is.numeric(start[i]) && is.numeric(end[i])) paste0(":", start[i], "-", end[i]) else ""))
    q = finemap_query(chr[i], start[i], end[i], strain1, strain2, consequence, impact)
    genehopper_request(q)
  })


  geno = as.data.frame(data.table::rbindlist(res))


  # Convert to respective data types
  geno[geno == "-" | geno == "."] = NA
  geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")] =
    vapply(geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")], as.numeric, rep(numeric(1), nrow(geno)))


  # Add comments
  comment(geno) = comment(res[[1]])


  # Calculate reduction factors
  message("Calculate reduction factors...")
  geno.add_strains = geno[!(tolower(names(geno)) %in% c("chr", "pos", "rsid", "ref", "alt", "consequences", tolower(strain1), tolower(strain2)))]
  rf = comb(geno.add_strains, min_strain_benef, max_set_size)


  rf = cbind(strain1 = rep(strain1, nrow(rf)), strain2 = rep(strain2, nrow(rf)), rf[rev(order(rf$min)),])
  row.names(rf) = seq_len(nrow(rf))


  return(list(genotypes = geno, reduction = rf))
}


#'Strain combination builder
#'@description Generate strain sets and calculate reduction factors
#'@param geno Data frame of genotypes for additional strains.
#'@param max_set_size Maximum set of strains. Default is 3.
#'@param min_strain_benef Minimum reduction factor (min) of a single strain. Default is 0.1.
#'@return Dataframe
#'@keywords internal
comb = function(geno, min_strain_benef = 0.1, max_set_size = 3){

  combination = NULL
  res.list = list()
  n_strains = ncol(geno)

  for(i in seq_len(max_set_size)){

    combs = as.data.frame(gtools::combinations(ncol(geno), i, colnames(geno)))
    message(paste0(nrow(combs), " combinations for set size ", i))
    red = reduction(combs, geno)

    if(i == 1 && max_set_size > 1 && min_strain_benef > 0){
      geno = geno[names(geno) %in% unlist(combs[red$min >= min_strain_benef,])]
      message(paste0(ncol(geno), " strains of ",  n_strains, " with min_strain_benef >= ", min_strain_benef))
    }


    # Extract and check for good combos
    df = cbind(tidyr::unite(combs, combination, sep=","), red)
    df$n = i
    res.list[[i]] = df
  }

  final = data.table::rbindlist(res.list)
  data.table::setDF(final)

  return(final)
}


#'Reduction factor calculation
#'@description Generate strain sets and calculate reduction factors
#'@param combs Data frame of strain sets.
#'@param geno Data frame of genotypes for additional strains.
#'@return Dataframe
#'@keywords internal
reduction = function(combs, geno){

  res = t(apply(combs, 1, function(x) {

    geno.subset = dplyr::as_tibble(geno[,x]) # For single column data frames
    rf = dplyr::mutate(dplyr::group_by_all(geno.subset), reduction_factor = 1-(dplyr::n()/nrow(geno.subset)))

    return(c(mean(rf$reduction_factor), min(rf$reduction_factor), max(rf$reduction_factor)))
  }))

  res = as.data.frame(res)
  colnames(res) = c("mean", "min", "max")

  return(res)
}


#'Best combinations
#'@description Get best strain combinations
#'@param rf Reduction factors data frame.
#'@param n_top Number if combinations to be returned.
#'@return Data frame
#'@examples l = mmusprio("chr1", start=5000000, end=6000000, strain1="C57BL_6J", strain2="AKR_J")
#'
#'get_top(l$reduction, 3)
#'@export
get_top = function(rf, n_top){

  rf = rf[rev(order(rf$min)),]
  top.n = rf[seq_len(min(nrow(rf), n_top)), ]

  return(top.n)
}


#'Visualize
#'@description Visualize reduction factors
#'@param geno Genotype data frame.
#'@param rf Reduction factor data frame.
#'@param n_top Number if combinations to be returned.
#'@return Dataframe
#'@examples l = mmusprio("chr1", start=5000000, end=6000000, strain1="C57BL_6J", strain2="AKR_J")
#'
#'plots = vis_reduction_factors(l$genotypes, l$reduction, 2)
#'
#'plots[[1]]
#'plots[[2]]
#'@export
vis_reduction_factors = function(geno, rf, n_top){

  pos = strain = allele = NULL

  chr = geno$chr[1]
  start = min(geno$pos)
  end = max(geno$pos)
  strain1 = rf$strain1[1]
  strain2 = rf$strain2[1]
  top.n = get_top(rf, n_top)

  geno = geno[order(geno$pos),]
  geno$pos = seq_len(nrow(geno))

  plots = list()

  for(i in seq_len(nrow(top.n))){

    strain.comb = unlist(strsplit(top.n$combination[i], ","))
    geno.filtered = geno[(tolower(names(geno)) %in% c("chr", "pos", tolower(c(strain1, strain2, strain.comb))))]

    geno.melt = reshape2::melt(geno.filtered, id.vars=c("chr", "pos"), variable.name="strain", value.name="allele")

    geno.melt$strain = ordered(geno.melt$strain, levels=rev(c(strain1, strain2, strain.comb)))

    p = ggplot2::ggplot(data=geno.melt, ggplot2::aes(x=pos, y=strain, fill=factor(allele))) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values=c("#3b5998", "#e9ebee","#f6f7f9"), guide = ggplot2::guide_legend(title = "Allele")) +
      ggplot2::theme_bw() +
      ggplot2::labs(y="Strain", x="Position", title=paste0("Region: chr", chr, ":", scales::comma(start), "-", scales::comma(end), " (GRCm38)"),
           subtitle=paste0("Additional strains: ", paste0(strain.comb, collapse=", "), " Reduction factor: ", round(top.n$min[i]), digits=3)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour="black"),
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0))

    plots[[paste0("top", i)]] = p
  }

  return(plots)
}

