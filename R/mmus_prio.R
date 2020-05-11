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
#   https://github.com/khodaygani/Thesis


source("R/send_request.R")
source("R/make_query.R")


#'MMUS Prio
#'@description Prioritize mouse strains
#'@param chr Chromosome name.
#'@param start Optional chromosomal start position of target region (GRCm38). NA by default.
#'@param end Optional chromosomal end position of target region (GRCm38). NA by default.
#'@param strain1 First strain.
#'@param strain2 Second strain.
#'@param consequence Vector containing consequence types. NA by default.
#'@param impact Vector containing impact types. NA by default.
#'@param max_set_size Maximum set of strains.
#'@param min_strain_benef Minimum reduction factor (min) of a single strain.
#'@return Dataframe
#'@examples df = mmusprio("chr1", start=5000000, end=6000000)
#'mmusprio("chr1", start=5000000, end=6000000, strain1="C57BL_6J", strain2="AKR_J")
#'@export
mmusprio = function(chr, start = NA, end = NA, strain1, strain2, consequence = NA, impact = NA, min_strain_benef = 0.1, max_set_size = 3){


  # Check if there is an internet connection
  if (!curl::has_internet())
    stop("No internet connection detected...")


  stopifnot(is.vector(strain1), is.vector(strain2), length(strain1) == 1, length(strain2) == 1)


  # Create URL
  message("Build query...")
  q = finemap_query(chr, start, end, strain1, strain2, consequence, impact)


  # Send query
  message("Retrieve data...")
  geno = genehopper_request(q)


  # Convert to respective data types
  geno[geno == "-" | geno == "."] = NA
  geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")] =
    sapply(geno[! names(geno) %in% c("rsid", "ref", "alt", "consequences")], as.numeric)


  # Calculate reduction factors
  message("Calculate reduction factors...")
  geno.add_strains = res[!(tolower(names(geno)) %in% c("chr", "pos", "rsid", "ref", "alt", "consequences", tolower(strain1), tolower(strain2)))]
  rf = comb(geno.add_strains, min_strain_benef, max_set_size)

  rf = cbind(strain1 = rep(strain1, nrow(rf)), strain2 = rep(strain2, nrow(rf)), rf[rev(order(rf$min)),])
  row.names(rf) <- 1:nrow(rf)

  return(list(genotypes = geno, reduction = rf))
}


#'@description Generate strain sets and calculate reduction factors
#'@param geno Data frame of genotypes for additional strains.
#'@param max_set_size Maximum set of strains.
#'@param min_strain_benef Minimum reduction factor (min) of a single strain.
#'@return Dataframe
#'@keywords internal
comb = function(geno, min_strain_benef = 0.1, max_set_size = 3){

  res.list = list()
  n_strains = ncol(geno)

  for(i in 1:max_set_size){

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


#'@description Generate strain sets and calculate reduction factors
#'@param combs Data frame of strain sets.
#'@param geno Data frame of genotypes for additional strains.
#'@return Dataframe
#'@keywords internal
reduction = function(combs, geno){

  res = t(apply(combs, 1, function(x) {

    geno.subset = dplyr::as_tibble(geno[,x]) # For single column data frames
    rf = dplyr::mutate(dplyr::group_by_all(geno.subset), reduction_factor = 1-(n()/nrow(geno.subset)))

    return(c(mean(rf$reduction_factor), min(rf$reduction_factor), max(rf$reduction_factor)))
  }))

  res = as.data.frame(res)
  colnames(res) = c("mean", "min", "max")

  return(res)
}


r = mmusprio("chr1", start=5000000, end=6000000, strain1="DBA_1J", strain2="AKR_J")
geno = r$genotypes
reduction_factors = r$reduction


#'@description Get best combinations
#'@param rf Reduction factors data frame.
#'@param n_top Number if combinations to be returned.
#'@return Dataframe
#'@examples r = mmusprio("chr1", start=5000000, end=6000000, strain1="C57BL_6J", strain2="AKR_J")
#'get_top(r.reduction_factors, 3)
#'@export
get_top = function(rf, n_top){

  reduction_factors = reduction_factors[rev(order(reduction_factors$min)),]
  top.n = reduction_factors[1:min(nrow(reduction_factors), n_top), ]

  return(top.n)
}


#'@description Get best combinations
#'@param geno Genotype data frame.
#'@param rf Reduction factor data frame.
#'@param n_top Number if combinations to be returned.
#'@return Dataframe
#'@examples r = mmusprio("chr1", start=5000000, end=6000000, strain1="C57BL_6J", strain2="AKR_J")
#'plots = vis_reduction_factors(r$genotypes, r$reduction, 2)
#'plots[[1]]
#'plots[[2]]
#'@export
vis_reduction_factors = function(geno, rf, n_top){

  chr = geno$chr[1]
  start = min(geno$pos)
  end = max(geno$pos)
  strain1 = reduction_factors$strain1[1]
  strain2 = reduction_factors$strain2[1]
  top.n = get_top(reduction_factors, n_top)

  geno = geno[order(geno$pos),]
  geno$pos = 1:nrow(geno)

  plots = list()

  for(i in 1:nrow(top.n)){

    strain.comb = unlist(strsplit(top.n$combination[i], ","))
    geno.filtered = geno[(tolower(names(geno)) %in% c("chr", "pos", tolower(c(strain1, strain2, strain.comb))))]

    geno.melt = reshape2::melt(geno.filtered, id.vars=c("chr", "pos"), variable.name="strain", value.name="allele")

    geno.melt$strain = ordered(geno.melt$strain, levels=rev(c(strain1, strain2, strain.comb)))

    p = ggplot2::ggplot(data=geno.melt, ggplot2::aes(x=pos, y=strain, fill=factor(allele))) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values=c("#3b5998", "#e9ebee","#f6f7f9"), guide = ggplot2::guide_legend(title = "Allele")) +
      ggplot2::theme_bw() +
      ggplot2::labs(y="Strain", x="Position", title=paste0("Region: chr", chr, ":", scales::comma(start), "-", scales::comma(end), " (GRCm38)"),
           subtitle=paste0("Additional strains: ", paste0(strain.comb, collapse=", "), " Reduction factor: ", top.n$min[i])) +
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





# if(i == 1 && max_set_size > 1){
#
#   strain2min_reduction = hash::hash(keys=unlist(combs), values=red$min)
#
#   best_strain = unlist(combs[red$min == max(red$min),])
#   message(best_strain)
#
#   bad_red_strains = combs[red$min < min_tot_benef - strain2min_reduction[[best_strain]],]
#   bad_red_strains = bad_red_strains[bad_red_strains != best_strain]
#
#   message(paste0(length(bad_red_strains), " bad strains"))
#
#   other_combs = as.data.frame(gtools::combinations(length(bad_red_strains), max_set_size-i, bad_red_strains))
#
#   message(paste0(nrow(other_combs), " combs"))
#
#   r = t(apply(other_combs, 1, function(x){
#     sapply(x, function(y)  strain2min_reduction[[y]])
#   }))
#   if(max_set_size-i > 1) r = rowSums(r)
#
#   other_combs.remain = other_combs[r >= min_tot_benef - strain2min_reduction[[best_strain]],]
#
#   bad_red_strains.remain = unique(unlist(other_combs.remain))
#
#   message(paste0(length(bad_red_strains.remain), " bad strains remain"))
# }
