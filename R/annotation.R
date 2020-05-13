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


#'Annotate consequences
#'@description Request variant consequences from Variant Effect Predictor (VEP) via Ensembl Rest Service.
#'@param snps Data frame including columns rsid, ref, alt.
#'@param species Species name, e.g. mouse (GRCm38) or human (GRCh38).
#'@return Data frame.
#'@examples geno = mmusfinemap("chr1", start=5000000, end=6000000,
#'strain1=c("C57BL_6J"), strain2=c("AKR_J", "A_J", "BALB_cJ"))
#'
#'df = annotate_consequences(geno, "mouse")
#'@export
annotate_consequences = function(snps, species){

  df.split = df_split(snps[!is.na(snps$rsid),], 200)

  res = lapply(df.split, function(x) {message(paste0("Query ", nrow(x), " variants..."))
               ensembl_rest_vep(x, "mouse")})

  return(data.table::rbindlist(res))
}


#'Request variant consequences from Variant Effect Predictor (VEP) via Ensembl Rest Service
#'@param snps Data frame including columns rsid, ref, alt.
#'@param species Species name, e.g. mouse or human.
#'@return Data frame.
#'@keywords internal
ensembl_rest_vep = function(snps, species){


  # Request
  server = "https://rest.ensembl.org"
  ext = paste0("/vep/", species, "/id")
  r = httr::POST(paste(server, ext, sep = ""), httr::content_type("application/json"), httr::accept("application/json"),
                 body = paste0('{ "ids" : ["', paste0(snps$rsid[1:200], collapse='","'), '" ] }'))


  httr::stop_for_status(r)


  res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))


  # Extract consequences
  cons.list = lapply(1:length(res$input), function(i){res$transcript_consequences[[i]]$snp = res$input[[i]]
  res$transcript_consequences[[i]]$consequence_terms = lapply(1:length(res$transcript_consequences[[i]]$consequence_terms),
                                                              function(j) paste0(res$transcript_consequences[[i]]$consequence_terms[[j]], collapse=","))
  return(res$transcript_consequences[[i]])})


  # Reformat
  final = data.table::rbindlist(cons.list, fill=TRUE)
  final[final == "NULL"] = NA
  cols = c("snp", "variant_allele", "consequence_terms", "impact", "transcript_id", "gene_id", "gene_symbol", "gene_symbol_source", "strand", "distance", "cdna_end", "cdna_start", "codons", "cds_end", "amino_acids", "protein_start", "cds_start", "protein_end", "sift_score", "sift_prediction")
  cols.missing = cols[! cols %in% colnames(final)]
  if(length(cols.missing) > 0) final[, cols.missing] = NA
  final = subset(final, select=cols)

  final = as.data.frame(matrix(unlist(final), ncol=ncol(final)))
  colnames(final) = cols


  # Filter for alleles
  final = final[final$snp %in% snps$rsid & (final$variant_allele %in% snps$ref | final$variant_allele %in% snps$alt),]


  return(final)
}


#'Splits data frame df into subsets with maximum n rows
#'@param df Data frame.
#'@param n Max number of rows per subset.
#'@return List of data frames.
#'@keywords internal
df_split = function(df, n) {

  l = nrow(df)
  s = ceiling(l/n)

  return(lapply(seq_len(s), function(i){
    to = min(l, i*200)
    from = i*200-n+1
    return(df[from:to,])
  }))
}
