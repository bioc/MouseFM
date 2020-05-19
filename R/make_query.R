# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#   roxygen2::roxygenise()
#   BiocCheck::BiocCheck(/path/to/project)


source("R/meta.R")


#'Filter query builder
#'@param chr Vector of chromosome names.
#'@param start Optional vector of chromosomal start positions of target regions (GRCm38).
#'@param end Optional vector of chromosomal end positions of target regions (GRCm38).
#'@param consequence Optional vector of consequence types.
#'@param impact Optional vector of impact types.
#'@return Query string.
#'@keywords internal
filter_query = function(chr,
                        start = NULL,
                        end = NULL,
                        consequence = NULL,
                        impact = NULL) {
    stopifnot(!missing(chr))
    stopifnot(is.null(start) ||
                  (is.numeric(start) &&
                       start %% 1 == 0 && start > 0))
    stopifnot(is.null(end) ||
                  (is.numeric(end) && end %% 1 == 0 &&
                       end > 0))
    stopifnot(is.null(consequence) || is.vector(consequence))
    stopifnot(is.null(impact) || is.vector(impact))


    # Remove scientific notation in printing
    options(scipen = 999)


    # Build URL
    q = paste0('http://mmusserv.genehopper.de/rest/mmusfilter/',
               chr)


    if (is.numeric(start) && is.numeric(end)) {
        q = paste0(q, ":", start, "-", end)
    }


    pars = NULL
    if (!is.null(consequence) && is.vector(consequence)) {
        consequence = tolower(unique(consequence))
        all_consequences = avail_consequences()$consequence

        if (!all(unlist(lapply(consequence, function(x)
            x %in% all_consequences)))) {
            stop(paste0(
                "Accepted values for argument 'consequences' are: ",
                paste(all_consequences, collapse = ", ")
            ))
        } else {
            pars = c(pars,
                     paste0("consequence=",
                            consequence,
                            collapse = "&"))
        }
    }


    if (!is.null(impact) && is.vector(impact)) {
        impacts = toupper(unique(impact))
        all_impacts = unique(avail_consequences()$impact)

        if (!all(unlist(lapply(impact, function(x)
            x %in% all_impacts)))) {
            stop(paste0(
                "Accepted values for argument 'impact' are: ",
                paste0(all_impacts, collapse = ", ")
            ))
        } else {
            pars = c(pars,
                     paste0("impact=", impact, collapse = "&"))
        }
    }


    if (length(pars) > 0) {
        q = paste0(q, "?", paste0(pars, collapse = "&"))
    }

    return(q)
}


#'Finemap query builder
#'@param chr Vector of chromosome names.
#'@param start Optional vector of chromosomal start positions of target regions (GRCm38).
#'@param end Optional vector of chromosomal end positions of target regions (GRCm38).
#'@param strain1 First strain(s).
#'@param strain2 Second strain(s).
#'@param consequence Optional vector of consequence types.
#'@param impact Optional vector of impact types.
#'@param thr1 Number discordant strains in strain1. Between 0 and length(strain1)-1. 0 by default.
#'@param thr2 Number discordant strains in strain2. Between 0 and length(strain2)-1. 0 by default.
#'@return Query string.
#'@keywords internal
finemap_query = function(chr,
                         start = NULL,
                         end = NULL,
                         strain1,
                         strain2,
                         consequence = NULL,
                         impact = NULL,
                         thr1 = 0,
                         thr2 = 0) {
    stopifnot(!missing(chr))
    stopifnot(is.null(start) ||
                  (is.numeric(start) &&
                       start %% 1 == 0 && start > 0))
    stopifnot(is.null(end) ||
                  (is.numeric(end) && end %% 1 == 0 &&
                       end > 0))
    stopifnot(is.null(consequence) || is.vector(consequence))
    stopifnot(is.null(impact) || is.vector(impact))
    stopifnot(is.vector(strain1), is.vector(strain2))
    stopifnot(is.numeric(thr1) &&
                  thr1 >= 0 && thr1 < length(strain1))
    stopifnot(is.numeric(thr2) &&
                  thr2 >= 0 && thr2 < length(strain2))


    # Remove scientific notation in printing
    options(scipen = 999)


    # Build URL
    q = paste0('http://mmusserv.genehopper.de/rest/mmusfinemap/',
               chr)


    if (is.numeric(start) && is.numeric(end)) {
        q = paste0(q, ":", start, "-", end)
    }


    pars = NULL

    pars = c(paste0("thr1=", thr1), paste0("thr2=", thr2))


    all_strains = tolower(avail_strains()$id)


    if (!all(unlist(lapply(tolower(unique(c(
        strain1, strain2
    ))), function(x)
        x %in% all_strains)))) {
        stop(paste0(
            "Accepted strains are are:\n",
            paste(avail_strains()$id, collapse = "\n")
        ))
    }


    if (length(intersect(strain1, strain2)) > 0) {
        stop(paste0("strain1 and strain2 are overlapping"))
    }


    pars = c(pars, paste0("strain1=", strain1, collapse = "&"))
    pars = c(pars, paste0("strain2=", strain2, collapse = "&"))


    if (!is.null(consequence) && is.vector(consequence)) {
        consequence = tolower(unique(consequence))
        all_consequences = avail_consequences()$consequence

        if (!all(unlist(lapply(consequence, function(x)
            x %in% all_consequences)))) {
            stop(paste0(
                "Accepted values for argument 'consequences' are: ",
                paste(all_consequences, collapse = "\n")
            ))
        } else {
            pars = c(pars,
                     paste0("consequence=",
                            consequence,
                            collapse = "&"))
        }
    }


    if (!is.null(impact) && is.vector(impact)) {
        impacts = toupper(unique(impact))
        all_impacts = unique(avail_consequences()$impact)


        if (!all(unlist(lapply(impact, function(x)
            x %in% all_impacts)))) {
            stop(paste0(
                "Accepted values for argument 'impact' are: ",
                paste0(all_impacts, collapse = ", ")
            ))
        } else {
            pars = c(pars,
                     paste0("impact=", impact, collapse = "&"))
        }
    }


    if (length(pars) > 0) {
        q = paste0(q, "?", paste0(pars, collapse = "&"))
    }


    return(q)
}
