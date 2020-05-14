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


#'Send HTTP request to MMUS Server
#'@param q Query string
#'@param n.tries Number of tries
#'@param method HTTP method to use
#'@return Data frame.
#'@keywords internal
#genehopper_request("http://mmusserv.genehopper.de/rest/mmusprio/chr1:5000000-6000000")
genehopper_request = function(q, n.tries=2, method="GET"){

  stopifnot(is.numeric(n.tries))
  stopifnot(method == "GET" || method == "POST")
  stopifnot(nchar(q) > 0)

  # Send HTTP request and retrieve response
  while(n.tries > 0){

    if(method == "GET"){
      response = httr::GET(q)
    }
    else if(method == "POST"){
      response = httr::POST(q)
    }

    if(!httr::http_error(response)){
      break
    }

    Sys.sleep(3)
    n.tries = n.tries-1
  }


  if(n.tries == 0){
    stop("Web server seems to be down. Try again later!")
  }


  result = httr::content(response)
  result = unlist(strsplit(result , "\n"))


  # Display error message from server
  if(length(grep("^#", result)) == 0){
    warning(result)
    return()
  }


  # Extract
  meta = grep("^#", result, value = TRUE)
  data = grep("^[^#]", result, value = TRUE)


  # If no results found
  if(length(data)-1 <= 0) {
    message("No results found")
    return()
  }


  # Convert to dataframe
  header = unlist(strsplit(data[1], "\t"))
  ncols = length(header)
  data = unlist(strsplit(data[-1], "\t"))

  m = matrix(data, ncol=ncols, byrow=TRUE)
  d = as.data.frame(m, stringsAsFactors=FALSE)
  colnames(d) = header
  comment(d) = meta


  return(d)
}

