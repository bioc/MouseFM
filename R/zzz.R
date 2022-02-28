.onLoad <- function(libname, pkgname) {
  options(url = "http://mousefm.genehopper.de/rest/finemap/")
  
  packageStartupMessage("
  ---------
  
  For example usage please run: vignette('MouseFM')
  
  Github Repo: https://github.com/matmu/MouseFM
  MouseFM Backend: https://github.com/matmu/MouseFM-Backend  

  Citation appreciated:
  Munz M, Khodaygani M, Aherrahrou Z, Busch H, Wohlers I (2021) In silico candidate variant and gene identification using inbred mouse strains. PeerJ. doi:10.7717/peerj.11017
                        
  ---------", domain = NULL, appendLF = TRUE)
}
