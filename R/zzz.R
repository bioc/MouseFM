.onLoad <- function(libname, pkgname) {
  
  packageStartupMessage("
  ---------
  For example usage please run: vignette('MouseFM')
  
  Support me: http://matthiasmunz.de/support_me

  Citation appreciated:
  Munz M, Khodaygani M, Aherrahrou Z, Busch H, Wohlers I (2021) In silico candidate variant and gene identification using inbred mouse strains. PeerJ. doi:10.7717/peerj.11017
  
  Github Repo: https://github.com/matmu/MouseFM
  MouseFM Backend: https://github.com/matmu/MouseFM-Backend
  ---------", domain = NULL, appendLF = TRUE)
}
