# In-silico methods for strain prioritization and finemapping of genetic regions in mice.
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?hashtags=asd&url=https://www.biorxiv.org/content/...)

&nbsp;+ [Introduction](#Introduction)\
&nbsp;+ [Installation](#Installation)\
&nbsp;+ [Usage](#Usage)\
&nbsp;&nbsp;&nbsp;&nbsp;|-- [Accepted query terms](#Accepted-query-terms)\
&nbsp;&nbsp;&nbsp;&nbsp;|-- [Optional parameters](#Optional-parameters)\
&nbsp;&nbsp;&nbsp;&nbsp;|-- [Meta information](#Meta-information)\
&nbsp;+ [Try out online](#Try-out-online)\
&nbsp;+ [Authors](#Authors)\
&nbsp;+ [Citation](#Citation)\
&nbsp;+ [License](#License)


## Introduction
This **R** package provides methods for **finemapping of inbred mice* by taking advantage of their **very high homozygosity rate** (>95% ). For one ore more chromosomal regions (GRCm38), the method "" extracts homozygous SNPs for which the allele differs between two sets of strains (e.g. case vs controls) and respective causal SNP/gene candidates.

Method "" allows to select additional strains which best resolve a QTL found by a crossing experiment of two inbred mouse strains. These additional strains can then be used to refine the QTL in further crossing experiments.


## Installation
```R
devtools::install_github('matmu/mmus', build_vignettes = TRUE)
```

**Please note**: A valid internet connection (HTTP port: 80) is required in order to install and use the package.

## Help pages
```R
browseVignettes("mmus")
help(package="mmus")
```

## Loading package
```{r}
library(mmus)
```

## Usage


### Meta information
Column descriptions of the received data frame can be accessed by calling:

```R
df = 
comment(df)
```

## Try out online
If you want to try out the R package online, there is an example **Google Colaboratory** project at

https://colab.research.google.com/drive/1i1sjQHCjaw2wYzVBnXQ9iaghnk-jSU95#scrollTo=5Hi6sCe7SPFb

To run the project, make a private copy or open the project in playground mode and sign in to Google. 


## Authors
Matthias Munz [![](https://img.shields.io/twitter/follow/_MatthiasMunz?label=Follow&style=social)](https://img.shields.io/twitter/follow/_MatthiasMunz?label=Follow&style=social)\
University of Lübeck, Germany


## Citation
Please cite the following article when using our tool:


## License
GNU General Public License v3.0


