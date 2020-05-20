# MouseFM: In-silico methods for finemapping of genetic regions in inbred mice
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?hashtags=asd&url=https://www.biorxiv.org/content/...)

&nbsp;+ [Introduction](#Introduction)\
&nbsp;+ [Available inbred strains](#Available-inbred-strains)\
&nbsp;+ [Installation](#Installation)\
&nbsp;+ [Help pages](#Help-pages)\
&nbsp;+ [Authors](#Authors)\
&nbsp;+ [License](#License)


## Introduction
This **R** package provides methods for **genetic finemapping** in **inbred mice** by taking advantage of their **very high homozygosity rate** (>95%). 

For one ore more chromosomal regions (**GRCm38**), method `finemap` extracts homozygous single nucleotide variants (SNVs) for which the allele differs between two sets of strains (e.g. case vs controls) and outputs respective causal SNV/gene candidates.

Method `prio` allows to select strain combinations which best refine a specified genetic region. E.g. if a crossing experiment with two inbred mouse strains 'strain1' and 'strain2' resulted in a QTL, the outputted strain combinations can be used to refine the respective region in further crossing experiments for selecting causal SNP/gene candidates.

Method `fetch` allows to fetch genotypes for a specific region of interest.


## Available inbred strains
| No | Strain             | No | Strain    | No | Strain      |
|----|--------------------|----|-----------|----|-------------|
| 1  | 129P2/OlaHsd       | 14 | C57BR/cdJ | 26 | MOLF/EiJ    |
| 2  | 129S1/SvImJ        | 15 | C57L/J    | 27 | NOD/ShiLtJ  |
| 3  | 129S5/SvEvBrd      | 16 | C58/J     | 28 | NZB/B1NJ    |
| 4  | A/J                | 17 | CAST/EiJ  | 29 | NZO/HlLtJ   |
| 5  | AKR/J              | 18 | CBA/J     | 30 | NZW/LacJ    |
| 6  | BALB/cJ            | 19 | DBA/1J    | 31 | PWK/PhJ     |
| 7  | BTBR               | 20 | DBA/2J    | 32 | RF/J        |
| 8  | BUB/BnJ            | 21 | FVB/NJ    | 33 | SEA/GnJ     |
| 9  | C3H/HeH            | 22 | I/LnJ     | 34 | SPRET/EiJ   |
| 10 | C3H/HeJ            | 23 | KK/HiJ    | 35 | ST/bJ       |
| 11 | C57BL/10J          | 24 | LEWES/EiJ | 36 | WSB/EiJ     |
| 12 | C57BL/6            | 25 | LP/J      | 37 | ZALENDE/EiJ |
| 13 | C57BL/6NJ          |    |           |    |             |


## Installation
```R
devtools::install_github('matmu/MouseFM', build_vignettes = TRUE)
```

Loading package
```{r}
library(MouseFM)
```

**Please note**: A valid internet connection (HTTP port: 80) is required in order to install and use the package.


## Help pages
Multiple vignettes exist that guide you through general functionality of **MouseFM**.
```R
browseVignettes("MouseFM")
```

To see the help pages for each specific funtion:
```R
help(package="MouseFM")
```


## Authors
Matthias Munz [![](https://img.shields.io/twitter/follow/_MatthiasMunz?label=Follow&style=social)](https://img.shields.io/twitter/follow/_MatthiasMunz?label=Follow&style=social)\
University of Lübeck, Germany


## License
GNU General Public License v3.0
