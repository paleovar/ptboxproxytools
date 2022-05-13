# PTBoxProxytools - Statistical tools for processing and analyzing paleo data time series

As part of the [PaleoToolBox (PTBox)](https://palmodapp.cloud.dkrz.de/), this package provides tools and helpers for processing and analyzing paleo data time series.

## Introduction

The package can be installed using `devtools` in [R](https://www.r-project.org/):

`require(devtools)` 

`devtools::install_git("https://github.com/paleovar/ptboxproxytools", build_vignettes = TRUE)`

`library(PTBoxProxytools)` 

If you'd like to contribute to the package, you might want to clone it (type `git clone https://github.com/paleovar/ptboxproxytools.git` in the terminal) and open the `.Rproj` in `RStudio`. Please don't hesitate to report issues, bugs, etc. to the authors.


After the installation run
```R
vignette('PTBoxProxytools_howto')
```
for a general introduction into the conventions of `PTBox`. Also check out the help pages with
```R
?PTBoxProxytools::PTBoxProxytools
```

## Responsible for this repository:
Developers: *[Nils Weitzel](https://github.com/nilsweitzel), [Beatrice Ellerhoff](https://github.com/bellerhoff) and [Moritz Adam](https://github.com/MoritzAdam)*

Contributors to the [PTBox](https://palmodapp.cloud.dkrz.de/) framework (alphabetically): *Moritz Adam, Jean-Philippe Baudouin, Oliver Bothe, Manuel Chevalier, Beatrice Ellerhoff, Patrizia Schoch, Kira Rehfeld, Nils Weitzel*

## Acknowledgements

We thank the [R Core team](https://www.R-project.org/) and the developers of all packages that `PTBoxProxytools` buids on, particularly the developers of the [`PaleoSpec`](https://github.com/EarthSystemDiagnostics/paleospec) and [`nest`](https://github.com/krehfeld/nest) package. Please see `citation()` for details on the R Core Team and `citation("packagename")` for details on the developers of individual packages.

The development of this package has been supported by the [PalMod](https://www.palmod.de/) project (subproject no. 01LP1926C). We thank all contributors and members of PalMod, the [Earth's climate and environmental dynamics](https://www.iup.uni-heidelberg.de/en/research/paleoclimate-dynamics) and [SPACY](https://uni-tuebingen.de/climatology/) group for discussion at different stages of this package. 

*The Authors, April 2022*
