
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DabomAsotinSthd

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KevinSee/DabomAsotinSthd/master?urlpath=rstudio)

This repository contains the data and code for running the **D**am
**A**dult **B**ranch **O**ccupancy **M**odel
([DABOM](https://github.com/KevinSee/DABOM)) for adult steelhead
returning past the Asotin weir in the Asotin River. This model estimates
escapement past various locations in the Asotin basin using detections
of PIT tagged fish. The full methods are described in the paper:

> Waterhouse, L., White, J., See, K.E., Murdoch, A. R. and Semmens,
> B.X., (2020). *A Bayesian nested patch occupancy model to estimate
> steelhead movement and abundance*. Ecological Applications
> <https://doi.org/10.1002/eap.2202>

The most recent estimates for various spawn years can be found in the
“outgoing/estimates” folder that is available after cloning or
downloading this compendium. There is a manual describing the
step-by-step instructions to generate these results, found in the
“analysis/paper” folder.

The user can find more information related to installation and use of
the [PITcleanr](https://kevinsee.github.io/PITcleanr/) and
[DABOM](https://kevinsee.github.io/DABOM/) packages on their package
websites.

## Contents

The **analysis** directory contains:

-   [:file_folder: R_scripts](/analysis/R_scripts): various R scripts.
    The filenames are numbered to show what order they should be run in.
-   [:file_folder: paper](/analysis/paper): R Markdown source document
    for a manual. Includes code to reproduce the analysis for example
    dataset. It also has a rendered version, `manual.html`, suitable for
    reading (the code can be toggled on/off).
-   [:file_folder: data](/analysis/data): Data used in the analysis.
-   [:file_folder: figures](/analysis/figures): Plots and other
    illustrations
-   [:file_folder:
    supplementary-materials](/analysis/supplementary-materials):
    Supplementary materials including notes and other documents prepared
    and collected during the analysis.

## How to run in your browser or download and run locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the compendium as a zip from from this URL:
[master.zip](/archive/master.zip). After unzipping: - open the `.Rproj`
file in RStudio - run `devtools::install()` to ensure you have the
packages this analysis depends on (also listed in the
[DESCRIPTION](/DESCRIPTION) file). - finally, open
`analysis/paper/manual.Rmd` and knit to produce the `manual.html`, or
run `rmarkdown::render("analysis/paper/manual.Rmd")` in the R console

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.
