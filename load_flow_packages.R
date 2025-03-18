## -----
## <<Required_Packages>>
## check for Bioconductor and install if not available
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")


## -----
## load packages or install if not available
## have to split these out by bioconductor vs. non-bioconductor
## non-bioconductor
package_installer <- function(x){
    if(!requireNamespace(x, quietly = TRUE))
        install.packages(x, dependencies = TRUE,
                         repos = "http://cran.wustl.edu/",
                         quiet = TRUE, INSTALL_opts = '--no-lock')}
packages <- c("colorspace", "lattice",
              "ggvis", "dygraphs",
              "DescTools", "viridis",
              "latticeExtra", "session",
              "grid", "lmerTest", "lme4",
              "RColorBrewer")
sapply(X = packages, FUN = package_installer)
sapply(X = packages, FUN = require, character.only = TRUE)
sapply(X = packages, FUN = library, character.only = TRUE)


## -----
## bioconductor
bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x, INSTALL_opts = '--no-lock')}
bioc_packages <-  c("flowCore", "flowViz",
                    "flowUtils", "flowFP",
                    "geneplotter", "Rgraphviz", "ggcyto")
sapply(X = bioc_packages, FUN = bioc_package_installer)
sapply(X = bioc_packages, FUN = require, character.only = TRUE)


## -----
## required for merging flowsets into a single flowframe
source(file = "https://raw.githubusercontent.com/mac230/flow_scripts/master/set2frame.R")
source(file = "~/emacs/R/functions/my_box_stats.R")
source("~/emacs/R/functions/color_setup_BY_RM.R")
