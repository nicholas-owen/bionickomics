
#' Functions to support installation and loading of CRAN and Bioconductor packages in R
#'
#' @param x a package_name
#'
#' @return
#' @export
#'
#' @examples
#'

install_package <- function(x){
  options(repos = list(CRAN="http://cran.rstudio.com/"))
  if (! isTRUE(x %in% .packages(all.available = TRUE)) && any(available.packages()[, 1] == x)) {
    # update.packages(ask=F) # update dependencies, if any.
    eval(parse(text = paste("install.packages('", x, "')", sep = "")))
  }

  ## CHECK IF ON BIOCONDUCTOR
  if (! isTRUE(x %in% .packages(all.available = TRUE))) {
    bcPackages <- BiocManager::available()
    if (any(bcPackages == x)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      eval(parse(text = paste("BiocManager::install('", x, "', ask=FALSE)", sep = "")))
    }
  }
}


#' Functions to support installation and loading of CRAN and Bioconductor packages in R
#'
#' @param x a package_name
#'
#' @return
#' @export
#'
#' @examples
#'


load_pack <- function(x, warn_conflicts=T){
  x <- as.character(substitute(x));
  install_package(x)

  ## load it using a library function so that load_pack errors if package is still not ins
  eval(parse(text = paste("base::library(", x, ",  quietly=T, warn.conflicts=", warn_conflicts, ")", sep = "")))
}


#' Functions to support installation and loading of CRAN and Bioconductor packages in R
#'
#' @param x a package_name
#'
#' @return
#' @export
#'
#' @examples
#'

check_version <- function(pkg_name, min_version) {
  cur_version <- packageVersion(pkg_name)
  if (cur_version < min_version) stop(sprintf("Package %s needs a newer version,
               found %s, need at least %s", pkg_name, cur_version, min_version))
}


##EXAMPLE : load_pack(stringr)
