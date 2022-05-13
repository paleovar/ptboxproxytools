#' @keywords internal
"_PACKAGE"


.onLoad <- function(libname, pkgname) {
    .PTBoxProxytools <<- new.env(parent = emptyenv())

    # check if PTBoxProxydata is installed (in particular for Proxyzoo support and zoo_apply)
    # first case for normal installation, second case for development with `devtools::load_all('/path/to/PtBoxProxydata')`
    if (!any("PTBoxProxydata" %in% installed.packages(), exists(".PTBoxProxydata"))) {
        stop("PTBoxProxydata not properly installed\n see github.com/paleovar/PTBoxProxydata for instructions")
    }
}
