#' document deprecated functions
#'
#' @title deprecated functions of mesheR
#' @name deprecated
#' @rdname mesheR-deprecated
#' @keywords internal
NULL

#' @rdname mesheR-deprecated
#' @export
plotDisplacementField <- function (...)
{
  .Deprecated("plotDisplacementField", package="mesheR")
  plot(...)
}
