#' Generate Timestamp String
#'
#' @description
#' Creates a formatted character string representing the current system time.
#' The format is "YYYY/MM/DD HH:MM:SS" (year/month/day hour:minute:second).
#'
#' @return Character string with current time in "YYYY/MM/DD HH:MM:SS" format.
#'
#' @examples
#' # Current time as formatted string
#' TimeStamp()
#' # Returns something like: "2025/08/19 10:30:45"
#'
#' @keywords SigBridgeR_internal
#'
TimeStamp <- function() {
    format(Sys.time(), "%Y/%m/%d %H:%M:%S")
}
