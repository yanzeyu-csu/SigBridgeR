#' @title Generate Timestamp String
#'
#' @description
#' Creates a formatted character string representing the current system time.
#' The format is "YYYY/MM/DD HH:MM:SS" (year/month/day hour:minute:second).
#'
#' @return Character string with current time in "YYYY/MM/DD HH:MM:SS" format.
#'
#' @examples
#' \dontrun{
#' # Current time as formatted string
#' TimeStamp()
#' # Returns something like: "2025/06/15 16:04:00"
#' }
#'
#' @keywords internal
#'
TimeStamp <- function() {
    format(Sys.time(), "%Y/%m/%d %H:%M:%S")
}
