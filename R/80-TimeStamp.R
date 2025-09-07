#' @title Generate Timestamp String
#'
#' @description
#' Creates a formatted character string representing the current system time.
#' The format is "YYYY/mm/DD HH:MM:SS" (year/month/day hour:minute:second).
#'
#' @return Character string with current time in "YYYY/mm/DD HH:MM:SS" format.
#'         If system time is unavailable, returns a fixed timestamp.
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
    time <- tryCatch(
        Sys.time(),
        error = function(e) as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
    )
    format(time, "%Y/%m/%d %H:%M:%S")
}
