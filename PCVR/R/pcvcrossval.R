#' Generate sequence of indices for cross-validation
#'
#' @description
#' Generates and returns sequence of object indices for each segment in random segmented
#' cross-validation
#'
#' @param cv
#' cross-validation settings, can be a number, a list or a vector with integers.
#' @param nobj
#' number of objects in a dataset
#' @param resp
#' vector or matrix with response values to use in case of venetian blinds
#'
#' @return
#' vector with object indices for each segment
#'
#' @details
#' Parameter `cv` defines how to split the rows of the training set. The split is similar
#' to cross-validation splits, as PCV is based on cross-validation. This parameter can have
#' the following values:
#'
#' 1. A number (e.g. `cv = 4`). In this case this number specifies number of segments for random
#' splits, except `cv = 1` which is a special case for leave-one-out (full cross-validation).
#'
#' 2. A list with 2 values: `list("name", nseg)`. In this case `"name"` defines the way to make
#' the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random
#' splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a
#' number of segments for splitting the rows into. For example, `cv = list("ven", 4)`, shown in
#' the code examples above, tells PCV to use Venetian blinds splits with 4 segments.
#'
#' 3. A vector with integer numbers, e.g. `cv = c(1, 2, 3, 1, 2, 3, 1, 2, 3)`. In this case number
#' of values in this vector must be the same as number of rows in the training set. The values
#' specify which segment a particular row will belong to. In case of the example shown here, it
#' is assumed that you have 9 rows in the calibration set, which will be split into 3 segments.
#' The first segment will consist of measurements from rows 1, 4 and 7.
#'
#' @export
pcvcrossval <- function(cv = 1, nobj = NULL, resp = NULL) {

   # get cross-validation parameters
   if (is.null(nobj)) nobj <- length(resp)

   # if user already provided matrix with values - return it
   if (is.numeric(cv) && length(cv) == nobj) return(as.matrix(cv))

   p <- getCrossvalParams(cv = cv, nobj = nobj)
   if (!(p$type %in% c("rand", "ven", "loo"))) {
      stop("Wrong name for cross-validation method.")
   }

   # check number of repetitions
   if (p$nrep < 1 || p$nrep > 100) {
      stop("Wrong value for cv repetitions (should be between 1 and 100).")
   }

   # check number of segments
   if (p$nseg < 2 || p$nseg > nobj) {
      stop("Wrong value for number of segments (should be between 2 and number of objects).")
   }

   if (p$type == "loo") {
      return(matrix(seq_len(nobj), ncol = 1))
   }

   if (p$type == "rand") {
      return(sapply(seq_len(p$nrep), function(i) rep(seq_len(p$nseg), length.out = nobj)[sample(nobj)]))
   }

   if (p$type == "ven") {
      ind <- if (is.null(resp)) seq_len(nobj) else order(resp)
      return(matrix(rep(seq_len(p$nseg), length.out = nobj)[ind], ncol = 1))
   }

   stop("Something went wrong.")
}


#' Returns parameters for cross-validation based on 'cv' value
#'
#' @param cv
#' settings for cross-validation provided by user
#' @param nobj
#' number of objects in calibration set
#'
getCrossvalParams <- function(cv, nobj) {

   nrep <- 1

   # random
   if (is.numeric(cv)) {
      return(
         list(
            type = "rand",
            nrep = 1,
            nseg = if (cv == 1) nobj else cv
         )
      )
   }

   # leave one out
   type <- cv[[1]]
   if (type == "loo") {
      return(
         list(
            type = "loo",
            nrep = 1,
            nseg = nobj
         )
      )
   }

   # venetian blinds
   nseg <- cv[[2]]
   if (type == "ven") {
      return(
         list(
            type = "ven",
            nrep = nrep,
            nseg = nseg
         )
      )
   }

   nrep <- if (length(cv) == 3) cv[[3]] else 1
   return(
      list(
         type = type,
         nrep = nrep,
         nseg = nseg
      )
   )
}
