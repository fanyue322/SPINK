#' @examples
#' \dontrun{
#' library(SPINK)
#' library(Seurat)
#' 
#' # Loading example data
#' data(object)
#' 
#' # 1. preprocess
#' object <- spink_preprocess(
#'   object = object,
#'   group.by = "domain",
#'   domain = "R3",
#'   refGenome = "hg38",
#'   num.core = 10
#' )
#' 
#' # 2. analysis
#' object <- spink_analysis(
#'   object = object)
#' 
#' # 3. Check Results
#' results <- GetLinkResult(object)
#' head(results)
#' }

NULL