#' @title
#' ds.metaFor
#'
#' @description
#' TODO
#'
#' @export
#'

ds.metaFor <- function(ests=NULL, SEs=NULL){
library(metagen)
library(metafor)
rma.REML<-rma(yi=ests, sei=SEs, method="REML")
rma.ML<-rma(yi=ests, sei=SEs, method="ML")
rma.FE<-rma(yi=ests, sei=SEs, method="FE")

return(list(rma.FE=rma.FE,rma.ML=rma.ML,rma.REML=rma.REML))

}
