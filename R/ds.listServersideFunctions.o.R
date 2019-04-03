#' @title ds.listServersideFunctions.o calling no server-side functions
#' @description Lists all current server-side functions
#' @details Uses dsadmin.get_methods function from opaladmin package to list all
#' assign and aggregate functions on the available opal servers. The only choice of
#' arguments is in datasources; i.e. which studies to interrogate. Once the studies have 
#' been selected ds.listServersideFunctions.o lists all assign functions for all 
#' of these studies and then all aggregate functions for all of them.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. If the <datasources>
#' argument is not specified the default set of connections will be used: see \link{datashield.connections_default}.
#' @return list containing all serverside functions by study. Firstly lists assign
#' and then aggregate functions.
#' @author Burton, PR.
#' @export
#' @import DSI
ds.listServersideFunctions.o<-function(datasources=NULL){
  .Deprecated("DSI::datashield.methods")
  
  assign.funs <- datashield.methods(datasources, 'assign')
  aggregate.funs <- datashield.methods(datasources, 'aggregate')
  rbind(assign.funs, aggregate.funs)
}

#ds.listServersideFunctions.o
