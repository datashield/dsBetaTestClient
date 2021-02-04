#'
#' @title Computes server-side roc curve
#' @description This function computes the roc curve
#'  of a given server-side pair of prediction and reference vector. 
#'
#' @details  This function is similar to the function \code{roc} in pROC.
#' 
#' Server function called: \code{rocDS}
#' @param reference a logical (boolean) vector refering to the "truth" on which case is positive/true or negative/false.
#' @param prediction a numeric vector (values between 0 and 1) holding the tendency to lean towards a positive or negative outcome.
#' @param checks logical. If TRUE  optional checks of model
#' components will be undertaken. Default is FALSE to save time. 
#' It is suggested that checks
#' should only be undertaken once the function call has failed. 
#' @param datasources a list of \code{\link{DSConnection-class}} 
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return \code{ds.roc} returns to the client-side plotting a roc curve.
#' 
#' @author Matthis A
#' @export
#' @examples
#' \dontrun{
#'
#'  ## Version 6, for version 5 see the Wiki
#'   
#'   # connecting to the Opal servers
#' 
#'   require('DSI')
#'   require('DSOpal')
#'   require('dsBaseClient')
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1", 
#'                  url = "http://192.168.56.100:8080/", 
#'                  user = "administrator", password = "datashield_test&", 
#'                  table = "CNSIM.CNSIM1", driver = "OpalDriver")
#'   builder$append(server = "study2", 
#'                  url = "http://192.168.56.100:8080/", 
#'                  user = "administrator", password = "datashield_test&", 
#'                  table = "CNSIM.CNSIM2", driver = "OpalDriver")
#'   builder$append(server = "study3",
#'                  url = "http://192.168.56.100:8080/", 
#'                  user = "administrator", password = "datashield_test&", 
#'                  table = "CNSIM.CNSIM3", driver = "OpalDriver")
#'   logindata <- builder$build()
#'   
#'   connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D") 
#'   
#'   #Calculate the roc curve of prediction and reference vectors in the server-side
#'   
#'   ds.roc(reference = "D$muac",
#'           prediction = "D$sex",
#'           checks = FALSE,
#'           datasources = connections)
#'              
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(connections)
#' }
#'
ds.roc <- function(reference=NULL, prediction=NULL, checks=FALSE, datasources=NULL){

#####################################################################################
#MODULE 1: IDENTIFY DEFAULT DS CONNECTIONS                                          #
  # look for DS connections                                                         #
  if(is.null(datasources)){                                                         #
    datasources <- datashield.connections_find()                                    #
  }                                                                                 #
#####################################################################################

Nstudies <- length(datasources)

checkarg(reference, "logical", "testing", checks, datasources)
checkarg(prediction, "numeric", "prediction", checks, datasources)

cally <- paste0("rocDS(", reference, ",", prediction, ")")

ss.obj <- DSI::datashield.aggregate(datasources, as.symbol(cally))

  

  ret = ss.obj[[1]]
  plot(1 - ret$Specificities, ret$Sensitivities, type="l")
  points(ret$Threshold, col="red")
  for(j in 2:Nstudies){
    ret <- ss.obj[[j]]
    lines(1 - ret$Specificities, ret$Sensitivities)
    points(ret$Threshold, col="red")
  }
  return(ss.obj)
}
#ds.roc

checkarg <- function(x, validtype, usage, checks, datasources){
#####################################################################################
#MODULE 2: SET UP KEY VARIABLES ALLOWING FOR DIFFERENT INPUT FORMATS                #
  if(is.null(x)){                                                                   #
    stop(paste0("Please provide the name of the input vector for ", usage, "!"), call.=FALSE)               #
  }                                                                                 #
  # the input variable might be given as a variable in a data frame (i.e. D$x)      #
  # or just as a vector not attached to a table (i.e. x)                            #
  # we have to make sure the function deals with each case                          #
  xnames <- extract(x)                                                              #
  varname <- xnames$elements                                                        #
  obj2lookfor <- xnames$holders                                                     #
#####################################################################################

###############################################################################################
#MODULE 3: GENERIC OPTIONAL CHECKS TO ENSURE CONSISTENT STRUCTURE OF KEY VARIABLES            #
#IN DIFFERENT SOURCES                                                                         #
  # beginning of optional checks - the process stops and reports as soon as one               #
  #check fails                                                                                #
                                                                                              #
  if(checks){                                                                                 #
    message(paste0(" -- Verifying the variables for ", usage, " in the model"))                                       #
                                                                                              #
    # check if the input object(s) is(are) defined in all the studies                           #
    if(is.na(obj2lookfor)){                                                                     #
      defined <- isDefined(datasources, varname)                                                #
    }else{                                                                                      #
      defined <- isDefined(datasources, obj2lookfor)                                            #
    }                                                                                           #
                                                                                                #
    # call the internal function that checks the input object is suitable in all studies        #
    varClass <- checkClass(datasources, x)                                                      #
    # the input object must be a numeric or an integer vector                                   #
    if(!(validtype %in% varClass)){                                          #
      stop(paste0("The input vector for ", usage, " must be of type ", validtype, "."), call.=FALSE)             #
    }                                                                                           #
  }                                                                                             #
###############################################################################################
}
