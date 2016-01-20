#' @title
#' ds.glm.REMA
#'
#' @description
#' TODO
#'
#' @export
#'

ds.glm.REMA <- function(formula=NULL, data=NULL, family=NULL, offset=NULL, weights=NULL,
               checks=FALSE, maxit=15, CI=0.95, viewIter=FALSE, datasources=NULL) {
  
 # if no opal login details are provided look for 'opal' objects in the environment
  if(is.null(datasources)){
    datasources <- findLoginObjects()
  }
  
  # verify that 'formula' was set
  if(is.null(formula)){
    stop(" Please provide a valid regression formula!", call.=FALSE)
  }

# check if user gave offset or weights directly in formula, if so the argument 'offset' or 'weights'
# to provide name of offset or weights variable
    if(sum(as.numeric(grepl('offset', formula, ignore.case=TRUE)))>0 ||
       sum(as.numeric(grepl('weights', formula, ignore.case=TRUE)))>0)
	{
       cat("\n\n WARNING: you may have specified an offset or regression weights")
       cat("\n as part of the model formula. In ds.glm (unlike the usual glm in R)")
       cat("\n you must specify an offset or weights separately from the formula")
       cat("\n using the offset or weights argument.\n\n")
	}

  formula <- as.formula(formula)

  
  # check that 'family' was set
  if(is.null(family)){
    stop(" Please provide a valid 'family' argument!", call.=FALSE)
  }
  
  # if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if(!(is.null(data))){
    defined <- isDefined(datasources, data)
  }
  
  # beginning of optional checks - the process stops if any of these checks fails #
  if(checks){
    message(" -- Verifying the variables in the model")
    # call the function that checks the variables in the formula are defined (exist) on the server site and are not missing at complete
 
   glmChecks(formula, data, offset, weights, datasources)
  }else{
    #message("WARNING:'checks' is set to FALSE; variables in the model are not checked and error messages may not be intelligible!")
  }

#MOVE ITERATION COUNT BEFORE ASSIGNMENT OF beta.vect.next  
#Iterations need to be counted. Start off with the count at 0
  #and increment by 1 at each new iteration
  iteration.count<-0

  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
 

#ARBITRARY LENGTH FOR START BETAs AT THIS STAGE BUT IN LEGAL TRANSMISSION FORMAT ("0,0,0,0,0")
 beta.vect.next <- rep(0,5)
 beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  

#IDENTIFY THE CORRECT DIMENSION FOR START BETAs VIA CALLING FIRST COMPONENT OF glmDS
 
   cally1 <- call('glmDS1.5s', formula, family, weights, data)
   
   study.summary <- datashield.aggregate(datasources, cally1)
#  num.par.glm<-study.summary$study1$dimX[2]
   num.par.glm<-study.summary[[1]][[1]][[2]]

coef.names<-study.summary[[1]][[2]]

y.invalid<-NULL
Xpar.invalid<-NULL
w.invalid<-NULL

	for(ss in 1:numstudies)
	{
   	y.invalid<-c(y.invalid,study.summary[[ss]][[3]])
	Xpar.invalid<-rbind(Xpar.invalid,study.summary[[ss]][[4]])
   	w.invalid<-c(w.invalid,study.summary[[ss]][[5]])
	}

y.invalid<-as.matrix(y.invalid)
sum.y.invalid<-sum(y.invalid)
dimnames(y.invalid)<-list(names(datasources),"y vector")

Xpar.invalid<-as.matrix(Xpar.invalid)
sum.Xpar.invalid<-sum(Xpar.invalid)
dimnames(Xpar.invalid)<-list(names(datasources),coef.names)


w.invalid<-as.matrix(w.invalid)
sum.w.invalid<-sum(w.invalid)
dimnames(w.invalid)<-list(names(datasources),"Weights")

output.blocked.information<-"OUTPUT BLOCKED: ANY VALUE OF 1 DENOTES A POTENTIAL DISCLOSURE RISK"

if(sum.y.invalid>0||sum.Xpar.invalid>0||sum.w.invalid>0){
    return(list(output.blocked.information,
		    y.data=y.invalid,
                X.data=Xpar.invalid,
                weights.data=w.invalid))
  }
########################################
   
#ADD
   study.summary.0<-study.summary   

   beta.vect.next <- rep(0,num.par.glm)
   beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")

   #NOW CALL SECOND COMPONENT OF glmDS TO GENERATE SCORE VECTORS AND INFORMATION MATRICES
    calltext <- call('glmREMADS', formula, family, beta.vect=beta.vect.temp, offset, weights, data)

    study.summary <- datashield.aggregate(datasources, calltext)

#####MODIFY TO REFLECT INTEGRATION IS NOW IN SERVER-SIDE FUNCTION

#####CHECK THIS CODE IS OK IF FUNCTION WONT WORK
  

#IF NEED TO INTEGRATE RETURNED OUTPUT
      .select <- function(l, field) {
      lapply(l, function(obj) {obj[[field]]})
    }

    disclosure.risk.total<-Reduce(f="+", .select(study.summary, 'disclosure.risk'))

    if(disclosure.risk.total>0){
	message("DISCLOSURE RISK IN y.vect, X.mat or w.vect IN AT LEAST ONE STUDY")
	output.blocked.information.1<-"POTENTIAL DISCLOSURE RISK AND CLIENTSIDE FUNCTION MODIFIED BY USER"
	output.blocked.information.2<-"SO SCORE VECTOR AND INFORMATION MATRIX DESTROYED IN AT-RISK STUDY"

   	return(list(output.blocked.information.1,
			output.blocked.information.2))
	}
   
    Nvalid.total<-Reduce(f="+", .select(study.summary, 'Nvalid'))
    Nmissing.total<-Reduce(f="+", .select(study.summary, 'Nmissing'))
    Ntotal.total<-Reduce(f="+", .select(study.summary, 'Ntotal'))
 
    

      nsubs.total<-Reduce(f="+", .select(study.summary, 'numsubs'))
      f <- study.summary[[1]]$family
	formulatext<-study.summary[[1]]$formulatext
   

#    if(viewIter){
#      #For ALL iterations summarise model state after current iteration
#      message("SUMMARY OF MODEL STATE after iteration ", iteration.count)
#      message("Current deviance ", dev.total," on ",(nsubs.total-length(beta.vect.next)), " degrees of freedom")
#      message("Convergence criterion ",converge.state," (", converge.value,")")
#      
#      message("\nbeta: ", paste(as.vector(beta.vect.next), collapse=" "))
#      
#      message("\nInformation matrix overall:")
#      message(paste(capture.output(info.matrix.total), collapse="\n"))
#      
#      message("\nScore vector overall:")
#      message(paste(capture.output(score.vect.total), collapse="\n"))
#      
#      message("\nCurrent deviance: ", dev.total, "\n")
#    }
#  }
#  if(!viewIter){
#    #For ALL iterations summarise model state after current iteration
#    message("SUMMARY OF MODEL STATE after iteration ", iteration.count)
#    message("Current deviance ", dev.total," on ",(nsubs.total-length(beta.vect.next)), " degrees of freedom")
#    message("Convergence criterion ",converge.state," (", converge.value,")")
#    
#    message("\nbeta: ", paste(as.vector(beta.vect.next), collapse=" "))
#    
#    message("\nInformation matrix overall:")
#    message(paste(capture.output(info.matrix.total), collapse="\n"))
#    
#    message("\nScore vector overall:")
#    message(paste(capture.output(score.vect.total), collapse="\n"))
#    
#    message("\nCurrent deviance: ", dev.total, "\n")
#  }
  
  #If convergence has been obtained, declare final (maximum likelihood) beta vector,
  #and calculate the corresponding standard errors, z scores and p values
  #(the latter two to be consistent with the output of a standard GLM analysis)
  #Then print out final model summary    

study.specific.output<-study.summary

      glm.REMA.output <- list(
			Nvalid=Nvalid.total,
			Nmissing=Nmissing.total,
			Ntotal=Ntotal.total,
			study.specific.output=study.specific.output,
     			formula=formulatext,
	  		family=f
    			)

				return(glm.REMA.output)
  
}

