
#' @title ds.lmerSLMA.o calling lmerDS2.o
#' @description Fits a linear mixed effects model (lme) on data from a single or multiple sources
#' @details  Fits a linear mixed effects model (lme) on data from a single source or from multiple sources.
#' In the latter case, the lme is fitted to convergence in each data source and the
#' estimates and standard errors
#' returned from each study separately. When these are then pooled using a function such as
#' ds.metafor, this is a form of study level meta-analysis (SLMA). The SLMA approach offers some advantages
#'  when there is marked heterogeneity
#' between sources that cannot simply be corrected with fixed effects each reflecting a study
#' or centre effect. In particular fixed effects cannot simply be used in this way when there
#' there is heterogeneity in the effect of scientific interest.
#' @param formula Denotes an object of class formula which is a character string which describes
#' the model to be fitted. Most shortcut notation allowed by Rlme4's lmer() function is
#' also allowed by ds.lmerSLMA. Many lmes can be fitted very simply using a formula like:
#' "y~a+b+(1|c)" which simply means fit a lme with y as the outcome variable with a and b as fixed effects, and c
#' as a random effect or grouping factor. This allows for a random slope between groups.
#' By default all such models also include an intercept (regression constant) term.
#' It is also possible to fit models with random slopes by specifying a model such as "y~a+b+(1+b|c)"
#' @param offset  A character string specifying the name of a variable to be used as
#' an offset (effectively
#' a component of the linear predictor of the lme which has a known coefficient a-priori and so does not need to be
#' estimated by the model).
#' @param weights A character string specifying the name of a variable containing
#' prior regression
#' weights for the fitting process.
#' @param dataName A character string specifying the name of an (optional) dataframe that contains
#' all of the variables in the lme formula. This avoids you having to specify the name of the
#' dataframe in front of each covariate in the formula e.g. if the dataframe is
#' called 'DataFrame' you
#' avoid having to write: "DataFrame$y~DataFrame$a+DataFrame$b+(1|DataFrame$c)"
#' Processing stops if a non existing data frame is indicated.
#' @param combine.with.metafor This argument is Boolean. If TRUE (the default) the
#' estimates and standard errors for each regression coefficient are pooled across
#' studies using random effects meta-analysis under maximum likelihood (ML),
#' restricted maximum likelihood (REML), or fixed effects meta-analysis (FE).
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. If the <datasources>
#' the default set of connections will be used: see \link{datashield.connections_default}.
#' @param REML A boolean indicating whether REstricted Maximum Likelihood (REML) should be used for
#' the parameter estimation criterion. This is useful for likelihood ratio tests when assessing the quality of the
#' fixed effects portion of the model. The idea is to compare the models with and without the fixed effects to see
#' if they are significantly different (e.g. using ANOVA). REML assumes that fixed effects structure is correct and
#' so for this type of comparison, it should not be used.
#' @return many of the elements of the output list returned by ds.lmerSLMA.o from
#' each study separately are
#' equivalent to those from lmer() in lme4 with potentially disclosive elements
#' such as individual-level residuals and linear predictors blocked.
#' The return results from each study appear first in the return list with one
#' block of results from each study in the order they appear in datasources.
#' As regards the elements within each study the most important
#' elements are included last in the return list because they then appear at the
#' bottom of a simple print out of the return object. In reverse order, these
#' key elements in reverse order are:
#' @return coefficients:- a matrix in which the first column contains the names of
#' all of the regression parameters (coefficients) in the model, the second column
#' contains the estimated values, the third their corresponding standard errors,
#' the fourth the ratio of estimate/standard error and the fifth the p-value
#' treating that as a standardised normal deviate
#' #' @return CorrMatrix:- the correlation matrix of parameter estimates
#' @return VarCovMatrix:- the variance covariance matrix of parameter estimates
#' @return weights:- the vector (if any) holding regression weights
#' @return offset:- the vector (if any) holding an offset (enters glm with a coefficient of 1.0)
#' @return cov.scaled:- equivalent to VarCovMatrix
#' @return Nmissing:- the number of missing observations in the given study
#' @return Nvalid:- the number of valid (non-missing) observations in the given study
#' @return Ntotal:- the total number of observations in the given study (Nvalid+Nmissing)
#' @return data:- - equivalent to input parameter dataName (above)
#' @return call:- - summary of key elements of the call to fit the model
#' @return there are a small number of more esoteric items of information returned
#' by ds.lmerSLMA.o. Additional information about these can be found in the help
#' file for the lmer() function in the lme4 package.
#' @return input.beta.matrix.for.SLMA:- a matrix containing the vector of coefficient
#' estimates from each study. In combination with the corresponding standard errors
#' (see input.se.matrix.for.SLMA) these can be imported directly into a study level
#' meta-analysis (SLMA) package such as metafor to generate estimates pooled via SLMA
#' @return input.se.matrix.for.SLMA:- a matrix containing the vector of standard error
#' estimates for coefficients from each study. In combination with the coefficients
#' (see input.beta.matrix.for.SLMA) these can be imported directly into a study level
#' meta-analysis (SLMA) package such as metafor to generate estimates pooled via SLMA
#' @return SLMA.pooled.estimates:- a matrix containing pooled estimates for each
#' regression coefficient across all studies with pooling under SLMA via
#' random effects meta-analysis under maximum likelihood (ML), restricted maximum
#' likelihood (REML) or via fixed effects meta-analysis (FE)
#' @author Tom Bishop
#' @export
ds.lmerSLMA.o<-function(formula=NULL, offset=NULL, weights=NULL, combine.with.metafor=TRUE,dataName=NULL,
                       checks=FALSE, maxit=15, datasources=NULL, REML=TRUE) {
  
  # details are provided look for 'opal' objects in the environment
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
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
  
  # set family to gaussian
  family <- 'gaussian'
  
  # if the argument 'dataName' is set, check that the data frame is defined (i.e. exists) on the server site
  if(!(is.null(dataName))){
    defined <- isDefined(datasources, dataName)
  }
  
  # beginning of optional checks - the process stops if any of these checks fails #
  if(checks){
    message(" -- Verifying the variables in the model")
    # call the function that checks the variables in the formula are defined (exist) on the server site and are not missing at complete
    
    glmChecks(formula, dataName, offset, weights, datasources)
  }else{
    #message("WARNING:'checks' is set to FALSE; variables in the model are not checked and error messages may not be intelligible!")
  }

  
  #NOW CALL SECOND COMPONENT OF glmDS TO GENERATE SCORE VECTORS AND INFORMATION MATRICES
  
  cally2 <- call('lmerSLMADS2.o', formula, offset, weights, dataName, REML)
  
  study.summary <- datashield.aggregate(datasources, cally2)
  
  numstudies<-length(datasources)
  
  #INTEGRATE RETURNED OUTPUT
  .select <- function(l, field) {
    lapply(l, function(obj) {obj[[field]]})
  }
  
  disclosure.risk.total<-Reduce(f="+", .select(study.summary, 'disclosure.risk'))
  
  disclosure.risk<-NULL
  errorMessage<-NULL
  
  for(ss in 1:numstudies){
    disclosure.risk<-c(disclosure.risk,study.summary[[ss]]$disclosure.risk)
    errorMessage<-c(errorMessage,study.summary[[ss]]$errorMessage)
  }
  
  disclosure.risk<-as.matrix(disclosure.risk)
  dimnames(disclosure.risk)<-list(names(datasources),"RISK OF DISCLOSURE")
  
  
  errorMessage<-as.matrix(errorMessage)
  dimnames(errorMessage)<-list(names(datasources),"ERROR MESSAGES")
  
  if(disclosure.risk.total>0){
    message("DISCLOSURE RISK IN y.vect, X.mat OR w.vect \nOR MODEL OVERPARAMETERIZED IN AT LEAST ONE STUDY")
    output.blocked.information.1<-"POTENTIAL DISCLOSURE RISK IN y.vect, X.mat OR w.vect"
    output.blocked.information.2<-"OR MODEL OVERPARAMETERIZED IN AT LEAST ONE STUDY."
    output.blocked.information.3<-"FURTHERMORE, CLIENTSIDE FUNCTION MAY HAVE BEEN MODIFIED"
    output.blocked.information.4<-"SO SCORE VECTOR AND INFORMATION MATRIX DESTROYED IN ALL AT-RISK STUDIES"
    
    return(list(output.blocked.information.1,
                output.blocked.information.2,
                output.blocked.information.3,
                output.blocked.information.4
    ))
  }
  
    for(j in 1:numstudies){
    inv.diag.se<-1/sqrt(Matrix::diag(study.summary[[j]]$cov.scaled))

    cor.matrix<-t(diag(inv.diag.se))%*%study.summary[[j]]$cov.scaled%*%(diag(inv.diag.se))
    cor.matrix@Dimnames = study.summary[[j]]$cov.scaled@Dimnames
    study.summary[[j]]$VarCovMatrix<-study.summary[[j]]$cov.scaled
    study.summary[[j]]$CorrMatrix<-cor.matrix
    study.summary[[j]]$RE<-as.data.frame(study.summary[[j]]$RE)
    
  }
  
  
  #ARRANGE betas AND ses AS RETURN OBJECTS TO FEED EASILY
  #INTO A RANDOM EFFECTS META-ANALYSIS FUNCTION SUCH AS metafor
  
  
  #MAKE SURE THAT IF SOME STUDIES HAVE MORE PARAMETERS IN THE
  #FITTED glm (eg BECAUSE OF ALIASING) THE FINAL RETURN MATRICES
  #HAVE ENOUGH ROWS TO FIT THE MAXIMUM LENGTH
  
  
   numcoefficients.max<-0
   
   for(g in 1:numstudies){
     if(length(study.summary[[g]]$coefficients[,1])>numcoefficients.max){
       numcoefficients.max<-length(study.summary[[g]]$coefficients[,1])
     }
   }
   
  numcoefficients<-numcoefficients.max

  betamatrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)
  sematrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)


  for(k in 1:numstudies){
    betamatrix[,k]<-study.summary[[k]]$coefficients[,1]
    sematrix[,k]<-study.summary[[k]]$coefficients[,2]
  }

  ################################################
  #ANNOTATE OUTPUT MATRICES WITH STUDY INDICATORS#
  ################################################

  study.names.list<-NULL
  betas.study.names.list<-NULL
  ses.study.names.list<-NULL
  
  for(v in 1:numstudies){

    study.names.list<-c(study.names.list,paste0("study",as.character(v)))
    betas.study.names.list<-c(betas.study.names.list,paste0("betas study ",as.character(v)))
    ses.study.names.list<-c(ses.study.names.list,paste0("ses study ",as.character(v)))
  }

  dimnames(betamatrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]], betas.study.names.list)
  dimnames(sematrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]], ses.study.names.list)

  output.summary.text<-paste0("list(")

  for(u in 1:numstudies){
    output.summary.text<-paste0(output.summary.text,"study",as.character(u),"=study.summary[[",as.character(u),"]],"," ")
  }

  output.summary.text.save<-output.summary.text
  output.summary.text<-paste0(output.summary.text,"input.beta.matrix.for.SLMA=as.matrix(betamatrix),input.se.matrix.for.SLMA=as.matrix(sematrix))")


  output.summary<-eval(parse(text=output.summary.text))
  #######################END OF ANNOTATION CODE

  SLMA.pooled.ests.matrix<-matrix(NA,nrow<-numcoefficients,ncol=6)


  if(!combine.with.metafor){
    return(output.summary)
  }


  # #IF combine.with.metafor == TRUE, FIRST CHECK THAT THE MODELS IN EACH STUDY MATCH
  # #IF THERE ARE DIFFERENT NUMBERS OF PARAMETERS THE ANALYST WILL
  # #HAVE TO USE THE RETURNED MATRICES FOR betas AND ses TO DETERMINE WHETHER
  # #COMBINATION ACROSS STUDIES IS POSSIBLE AND IF SO, WHICH PARAMETERS GO WITH WHICH
  # if(combine.with.metafor){
  #   numstudies<-length(datasources)
  #   coefficient.vectors.match<-TRUE
  #   for(j in 1:(numstudies-1)){
  #     if(dim(study.summary[[j]]$coefficients)[1]!=dim(study.summary[[(j+1)]]$coefficients)[1])coefficient.vectors.match<-FALSE
  #   }
  #   if(!coefficient.vectors.match){
  #     cat("\n\nModels in different sources vary in structure\nplease match coefficients for meta-analysis individually\n\n")
  #     return(output.summary)
  #   }
  #   
  #   #IF combine.with.metafor == TRUE AND MODEL STRUCTURES MATCH ACROSS ALL STUDIES
  #   #CREATE STUDY LEVEL META-ANALYSIS (SLMA) ESTIMATES FOR ALL PARAMETERS
  #   #USING metafor() AND THREE APPROACHES TO SLMA: SLMA UNDER MAXIMUM LIKELIHOOD (SMLA-ML)
  #   #SLMA UNDER RESTRICTED MAXIMUM LIKELIHOOD (SMLA-REML) AND USING FIXED EFFECTS (SLMA-FE)
  #   
  #   dimnames(SLMA.pooled.ests.matrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]],
  #                                           c("pooled.ML","se.ML","pooled.REML","se.REML","pooled.FE","se.FE"))
  #   
  #   for(p in 1:numcoefficients){
  #     rma.ML<-rma(yi=betamatrix[p,], sei=sematrix[p,], method="ML")
  #     rma.REML<-rma(yi=betamatrix[p,], sei=sematrix[p,], method="REML")
  #     rma.FE<-rma(yi=betamatrix[p,], sei=sematrix[p,], method="FE")
  #     
  #     SLMA.pooled.ests.matrix[p,1]<-rma.ML$beta
  #     SLMA.pooled.ests.matrix[p,2]<-rma.ML$se
  #     
  #     SLMA.pooled.ests.matrix[p,3]<-rma.REML$beta
  #     SLMA.pooled.ests.matrix[p,4]<-rma.REML$se
  #     
  #     SLMA.pooled.ests.matrix[p,5]<-rma.FE$beta
  #     SLMA.pooled.ests.matrix[p,6]<-rma.FE$se
  #     
  #   }
  #   
  # }
  # 
  # output.summary.plus.pooled.SLMA.text<-paste0(output.summary.text.save,
  #                                              "input.beta.matrix.for.SLMA=as.matrix(betamatrix),input.se.matrix.for.SLMA=as.matrix(sematrix),SLMA.pooled.estimates=SLMA.pooled.ests.matrix)")
  # 
  # 
  # output.summary.plus.pooled.SLMA<-eval(parse(text=output.summary.plus.pooled.SLMA.text))
  # 
  # 
  #  return(output.summary.plus.pooled.SLMA)
  #return(study.summary)
}

# ds.lmerSLMA.o
