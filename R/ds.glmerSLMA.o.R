
#' @title ds.glmSLMA.o calling glmDS1.o, glmSLMADS2.o
#' @description Fits a generalized linear model (glm) on data from a single or multiple sources
#' @details Fits a thing

#' @export
ds.glmerSLMA.o<-function(formula=NULL, family=NULL, offset=NULL, weights=NULL, combine.with.metafor=TRUE,dataName=NULL,
                       checks=FALSE, maxit=15, datasources=NULL) {
  
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
  
  
  # check that 'family' was set
  if(is.null(family)){
    stop(" Please provide a valid 'family' argument!", call.=FALSE)
  }
  
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
  
  #MOVE ITERATION COUNT BEFORE ASSIGNMENT OF beta.vect.next
  #Iterations need to be counted. Start off with the count at 0
  #and increment by 1 at each new iteration
  #iteration.count<-0
  
  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
  
  
  #ARBITRARY LENGTH FOR START BETAs AT THIS STAGE BUT IN LEGAL TRANSMISSION FORMAT ("0,0,0,0,0")
  #beta.vect.next <- rep(0,5)
  #beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  
  
  #IDENTIFY THE CORRECT DIMENSION FOR START BETAs VIA CALLING FIRST COMPONENT OF glmDS
  
  #cally1 <- call('glmDS1.o', formula, family, weights, dataName)
  
  #study.summary.0 <- datashield.aggregate(datasources, cally1)
  
  
  # at.least.one.study.data.error<-0
  # 
  # for(hh in 1:numstudies) {
  #   if(study.summary.0[[hh]]$errorMessage!="No errors"){
  #     at.least.one.study.data.error<-1
  #   }
  # }
  # 
  # 
  # num.par.glm<-NULL
  # coef.names<-NULL
  # 
  # if(at.least.one.study.data.error==0){
  #   num.par.glm<-study.summary.0[[1]][[1]][[2]]
  #   coef.names<-study.summary.0[[1]][[2]]
  # }
  # 
  # y.invalid<-NULL
  # Xpar.invalid<-NULL
  # w.invalid<-NULL
  # glm.saturation.invalid<-NULL
  # errorMessage<-NULL
  # 
  # for(ss in 1:numstudies)
  # {
  #   y.invalid<-c(y.invalid,study.summary.0[[ss]][[3]])
  #   Xpar.invalid<-rbind(Xpar.invalid,study.summary.0[[ss]][[4]])
  #   w.invalid<-c(w.invalid,study.summary.0[[ss]][[5]])
  #   glm.saturation.invalid <-c(glm.saturation.invalid,study.summary.0[[ss]][[6]])
  #   errorMessage<-c(errorMessage,study.summary.0[[ss]][[7]])
  # }
  # 
  # y.invalid<-as.matrix(y.invalid)
  # sum.y.invalid<-sum(y.invalid)
  # dimnames(y.invalid)<-list(names(datasources),"Y VECTOR")
  # 
  # Xpar.invalid<-as.matrix(Xpar.invalid)
  # sum.Xpar.invalid<-sum(Xpar.invalid)
  # dimnames(Xpar.invalid)<-list(names(datasources),coef.names)
  # 
  # 
  # w.invalid<-as.matrix(w.invalid)
  # sum.w.invalid<-sum(w.invalid)
  # dimnames(w.invalid)<-list(names(datasources),"WEIGHT VECTOR")
  # 
  # glm.saturation.invalid<-as.matrix(glm.saturation.invalid)
  # sum.glm.saturation.invalid<-sum(glm.saturation.invalid)
  # dimnames(glm.saturation.invalid)<-list(names(datasources),"MODEL OVERPARAMETERIZED")
  # 
  # errorMessage<-as.matrix(errorMessage)
  # dimnames(errorMessage)<-list(names(datasources),"ERROR MESSAGES")
  # 
  # 
  # 
  # output.blocked.information.1<-"MODEL FITTING TERMINATED AT FIRST ITERATION:"
  # output.blocked.information.2<-"ANY VALUES OF 1 IN THE FOLLOWING TABLES DENOTE"
  # output.blocked.information.3<-"POTENTIAL DISCLOSURE RISKS. PLEASE USE THE ARGUMENT"
  # output.blocked.information.4<-"[datasources=] TO EXCLUDE STUDIES WITH DATA ERRORS"
  # 
  # 
  # 
  # 
  # if(sum.y.invalid>0||sum.Xpar.invalid>0||sum.w.invalid>0||sum.glm.saturation.invalid>0||at.least.one.study.data.error==1){
  #   return(list(
  #     output.blocked.information.1,
  #     output.blocked.information.2,
  #     output.blocked.information.3,
  #     output.blocked.information.4,
  #     y.vector.error=y.invalid,
  #     X.matrix.error=Xpar.invalid,
  #     weight.vector.error=w.invalid,
  #     glm.overparameterized=glm.saturation.invalid,
  #     errorMessage=errorMessage
  #   ))
  #   stop("DATA ERROR")
  # }
  # 
  # 
  # 
  # beta.vect.next <- rep(0,num.par.glm)
  # beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  
  
  #Provide arbitrary starting value for deviance to enable subsequent calculation of the
  #change in deviance between iterations
  #dev.old<-9.99e+99
  
  #Convergence state needs to be monitored.
  #converge.state<-FALSE
  
  #Define a convergence criterion. This value of epsilon corresponds to that used
  #by default for GLMs in R (see section S3 for details)
  #epsilon<-1.0e-08
  
  #f<-NULL
  
  #NOW CALL SECOND COMPONENT OF glmDS TO GENERATE SCORE VECTORS AND INFORMATION MATRICES
  
  cally2 <- call('glmerSLMADS2.o', formula, family, offset, weights, dataName)
  
  study.summary <- datashield.aggregate(datasources, cally2)
  
  numstudies<-length(datasources)
  
  # for(j in 1:numstudies){
  #   inv.diag.se<-1/sqrt(diag(study.summary[[j]]$cov.scaled))
  #   
  #   cor.matrix<-t(diag(inv.diag.se))%*%study.summary[[j]]$cov.scaled%*%(diag(inv.diag.se))
  #   study.summary[[j]]$VarCovMatrix<-study.summary[[j]]$cov.scaled
  #   study.summary[[j]]$CorrMatrix<-cor.matrix
  # }
  
  
  #ARRANGE betas AND ses AS RETURN OBJECTS TO FEED EASILY
  #INTO A RANDOM EFFECTS META-ANALYSIS FUNCTION SUCH AS metafor
  
  
  #MAKE SURE THAT IF SOME STUDIES HAVE MORE PARAMETERS IN THE
  #FITTED glm (eg BECAUSE OF ALIASING) THE FINAL RETURN MATRICES
  #HAVE ENOUGH ROWS TO FIT THE MAXIMUM LENGTH
  
  
  # numcoefficients.max<-0
  # 
  # for(g in 1:numstudies){
  #   if(length(study.summary[[g]]$coefficients[,1])>numcoefficients.max){
  #     numcoefficients.max<-length(study.summary[[g]]$coefficients[,1])
  #   }
  # }
  # 
  # numcoefficients<-numcoefficients.max
  # 
  # betamatrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)
  # sematrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)
  # 
  # 
  # for(k in 1:numstudies){
  #   betamatrix[,k]<-study.summary[[k]]$coefficients[,1]
  #   sematrix[,k]<-study.summary[[k]]$coefficients[,2]
  # }
  # 
  # ################################################
  # #ANNOTATE OUTPUT MATRICES WITH STUDY INDICATORS#
  # ################################################
  # 
  # study.names.list<-NULL
  # betas.study.names.list<-NULL
  # ses.study.names.list<-NULL
  # 
  # for(v in 1:numstudies){
  #   
  #   study.names.list<-c(study.names.list,paste0("study",as.character(v)))
  #   betas.study.names.list<-c(betas.study.names.list,paste0("betas study ",as.character(v)))
  #   ses.study.names.list<-c(ses.study.names.list,paste0("ses study ",as.character(v)))
  # }
  # 
  # dimnames(betamatrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]], betas.study.names.list)
  # dimnames(sematrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]], ses.study.names.list)
  # 
  # output.summary.text<-paste0("list(")
  # 
  # for(u in 1:numstudies){
  #   output.summary.text<-paste0(output.summary.text,"study",as.character(u),"=study.summary[[",as.character(u),"]],"," ")
  # }
  # 
  # output.summary.text.save<-output.summary.text
  # output.summary.text<-paste0(output.summary.text,"input.beta.matrix.for.SLMA=as.matrix(betamatrix),input.se.matrix.for.SLMA=as.matrix(sematrix))")
  # 
  # 
  # output.summary<-eval(parse(text=output.summary.text))
  # #######################END OF ANNOTATION CODE
  # 
  # SLMA.pooled.ests.matrix<-matrix(NA,nrow<-numcoefficients,ncol=6)
  # 
  # 
  # if(!combine.with.metafor){
  #   return(output.summary)
  # }
  # 
  # 
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
  # return(output.summary.plus.pooled.SLMA)
  return(study.summary)
}

# ds.glmerSLMA.o
