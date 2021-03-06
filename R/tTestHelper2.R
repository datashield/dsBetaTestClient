#' 
#' @title Uses glm to compute means of numeric vector across factor vector
#' @description This is an internal function.
#' @details internal function
#' @param formula a character form of a formula (e.g. 'a~b')
#' @param CI a numeric, the confidence interval
#' @param family a description of the error distribution function to use in the model
#' @param datasources see the calling function
#' @keywords internal see the calling function
#' @return results similar to that of a t.test 
#'
tTestHelper2 <- function(formula, CI, datasources) {
  
  # The variables names 
  a <- unlist(strsplit(formula, split='~'))[1]
  b <- unlist(strsplit(formula, split='~'))[2]
  
  # check if the variables in the formula are defined in all the studies, if not defined a message is thrown and the process stops
  defined <- dsBaseClient:::isDefined(datasources, a)
  defined <- dsBaseClient:::isDefined(datasources, b)
  
  # Check that the variables in the formula are of the right type: 'numeric'/'integer' for the outcome and 'factor' for the covariate
  typ1 <- dsBaseClient:::checkClass(datasources, a)
  if(typ1 != "numeric" & typ1 != "integer"){
    stop(paste0(" ", a, " must be a numeric vector!"), call.=FALSE)
  }
  typ2 <- dsBaseClient:::checkClass(datasources, b)
  if(typ2 != "factor"){
    stop(paste0(" ", b, " must be a factor vector!"), call.=FALSE)
  }else{
    # the covariate must have only two categories (as for t.test in R), so stop and throw a message otherwise    
    cally <- paste0("levels(", b, ")")
    levels_all <- opal::datashield.aggregate(datasources, as.symbol(cally))
    classes <- unique(unlist(levels_all))
    if(length(classes) != 2){
      stop(paste0(" ", b, " must two and only two categories!"), call.=FALSE)
    }
  }
  
  # turn the formula provided as a character into a formula object
  formula <- stats::as.formula(formula)
  family <- 'gaussian'
  
  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
  
  # start beta values
  beta.vect.next <- c(0,0)
  beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  
  # Iterations need to be counted. Start off with the count at 0
  # and increment by 1 at each new iteration
  iteration.count <- 0
  
  # identify the correct dimension for start beta coeffs by calling the 1st component of glmDS
  cally1 <- call('glmDS1', formula, family, beta.vect=beta.vect.temp, NULL)
  
  study.summary <- opal::datashield.aggregate(datasources, cally1)
  num.par.glm <- study.summary[[1]][[1]][[2]]
  
  beta.vect.next <- rep(0, num.par.glm)
  beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  
  
  # Provide arbitrary starting value for deviance to enable subsequent calculation of the
  # change in deviance between iterations
  dev.old <- 9.99e+99
  
  #Convergence state needs to be monitored.
  converge.state <- FALSE
  
  # Define a convergence criterion. This value of epsilon corresponds to that used
  # by default for GLMs in R (see section S3 for details)
  epsilon <- 1.0e-08
  
  f <- NULL
  
  while(!converge.state && iteration.count < 20) {
    
    iteration.count <- iteration.count+1
    
    # now call second component of glmDS to generate score vectors and informations matrices
    cally2 <- call('glmDS2', formula, family, beta.vect=beta.vect.temp, NULL, NULL, NULL)
    
    study.summary <- opal::datashield.aggregate(datasources, cally2)
    
    
    .select <- function(l, field) {
      lapply(l, function(obj) {obj[[field]]})
    }
    
    info.matrix.total<-Reduce(f="+", .select(study.summary, 'info.matrix'))
    score.vect.total<-Reduce(f="+", .select(study.summary, 'score.vect'))
    dev.total<-Reduce(f="+", .select(study.summary, 'dev'))
    
    if(iteration.count==1) {
      # Sum participants only during first iteration.
      nsubs.total <- Reduce(f="+", .select(study.summary, 'numsubs'))
      # Save family
      f <- study.summary[[1]]$family
    }
    
    # create variance covariance matrix as inverse of information matrix
    variance.covariance.matrix.total<-solve(info.matrix.total)
    
    # create beta vector update terms
    beta.update.vect<-variance.covariance.matrix.total %*% score.vect.total
    
    # add update terms to current beta vector to obtain new beta vector for next iteration
    if(iteration.count==1){
      beta.vect.next<-rep(0,length(beta.update.vect))
    }
    
    beta.vect.next<-beta.vect.next+beta.update.vect
    beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  
    # calculate value of convergence statistic and test whether meets convergence criterion
    converge.value<-abs(dev.total-dev.old)/(abs(dev.total)+0.1)
    if(converge.value<=epsilon)converge.state<-TRUE
    if(converge.value>epsilon)dev.old<-dev.total
    
  }
  
  # if convergence has been obtained, declare final (maximum likelihood) beta vector,
  # and calculate the corresponding standard errors, t scores and p values
  # (the latter two to be consistent with the output of a standard GLM analysis)
  # then print out final model summary
  if(converge.state){
    
    beta.vect.final<-beta.vect.next
    
    scale.par <- 1
    scale.par <- dev.total / (nsubs.total-length(beta.vect.next))
    
    t.df <- (nsubs.total-length(beta.vect.next))
    se.vect.final <- sqrt(diag(variance.covariance.matrix.total)) * sqrt(scale.par)
    z.vect.final <- beta.vect.final/se.vect.final
    pval.vect.final <- 2*stats::pt(q=-abs(z.vect.final),df=t.df)
    parameter.names <- names(score.vect.total[,1])
    model.parameters <- cbind(beta.vect.final,se.vect.final,z.vect.final,pval.vect.final)
    dimnames(model.parameters) <- list(parameter.names,c("Estimate","Std. Error","t-value","p-value"))
    
    if(CI > 0){
      
      ci.mult <- stats::qt(p=(1-(1-CI)/2),df=t.df)
      low.ci.lp <- model.parameters[,1]-ci.mult*model.parameters[,2]
      hi.ci.lp <- model.parameters[,1]+ci.mult*model.parameters[,2]
      estimate.lp <- model.parameters[,1]
      
      estimate.natural <- estimate.lp
      low.ci.natural <- low.ci.lp
      hi.ci.natural <- hi.ci.lp
      name1 <- paste0("low",CI,"CI")
      name2 <- paste0("high",CI,"CI")
      ci.mat <- cbind(low.ci.lp,hi.ci.lp)
      dimnames(ci.mat) <- list(NULL,c(name1,name2))   
   }
    
    model.parameters <- cbind(model.parameters,ci.mat)
    
    formulatext <- Reduce(paste, deparse(formula))
    
    glmds <- list(
      formula = formulatext,
      coefficients = model.parameters,
      dev = dev.total,
      nsubs = nsubs.total,
      df = (nsubs.total-length(beta.vect.next)),
      iter = iteration.count
    )
    
    class(glmds) <- 'glmds'
   
    # get the values that one would expect from a t.test: means for each category, p.value etc...
    meanA <- round(glmds$coefficients[1,1],4)
    meanB <- round(meanA + glmds$coefficients[2,1],4)
    meanLabel <- paste0("mean in group ", classes[1], " & ",  "mean in group ", classes[2])
    output <- list(
      statistic = glmds$coefficients[2,3],
      parameter = glmds$df,
      p.value = glmds$coefficients[2,4],
      conf.int = c(glmds$coefficients[2,5], glmds$coefficients[2,6]),
      estimate = c(meanA, meanB),
      alternative = "two.sided",
      method = "GLM",
      data.name = paste0(a, " by ", b)
    )
    # print results to screen in a 't.test' way
    mylabels <- c("t ", "df ", "p.value ", paste0(CI*100, " percent confidence interval"), "sample estimates:\nmeanLabel")
    message("GLM to assess difference of statistical means across two categories")
    message(paste0("data: ", a, " by ", b))
    message(paste0("t = ", round(glmds$coefficients[2,3],4), ", df = ", glmds$df, ", p.value = ", signif(glmds$coefficients[2,4]),4))
    message(paste0(CI*100, " percent confidence interval:"))
    message(signif(output$conf.int),5)
    message(paste0("sample estimates: ", meanLabel))
    message(paste0(output$estimate[1], " ", output$estimate[2]))
   
    return(output)
  } else {
    warning(paste("GLM did not converge after", 20, "iterations."))
    return(NULL)
  }
}
