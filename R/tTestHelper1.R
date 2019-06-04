#' 
#' @title runs a t-test for two continuous variables
#' @description This is an internal function.
#' @details internal function
#' @param type see the calling function
#' @param alternative see the calling function
#' @param mu see the calling function 
#' @param paired see the calling function
#' @param var.equal see the calling function
#' @param conf.level see the calling function
#' @param datasources see the calling function
#' @keywords internal see the calling function
#' @return the results of the t-test 
#'
tTestHelper1 <- function(x, y, type, alternative, mu, paired, var.equal, conf.level, datasources){
  # get the names of the variables used for the analysis
  # the input variable might be given as column table (i.e. D$x)
  # or just as a vector not attached to a table (i.e. x)
  # we have to make sure the function deals with each case
  if(is.null(y)){
    xname <- extract(x)
    dname = xname$elements
    variables <- dname
  }else{
    xname <- extract(x)
    yname <- extract(y)
    dname1 = xname$elements
    dname2 = yname$elements
    variables <- c(dname1, dname2)
    dname = paste(dname1, 'and', dname2)
  }
  
  # call the function that checks theinput variables are defined in all the studies
  if(is.null(y)){
    obj2lookfor <- xname$holders
    if(is.na(obj2lookfor)){
      defined <- isDefined(datasources, variables[1])
    }else{
      defined <- isDefined(datasources, obj2lookfor)
    }
  }else{
    obj2lookfor1 <- xname$holders
    obj2lookfor2 <- yname$holders
    if(is.na(obj2lookfor1)){
      defined <- isDefined(datasources, variables[1])
    }else{
      defined <- isDefined(datasources, obj2lookfor1)
    }
    if(is.na(obj2lookfor2)){
      defined <- isDefined(datasources, variables[2])
    }else{
      defined <- isDefined(datasources, obj2lookfor2)
    }
  }
  
  # call the internal function that checks an input object is of the same class in all studies.
  if(is.null(y)){
    typ1 <- checkClass(datasources, x)
  }else{
    typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, y)
  }
  
  # number of studies
  num.sources = length(datasources)
  
  if(type == "combine"){
    
    # Performs t-test on merged data sets
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
      stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                   conf.level < 0 || conf.level > 1)) 
      stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
      if (paired) {
        cally = paste0("complete.cases(",  x, ",", y, ")")
        opal::datashield.assign(datasources, 'pair.compl.obs', as.symbol(cally))
        cally = paste0("subsetDS('",  x, "', pair.compl.obs)")
        opal::datashield.assign(datasources, 'xok', as.symbol(cally))
        cally = paste0("subsetDS('",  y, "', pair.compl.obs)")
        opal::datashield.assign(datasources, 'yok', as.symbol(cally))
      } else {
        cally = paste0("complete.cases(",  x, ")")
        opal::datashield.assign(datasources, 'not.na.x', as.symbol(cally))
        cally = paste0("subsetDS('",  x, "',not.na.x)")
        opal::datashield.assign(datasources, 'xok', as.symbol(cally))
        
        cally = paste0("complete.cases(",  y, ")")
        opal::datashield.assign(datasources, 'not.na.y', as.symbol(cally))
        cally = paste0("subsetDS('",  y, "',not.na.y)")
        opal::datashield.assign(datasources, 'yok', as.symbol(cally))
      }
    } else {
      if (paired) 
        stop("'y' is missing for paired test")
      cally = paste0("complete.cases(",  x, ")")
      opal::datashield.assign(datasources, 'not.na.x', as.symbol(cally))
      cally = paste0("subsetDS('",  x, "', not.na.x)")
      opal::datashield.assign(datasources, 'xok', as.symbol(cally))
      cally = paste0("as.null(",  x, ")")
      opal::datashield.assign(datasources, 'yok', as.symbol(cally)) # does not matter that as.null(x) since we just want to set y to NULL
    }    
    
    if (paired) {
      cally = paste0("(yok)","*(",-1,")")
      opal::datashield.assign(datasources, 'minus_y', as.symbol(cally))
      # opal::datashield.assign(datasources, 'dummy', quote(cbind(xok, minus_y)))
      opal::datashield.assign(datasources, 'xok', as.symbol("xok+minus_y"))
      opal::datashield.assign(datasources, 'yok', as.symbol("as.null(yok)"))
    }
    
    length.local.x = opal::datashield.aggregate(datasources, as.symbol("NROW(xok)"))
    mean.local.x = opal::datashield.aggregate(datasources, as.symbol("meanDS(xok)"))
    var.local.x = opal::datashield.aggregate(datasources, as.symbol("varDS(xok)"))
    
    length.total.x = 0
    sum.weighted.x = 0
    
    for (i in 1:num.sources) 
      if (!is.null(length.local.x[[i]]))
        length.total.x = length.total.x + length.local.x[[i]]
    
    for (i in 1:num.sources) 
      if (!is.null(length.local.x[[i]]))
        sum.weighted.x = sum.weighted.x + length.local.x[[i]]*mean.local.x[[i]]
    
    if (!is.na(sum.weighted.x))
      mean.global.x = sum.weighted.x/length.total.x else
        stop(paste("Check the data supplied: global ", variables[1], " mean is NA", sep=""))
    estimate = mean.global.x
    
    nrows_var.x = NROW(var.local.x[[1]])
    ncols_var.x = NCOL(var.local.x[[1]])
    dummy.sum.x = matrix(0, nrows_var.x, ncols_var.x)
    
    for (i in 1:num.sources) {
      if (!is.null(var.local.x[[i]]) & !is.null(mean.local.x[[i]]))
        if (!is.na(var.local.x[[i]]) & !is.na(mean.local.x[[i]])) {
          var.weight.x = (length.local.x[[i]]-1)*var.local.x[[i]]
          add.elem.x = length.local.x[[i]]*(mean.local.x[[i]]%x%t(mean.local.x[[i]]))
          dummy.sum.x = dummy.sum.x +var.weight.x+add.elem.x
        }      
    }
    mean.global.products.x = length.total.x*(mean.global.x%x%t(mean.global.x))
    var.global.x = 1/(length.total.x-1)*(dummy.sum.x-mean.global.products.x)
    
    null.y = opal::datashield.aggregate(datasources, as.symbol("is.null(yok)"))
    null.y = unlist(null.y)
    
    if (all(null.y)) {
      if (length.total.x < 2) 
        stop("not enough 'x' observations")
      df <- length.total.x - 1
      stderr <- sqrt(var.global.x/length.total.x)
      if (stderr < 10 * .Machine$double.eps * abs(mean.global.x)) 
        stop("data are essentially constant")
      tstat <- (mean.global.x - mu)/stderr
      method <- ifelse(paired, "Paired t-test", "One Sample t-test")
      names(estimate) <- ifelse(paired, "mean of the differences", paste("mean of", variables[1], sep=""))
    } else {
      length.local.y = opal::datashield.aggregate(datasources, as.symbol("NROW(yok)"))
      
      length.total.y = 0
      sum.weighted.y = 0
      
      for (i in 1:num.sources) 
        if (!is.null(length.local.y[[i]]))
          length.total.y = length.total.y + length.local.y[[i]]
      
      if (length.total.x < 1 || (!var.equal && length.total.x < 2)) 
        stop(paste("not enough ", variables[1], "observations", sep=""))
      if (length.total.y < 1 || (!var.equal && length.total.y < 2)) 
        stop(paste("not enough ", variables[2], "observations", sep=""))
      if (var.equal && length.total.x + length.total.y < 3) 
        stop("not enough observations")
      
      mean.local.y = opal::datashield.aggregate(datasources, as.symbol("meanDS(yok)"))
      var.local.y = opal::datashield.aggregate(datasources, as.symbol("varDS(yok)"))
      method <- paste(if (!var.equal) 
        "Welch", "Two Sample t-test")
      
      length.total.y = 0
      sum.weighted.y = 0
      
      for (i in 1:num.sources) 
        if (!is.null(length.local.y[[i]])) {
          length.total.y = length.total.y + length.local.y[[i]]
          sum.weighted.y = sum.weighted.y + length.local.y[[i]]*mean.local.y[[i]]
        }
      if (!is.na(sum.weighted.y))
        mean.global.y = sum.weighted.y/length.total.y else
          stop(paste("Check the data supplied: global ", variables[2], " mean is NA", sep=""))
      
      estimate <- c(mean.global.x, mean.global.y)
      names(estimate) <- c(paste("mean of ", variables[1], sep=""), paste("mean of ", variables[2], sep=""))
      
      nrows_var.y = NROW(var.local.y[[1]])
      ncols_var.y = NCOL(var.local.y[[1]])
      dummy.sum.y = matrix(0, nrows_var.y, ncols_var.y)
      
      for (i in 1:num.sources) {
        if (!is.null(var.local.y[[i]]) & !is.null(mean.local.y[[i]]))
          if (!is.na(var.local.y[[i]]) & !is.na(mean.local.y[[i]])) {
            var.weight.y = (length.local.y[[i]]-1)*var.local.y[[i]]
            add.elem.y = length.local.y[[i]]*(mean.local.y[[i]]%x%t(mean.local.y[[i]]))
            dummy.sum.y = dummy.sum.y +var.weight.y+add.elem.y
          }      
      }
      mean.global.products.y = length.total.y*(mean.global.y%x%t(mean.global.y))
      var.global.y = 1/(length.total.y-1)*(dummy.sum.y-mean.global.products.y)
      
      if (var.equal) {
        df <- length.total.x + length.total.x - 2
        v <- 0
        if (length.total.x > 1) 
          v <- v + (length.total.x - 1) * var.global.x
        if (length.total.y > 1) 
          v <- v + (length.total.y - 1) * var.global.y
        v <- v/df
        stderr <- sqrt(v * (1/length.total.x + 1/length.total.y))
      } else {
        stderrx <- sqrt(var.global.x/length.total.x)
        stderry <- sqrt(var.global.y/length.total.y)
        stderr <- sqrt(stderrx^2 + stderry^2)
        df <- stderr^4/(stderrx^4/(length.total.x - 1) + stderry^4/(length.total.y - 1))
      }
      if (stderr < 10 * .Machine$double.eps * max(abs(mean.global.x), 
                                                  abs(mean.global.y))) 
        stop("data are essentially constant")
      tstat <- (mean.global.x - mean.global.y - mu)/stderr
    }
    
    
    if (alternative == "less") {
      pval <- stats::pt(tstat, df)
      cint <- c(-Inf, tstat + stats::qt(conf.level, df))
    } else if (alternative == "greater") {
      pval <- stats::pt(tstat, df, lower.tail = FALSE)
      cint <- c(tstat - stats::qt(conf.level, df), Inf)
    } else {
      pval <- 2 * stats::pt(-abs(tstat), df)
      alpha <- 1 - conf.level
      cint <- stats::qt(1 - alpha/2, df)
      cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y)) 
      "difference in means" else
        "mean"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
                 conf.int = cint, estimate = estimate, null.value = mu, 
                 alternative = alternative, method = method, data.name=dname)
    class(rval) <- "htest"
    
    # delete files that are no more required
    opal::datashield.rm(datasources, 'pair.compl.obs')
    opal::datashield.rm(datasources, 'xok')
    opal::datashield.rm(datasources, 'yok')
    opal::datashield.rm(datasources, 'not.na.x')
    opal::datashield.rm(datasources, 'minus_y')  
    
    return(rval)
    
  }else{
    if(type == "split"){
      cally <- paste0("t.test(", x, ",", y, ",alternative='",alternative, "',mu=",mu, ",paired=",paired, ",var.equal=",var.equal, ",conf.level=",conf.level,")")
      results <- opal::datashield.aggregate(datasources, as.symbol(cally))
      return(results)
    }else{
      stop('Function argument "type" has to be either "combine" or "split"')
    }
  }
}