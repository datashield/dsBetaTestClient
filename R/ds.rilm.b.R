#' 
#' @title Random Intercept Linear Model (RILM)
#' @description Runs a Random Intercept Linear Model on horizontally-partinioned data
#' @details  see description
#' @param ... the explanatory variables, could be factor or numeric vectors
#' @param y the dependent variable
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal datasources.
#' @return regression outcome
#' @author Avraam D; Jones EM; van den Heuvel ER
#' @export
#' 
ds.rilm.b <- function(..., y=NULL, datasources=NULL){

##############################
is.whole <- function(x){
    start <- 1
    end <- length(x) + 1
    while (start < end) {
        y <- x[start]
        test <- floor(y)
        if (y == test) {
            if (start == 1) {
                result = TRUE
            }
            else {
                result <- c(result, TRUE)
            }
        }
        else {
            if (start == 1) {
                result = FALSE
            }
            else {
                result <- c(result, FALSE)
            }
        }
        start <- start + 1
    }
    return(result)
}
###############################



	# Set number of decimals.
	options(digits = 20)

	# INPUT VARIABLES
	# Max number of iterations
	Max <- 1 * (10^8)
	# Convergence Criteria
	Conf <- (1 / (10^8))
	# Counting Variable
	Count <- 1
	# Print a check every [Milestone] iterations
	Milestone <- 1000000

	xvars <- cbind(...)
	N.xvar <- length(xvars)
		
	if(is.null(y)){
		stop("Please provide the y-variable", call.=FALSE)
	}
	
	if(is.null(datasources)){
		datasources <- findLoginObjects()
	}

	# Find the number of centers.
	N.centres <- length(datasources)
			
	call <- paste0("rilmDS.b(", paste(..., sep=", "), ", ", "yvar=", y, ")")
	output <- opal::datashield.aggregate(datasources, call)

	N.subs <- list()
	XX_Matrixes <- list()
	XY_Matrixes <- list()
	YY_Matrixes <- list()

	for (i in 1:N.centres){
		N.subs[[i]] <- output[[i]][[1]]
		XX_Matrixes[[i]] <- output[[i]][[2]]
		XY_Matrixes[[i]] <- output[[i]][[3]]
		YY_Matrixes[[i]] <- output[[i]][[4]]
	}

	# All the functions used in the procedure are here:
	  
	# Function for calculating R(k+1): ZigmaS2=ZigmaS2(k), ZigmaR2=ZigmaR2(k)
	fR <- function(ZigmaS2, ZigmaR2){
		(ZigmaS2) / (ZigmaS2 + ZigmaR2)
	}
	  
	# Function for calculating Rho(k+1): i=group (center), R=R(k+1), A=X'X list
	fRho <- function(i, R, A){
		(A[[i]][1,1] * R) / (1 + ((A[[i]][1,1] - 1) * R))
	}
	  
	# Function for Tau^2(k+1): i=group (center), ZigmaS2=ZigmaS2(k), Rho=Rho(k+1)
	fTau <- function(i, ZigmaS2, Rho){
		ZigmaS2 * (1 - Rho[i,1])
	}

	# Function for Delta(k+1): i=group (center), Beta=Beta(k), Rho=Rho(k+1), A=X'X list
	fDelta <- function(i, Beta, Rho, A, XY){
		# Calculate Rho
		Rho[i,1] * (1/A[[i]][1,1]) * (XY[[i]][1,1] - (t(Beta) %*% A[[i]][,1]))
	}

	# Function for Beta(k+1): m=number of centers, Delta=Delta(k+1), A=X'X List, XY=X'Y List
	fBeta <- function(m, Delta, A, XY){
		# Create Hulp List.
		Hulp <- list()
		Asum <- Reduce('+', A)
		XYsum <- Reduce('+', XY)
		for (i in 1:m){
			Hulp[[i]] <- Delta[i,1] * A[[i]][,1]
		}
		DeltaXsum <- Reduce('+', Hulp)
		# Calculate the Beta vector.
		(solve(Asum)) %*% (XYsum - DeltaXsum) 
		#(chol2inv(Asum,size=ncol(Asum),LINPACK=FALSE)) %*% (XYsum - DeltaXsum)
		#qr.solve(Asum,XYsum - DeltaXsum,tol=1e-10)
	}

	# Function for ZigmaS2(k+1): m=number of centers, Delta2=Delta2(k+1), Tau2=Tau2(k+1).
	fZigmaS2 <- function(m, Delta2, Tau2){
		(1/m) * sum(Delta2 + Tau2)
	}

	# Function for ZigmaR(k+1): m=number of centers,Beta=Beta(k+1),Delta=Delta(k+1),Delta2=Delta2(k+1),Tau2=Tau2(k+1),
	# XY=X'Y List,Y=I'Y List,Y2=Y'Y List,A=X'X List
	fZigmaR2 <- function(m, Beta, Delta, Delta2, Tau2, XY, YY, A){
		Term1 <- list()
		Term2 <- list()
		Term3 <- list()
		Term4 <- list()
		Term5 <- list()
		Term6 <- list()
		Term7 <- list()
		for (i in 1:m){
			Term1[[i]] <- A[[i]][1,1] * Tau2[i,1]
			Term2[[i]] <- YY[[i]]
			Term3[[i]] <- 2 * t(XY[[i]]) %*% Beta
			Term4[[i]] <- 2 * Delta[i,1] * XY[[i]][1,1]
			Term5[[i]] <- t(Beta) %*% A[[i]] %*% Beta
			Term6[[i]] <- 2 * Delta[i,1] * A[[i]][,1] %*% Beta
			Term7[[i]] <- A[[i]][1,1] * Delta2[i,1]
		}
		(1/Reduce('+', A)[1,1]) * (Reduce('+', Term1) + Reduce('+', Term2) - Reduce('+', Term3) - Reduce('+', Term4) + Reduce('+', Term5) + Reduce('+', Term6) + Reduce('+', Term7))
	}

	# Function for checking if the absolute difference between (k) and (k+1) is smaller than the convergence criteria.
	# New_Beta=Beta(k+1), New_ZigmaS2=ZigmaS2(k+1), New_ZigmaR2=ZigmaR2(k+1)
	# Old_Beta=Beta(k), Old_ZigmaS2=ZigmaS2(k), Old_ZigmaR2=ZigmaR2(k), Conf=Confergence Criteria
	fControle <- function(New_Beta, Old_Beta, New_ZigmaS2, Old_ZigmaS2, New_ZigmaR2, Old_ZigmaR2, Conf){
		Dummy1 <- ifelse(abs(New_Beta - Old_Beta) >= Conf, 1, 0)
		Dummy2 <- ifelse(abs(New_ZigmaS2 - Old_ZigmaS2) >= Conf, 1, 0)
		Dummy3 <- ifelse(abs(New_ZigmaR2 - Old_ZigmaR2) >= Conf, 1, 0)
		sum(Dummy1) + Dummy2 + Dummy3
	}


	# CREATE MATRIX.
	# Rho.
	Rho <- matrix(nrow=N.centres, ncol=1)
	rownames(Rho) <- c(1:N.centres)
	colnames(Rho) <- "Rho"
	# Tau2.
	Tau2 <- matrix(nrow=N.centres, ncol=1)
	rownames(Tau2) <- c(1:N.centres)
	colnames(Tau2) <- "Tau2"
	# Delta.
	Delta <- matrix(nrow=N.centres, ncol=1)
	rownames(Delta) <- c(1:N.centres)
	colnames(Delta) <- "Delta"
	# Delta2.
	Delta2 <- matrix(nrow=N.centres, ncol=1)
	rownames(Delta2) <- c(1:N.centres)
	colnames(Delta2) <- "Delta2"
	# Beta.
	Beta <- matrix(nrow=(N.xvar+1), ncol=1)
	rownames(Beta) <- c(0:(N.xvar))
	colnames(Beta) <- "Beta"

	# INITIAL VALUES
	Beta[,1]   <- rep(1, (N.xvar+1))
	ZigmaS <- 0
	ZigmaR <- 1
	Test <- 1

	# Save Initial Values.
	Start_Beta <- Beta
	Start_ZigmaS <- ZigmaS
	Start_ZigmaR <- ZigmaR

	  
	# Set the loop to continue till either max iterations are reached or all differences are lower than the Conf Crit.
	while (Count <= Max && Test > 0){
	   
		# Calculate R(k+1).
		R <- fR(ZigmaS, ZigmaR)
	   
		# Calculate Rho(k+1), Tau2(k+1) and Delta(k+1) for each center.
		for (i in 1:N.centres){
			Rho[i,1] <- fRho(i, R, XX_Matrixes)
			Tau2[i,1] <- fTau(i, ZigmaS, Rho)
			Delta[i,1] <- fDelta(i, Beta, Rho, XX_Matrixes, XY_Matrixes)
		}
		
		# Calculate Delta2(k+1) and name the rows and columns.
		Delta2 <- Delta*Delta
		rownames(Delta2) <- c(1:N.centres)
		colnames(Delta2) <- "Delta2"
	  
		# Calculate Beta(k+1),ZigmaS2(k+1) and ZigmaR2(k+1).
		New_Beta <- fBeta(N.centres, Delta, XX_Matrixes, XY_Matrixes)
		New_ZigmaS <- fZigmaS2(N.centres, Delta2, Tau2)
		New_ZigmaR <- fZigmaR2(N.centres, Beta, Delta, Delta2, Tau2, XY_Matrixes, YY_Matrixes, XX_Matrixes)
		
		# Count the number of differences that are higher than the confergence criteria.
		Test <- fControle(New_Beta, Beta, New_ZigmaS, ZigmaS, New_ZigmaR, ZigmaR, Conf)
		Diff_Beta <- abs(New_Beta - Beta) <= Conf
		Diff_ZigmaS2 <- abs(New_ZigmaS - ZigmaS) <= Conf
		Diff_ZigmaR2 <- abs(New_ZigmaR - ZigmaR) <= Conf
	  
		# Prepare ZigmaS2,ZigmaR2 and Beta for the next iteration by replacing the (k) value by the (k+1) value.
		ZigmaS <- New_ZigmaS
		ZigmaR <- New_ZigmaR
		Beta <- New_Beta
		# Name the rows and columns of Beta.
		rownames(Beta) <- c(0:(N.xvar))
		colnames(Beta) <- "Beta"
	  
		# Print something to see that the program is running / working.
		if(is.whole(Count / Milestone)){
			options(digits = 20)
			cat("The program is at iteration", Count, "\n")
		}
		
		Count <- Count + 1
	  
	}

	cat("Is the difference between k+1 and k for Beta0 lower than the convergence criteria? \n That's", Diff_Beta[1,1],"\n")
	cat("Is the difference between k+1 and k for Beta1 lower than the convergence criteria? \n That's", Diff_Beta[2,1],"\n")
	cat("Is the difference between k+1 and k for ZigmaS lower than the convergence criteria? \n That's", Diff_ZigmaS2,"\n")
	cat("Is the difference between k+1 and k for ZigmaR lower than the convergence criteria? \n That's", Diff_ZigmaR2,"\n")
	cat("Has the program reached the max number of iterations allowed? \n That's", Count>=Max,"\n")
	cat("The program stopped at iteration",Count-1,"\n")
	cat("The convergence criteria value is",Conf,"\n")
	cat("\n")
	cat("The start values are:","\n")
	for (i in 1:(N.xvar+1)){
		cat("Beta", i-1, "is", Start_Beta[i,1], "\n")
	}
	cat("ZigmaS is",Start_ZigmaS,"\n")
	cat("ZigmaR is",Start_ZigmaR,"\n")
	cat("\n")
	cat("The end values are:","\n")
	for (i in 1:(N.xvar+1)){
		cat("Beta",i-1,"is",Beta[i,1],"\n")
	}
	cat("ZigmaS is",ZigmaS,"\n")
	cat("ZigmaR is",ZigmaR,"\n")

	# Compute standard errors for beta's:

	# Create the XVX_i, i=1,..,N.centres matrices (p+1 x p+1) by defining list:
	# Nb: In the calculation-notes: "sigma^2"=ZigmaR (se residuals), "tau^2"=ZigmaS (se random intercept/study effect)
	XVX_matrix <- list()
	for (i in 1:N.centres){
		# XVX_i matrix for each study i=1,...,N.centres
		ni <- N.subs[[i]]                # number participants in the study
		# cat("ni is",ni,"\n")
		vi <- ZigmaR+ni*ZigmaS           # for easy notation
		# cat("vi is",vi,"\n")
		XVX_matrix[[i]] <- matrix(nrow = (N.xvar + 1), ncol = (N.xvar + 1))   # XVX_i is a (p+1) by (p+1) matrix
		XVX_matrix[[i]][1,1] <- ni * (vi - ni * ZigmaS) / (ZigmaR * vi)
		# cat("XVX_matrix [[",i,"]] is ",XVX_matrix[[i]][1,1],"\n")
	    for (j in 2:(N.xvar + 1)){ 	     # filling the remaining part of the first row and first column of matrix XVX_i
			XVX_matrix[[i]][1,j] <- (vi-ni * ZigmaS) * XX_Matrixes[[i]][1,j] / (ZigmaR * vi)
			XVX_matrix[[i]][j,1] <- XVX_matrix[[i]][1,j]  
			# cat("XVX_matrix [[1",j,"]] is ",XVX_matrix[[i]][1,j],"\n")
			# cat("XVX_matrix [[",j,",1]] is ",XVX_matrix[[i]][j,1],"\n")
		}
		for (j in 2:(N.xvar + 1)){       # filling remaining part of the XVX_i matrix
			for (k in j:(N.xvar + 1)){
				XVX_matrix[[i]][j,k] <- (vi * XX_Matrixes[[i]][j,k] - ZigmaS * XX_Matrixes[[i]][1,j] * XX_Matrixes[[i]][1,k]) / (ZigmaR * vi)
				XVX_matrix[[i]][k,j] <- XVX_matrix[[i]][j,k]  
				# cat("XVX_matrix [[",j, k,"]] is ",XVX_matrix[[i]][j,k],"\n")
				# cat("XVX_matrix [[",k ,j,"]] is ",XVX_matrix[[i]][k,j],"\n")
			} 
		}
	}

	# summing of the N.centres times different (N.xvar+1) times (N.xvar+1) XVX_matrices:
	XVX_total <- matrix(c(rep(0, (N.xvar + 1) * (N.xvar + 1))), nrow = (N.xvar + 1), ncol = (N.xvar + 1))
	for (i in 1:N.centres){
		XVX_total <- XVX_total + XVX_matrix[[i]]
	}

	# inverting this matrix results in the covariance matrix of Betas: 
	covbetas <- solve(XVX_total)
	sesbetas <- c(rep(0, (N.xvar + 1)))
	for (i in 1:(N.xvar + 1)){
		sesbetas[i] <- sqrt(covbetas[i,i]) # square root of the diagonal of the covariance matrix
		cat("The standard error for beta ", i-1, "is ", sesbetas[i], "\n")
	}

}



