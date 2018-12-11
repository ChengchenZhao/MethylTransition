################################################################################################
###################################    Built-in functions    ###################################
################################################################################################
NewtonsMethod <- function(func = objfun, x0, tol = 1e-5, n.max = 50,...){
	suppressMessages(library(numDeriv))
	x <- x0
	g <- grad(func, x, ...)
	h <- hessian(func, x, ...)

	n <- 0
	while( max(abs(g))>tol && n<n.max ){
		x <- x-solve(h,g)
		g <- grad(func, x, ...)
		h <- hessian(func, x, ...)
		n <- n+1
	}
	if(n == n.max){
		# cat('NewtonsMethod failed to converge\n')
		return(x)
	}
	return(x)
}

CostFunction <- function(x,observe){
	# x[1] : u
	# x[2] : d
	# x[3] : p
	predict <- c((-x[1] + 1)**4,
				 0.5*(-x[1] + 1)**2*(0.5*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.5*(-x[1] + 1)**2) + 0.25*(-x[1] + 1)**2*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2),
				 (1/3)*x[2]*(-x[3] + 1)*(-x[1] + 1)**3 + (2/3)*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2),
				 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2),
				 1.0*x[2]**2*(-x[3] + 1)**2*(-x[1] + 1)**2,
				 4*x[1]*(-x[1] + 1)**3,
				 1.0*x[1]*(-x[1] + 1)*(0.5*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.5*(-x[1] + 1)**2) + 0.5*x[1]*(-x[1] + 1)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + 0.5*(-x[1] + 1)**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*x[1]*(-x[1] + 1)) + 0.25*(-x[1] + 1)**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 0.5*(-x[1] + 1)**2*(0.5*x[1]*(-x[1] + 1) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.25*(-x[1] + 1)**2*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 (2/3)*x[2]*x[1]*(-x[3] + 1)*(-x[1] + 1)**2 + (1/3)*(-x[1] + 1)**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/6)*(-x[1] + 1)**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + (2/3)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2) + (2/3)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + (2/3)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2),
				 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1)) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + 1.0*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2),
				 4*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 6*x[1]**2*(-x[1] + 1)**2,
				 0.5*x[1]**2*(0.5*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.5*(-x[1] + 1)**2) + 0.25*x[1]**2*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + 1.0*x[1]*(-x[1] + 1)*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*x[1]*(-x[1] + 1)) + 0.5*x[1]*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 1.0*x[1]*(-x[1] + 1)*(0.5*x[1]*(-x[1] + 1) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.5*x[1]*(-x[1] + 1)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.5*(-x[1] + 1)**2*(0.5*x[1]**2 + 0.5*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.25*(-x[1] + 1)**2*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
				 (1/3)*x[2]*x[1]**2*(-x[3] + 1)*(-x[1] + 1) + (2/3)*x[1]*(-x[1] + 1)*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/3)*x[1]*(-x[1] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/3)*(-x[2] + 1)*(-x[1] + 1)**2*(-x[3]*x[1] + x[3] + x[1]) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2) + (2/3)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2) + (2/3)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + (2/3)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1))*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 1.0*x[2]*(-x[3] + 1)*(-x[1] + 1)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.25*x[2]*(-x[3] + 1)*(-x[1] + 1)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 1.0*(-x[2] + 1)*(0.25*x[2]*(-x[3] + 1)*(-x[1] + 1) + 0.25*(-x[1] + 1)**2)*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[2]*(-x[3] + 1)*(-x[1] + 1) + (-x[1] + 1)**2)*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 1.0*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 2.0*x[2]*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)*(-x[3]*x[1] + x[3] + x[1]) + 4*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 4*x[1]**3*(-x[1] + 1),
				 0.5*x[1]**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*x[1]*(-x[1] + 1)) + 0.25*x[1]**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + 0.5*x[1]**2*(0.5*x[1]*(-x[1] + 1) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 0.25*x[1]**2*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*x[1]*(-x[1] + 1)*(0.5*x[1]**2 + 0.5*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.5*x[1]*(-x[1] + 1)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
				 (1/3)*x[1]**2*(0.5*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.5*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (1/6)*x[1]**2*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*x[1]*(-x[2] + 1)*(-x[1] + 1)*(-x[3]*x[1] + x[3] + x[1]) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1)) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + (2/3)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1)) + (2/3)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 1.0*(-x[2] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*x[1]*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + x[1]*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(-x[2] + 1)*(0.25*x[1]*(-x[1] + 1) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[1]*(-x[1] + 1) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)) + 1.0*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1)),
				 4*(-x[2] + 1)*(0.25*x[2]*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]) + 1.0*(-x[2] + 1)*(x[2]*(-x[3]*x[1] + x[3] + x[1]) + (-x[2] + 1)*(-x[3] + 1)*(-x[1] + 1))*(-x[3]*x[1] + x[3] + x[1]),
				 x[1]**4,
				 0.5*x[1]**2*(0.5*x[1]**2 + 0.5*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])) + 0.25*x[1]**2*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
				 (1/3)*x[1]**2*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]) + (2/3)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1])),
				 1.0*(-x[2] + 1)*(0.25*x[1]**2 + 0.25*(-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(-x[3]*x[1] + x[3] + x[1]) + 0.25*(-x[2] + 1)*(x[1]**2 + (-x[2] + 1)*(-x[3]*x[1] + x[3] + x[1]))*(-x[3]*x[1] + x[3] + x[1]),
				 1.0*(-x[2] + 1)**2*(-x[3]*x[1] + x[3] + x[1])**2)
	y <- 0
	for (each in seq(length(observe))){
		y <- y + (predict[each]-observe[each])^2

	}
	return (y)
}

TransitionMatrixCellCycle <- function(observation_matrix,cell_cycle_times){
	if (cell_cycle_times==1){
		return(as.vector(t(observation_matrix)))
	}else{
		suppressMessages(library(pracma))
		tmp <- expm(1/cell_cycle_times*logm(observation_matrix))
		return(as.vector(t(tmp)))
	}
}

################################################################################################
###################################     Library functions    ###################################
################################################################################################

#' @title "ParameterEstimation"
#' @description Estimated the parameters that represent the probabilities of three active DNA methylation change types during n cell cycle(s).
#' @param observation_matrix The transition matrix(5X5) from the original state to the terminational state. \cr
#' 	\tabular{cccccc}{
#' 		\tab original_class1 \tab original_class2 \tab original_class3 \tab original_class4 \tab original_class5\cr
#' 		terminational_class1 \tab a1 \tab b1 \tab c1 \tab d1 \tab e1\cr
#' 		terminational_class2 \tab a2 \tab b2 \tab c2 \tab d2 \tab e2\cr
#' 		terminational_class3 \tab a3 \tab b3 \tab c3 \tab d3 \tab e3\cr
#' 		terminational_class4 \tab a4 \tab b4 \tab c4 \tab d4 \tab e4\cr
#' 		terminational_class5 \tab a5 \tab b5 \tab c5 \tab d5 \tab e5
#' 	}\cr
#' 	# observation_matrix[1,1](a1) is the ratio of original_class1 to terminational_class1,observation_matrix[1,2](b1) is the ratio of original_class2 to terminational_class1 and so on\cr
#' 	# observation_matrix[2,1](a2) is the ratio of original_class1 to terminational_class2,observation_matrix[2,2](b2) is the ratio of original_class2 to terminational_class2 and so on\cr
#' 	# The sum of the ratio that original_class1 change to 5 classes should be 1. That is a1+a2+a3+a4+a5=1.\cr
#' @param iter The iteration times of the parameter estimation using the Newton-Raphson method with different initial guesses.
#' @param cell_cycle The cell cycle times from the original state to the terminational state.
#' @return estimated_parameters The estimated parameters using the maximum likelihood estimation and the Newton-Raphson method.
#' @return predicted_matrix The calculated transition matrix using the estimated parameters.
#' @details The transition matrix of this model describes the changes of DNA methylation during one cell cycle in three steps:
#' passive demethylation by DNA replication, active DNA methylation changes affected by DNA methylation-modifying enzymes
#' and DNA methylation combinations during homologous recombination.
#' @references .
#' @examples # Let's start from a transtion matrix
#' observation_matrix <- matrix(c(0.9388,0.0952,0.0377,0,0,0.0497,0.5873,0.1887,0.0344,0.0149,0.0096,0.2381,0.4151,0.0653,0.0876,0.0019,0.0635,0.3396,0.6151,0.253,0,0.0159,0.0189,0.2852,0.6444),5,5,byrow=T)
#'
#' # The DNA methylation states change from original state to the terminational state after 1 cell cycle.
#' ParameterEstimation(observation_matrix,iter=30,cell_cycle=1)
#'
#' # The DNA methylation states change from original state to the terminational state after 2 cell cycle.
#' # ParameterEstimation(observation_matrix,iter=1,cell_cycle=2)
#' # if this function was not successful to estimated the function, you may try more iterations with different initial guesses
#' ParameterEstimation(observation_matrix,iter=50,cell_cycle=2)
#'
#' # The DNA methylation states change from original state to the terminational state after 30 cell cycle.
#' ParameterEstimation(observation_matrix,iter=50,cell_cycle=30)
#' @export ParameterEstimation

ParameterEstimation <- function(observation_matrix,iter=50,cell_cycle=1){
	stopifnot(is.numeric(observation_matrix), is.matrix(observation_matrix))
	stopifnot(is.numeric(iter), is.numeric(cell_cycle))
	observe <- TransitionMatrixCellCycle(observation_matrix,cell_cycle)
	para_maxtrix <- c()
	cat("\n")
	cat("\tParameter Estimation\n")
	cat("\n")
	cat(paste("ParameterEstimation for ",cell_cycle," cell cycle(s) is running ...\n",sep=""))
	for (i in seq(iter)){
		if (i%%10==0){
			cat(paste("\t",i,"/",iter," iterations...\n",sep=""))
		}
		x0 <- runif(3, min=0, max=1)
		tryCatch({y <- NewtonsMethod(func=CostFunction, x0=x0,observe=observe);
		para_maxtrix <- rbind(para_maxtrix,y);},
		error=function(e){i})
	}
	selected_parameter_matrix <- para_maxtrix[which(para_maxtrix[,1]<=1 & para_maxtrix[,1]>=0 & para_maxtrix[,2]<=1 & para_maxtrix[,2]>=0 & para_maxtrix[,3]<=1 & para_maxtrix[,3]>=0),]
	if(is.null(nrow(selected_parameter_matrix))){
		selected_parameter <- round(selected_parameter_matrix,6)
	}else if(nrow(selected_parameter_matrix)==0){
		stop(paste("ParameterEstimation failed to estimated the parameters in ",iter," iterations. Please try more iterations with different initial guesses by choose a bigger 'iter'",sep=""))
	}else{
		cost_value <- round(apply(selected_parameter_matrix,1,CostFunction,observe),6)
		optimal_para <- selected_parameter_matrix[which(cost_value==min(cost_value)),]
		if (!is.null(nrow(optimal_para))){
			selected_parameter <- round(apply(optimal_para,2,mean,na.rm=T),6)
		}else{
			selected_parameter <- round(optimal_para,6)
		}
	}
	names(selected_parameter) <- c("DeNovoMethylation_u","DeMethylation_d","SemiMethylation_p")
	predicted_matrix <- MethylationTransMatrix(selected_parameter[1],selected_parameter[2],selected_parameter[3])
	colnames(predicted_matrix) <- c("original_class1","original_class2","original_class3","original_class4","original_class5")
	rownames(predicted_matrix) <- c("terminational_class1","terminational_class2","terminational_class3","terminational_class4","terminational_class5")
	cat("\n")
	return(list("estimated_parameters"=selected_parameter,"predicted_matrix"=predicted_matrix))
}
