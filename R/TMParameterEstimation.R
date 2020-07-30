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
		# suppressMessages(library(expm))
		# tmp <- expm(1/cell_cycle_times*logm(observation_matrix,method="Eigen"))
		return(as.vector(t(tmp)))
	}
}

################################################################################################
###################################     Library functions    ###################################
################################################################################################

#' @title "TMParameterEstimation"
#' @description Estimated the parameters that represent the probabilities of three active DNA methylation change types during n cell cycle(s).
#' @param observation_matrix The transition matrix(\eqn{5 \times 5}) from the original state to the terminal state. \cr
#' 	\tabular{cccccc}{
#' 		\tab original_class1 \tab original_class2 \tab original_class3 \tab original_class4 \tab original_class5\cr
#' 		terminal_class1 \tab a1 \tab b1 \tab c1 \tab d1 \tab e1\cr
#' 		terminal_class2 \tab a2 \tab b2 \tab c2 \tab d2 \tab e2\cr
#' 		terminal_class3 \tab a3 \tab b3 \tab c3 \tab d3 \tab e3\cr
#' 		terminal_class4 \tab a4 \tab b4 \tab c4 \tab d4 \tab e4\cr
#' 		terminal_class5 \tab a5 \tab b5 \tab c5 \tab d5 \tab e5\cr
#' 	}
#' 	# observation_matrix[1,1](a1) is the ratio of original_class1 to terminal_class1,observation_matrix[1,2](b1) is the ratio of original_class2 to terminal_class1 and so on\cr
#' 	# observation_matrix[2,1](a2) is the ratio of original_class1 to terminal_class2,observation_matrix[2,2](b2) is the ratio of original_class2 to terminal_class2 and so on\cr
#' 	# The sum of the ratio that original_class1 change to 5 classes should be 1. That is a1+a2+a3+a4+a5=1.\cr
#' @param iter The iteration times of the parameter estimation using the Newton-Raphson method with different initial guesses.
#' @param cell_cycle The cell cycle times from the original state to the terminal state.
#' @return estimated_parameters The estimated parameters using the maximum likelihood estimation and the Newton-Raphson method.
#' @return predicted_matrix The calculated transition matrix using the estimated parameters.
#' @details The transition matrix of this model describes the changes of DNA methylation during one cell cycle in three steps:
#' passive demethylation by DNA replication, active DNA methylation changes affected by DNA methylation-modifying enzymes
#' and DNA methylation combinations during homologous recombination.
#' For each CpG site in a chromsome, the methylation states are one of these four types: 0-0, 0-1, 1-0, 1-1.
#' The transition matrix after DNA replication would be :
#' 	\tabular{cccccc}{
#' 		\tab original(0-0) \tab original(0-1) \tab original(1-0) \tab original(1-1) \cr
#' 		after_replication(0-0) \tab \eqn{1} \tab \eqn{a} \tab \eqn{1-a} \tab \eqn{0} \cr
#' 		after_replication(0-1) \tab \eqn{0} \tab \eqn{1-a} \tab \eqn{0} \tab \eqn{1-a} \cr
#' 		after_replication(1-0) \tab \eqn{0} \tab \eqn{0} \tab \eqn{0} \tab \eqn{a} \cr
#' 		after_replication(1-1) \tab \eqn{0} \tab \eqn{0} \tab \eqn{0} \tab \eqn{0} \cr
#' 	}
#' 	among this matrix ,\eqn{a} is the methylation change probability with DNA replication and equal to 0.5.
#' 	Then the transition matrix after active DNA methylation changes would be:
#' 	\tabular{cccccc}{
#' 		\tab after_replication(0-0) \tab after_replication(0-1) \tab after_replication(1-0) \tab after_replication(1-1) \cr
#' 		after_enzymemodifying(0-0) \tab \eqn{(1-u) \times (1-u)} \tab \eqn{(1-u-p+u \times p) \times d} \tab \eqn{d \times (1-u-p+u \times p)} \tab \eqn{0} \cr
#' 		after_enzymemodifying(0-1) \tab \eqn{u \times (1-u)} \tab \eqn{(1-u-p+u \times p) \times (1-d)} \tab \eqn{d \times (u+p-u \times p)} \tab \eqn{0} \cr
#' 		after_enzymemodifying(1-0) \tab \eqn{u \times (1-u)} \tab \eqn{(u+p-u \times p) \times d} \tab \eqn{(1-d) \times (1-u-p+u \times p)} \tab \eqn{0} \cr
#' 		after_enzymemodifying(1-1) \tab \eqn{u \times u} \tab \eqn{(u+p-u \times p) \times (1-d)} \tab \eqn{(1-d) \times (u+p-u \times p)} \tab \eqn{1} \cr
#' 	}
#' 	The paramater \eqn{u} described the methylation probablity on CpG site.
#' 	The paramater \eqn{d} described the de-methylation probablity on 5mCpG site.
#' 	The paramater \eqn{p} described the methylation probablity on semi-CpG site.
#' 	Thus the mathlation state change is \eqn{P = Penzymemodifying \cdot Preplication \cdot original_classes}
#' 	We using \eqn{t_{i,j}} represent the \eqn{i} and \eqn{j} vector of the matrix \eqn{Penzymemodifying \cdot Preplication}
#' 	Then two chromsomes are combined during homologous recombination. The observed DNA methylation of each CpG site is the combination of the methylation types in both chromsomes.
#'  The observed transition matrix would be :
#' 	\tabular{cccccc}{
#' 		\tab original_class1(0) \tab original_class2(1/4) \tab original_class3(1/2) \tab original_class4(3/4) \tab original_class5(1)\cr
#' 		terminal_class1(0) \tab \eqn{x_{1,1}} \tab \eqn{x_{1,2}} \tab \eqn{x_{1,3}} \tab \eqn{x_{1,4}} \tab \eqn{x_{1,5}}\cr
#' 		terminal_class2(1/4) \tab \eqn{x_{2,1}} \tab \eqn{x_{2,2}} \tab \eqn{x_{2,3}} \tab \eqn{x_{2,4}} \tab \eqn{x_{2,5}}\cr
#' 		terminal_class3(1/2) \tab \eqn{x_{3,1}} \tab \eqn{x_{3,2}} \tab \eqn{x_{3,3}} \tab \eqn{x_{3,4}} \tab \eqn{x_{3,5}}\cr
#' 		terminal_class4(3/4) \tab \eqn{x_{4,1}} \tab \eqn{x_{4,2}} \tab \eqn{x_{4,3}} \tab \eqn{x_{4,4}} \tab \eqn{x_{4,5}}\cr
#' 		terminal_class5(1) \tab \eqn{x_{5,1}} \tab \eqn{x_{5,2}} \tab \eqn{x_{5,3}} \tab \eqn{x_{5,4}} \tab \eqn{x_{5,5}}\cr
#' 	}
#' 	and \deqn{x_{1,1}=t_{1,1} \times t_{1,1}}
#' \deqn{x_{1,1}=t_{{,1},1} \times t_{{,1},1}}
#' \deqn{x_{1,2}=1/4 \times (t_{1,1} \times t_{1,2}+t_{1,1} \times t_{1,3}+t_{1,2} \times t_{1,1}+t_{1,3} \times t_{1,1})}
#' \deqn{x_{1,3}=1/6 \times (t_{1,1} \times t_{1,4}+t_{1,2} \times t_{1,2}+t_{1,2} \times t_{1,3}+t_{1,3} \times t_{1,2}+t_{1,3} \times t_{1,3}+t_{1,4} \times t_{1,1})}
#' \deqn{x_{1,4}=1/4 \times (t_{1,2} \times t_{1,4}+t_{1,3} \times t_{1,4}+t_{1,4} \times t_{1,2}+t_{1,4} \times t_{1,3})}
#' \deqn{x_{1,5}=t_{1,4} \times t_{1,4}}
#' \deqn{x_{2,1}=t_{1,1} \times t_{2,1}+t_{1,1} \times t_{3,1}+t_{2,1} \times t_{1,1}+t_{3,1} \times t_{1,1}}
#' \deqn{x_{2,2}=1/4 \times (t_{1,1} \times t_{2,2}+t_{1,1} \times t_{2,3}+t_{1,2} \times t_{2,1}+t_{1,3} \times t_{2,1}+t_{1,1} \times t_{3,2}+t_{1,1} \times t_{3,3}+t_{1,2} \times t_{3,1}+t_{1,3} \times t_{3,1}+t_{2,1} \times t_{1,2}+t_{2,1} \times t_{1,3}+t_{2,2} \times t_{1,1}+t_{2,3} \times t_{1,1}+t_{3,1} \times t_{1,2}+t_{3,1} \times t_{1,3}+t_{3,2} \times t_{1,1}+t_{3,3} \times t_{1,1})}
#' \deqn{x_{2,3}=1/6 \times (t_{1,1} \times t_{2,4}+t_{1,2} \times t_{2,2}+t_{1,2} \times t_{2,3}+t_{1,3} \times t_{2,2}+t_{1,3} \times t_{2,3}+t_{1,4} \times t_{2,1}+t_{1,1} \times t_{3,4}+t_{1,2} \times t_{3,2}+t_{1,2} \times t_{3,3}+t_{1,3} \times t_{3,2}+t_{1,3} \times t_{3,3}+t_{1,4} \times t_{3,1}+t_{2,1} \times t_{1,4}+t_{2,2} \times t_{1,2}+t_{2,2} \times t_{1,3}+t_{2,3} \times t_{1,2}+t_{2,3} \times t_{1,3}+t_{2,4} \times t_{1,1}+t_{3,1} \times t_{1,4}+t_{3,2} \times t_{1,2}+t_{3,2} \times t_{1,3}+t_{3,3} \times t_{1,2}+t_{3,3} \times t_{1,3}+t_{3,4} \times t_{1,1})}
#' \deqn{x_{2,4}=1/4 \times (t_{1,2} \times t_{2,4}+t_{1,3} \times t_{2,4}+t_{1,4} \times t_{2,2}+t_{1,4} \times t_{2,3}+t_{1,2} \times t_{3,4}+t_{1,3} \times t_{3,4}+t_{1,4} \times t_{3,2}+t_{1,4} \times t_{3,3}+t_{2,2} \times t_{1,4}+t_{2,3} \times t_{1,4}+t_{2,4} \times t_{1,2}+t_{2,4} \times t_{1,3}+t_{3,2} \times t_{1,4}+t_{3,3} \times t_{1,4}+t_{3,4} \times t_{1,2}+t_{3,4} \times t_{1,3})}
#' \deqn{x_{2,5}=t_{1,4} \times t_{2,4}+t_{1,4} \times t_{3,4}+t_{2,4} \times t_{1,4}+t_{3,4} \times t_{1,4}}
#' \deqn{x_{3,1}=t_{1,1} \times t_{4,1}+t_{2,1} \times t_{2,1}+t_{2,1} \times t_{3,1}+t_{3,1} \times t_{2,1}+t_{3,1} \times t_{3,1}+t_{4,1} \times t_{1,1}}
#' \deqn{x_{3,2}=1/4 \times (t_{1,1} \times t_{4,2}+t_{1,1} \times t_{4,3}+t_{1,2} \times t_{4,1}+t_{1,3} \times t_{4,1}+t_{2,1} \times t_{2,2}+t_{2,1} \times t_{2,3}+t_{2,2} \times t_{2,1}+t_{2,3} \times t_{2,1}+t_{2,1} \times t_{3,2}+t_{2,1} \times t_{3,3}+t_{2,2} \times t_{3,1}+t_{2,3} \times t_{3,1}+t_{3,1} \times t_{2,2}+t_{3,1} \times t_{2,3}+t_{3,2} \times t_{2,1}+t_{3,3} \times t_{2,1}+t_{3,1} \times t_{3,2}+t_{3,1} \times t_{3,3}+t_{3,2} \times t_{3,1}+t_{3,3} \times t_{3,1}+t_{4,1} \times t_{1,2}+t_{4,1} \times t_{1,3}+t_{4,2} \times t_{1,1}+t_{4,3} \times t_{1,1})}
#' \deqn{x_{3,3}=1/6 \times (t_{1,1} \times t_{4,4}+t_{1,2} \times t_{4,2}+t_{1,2} \times t_{4,3}+t_{1,3} \times t_{4,2}+t_{1,3} \times t_{4,3}+t_{1,4} \times t_{4,1}+t_{2,1} \times t_{2,4}+t_{2,2} \times t_{2,2}+t_{2,2} \times t_{2,3}+t_{2,3} \times t_{2,2}+t_{2,3} \times t_{2,3}+t_{2,4} \times t_{2,1}+t_{2,1} \times t_{3,4}+t_{2,2} \times t_{3,2}+t_{2,2} \times t_{3,3}+t_{2,3} \times t_{3,2}+t_{2,3} \times t_{3,3}+t_{2,4} \times t_{3,1}+t_{3,1} \times t_{2,4}+t_{3,2} \times t_{2,2}+t_{3,2} \times t_{2,3}+t_{3,3} \times t_{2,2}+t_{3,3} \times t_{2,3}+t_{3,4} \times t_{2,1}+t_{3,1} \times t_{3,4}+t_{3,2} \times t_{3,2}+t_{3,2} \times t_{3,3}+t_{3,3} \times t_{3,2}+t_{3,3} \times t_{3,3}+t_{3,4} \times t_{3,1}+t_{4,1} \times t_{1,4}+t_{4,2} \times t_{1,2}+t_{4,2} \times t_{1,3}+t_{4,3} \times t_{1,2}+t_{4,3} \times t_{1,3}+t_{4,4} \times t_{1,1})}
#' \deqn{x_{3,4}=1/4 \times (t_{1,2} \times t_{4,4}+t_{1,3} \times t_{4,4}+t_{1,4} \times t_{4,2}+t_{1,4} \times t_{4,3}+t_{2,2} \times t_{2,4}+t_{2,3} \times t_{2,4}+t_{2,4} \times t_{2,2}+t_{2,4} \times t_{2,3}+t_{2,2} \times t_{3,4}+t_{2,3} \times t_{3,4}+t_{2,4} \times t_{3,2}+t_{2,4} \times t_{3,3}+t_{3,2} \times t_{2,4}+t_{3,3} \times t_{2,4}+t_{3,4} \times t_{2,2}+t_{3,4} \times t_{2,3}+t_{3,2} \times t_{3,4}+t_{3,3} \times t_{3,4}+t_{3,4} \times t_{3,2}+t_{3,4} \times t_{3,3}+t_{4,2} \times t_{1,4}+t_{4,3} \times t_{1,4}+t_{4,4} \times t_{1,2}+t_{4,4} \times t_{1,3})}
#' \deqn{x_{3,5}=t_{1,4} \times t_{4,4}+t_{2,4} \times t_{2,4}+t_{2,4} \times t_{3,4}+t_{3,4} \times t_{2,4}+t_{3,4} \times t_{3,4}+t_{4,4} \times t_{1,4}}
#' \deqn{x_{4,1}=t_{2,1} \times t_{4,1}+t_{3,1} \times t_{4,1}+t_{4,1} \times t_{2,1}+t_{4,1} \times t_{3,1}}
#' \deqn{x_{4,2}=1/4 \times (t_{2,1} \times t_{4,2}+t_{2,1} \times t_{4,3}+t_{2,2} \times t_{4,1}+t_{2,3} \times t_{4,1}+t_{3,1} \times t_{4,2}+t_{3,1} \times t_{4,3}+t_{3,2} \times t_{4,1}+t_{3,3} \times t_{4,1}+t_{4,1} \times t_{2,2}+t_{4,1} \times t_{2,3}+t_{4,2} \times t_{2,1}+t_{4,3} \times t_{2,1}+t_{4,1} \times t_{3,2}+t_{4,1} \times t_{3,3}+t_{4,2} \times t_{3,1}+t_{4,3} \times t_{3,1})}
#' \deqn{x_{4,3}=1/6 \times (t_{2,1} \times t_{4,4}+t_{2,2} \times t_{4,2}+t_{2,2} \times t_{4,3}+t_{2,3} \times t_{4,2}+t_{2,3} \times t_{4,3}+t_{2,4} \times t_{4,1}+t_{3,1} \times t_{4,4}+t_{3,2} \times t_{4,2}+t_{3,2} \times t_{4,3}+t_{3,3} \times t_{4,2}+t_{3,3} \times t_{4,3}+t_{3,4} \times t_{4,1}+t_{4,1} \times t_{2,4}+t_{4,2} \times t_{2,2}+t_{4,2} \times t_{2,3}+t_{4,3} \times t_{2,2}+t_{4,3} \times t_{2,3}+t_{4,4} \times t_{2,1}+t_{4,1} \times t_{3,4}+t_{4,2} \times t_{3,2}+t_{4,2} \times t_{3,3}+t_{4,3} \times t_{3,2}+t_{4,3} \times t_{3,3}+t_{4,4} \times t_{3,1})}
#' \deqn{x_{4,4}=1/4 \times (t_{2,2} \times t_{4,4}+t_{2,3} \times t_{4,4}+t_{2,4} \times t_{4,2}+t_{2,4} \times t_{4,3}+t_{3,2} \times t_{4,4}+t_{3,3} \times t_{4,4}+t_{3,4} \times t_{4,2}+t_{3,4} \times t_{4,3}+t_{4,2} \times t_{2,4}+t_{4,3} \times t_{2,4}+t_{4,4} \times t_{2,2}+t_{4,4} \times t_{2,3}+t_{4,2} \times t_{3,4}+t_{4,3} \times t_{3,4}+t_{4,4} \times t_{3,2}+t_{4,4} \times t_{3,3})}
#' \deqn{x_{4,5}=t_{2,4} \times t_{4,4}+t_{3,4} \times t_{4,4}+t_{4,4} \times t_{2,4}+t_{4,4} \times t_{3,4}}
#' \deqn{x_{5,1}=t_{4,1} \times t_{4,1}}
#' \deqn{x_{5,2}=1/4 \times (t_{4,1} \times t_{4,2}+t_{4,1} \times t_{4,3}+t_{4,2} \times t_{4,1}+t_{4,3} \times t_{4,1})}
#' \deqn{x_{5,3}=1/6 \times (t_{4,1} \times t_{4,4}+t_{4,2} \times t_{4,2}+t_{4,2} \times t_{4,3}+t_{4,3} \times t_{4,2}+t_{4,3} \times t_{4,3}+t_{4,4} \times t_{4,1})}
#' \deqn{x_{5,4}=1/4 \times (t_{4,2} \times t_{4,4}+t_{4,3} \times t_{4,4}+t_{4,4} \times t_{4,2}+t_{4,4} \times t_{4,3})}
#' \deqn{x_{5,5}=t_{4,4} \times t_{4,4}}
#' The cost function was defined by \deqn{f_{cost}=\sum_{i=1,j=1}^{n=5} (o_{i,j}-x_{i,j})} and minimized using the Newton-Raphson method.
#' @references \cite{Zhao, C. et.al.(2018). A DNA methylation state transition model reveals the programmed epigenetic heterogeneity in pre-implantation embryos. Under revision.}
#' @examples # Let's start from a transtion matrix
#' observation_matrix <- matrix(c(0.9388,0.0952,0.0377,0,0.0001,
#'                                0.0497,0.5873,0.1887,0.0344,0.0149,
#'                                0.0096,0.2381,0.4151,0.0653,0.0876,
#'                                0.0019,0.0635,0.3396,0.6151,0.253,
#'                                0,0.0159,0.0189,0.2852,0.6444),5,5,byrow=T)
#'
#' # The DNA methylation states change from original state to the terminal state after 1 cell cycle.
#' TMParameterEstimation(observation_matrix,iter=30,cell_cycle=1)
#'
#' # The DNA methylation states change from original state to the terminal state after 2 cell cycle.
#' # TMParameterEstimation(observation_matrix,iter=1,cell_cycle=2)
#' # if this function was not successful to estimated the function, you may try more iterations with different initial guesses
#' TMParameterEstimation(observation_matrix,iter=50,cell_cycle=2)
#'
#' # The DNA methylation states change from original state to the terminal state after 30 cell cycle.
#' TMParameterEstimation(observation_matrix,iter=50,cell_cycle=30)
#' @concept MethylTransition
#' @export TMParameterEstimation

TMParameterEstimation <- function(observation_matrix,iter=50,cell_cycle=1){
	stopifnot(is.numeric(iter), is.numeric(cell_cycle))
	if (!is.matrix(observation_matrix)){
		stop("The input observation transition ratio should be a matrix.")
	}
	if (!is.numeric(observation_matrix)){
		stop("The input observation transition ratio should be numeric number.")
	}
	if (ncol(observation_matrix) != 5 | nrow(observation_matrix) != 5){
		stop("The input transition matrix should be a \"5 \times 5\" matrix.\n\n")
	}
	if (cell_cycle <= 0){
		stop("Please selected a right number of cell cycles!\n\n")
	}
	if (iter <= 0){
		stop("Please selected a right number of iteration!\n\n")
	}
	if(abs(sum(observation_matrix[,1])-1) >= 1e5){
		warning("The total ratio of the 1st class(the coloum1) is not equal to 1, the estimation result may not accurate.")
	}
	if(abs(sum(observation_matrix[,2])-1) >= 1e5){
		warning("The total ratio of the 2nd class(the coloum2) is not equal to 1, the estimation result may not accurate.")
	}
	if(abs(sum(observation_matrix[,3])-1) >= 1e5){
		warning("The total ratio of the 3rd class(the coloum3) is not equal to 1, the estimation result may not accurate.")
	}
	if(abs(sum(observation_matrix[,4])-1) >= 1e5){
		warning("The total ratio of the 4th class(the coloum4) is not equal to 1, the estimation result may not accurate.")
	}
	if(abs(sum(observation_matrix[,5])-1) >= 1e5){
		warning("The total ratio of the 5th class(the coloum5) is not equal to 1, the estimation result may not accurate.")
	}
	if(sum(observation_matrix>1 | observation_matrix<0)>0){
		warning("There are ratios that are not in the range of zero to one, the estimation result may not accurate.")
	}
	observe <- TransitionMatrixCellCycle(observation_matrix,cell_cycle)
	para_maxtrix <- c()
	cat("\n")
	cat(paste("Parameter estimation for ",round(cell_cycle)," cell cycle(s) is running ...\n",sep=""))
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
		stop(paste("TMParameterEstimation failed to estimated the parameters in ",iter," iterations. Please try more iterations with different initial guesses by choose a bigger 'iter'",sep=""))
	}else{
		cost_value <- round(apply(selected_parameter_matrix,1,CostFunction,observe),12)
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
	rownames(predicted_matrix) <- c("terminal_class1","terminal_class2","terminal_class3","terminal_class4","terminal_class5")
	cat("\n")
	return(list("estimated_parameters"=selected_parameter,"predicted_matrix"=predicted_matrix))
}
