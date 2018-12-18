################################################################################################
###################################    Built-in functions    ###################################
################################################################################################
MethylationClassSummary <- function(x){
	cutoff_1 <- 1/8
	cutoff_2 <- 3/8
	cutoff_3 <- 5/8
	cutoff_4 <- 7/8
	a <- sum(x<=cutoff_1,na.rm=T)
	b <- sum(x>cutoff_1 & x<=cutoff_2,na.rm=T)
	c <- sum(x>cutoff_2 & x<=cutoff_3,na.rm=T)
	d <- sum(x>cutoff_3 & x<=cutoff_4,na.rm=T)
	e <- sum(x>cutoff_4,na.rm=T)
	total <- a+b+c+d+e
	return(round(c(a/total,b/total,c/total,d/total,e/total),4))
	# return(c(a/total,b/total,c/total,d/total,e/total))
}

MethylationTransMatrix <- function(u,d,p){
	#	0-0 0-1 1-0 1-1
	# 0-0 (1-u)(1-u)  (1-u)(1-p)d d(1-u)(1-p) 0
	# 0-1 u(1-u)  (1-p)(1-u)(1-d) d(u+p-up) 0
	# 1-0 u(1-u)  (u+p-up)d (1-d)(1-u)(1-p) 0
	# 1-1 uu  (u+p-up)(1-d) (1-d)(u+p-up) 1

	Q00_00 <- (1-u)*(1-u); Q01_00 <- (1-u)*(1-p)*d; Q10_00 <- d*(1-u)*(1-p) ; Q11_00 <- 0
	Q00_01 <- u*(1-u); Q01_01 <- (1-p)*(1-u)*(1-d); Q10_01 <- d*(u+p-u*p) ; Q11_01 <- 0
	Q00_10 <- u*(1-u); Q01_10 <- (u+p-u*p)*d; Q10_10 <- (1-d)*(1-u)*(1-p) ; Q11_10 <- 0
	Q00_11 <- u*u; Q01_11 <- (u+p-u*p)*(1-d); Q10_11 <- (1-d)*(u+p-u*p) ; Q11_11 <- 1
	# Q <- matrix(c(Q00_00,Q01_00,Q10_00,Q11_00,Q00_01,Q01_01,Q10_01,Q11_01,Q00_10,Q01_10,Q10_10,Q11_10,Q00_11,Q01_11,Q10_11,Q11_11),c(4,4),byrow = T)
	# print(Q)
	x00_00 <- Q00_00; x01_00 <- 1/2*(Q00_00+Q01_00); x10_00 <- 1/2*(Q00_00+Q10_00); x11_00 <- 1/2*(Q01_00+Q10_00)
	x00_01 <- Q00_01; x01_01 <- 1/2*(Q00_01+Q01_01); x10_01 <- 1/2*(Q00_01+Q10_01); x11_01 <- 1/2*(Q01_01+Q10_01)
	x00_10 <- Q00_10; x01_10 <- 1/2*(Q00_10+Q01_10); x10_10 <- 1/2*(Q00_10+Q10_10); x11_10 <- 1/2*(Q01_10+Q10_10)
	x00_11 <- Q00_11; x01_11 <- 1/2*(Q00_11+Q01_11); x10_11 <- 1/2*(Q00_11+Q10_11); x11_11 <- 1/2*(Q01_11+Q10_11)
	# x <- matrix(c(x00_00,x01_00,x10_00,x11_00,x00_01,x01_01,x10_01,x11_01,x00_10,x01_10,x10_10,x11_10,x00_11,x01_11,x10_11,x11_11),c(4,4),byrow = T)
	# print(x)
	tmp_matrix <- c(x00_00*x00_00,1/4*(x00_00*x01_00+x00_00*x10_00+x01_00*x00_00+x10_00*x00_00),1/6*(x00_00*x11_00+x01_00*x01_00+x01_00*x10_00+x10_00*x01_00+x10_00*x10_00+x11_00*x00_00),1/4*(x01_00*x11_00+x10_00*x11_00+x11_00*x01_00+x11_00*x10_00),x11_00*x11_00,x00_00*x00_01+x00_00*x00_10+x00_01*x00_00+x00_10*x00_00,1/4*(x00_00*x01_01+x00_00*x10_01+x01_00*x00_01+x10_00*x00_01+x00_00*x01_10+x00_00*x10_10+x01_00*x00_10+x10_00*x00_10+x00_01*x01_00+x00_01*x10_00+x01_01*x00_00+x10_01*x00_00+x00_10*x01_00+x00_10*x10_00+x01_10*x00_00+x10_10*x00_00),1/6*(x00_00*x11_01+x01_00*x01_01+x01_00*x10_01+x10_00*x01_01+x10_00*x10_01+x11_00*x00_01+x00_00*x11_10+x01_00*x01_10+x01_00*x10_10+x10_00*x01_10+x10_00*x10_10+x11_00*x00_10+x00_01*x11_00+x01_01*x01_00+x01_01*x10_00+x10_01*x01_00+x10_01*x10_00+x11_01*x00_00+x00_10*x11_00+x01_10*x01_00+x01_10*x10_00+x10_10*x01_00+x10_10*x10_00+x11_10*x00_00),1/4*(x01_00*x11_01+x10_00*x11_01+x11_00*x01_01+x11_00*x10_01+x01_00*x11_10+x10_00*x11_10+x11_00*x01_10+x11_00*x10_10+x01_01*x11_00+x10_01*x11_00+x11_01*x01_00+x11_01*x10_00+x01_10*x11_00+x10_10*x11_00+x11_10*x01_00+x11_10*x10_00),x11_00*x11_01+x11_00*x11_10+x11_01*x11_00+x11_10*x11_00,x00_00*x00_11+x00_01*x00_01+x00_01*x00_10+x00_10*x00_01+x00_10*x00_10+x00_11*x00_00,1/4*(x00_00*x01_11+x00_00*x10_11+x01_00*x00_11+x10_00*x00_11+x00_01*x01_01+x00_01*x10_01+x01_01*x00_01+x10_01*x00_01+x00_01*x01_10+x00_01*x10_10+x01_01*x00_10+x10_01*x00_10+x00_10*x01_01+x00_10*x10_01+x01_10*x00_01+x10_10*x00_01+x00_10*x01_10+x00_10*x10_10+x01_10*x00_10+x10_10*x00_10+x00_11*x01_00+x00_11*x10_00+x01_11*x00_00+x10_11*x00_00),1/6*(x00_00*x11_11+x01_00*x01_11+x01_00*x10_11+x10_00*x01_11+x10_00*x10_11+x11_00*x00_11+x00_01*x11_01+x01_01*x01_01+x01_01*x10_01+x10_01*x01_01+x10_01*x10_01+x11_01*x00_01+x00_01*x11_10+x01_01*x01_10+x01_01*x10_10+x10_01*x01_10+x10_01*x10_10+x11_01*x00_10+x00_10*x11_01+x01_10*x01_01+x01_10*x10_01+x10_10*x01_01+x10_10*x10_01+x11_10*x00_01+x00_10*x11_10+x01_10*x01_10+x01_10*x10_10+x10_10*x01_10+x10_10*x10_10+x11_10*x00_10+x00_11*x11_00+x01_11*x01_00+x01_11*x10_00+x10_11*x01_00+x10_11*x10_00+x11_11*x00_00),1/4*(x01_00*x11_11+x10_00*x11_11+x11_00*x01_11+x11_00*x10_11+x01_01*x11_01+x10_01*x11_01+x11_01*x01_01+x11_01*x10_01+x01_01*x11_10+x10_01*x11_10+x11_01*x01_10+x11_01*x10_10+x01_10*x11_01+x10_10*x11_01+x11_10*x01_01+x11_10*x10_01+x01_10*x11_10+x10_10*x11_10+x11_10*x01_10+x11_10*x10_10+x01_11*x11_00+x10_11*x11_00+x11_11*x01_00+x11_11*x10_00),x11_00*x11_11+x11_01*x11_01+x11_01*x11_10+x11_10*x11_01+x11_10*x11_10+x11_11*x11_00,x00_01*x00_11+x00_10*x00_11+x00_11*x00_01+x00_11*x00_10,1/4*(x00_01*x01_11+x00_01*x10_11+x01_01*x00_11+x10_01*x00_11+x00_10*x01_11+x00_10*x10_11+x01_10*x00_11+x10_10*x00_11+x00_11*x01_01+x00_11*x10_01+x01_11*x00_01+x10_11*x00_01+x00_11*x01_10+x00_11*x10_10+x01_11*x00_10+x10_11*x00_10),1/6*(x00_01*x11_11+x01_01*x01_11+x01_01*x10_11+x10_01*x01_11+x10_01*x10_11+x11_01*x00_11+x00_10*x11_11+x01_10*x01_11+x01_10*x10_11+x10_10*x01_11+x10_10*x10_11+x11_10*x00_11+x00_11*x11_01+x01_11*x01_01+x01_11*x10_01+x10_11*x01_01+x10_11*x10_01+x11_11*x00_01+x00_11*x11_10+x01_11*x01_10+x01_11*x10_10+x10_11*x01_10+x10_11*x10_10+x11_11*x00_10),1/4*(x01_01*x11_11+x10_01*x11_11+x11_01*x01_11+x11_01*x10_11+x01_10*x11_11+x10_10*x11_11+x11_10*x01_11+x11_10*x10_11+x01_11*x11_01+x10_11*x11_01+x11_11*x01_01+x11_11*x10_01+x01_11*x11_10+x10_11*x11_10+x11_11*x01_10+x11_11*x10_10),x11_01*x11_11+x11_10*x11_11+x11_11*x11_01+x11_11*x11_10,x00_11*x00_11,1/4*(x00_11*x01_11+x00_11*x10_11+x01_11*x00_11+x10_11*x00_11),1/6*(x00_11*x11_11+x01_11*x01_11+x01_11*x10_11+x10_11*x01_11+x10_11*x10_11+x11_11*x00_11),1/4*(x01_11*x11_11+x10_11*x11_11+x11_11*x01_11+x11_11*x10_11),x11_11*x11_11)
	return(matrix(tmp_matrix,c(5,5),byrow = T))
}

PredictionMethylationClass <- function(start,u,d,p){
	t <- MethylationTransMatrix(u,d,p)
	predict <- t%*%start
	return(as.vector(predict))
}

# PredictionMethylationClass <- function(start,u,d,p,total,loops){
# 	start_number <- start*total
# 	t <- MethylationTransMatrix(u,d,p)
# 	t <- t*total
# 	t1 <- c(rep(0,t[1,1]),rep(0.25,t[2,1]),rep(0.5,t[3,1]),rep(0.75,t[4,1]),rep(1,t[5,1]))
# 	t2 <- c(rep(0,t[1,2]),rep(0.25,t[2,2]),rep(0.5,t[3,2]),rep(0.75,t[4,2]),rep(1,t[5,2]))
# 	t3 <- c(rep(0,t[1,3]),rep(0.25,t[2,3]),rep(0.5,t[3,3]),rep(0.75,t[4,3]),rep(1,t[5,3]))
# 	t4 <- c(rep(0,t[1,4]),rep(0.25,t[2,4]),rep(0.5,t[3,4]),rep(0.75,t[4,4]),rep(1,t[5,4]))
# 	t5 <- c(rep(0,t[1,5]),rep(0.25,t[2,5]),rep(0.5,t[3,5]),rep(0.75,t[4,5]),rep(1,t[5,5]))
# 	predict_loop_result <- c()
# 	for (i in seq(loops)){
# 		tmp_predict <- c()
# 		tmp_predict <- c(tmp_predict,sample(t1,start_number[1],replace=T))
# 		tmp_predict <- c(tmp_predict,sample(t2,start_number[2],replace=T))
# 		tmp_predict <- c(tmp_predict,sample(t3,start_number[3],replace=T))
# 		tmp_predict <- c(tmp_predict,sample(t4,start_number[4],replace=T))
# 		tmp_predict <- c(tmp_predict,sample(t5,start_number[5],replace=T))
# 		predict_loop_result <- rbind(predict_loop_result,MethylationClassSummary(tmp_predict))
# 	}
# 	if (loops>=2){
# 		predict <- apply(predict_loop_result,2,mean)
# 	}else{
# 		predict <- predict_loop_result
# 	}
# 	return(as.vector(round(predict,4)))
# }

MethyClass2Mean <- function(x){
	if (length(x)!=5){
		print("The methylation level should be classified to FIVE classes(0,1/4,1/2,3/4,1).")
	}else{
		return(x[1]*0+x[2]*0.25+x[3]*0.5+x[4]*0.75+x[5]*1)
	}
}

################################################################################################
###################################     Library functions    ###################################
################################################################################################

#' @title "MethylCalculation"
#' @description Calculate the ratio of each methylation class after n cell cycle(s).
#' @param original_classes The original methylation classes
#' @param u The paramater that describing the methylation probablity on CpG site
#' @param d The paramater that describing the de-methylation probablity on 5mCpG site
#' @param p The paramater that describing the methylation probablity on semi-CpG site
#' @param total The total genes when do the calculation
#' @param loops The loop times for the calculation
#' @param cell_cycle The cell cycle times
#' @return terminational_classes The terminational methylation classes
#' @return average_methylation_level The average methylation level after [n] cell cycle(s).
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
#' 	}\cr
#' 	among this matrix ,\eqn{a} is the methylation change probability with DNA replication and equal to 0.5.
#' 	Then the transition matrix after active DNA methylation changes would be:
#' 	\tabular{cccccc}{
#' 		\tab after_replication(0-0) \tab after_replication(0-1) \tab after_replication(1-0) \tab after_replication(1-1) \cr
#' 		after_enzymemodifying(0-0) \tab \eqn{(1-u)×(1-u)} \tab \eqn{(1-u-p+u×p)×d} \tab \eqn{d×(1-u-p+u×p)} \tab \eqn{0} \cr
#' 		after_enzymemodifying(0-1) \tab \eqn{u×(1-u)} \tab \eqn{(1-u-p+u×p)×(1-d)} \tab \eqn{d×(u+p-u×p)} \tab \eqn{0} \cr
#' 		after_enzymemodifying(1-0) \tab \eqn{u×(1-u)} \tab \eqn{(u+p-u×p)×d} \tab \eqn{(1-d)×(1-u-p+u×p)} \tab \eqn{0} \cr
#' 		after_enzymemodifying(1-1) \tab \eqn{u×u} \tab \eqn{(u+p-u×p)×(1-d)} \tab \eqn{(1-d)×(u+p-u×p)} \tab \eqn{1} \cr
#' 	}\cr
#' 	The paramater \eqn{u} described the methylation probablity on CpG site.
#' 	The paramater \eqn{d} described the de-methylation probablity on 5mCpG site.
#' 	The paramater \eqn{p} described the methylation probablity on semi-CpG site.
#' 	Thus the mathlation state change is \eqn{P = Penzymemodifying \cdot Preplication \cdot original_classes}
#' 	We using \eqn{t_{i,j}} represent the \eqn{i} and \eqn{j} vector of the matrix \eqn{Penzymemodifying \cdot Preplication}
#' 	Then two chromsomes are combined during homologous recombination. The observed DNA methylation of each CpG site is the combination of the methylation types in both chromsomes.
#'  The observed transition matrix would be :
#' 	\tabular{cccccc}{
#' 		\tab original_class1(0) \tab original_class2(1/4) \tab original_class3(1/2) \tab original_class4(3/4) \tab original_class5(1)\cr
#' 		terminational_class1(0) \tab \eqn{x_{1,1}} \tab \eqn{x_{1,2}} \tab \eqn{x_{1,3}} \tab \eqn{x_{1,4}} \tab \eqn{x_{1,5}}\cr
#' 		terminational_class2(1/4) \tab \eqn{x_{2,1}} \tab \eqn{x_{2,2}} \tab \eqn{x_{2,3}} \tab \eqn{x_{2,4}} \tab \eqn{x_{2,5}}\cr
#' 		terminational_class3(1/2) \tab \eqn{x_{3,1}} \tab \eqn{x_{3,2}} \tab \eqn{x_{3,3}} \tab \eqn{x_{3,4}} \tab \eqn{x_{3,5}}\cr
#' 		terminational_class4(3/4) \tab \eqn{x_{4,1}} \tab \eqn{x_{4,2}} \tab \eqn{x_{4,3}} \tab \eqn{x_{4,4}} \tab \eqn{x_{4,5}}\cr
#' 		terminational_class5(1) \tab \eqn{x_{5,1}} \tab \eqn{x_{5,2}} \tab \eqn{x_{5,3}} \tab \eqn{x_{5,4}} \tab \eqn{x_{5,5}}
#' 	}\cr
#' 	and \deqn{x_{1,1}=t_{1,1}×t_{1,1}}
#' \deqn{x_{1,1}=t_{{,1},1}×t_{{,1},1}}
#' \deqn{x_{1,2}=1/4×(t_{1,1}×t_{1,2}+t_{1,1}×t_{1,3}+t_{1,2}×t_{1,1}+t_{1,3}×t_{1,1})}
#' \deqn{x_{1,3}=1/6×(t_{1,1}×t_{1,4}+t_{1,2}×t_{1,2}+t_{1,2}×t_{1,3}+t_{1,3}×t_{1,2}+t_{1,3}×t_{1,3}+t_{1,4}×t_{1,1})}
#' \deqn{x_{1,4}=1/4×(t_{1,2}×t_{1,4}+t_{1,3}×t_{1,4}+t_{1,4}×t_{1,2}+t_{1,4}×t_{1,3})}
#' \deqn{x_{1,5}=t_{1,4}×t_{1,4}}
#' \deqn{x_{2,1}=t_{1,1}×t_{2,1}+t_{1,1}×t_{3,1}+t_{2,1}×t_{1,1}+t_{3,1}×t_{1,1}}
#' \deqn{x_{2,2}=1/4×(t_{1,1}×t_{2,2}+t_{1,1}×t_{2,3}+t_{1,2}×t_{2,1}+t_{1,3}×t_{2,1}+t_{1,1}×t_{3,2}+t_{1,1}×t_{3,3}+t_{1,2}×t_{3,1}+t_{1,3}×t_{3,1}+t_{2,1}×t_{1,2}+t_{2,1}×t_{1,3}+t_{2,2}×t_{1,1}+t_{2,3}×t_{1,1}+t_{3,1}×t_{1,2}+t_{3,1}×t_{1,3}+t_{3,2}×t_{1,1}+t_{3,3}×t_{1,1})}
#' \deqn{x_{2,3}=1/6×(t_{1,1}×t_{2,4}+t_{1,2}×t_{2,2}+t_{1,2}×t_{2,3}+t_{1,3}×t_{2,2}+t_{1,3}×t_{2,3}+t_{1,4}×t_{2,1}+t_{1,1}×t_{3,4}+t_{1,2}×t_{3,2}+t_{1,2}×t_{3,3}+t_{1,3}×t_{3,2}+t_{1,3}×t_{3,3}+t_{1,4}×t_{3,1}+t_{2,1}×t_{1,4}+t_{2,2}×t_{1,2}+t_{2,2}×t_{1,3}+t_{2,3}×t_{1,2}+t_{2,3}×t_{1,3}+t_{2,4}×t_{1,1}+t_{3,1}×t_{1,4}+t_{3,2}×t_{1,2}+t_{3,2}×t_{1,3}+t_{3,3}×t_{1,2}+t_{3,3}×t_{1,3}+t_{3,4}×t_{1,1})}
#' \deqn{x_{2,4}=1/4×(t_{1,2}×t_{2,4}+t_{1,3}×t_{2,4}+t_{1,4}×t_{2,2}+t_{1,4}×t_{2,3}+t_{1,2}×t_{3,4}+t_{1,3}×t_{3,4}+t_{1,4}×t_{3,2}+t_{1,4}×t_{3,3}+t_{2,2}×t_{1,4}+t_{2,3}×t_{1,4}+t_{2,4}×t_{1,2}+t_{2,4}×t_{1,3}+t_{3,2}×t_{1,4}+t_{3,3}×t_{1,4}+t_{3,4}×t_{1,2}+t_{3,4}×t_{1,3})}
#' \deqn{x_{2,5}=t_{1,4}×t_{2,4}+t_{1,4}×t_{3,4}+t_{2,4}×t_{1,4}+t_{3,4}×t_{1,4}}
#' \deqn{x_{3,1}=t_{1,1}×t_{4,1}+t_{2,1}×t_{2,1}+t_{2,1}×t_{3,1}+t_{3,1}×t_{2,1}+t_{3,1}×t_{3,1}+t_{4,1}×t_{1,1}}
#' \deqn{x_{3,2}=1/4×(t_{1,1}×t_{4,2}+t_{1,1}×t_{4,3}+t_{1,2}×t_{4,1}+t_{1,3}×t_{4,1}+t_{2,1}×t_{2,2}+t_{2,1}×t_{2,3}+t_{2,2}×t_{2,1}+t_{2,3}×t_{2,1}+t_{2,1}×t_{3,2}+t_{2,1}×t_{3,3}+t_{2,2}×t_{3,1}+t_{2,3}×t_{3,1}+t_{3,1}×t_{2,2}+t_{3,1}×t_{2,3}+t_{3,2}×t_{2,1}+t_{3,3}×t_{2,1}+t_{3,1}×t_{3,2}+t_{3,1}×t_{3,3}+t_{3,2}×t_{3,1}+t_{3,3}×t_{3,1}+t_{4,1}×t_{1,2}+t_{4,1}×t_{1,3}+t_{4,2}×t_{1,1}+t_{4,3}×t_{1,1})}
#' \deqn{x_{3,3}=1/6×(t_{1,1}×t_{4,4}+t_{1,2}×t_{4,2}+t_{1,2}×t_{4,3}+t_{1,3}×t_{4,2}+t_{1,3}×t_{4,3}+t_{1,4}×t_{4,1}+t_{2,1}×t_{2,4}+t_{2,2}×t_{2,2}+t_{2,2}×t_{2,3}+t_{2,3}×t_{2,2}+t_{2,3}×t_{2,3}+t_{2,4}×t_{2,1}+t_{2,1}×t_{3,4}+t_{2,2}×t_{3,2}+t_{2,2}×t_{3,3}+t_{2,3}×t_{3,2}+t_{2,3}×t_{3,3}+t_{2,4}×t_{3,1}+t_{3,1}×t_{2,4}+t_{3,2}×t_{2,2}+t_{3,2}×t_{2,3}+t_{3,3}×t_{2,2}+t_{3,3}×t_{2,3}+t_{3,4}×t_{2,1}+t_{3,1}×t_{3,4}+t_{3,2}×t_{3,2}+t_{3,2}×t_{3,3}+t_{3,3}×t_{3,2}+t_{3,3}×t_{3,3}+t_{3,4}×t_{3,1}+t_{4,1}×t_{1,4}+t_{4,2}×t_{1,2}+t_{4,2}×t_{1,3}+t_{4,3}×t_{1,2}+t_{4,3}×t_{1,3}+t_{4,4}×t_{1,1})}
#' \deqn{x_{3,4}=1/4×(t_{1,2}×t_{4,4}+t_{1,3}×t_{4,4}+t_{1,4}×t_{4,2}+t_{1,4}×t_{4,3}+t_{2,2}×t_{2,4}+t_{2,3}×t_{2,4}+t_{2,4}×t_{2,2}+t_{2,4}×t_{2,3}+t_{2,2}×t_{3,4}+t_{2,3}×t_{3,4}+t_{2,4}×t_{3,2}+t_{2,4}×t_{3,3}+t_{3,2}×t_{2,4}+t_{3,3}×t_{2,4}+t_{3,4}×t_{2,2}+t_{3,4}×t_{2,3}+t_{3,2}×t_{3,4}+t_{3,3}×t_{3,4}+t_{3,4}×t_{3,2}+t_{3,4}×t_{3,3}+t_{4,2}×t_{1,4}+t_{4,3}×t_{1,4}+t_{4,4}×t_{1,2}+t_{4,4}×t_{1,3})}
#' \deqn{x_{3,5}=t_{1,4}×t_{4,4}+t_{2,4}×t_{2,4}+t_{2,4}×t_{3,4}+t_{3,4}×t_{2,4}+t_{3,4}×t_{3,4}+t_{4,4}×t_{1,4}}
#' \deqn{x_{4,1}=t_{2,1}×t_{4,1}+t_{3,1}×t_{4,1}+t_{4,1}×t_{2,1}+t_{4,1}×t_{3,1}}
#' \deqn{x_{4,2}=1/4×(t_{2,1}×t_{4,2}+t_{2,1}×t_{4,3}+t_{2,2}×t_{4,1}+t_{2,3}×t_{4,1}+t_{3,1}×t_{4,2}+t_{3,1}×t_{4,3}+t_{3,2}×t_{4,1}+t_{3,3}×t_{4,1}+t_{4,1}×t_{2,2}+t_{4,1}×t_{2,3}+t_{4,2}×t_{2,1}+t_{4,3}×t_{2,1}+t_{4,1}×t_{3,2}+t_{4,1}×t_{3,3}+t_{4,2}×t_{3,1}+t_{4,3}×t_{3,1})}
#' \deqn{x_{4,3}=1/6×(t_{2,1}×t_{4,4}+t_{2,2}×t_{4,2}+t_{2,2}×t_{4,3}+t_{2,3}×t_{4,2}+t_{2,3}×t_{4,3}+t_{2,4}×t_{4,1}+t_{3,1}×t_{4,4}+t_{3,2}×t_{4,2}+t_{3,2}×t_{4,3}+t_{3,3}×t_{4,2}+t_{3,3}×t_{4,3}+t_{3,4}×t_{4,1}+t_{4,1}×t_{2,4}+t_{4,2}×t_{2,2}+t_{4,2}×t_{2,3}+t_{4,3}×t_{2,2}+t_{4,3}×t_{2,3}+t_{4,4}×t_{2,1}+t_{4,1}×t_{3,4}+t_{4,2}×t_{3,2}+t_{4,2}×t_{3,3}+t_{4,3}×t_{3,2}+t_{4,3}×t_{3,3}+t_{4,4}×t_{3,1})}
#' \deqn{x_{4,4}=1/4×(t_{2,2}×t_{4,4}+t_{2,3}×t_{4,4}+t_{2,4}×t_{4,2}+t_{2,4}×t_{4,3}+t_{3,2}×t_{4,4}+t_{3,3}×t_{4,4}+t_{3,4}×t_{4,2}+t_{3,4}×t_{4,3}+t_{4,2}×t_{2,4}+t_{4,3}×t_{2,4}+t_{4,4}×t_{2,2}+t_{4,4}×t_{2,3}+t_{4,2}×t_{3,4}+t_{4,3}×t_{3,4}+t_{4,4}×t_{3,2}+t_{4,4}×t_{3,3})}
#' \deqn{x_{4,5}=t_{2,4}×t_{4,4}+t_{3,4}×t_{4,4}+t_{4,4}×t_{2,4}+t_{4,4}×t_{3,4}}
#' \deqn{x_{5,1}=t_{4,1}×t_{4,1}}
#' \deqn{x_{5,2}=1/4×(t_{4,1}×t_{4,2}+t_{4,1}×t_{4,3}+t_{4,2}×t_{4,1}+t_{4,3}×t_{4,1})}
#' \deqn{x_{5,3}=1/6×(t_{4,1}×t_{4,4}+t_{4,2}×t_{4,2}+t_{4,2}×t_{4,3}+t_{4,3}×t_{4,2}+t_{4,3}×t_{4,3}+t_{4,4}×t_{4,1})}
#' \deqn{x_{5,4}=1/4×(t_{4,2}×t_{4,4}+t_{4,3}×t_{4,4}+t_{4,4}×t_{4,2}+t_{4,4}×t_{4,3})}
#' \deqn{x_{5,5}=t_{4,4}×t_{4,4}}
#' @references .
#' @examples MethylCalculation(original_classes,u,d,p)
#' MethylCalculation(c(0.1,0.2,0.3,0.1,0.1),u=0.01,d=0.2,p=0.8,cell_cycle=1)
#' MethylCalculation(c(0.1,0.2,0.3,0.1,0.1),u=0.01,d=0.2,p=0.8,cell_cycle=10)
#' @export MethylCalculation

MethylCalculation <- function(original_classes,u,d,p,cell_cycle=1){
	stopifnot(is.numeric(original_classes), is.vector(original_classes))
	stopifnot(is.numeric(u), is.numeric(d), is.numeric(p), is.integer(cell_cycle))
	if (u > 1||u < -1||d > 1||d < -1||p > 1||p < -1){
		stop("The probablities shoud be a number in c(0,1)!\n\n")
	}
	for (i in seq(cell_cycle)){
		original_classes <- PredictionMethylationClass(original_classes,u,d,p)
	}
	end_mean <- MethyClass2Mean(original_classes)
	cat("\n")
	cat("\tMethylation Calculation\n")
	cat("\n")
	cat(paste("It calculated the ratio of each methylation states after ",cell_cycle ," cell cycle(s).\n"))
	cat("\n")
	return(list("terminational_classes"=original_classes,"average_methylation_level"=end_mean))
}
