################################################################################################
###################################    Built-in functions    ###################################
################################################################################################
MethylationClassClass2ClassSummary <- function(x,y){
	cutoff_1 <- 1/8
	cutoff_2 <- 3/8
	cutoff_3 <- 5/8
	cutoff_4 <- 7/8
	s1e1 <- sum(x<=cutoff_1 & y<=cutoff_1,na.rm=T)
	s1e2 <- sum(x<=cutoff_1 & y>cutoff_1&y<=cutoff_2,na.rm=T) # s1e2 means startClass1 -> endClass2 (not the row1col2). actually is col1row2.
	s1e3 <- sum(x<=cutoff_1 & y>cutoff_2&y<=cutoff_3,na.rm=T)
	s1e4 <- sum(x<=cutoff_1 & y>cutoff_3&y<=cutoff_4,na.rm=T)
	s1e5 <- sum(x<=cutoff_1 & y>cutoff_4,na.rm=T)
	s2e1 <- sum(x>cutoff_1&x<=cutoff_2 & y<=cutoff_1,na.rm=T)
	s2e2 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_1&y<=cutoff_2,na.rm=T)
	s2e3 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_2&y<=cutoff_3,na.rm=T)
	s2e4 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_3&y<=cutoff_4,na.rm=T)
	s2e5 <- sum(x>cutoff_1&x<=cutoff_2 & y>cutoff_4,na.rm=T)
	s3e1 <- sum(x>cutoff_2&x<=cutoff_3 & y<=cutoff_1,na.rm=T)
	s3e2 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_1&y<=cutoff_2,na.rm=T)
	s3e3 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_2&y<=cutoff_3,na.rm=T)
	s3e4 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_3&y<=cutoff_4,na.rm=T)
	s3e5 <- sum(x>cutoff_2&x<=cutoff_3 & y>cutoff_4,na.rm=T)
	s4e1 <- sum(x>cutoff_3&x<=cutoff_4 & y<=cutoff_1,na.rm=T)
	s4e2 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_1&y<=cutoff_2,na.rm=T)
	s4e3 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_2&y<=cutoff_3,na.rm=T)
	s4e4 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_3&y<=cutoff_4,na.rm=T)
	s4e5 <- sum(x>cutoff_3&x<=cutoff_4 & y>cutoff_4,na.rm=T)
	s5e1 <- sum(x>cutoff_4 & y<=cutoff_1,na.rm=T)
	s5e2 <- sum(x>cutoff_4 & y>cutoff_1&y<=cutoff_2,na.rm=T)
	s5e3 <- sum(x>cutoff_4 & y>cutoff_2&y<=cutoff_3,na.rm=T)
	s5e4 <- sum(x>cutoff_4 & y>cutoff_3&y<=cutoff_4,na.rm=T)
	s5e5 <- sum(x>cutoff_4 & y>cutoff_4,na.rm=T)
	s1_total <- s1e1+s1e2+s1e3+s1e4+s1e5
	s2_total <- s2e1+s2e2+s2e3+s2e4+s2e5
	s3_total <- s3e1+s3e2+s3e3+s3e4+s3e5
	s4_total <- s4e1+s4e2+s4e3+s4e4+s4e5
	s5_total <- s5e1+s5e2+s5e3+s5e4+s5e5
	s1e1_ratio <- s1e1/s1_total;s1e2_ratio <- s1e2/s1_total;s1e3_ratio <- s1e3/s1_total;s1e4_ratio <- s1e4/s1_total;s1e5_ratio <- s1e5/s1_total
	s2e1_ratio <- s2e1/s2_total;s2e2_ratio <- s2e2/s2_total;s2e3_ratio <- s2e3/s2_total;s2e4_ratio <- s2e4/s2_total;s2e5_ratio <- s2e5/s2_total
	s3e1_ratio <- s3e1/s3_total;s3e2_ratio <- s3e2/s3_total;s3e3_ratio <- s3e3/s3_total;s3e4_ratio <- s3e4/s3_total;s3e5_ratio <- s3e5/s3_total
	s4e1_ratio <- s4e1/s4_total;s4e2_ratio <- s4e2/s4_total;s4e3_ratio <- s4e3/s4_total;s4e4_ratio <- s4e4/s4_total;s4e5_ratio <- s4e5/s4_total
	s5e1_ratio <- s5e1/s5_total;s5e2_ratio <- s5e2/s5_total;s5e3_ratio <- s5e3/s5_total;s5e4_ratio <- s5e4/s5_total;s5e5_ratio <- s5e5/s5_total
	if (s1_total == 0){
		s1e1_ratio <- 0;s1e2_ratio <- 0;s1e3_ratio <- 0;s1e4_ratio <- 0;s1e5_ratio <- 0
	}
	if (s2_total == 0){
		s2e1_ratio <- 0;s2e2_ratio <- 0;s2e3_ratio <- 0;s2e4_ratio <- 0;s2e5_ratio <- 0
	}
	if (s3_total == 0){
		s3e1_ratio <- 0;s3e2_ratio <- 0;s3e3_ratio <- 0;s3e4_ratio <- 0;s3e5_ratio <- 0
	}
	if (s4_total == 0){
		s4e1_ratio <- 0;s4e2_ratio <- 0;s4e3_ratio <- 0;s4e4_ratio <- 0;s4e5_ratio <- 0
	}
	if (s5_total == 0){
		s5e1_ratio <- 0;s5e2_ratio <- 0;s5e3_ratio <- 0;s5e4_ratio <- 0;s5e5_ratio <- 0
	}
	tmp_matrix <- round(matrix(c(s1e1/s1_total,s1e2/s1_total,s1e3/s1_total,s1e4/s1_total,s1e5/s1_total,s2e1/s2_total,s2e2/s2_total,s2e3/s2_total,s2e4/s2_total,s2e5/s2_total,s3e1/s3_total,s3e2/s3_total,s3e3/s3_total,s3e4/s3_total,s3e5/s3_total,s4e1/s4_total,s4e2/s4_total,s4e3/s4_total,s4e4/s4_total,s4e5/s4_total,s5e1/s5_total,s5e2/s5_total,s5e3/s5_total,s5e4/s5_total,s5e5/s5_total),c(5,5),byrow = F),4)
	colnames(tmp_matrix) <- c("original_class1","original_class2","original_class3","original_class4","original_class5")
	rownames(tmp_matrix) <- c("terminal_class1","terminal_class2","terminal_class3","terminal_class4","terminal_class5")
	return(tmp_matrix)
}

################################################################################################
###################################     Library functions    ###################################
################################################################################################

#' @title "TransitionMatrixGeneration"
#' @description Generate the transition matrix from the original DNA methylation level to the terminal DNA methylation level which is suitable for TMParameterEstimation().
#' @param original_methyl The original methylation level of each CpG/gene/region.
#' @param terminal_methyl The paired terminal methylation level of each CpG/gene/region.
#' @return observation_matrix The transition ratio matrix(\eqn{5 \times 5}) from the original state to the terminal state.
#' DNA methylation level was classed to five classes (class1:0;class1:1/4;class1:1/2;class1:3/4;class5:1).
#' The ratio of each class change to others are described as blows:
#' 	\tabular{cccccc}{
#' 		\tab original_class1 \tab original_class2 \tab original_class3 \tab original_class4 \tab original_class5\cr
#' 		terminal_class1 \tab a1 \tab b1 \tab c1 \tab d1 \tab e1\cr
#' 		terminal_class2 \tab a2 \tab b2 \tab c2 \tab d2 \tab e2\cr
#' 		terminal_class3 \tab a3 \tab b3 \tab c3 \tab d3 \tab e3\cr
#' 		terminal_class4 \tab a4 \tab b4 \tab c4 \tab d4 \tab e4\cr
#' 		terminal_class5 \tab a5 \tab b5 \tab c5 \tab d5 \tab e5\cr
#' 	}
#' 	# observation_matrix[1,1](a1) is the ratio of original_class1 to terminal_class1,observation_matrix[1,2](b1) is the ratio of original_class2 to terminal_class1 and so on
#' 	# observation_matrix[2,1](a2) is the ratio of original_class1 to terminal_class2,observation_matrix[2,2](b2) is the ratio of original_class2 to terminal_class2 and so on
#' 	# The sum of the ratio that original_class1 change to 5 classes should be 1. That is a1+a2+a3+a4+a5=1.
#' @details .
#' @references \cite{Zhao, C. et.al.(2018). A DNA methylation state transition model reveals the programmed epigenetic heterogeneity in pre-implantation embryos. Under revision.}
#' @examples # Let's simulate the original DNA methylation level vector and the terminal one.
#' set.seed(0)
#' original_methyl <- runif(10000, min = 0, max = 1)
#' set.seed(1)
#' terminal_methyl <- runif(10000, min = 0, max = 1)
#' TransitionMatrixGeneration(original_methyl,terminal_methyl)
#' @concept MethylTransition
#' @export TransitionMatrixGeneration

TransitionMatrixGeneration <- function(original_methyl,terminal_methyl){
	stopifnot(is.vector(original_methyl))
	stopifnot(is.vector(terminal_methyl))
	if (length(original_methyl) != length(original_methyl)){
		stop("The number of methylation values between original and terminal are not match.")
	}
	cat("\n")
	cat("\tThe transtion matrix is :\n")
	return(MethylationClassClass2ClassSummary(original_methyl,terminal_methyl))
}
