% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethylCalculation.R
\name{MethylCalculation}
\alias{MethylCalculation}
\title{"MethylCalculation"}
\usage{
MethylCalculation(original_classes, u, d, p, cell_cycle = 1)
}
\arguments{
\item{original_classes}{The original methylation classes}

\item{u}{The paramater that describing the methylation probablity on CpG site}

\item{d}{The paramater that describing the de-methylation probablity on 5mCpG site}

\item{p}{The paramater that describing the methylation probablity on semi-CpG site}

\item{cell_cycle}{The cell cycle times}
}
\value{
terminal_classes The terminal methylation classes

average_methylation_level The average methylation level after [n] cell cycle(s).
}
\description{
Calculate the ratio of each methylation class after n cell cycle(s).
}
\details{
The transition matrix of this model describes the changes of DNA methylation during one cell cycle in three steps:
passive demethylation by DNA replication, active DNA methylation changes affected by DNA methylation-modifying enzymes
and DNA methylation combinations during homologous recombination.
For each CpG site in a chromsome, the methylation states are one of these four types: 0-0, 0-1, 1-0, 1-1.
The transition matrix after DNA replication would be :
	\tabular{cccccc}{
		\tab original(0-0) \tab original(0-1) \tab original(1-0) \tab original(1-1) \cr
		after_replication(0-0) \tab \eqn{1} \tab \eqn{a} \tab \eqn{1-a} \tab \eqn{0} \cr
		after_replication(0-1) \tab \eqn{0} \tab \eqn{1-a} \tab \eqn{0} \tab \eqn{1-a} \cr
		after_replication(1-0) \tab \eqn{0} \tab \eqn{0} \tab \eqn{0} \tab \eqn{a} \cr
		after_replication(1-1) \tab \eqn{0} \tab \eqn{0} \tab \eqn{0} \tab \eqn{0} \cr
	}
	among this matrix ,\eqn{a} is the methylation change probability with DNA replication and equal to 0.5.
	Then the transition matrix after active DNA methylation changes would be:
	\tabular{cccccc}{
		\tab after_replication(0-0) \tab after_replication(0-1) \tab after_replication(1-0) \tab after_replication(1-1) \cr
		after_enzymemodifying(0-0) \tab \eqn{(1-u )\times (1-u)} \tab \eqn{(1-u-p+u \times p )\times d} \tab \eqn{d \times (1-u-p+u \times p)} \tab \eqn{0} \cr
		after_enzymemodifying(0-1) \tab \eqn{u \times (1-u)} \tab \eqn{(1-u-p+u \times p )\times (1-d)} \tab \eqn{d \times (u+p-u \times p)} \tab \eqn{0} \cr
		after_enzymemodifying(1-0) \tab \eqn{u \times (1-u)} \tab \eqn{(u+p-u \times p )\times d} \tab \eqn{(1-d )\times (1-u-p+u \times p)} \tab \eqn{0} \cr
		after_enzymemodifying(1-1) \tab \eqn{u \times u} \tab \eqn{(u+p-u \times p )\times (1-d)} \tab \eqn{(1-d )\times (u+p-u \times p)} \tab \eqn{1} \cr
	}
	The paramater \eqn{u} described the methylation probablity on CpG site.
	The paramater \eqn{d} described the de-methylation probablity on 5mCpG site.
	The paramater \eqn{p} described the methylation probablity on semi-CpG site.
	Thus the mathlation state change is \eqn{P = Penzymemodifying \cdot Preplication \cdot original_classes}
	We using \eqn{t_{i,j}} represent the \eqn{i} and \eqn{j} vector of the matrix \eqn{Penzymemodifying \cdot Preplication}
	Then two chromsomes are combined during homologous recombination. The observed DNA methylation of each CpG site is the combination of the methylation types in both chromsomes.
 The observed transition matrix would be :
	\tabular{cccccc}{
		\tab original_class1(0) \tab original_class2(1/4) \tab original_class3(1/2) \tab original_class4(3/4) \tab original_class5(1)\cr
		terminal_class1(0) \tab \eqn{x_{1,1}} \tab \eqn{x_{1,2}} \tab \eqn{x_{1,3}} \tab \eqn{x_{1,4}} \tab \eqn{x_{1,5}}\cr
		terminal_class2(1/4) \tab \eqn{x_{2,1}} \tab \eqn{x_{2,2}} \tab \eqn{x_{2,3}} \tab \eqn{x_{2,4}} \tab \eqn{x_{2,5}}\cr
		terminal_class3(1/2) \tab \eqn{x_{3,1}} \tab \eqn{x_{3,2}} \tab \eqn{x_{3,3}} \tab \eqn{x_{3,4}} \tab \eqn{x_{3,5}}\cr
		terminal_class4(3/4) \tab \eqn{x_{4,1}} \tab \eqn{x_{4,2}} \tab \eqn{x_{4,3}} \tab \eqn{x_{4,4}} \tab \eqn{x_{4,5}}\cr
		terminal_class5(1) \tab \eqn{x_{5,1}} \tab \eqn{x_{5,2}} \tab \eqn{x_{5,3}} \tab \eqn{x_{5,4}} \tab \eqn{x_{5,5}}\cr
	}
	and \deqn{x_{1,1}=t_{1,1 }\times t_{1,1}}
\deqn{x_{1,1}=t_{{,1},1 }\times t_{{,1},1}}
\deqn{x_{1,2}=1/4 \times (t_{1,1 }\times t_{1,2}+t_{1,1 }\times t_{1,3}+t_{1,2 }\times t_{1,1}+t_{1,3 }\times t_{1,1})}
\deqn{x_{1,3}=1/6 \times (t_{1,1 }\times t_{1,4}+t_{1,2 }\times t_{1,2}+t_{1,2 }\times t_{1,3}+t_{1,3 }\times t_{1,2}+t_{1,3 }\times t_{1,3}+t_{1,4 }\times t_{1,1})}
\deqn{x_{1,4}=1/4 \times (t_{1,2 }\times t_{1,4}+t_{1,3 }\times t_{1,4}+t_{1,4 }\times t_{1,2}+t_{1,4 }\times t_{1,3})}
\deqn{x_{1,5}=t_{1,4 }\times t_{1,4}}
\deqn{x_{2,1}=t_{1,1 }\times t_{2,1}+t_{1,1 }\times t_{3,1}+t_{2,1 }\times t_{1,1}+t_{3,1 }\times t_{1,1}}
\deqn{x_{2,2}=1/4 \times (t_{1,1 }\times t_{2,2}+t_{1,1 }\times t_{2,3}+t_{1,2 }\times t_{2,1}+t_{1,3 }\times t_{2,1}+t_{1,1 }\times t_{3,2}+t_{1,1 }\times t_{3,3}+t_{1,2 }\times t_{3,1}+t_{1,3 }\times t_{3,1}+t_{2,1 }\times t_{1,2}+t_{2,1 }\times t_{1,3}+t_{2,2 }\times t_{1,1}+t_{2,3 }\times t_{1,1}+t_{3,1 }\times t_{1,2}+t_{3,1 }\times t_{1,3}+t_{3,2 }\times t_{1,1}+t_{3,3 }\times t_{1,1})}
\deqn{x_{2,3}=1/6 \times (t_{1,1 }\times t_{2,4}+t_{1,2 }\times t_{2,2}+t_{1,2 }\times t_{2,3}+t_{1,3 }\times t_{2,2}+t_{1,3 }\times t_{2,3}+t_{1,4 }\times t_{2,1}+t_{1,1 }\times t_{3,4}+t_{1,2 }\times t_{3,2}+t_{1,2 }\times t_{3,3}+t_{1,3 }\times t_{3,2}+t_{1,3 }\times t_{3,3}+t_{1,4 }\times t_{3,1}+t_{2,1 }\times t_{1,4}+t_{2,2 }\times t_{1,2}+t_{2,2 }\times t_{1,3}+t_{2,3 }\times t_{1,2}+t_{2,3 }\times t_{1,3}+t_{2,4 }\times t_{1,1}+t_{3,1 }\times t_{1,4}+t_{3,2 }\times t_{1,2}+t_{3,2 }\times t_{1,3}+t_{3,3 }\times t_{1,2}+t_{3,3 }\times t_{1,3}+t_{3,4 }\times t_{1,1})}
\deqn{x_{2,4}=1/4 \times (t_{1,2 }\times t_{2,4}+t_{1,3 }\times t_{2,4}+t_{1,4 }\times t_{2,2}+t_{1,4 }\times t_{2,3}+t_{1,2 }\times t_{3,4}+t_{1,3 }\times t_{3,4}+t_{1,4 }\times t_{3,2}+t_{1,4 }\times t_{3,3}+t_{2,2 }\times t_{1,4}+t_{2,3 }\times t_{1,4}+t_{2,4 }\times t_{1,2}+t_{2,4 }\times t_{1,3}+t_{3,2 }\times t_{1,4}+t_{3,3 }\times t_{1,4}+t_{3,4 }\times t_{1,2}+t_{3,4 }\times t_{1,3})}
\deqn{x_{2,5}=t_{1,4 }\times t_{2,4}+t_{1,4 }\times t_{3,4}+t_{2,4 }\times t_{1,4}+t_{3,4 }\times t_{1,4}}
\deqn{x_{3,1}=t_{1,1 }\times t_{4,1}+t_{2,1 }\times t_{2,1}+t_{2,1 }\times t_{3,1}+t_{3,1 }\times t_{2,1}+t_{3,1 }\times t_{3,1}+t_{4,1 }\times t_{1,1}}
\deqn{x_{3,2}=1/4 \times (t_{1,1 }\times t_{4,2}+t_{1,1 }\times t_{4,3}+t_{1,2 }\times t_{4,1}+t_{1,3 }\times t_{4,1}+t_{2,1 }\times t_{2,2}+t_{2,1 }\times t_{2,3}+t_{2,2 }\times t_{2,1}+t_{2,3 }\times t_{2,1}+t_{2,1 }\times t_{3,2}+t_{2,1 }\times t_{3,3}+t_{2,2 }\times t_{3,1}+t_{2,3 }\times t_{3,1}+t_{3,1 }\times t_{2,2}+t_{3,1 }\times t_{2,3}+t_{3,2 }\times t_{2,1}+t_{3,3 }\times t_{2,1}+t_{3,1 }\times t_{3,2}+t_{3,1 }\times t_{3,3}+t_{3,2 }\times t_{3,1}+t_{3,3 }\times t_{3,1}+t_{4,1 }\times t_{1,2}+t_{4,1 }\times t_{1,3}+t_{4,2 }\times t_{1,1}+t_{4,3 }\times t_{1,1})}
\deqn{x_{3,3}=1/6 \times (t_{1,1 }\times t_{4,4}+t_{1,2 }\times t_{4,2}+t_{1,2 }\times t_{4,3}+t_{1,3 }\times t_{4,2}+t_{1,3 }\times t_{4,3}+t_{1,4 }\times t_{4,1}+t_{2,1 }\times t_{2,4}+t_{2,2 }\times t_{2,2}+t_{2,2 }\times t_{2,3}+t_{2,3 }\times t_{2,2}+t_{2,3 }\times t_{2,3}+t_{2,4 }\times t_{2,1}+t_{2,1 }\times t_{3,4}+t_{2,2 }\times t_{3,2}+t_{2,2 }\times t_{3,3}+t_{2,3 }\times t_{3,2}+t_{2,3 }\times t_{3,3}+t_{2,4 }\times t_{3,1}+t_{3,1 }\times t_{2,4}+t_{3,2 }\times t_{2,2}+t_{3,2 }\times t_{2,3}+t_{3,3 }\times t_{2,2}+t_{3,3 }\times t_{2,3}+t_{3,4 }\times t_{2,1}+t_{3,1 }\times t_{3,4}+t_{3,2 }\times t_{3,2}+t_{3,2 }\times t_{3,3}+t_{3,3 }\times t_{3,2}+t_{3,3 }\times t_{3,3}+t_{3,4 }\times t_{3,1}+t_{4,1 }\times t_{1,4}+t_{4,2 }\times t_{1,2}+t_{4,2 }\times t_{1,3}+t_{4,3 }\times t_{1,2}+t_{4,3 }\times t_{1,3}+t_{4,4 }\times t_{1,1})}
\deqn{x_{3,4}=1/4 \times (t_{1,2 }\times t_{4,4}+t_{1,3 }\times t_{4,4}+t_{1,4 }\times t_{4,2}+t_{1,4 }\times t_{4,3}+t_{2,2 }\times t_{2,4}+t_{2,3 }\times t_{2,4}+t_{2,4 }\times t_{2,2}+t_{2,4 }\times t_{2,3}+t_{2,2 }\times t_{3,4}+t_{2,3 }\times t_{3,4}+t_{2,4 }\times t_{3,2}+t_{2,4 }\times t_{3,3}+t_{3,2 }\times t_{2,4}+t_{3,3 }\times t_{2,4}+t_{3,4 }\times t_{2,2}+t_{3,4 }\times t_{2,3}+t_{3,2 }\times t_{3,4}+t_{3,3 }\times t_{3,4}+t_{3,4 }\times t_{3,2}+t_{3,4 }\times t_{3,3}+t_{4,2 }\times t_{1,4}+t_{4,3 }\times t_{1,4}+t_{4,4 }\times t_{1,2}+t_{4,4 }\times t_{1,3})}
\deqn{x_{3,5}=t_{1,4 }\times t_{4,4}+t_{2,4 }\times t_{2,4}+t_{2,4 }\times t_{3,4}+t_{3,4 }\times t_{2,4}+t_{3,4 }\times t_{3,4}+t_{4,4 }\times t_{1,4}}
\deqn{x_{4,1}=t_{2,1 }\times t_{4,1}+t_{3,1 }\times t_{4,1}+t_{4,1 }\times t_{2,1}+t_{4,1 }\times t_{3,1}}
\deqn{x_{4,2}=1/4 \times (t_{2,1 }\times t_{4,2}+t_{2,1 }\times t_{4,3}+t_{2,2 }\times t_{4,1}+t_{2,3 }\times t_{4,1}+t_{3,1 }\times t_{4,2}+t_{3,1 }\times t_{4,3}+t_{3,2 }\times t_{4,1}+t_{3,3 }\times t_{4,1}+t_{4,1 }\times t_{2,2}+t_{4,1 }\times t_{2,3}+t_{4,2 }\times t_{2,1}+t_{4,3 }\times t_{2,1}+t_{4,1 }\times t_{3,2}+t_{4,1 }\times t_{3,3}+t_{4,2 }\times t_{3,1}+t_{4,3 }\times t_{3,1})}
\deqn{x_{4,3}=1/6 \times (t_{2,1 }\times t_{4,4}+t_{2,2 }\times t_{4,2}+t_{2,2 }\times t_{4,3}+t_{2,3 }\times t_{4,2}+t_{2,3 }\times t_{4,3}+t_{2,4 }\times t_{4,1}+t_{3,1 }\times t_{4,4}+t_{3,2 }\times t_{4,2}+t_{3,2 }\times t_{4,3}+t_{3,3 }\times t_{4,2}+t_{3,3 }\times t_{4,3}+t_{3,4 }\times t_{4,1}+t_{4,1 }\times t_{2,4}+t_{4,2 }\times t_{2,2}+t_{4,2 }\times t_{2,3}+t_{4,3 }\times t_{2,2}+t_{4,3 }\times t_{2,3}+t_{4,4 }\times t_{2,1}+t_{4,1 }\times t_{3,4}+t_{4,2 }\times t_{3,2}+t_{4,2 }\times t_{3,3}+t_{4,3 }\times t_{3,2}+t_{4,3 }\times t_{3,3}+t_{4,4 }\times t_{3,1})}
\deqn{x_{4,4}=1/4 \times (t_{2,2 }\times t_{4,4}+t_{2,3 }\times t_{4,4}+t_{2,4 }\times t_{4,2}+t_{2,4 }\times t_{4,3}+t_{3,2 }\times t_{4,4}+t_{3,3 }\times t_{4,4}+t_{3,4 }\times t_{4,2}+t_{3,4 }\times t_{4,3}+t_{4,2 }\times t_{2,4}+t_{4,3 }\times t_{2,4}+t_{4,4 }\times t_{2,2}+t_{4,4 }\times t_{2,3}+t_{4,2 }\times t_{3,4}+t_{4,3 }\times t_{3,4}+t_{4,4 }\times t_{3,2}+t_{4,4 }\times t_{3,3})}
\deqn{x_{4,5}=t_{2,4 }\times t_{4,4}+t_{3,4 }\times t_{4,4}+t_{4,4 }\times t_{2,4}+t_{4,4 }\times t_{3,4}}
\deqn{x_{5,1}=t_{4,1 }\times t_{4,1}}
\deqn{x_{5,2}=1/4 \times (t_{4,1 }\times t_{4,2}+t_{4,1 }\times t_{4,3}+t_{4,2 }\times t_{4,1}+t_{4,3 }\times t_{4,1})}
\deqn{x_{5,3}=1/6 \times (t_{4,1 }\times t_{4,4}+t_{4,2 }\times t_{4,2}+t_{4,2 }\times t_{4,3}+t_{4,3 }\times t_{4,2}+t_{4,3 }\times t_{4,3}+t_{4,4 }\times t_{4,1})}
\deqn{x_{5,4}=1/4 \times (t_{4,2 }\times t_{4,4}+t_{4,3 }\times t_{4,4}+t_{4,4 }\times t_{4,2}+t_{4,4 }\times t_{4,3})}
\deqn{x_{5,5}=t_{4,4 }\times t_{4,4}}
}
\examples{
# MethylCalculation(original_classes,u,d,p)
MethylCalculation(c(0.1,0.2,0.3,0.1,0.1),u=0.01,d=0.2,p=0.8,cell_cycle=1)
MethylCalculation(c(0.1,0.2,0.3,0.1,0.1),u=0.01,d=0.2,p=0.8,cell_cycle=10)
}
\references{
\cite{Zhao, C. et.al.(2019). A DNA methylation state transition model reveals the programmed epigenetic heterogeneity in pre-implantation embryos. Under revision.}
}
\concept{MethylTransition}
