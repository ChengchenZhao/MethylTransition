# MethylTransition overview

MethylTransition is a R package for the parameter estimation of DNA methylation transition in cell mitosis. The transition matrix of this model describes the changes of DNA methylation during one cell cycle in three steps: passive demethylation by DNA replication, active DNA methylation changes affected by DNA methylation-modifying enzymes and DNA methylation combinations during homologous recombination. This model includes three parameters, *u*, *d* and *p*, to represent the probabilities of three active DNA methylation change types, and those parameters are estimated by minimizing a cost function using the Newton-Raphson method. For now, if you use this code, please cite [the manuscript](https://zhanglab.tongji.edu.cn).

# Install

Install from GitHub :
```
# install.packages("devtools") # if devtools not be installed
devtools::install_github("ChengchenZhao/MethylTransition")
```

# Step by step guide

MethylTransition has two main functions,```ParameterEstimation()``` and ```MethylCalculation()```. ```ParameterEstimation()``` can be used to estimate the methylation change parameters using the initial and terminal DNA methylation states. ```MethylCalculation()``` can be used to calculate calculation the terminal proportion of each DNA methylation state using the initial DNA methylation states and a group of parameters. Here we show two examples of using these two functions.

### ```ParameterEstimation()```

You can use this function by running
```
ParameterEstimation(original_methyl, terminal_methyl, iter = 50, cell_cycle = 1)
```
-**original_methyl** The original methylation level of each CpG/gene/region.
-**terminal_methyl** The paired terminal methylation level of each CpG/gene/region.
-**iter** The iteration times of the parameter estimation using the Newton-Raphson method with different initial guesses.
-**cell_cycle** The cell cycle times from the original state to the terminal state.