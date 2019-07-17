# MethylTransition overview

MethylTransition is a R package for the parameter estimation of DNA methylation transition in cell mitosis. The transition matrix of this model describes the changes of DNA methylation during one cell cycle in three steps: passive demethylation by DNA replication, active DNA methylation changes affected by DNA methylation-modifying enzymes and DNA methylation combinations during homologous recombination. This model includes three parameters, *u*, *d* and *p*, to represent the probabilities of three active DNA methylation change types, and those parameters are estimated by minimizing a cost function using the Newton-Raphson method. For now, if you use this code, please cite Zhao, C. et.al.(2019). *A DNA methylation state transition model reveals the programmed epigenetic heterogeneity in pre-implantation embryos.* Under revision. 
<!-- # [the manuscript](https://zhanglab.tongji.edu.cn).-->

# Install

Install from GitHub :
```R
# install.packages("devtools") # if devtools not be installed
devtools::install_github("ChengchenZhao/MethylTransition")
```

# Step by step guide

MethylTransition has two main functions,`ParameterEstimation()` and `MethylCalculation()`. `ParameterEstimation()` can be used to estimate the methylation change parameters using the initial and terminal DNA methylation states. `MethylCalculation()` can be used to calculate calculation the terminal proportion of each DNA methylation state using the initial DNA methylation states and a group of parameters. Here we show two examples of using these two functions.

### ```ParameterEstimation()```

You can use this function by running
```R
ParameterEstimation(original_methyl, terminal_methyl, iter = 50, cell_cycle = 1)
```
- **original_methyl** - The original methylation level of each CpG/gene/region.
- **terminal_methyl** - The paired terminal methylation level of each CpG/gene/region.
- **iter** - The iteration times of the parameter estimation using the Newton-Raphson method with different initial guesses.
- **cell_cycle** - The cell cycle times from the original state to the terminal state.

Next is an example using simulation DNA methylation data.
Let's first simulate the original DNA methylation level vector and the terminal one.
```R
set.seed(0)
original_methyl <- runif(10000, min = 0, max = 1)
set.seed(1)
terminal_methyl <- runif(10000, min = 0, max = 1)
```

If this process goes through one cell cycle, we can set that "cell_cycle=1"
```R
ParameterEstimation(original_methyl, terminal_methyl, iter = 30, cell_cycle = 1)
```

If this process goes through two cell cycle, we then set that "cell_cycle=2"
```R
ParameterEstimation(original_methyl, terminal_methyl, iter = 1, cell_cycle = 2)
```

If this function was not successful to estimated the function, you may try more iterations with different initial guesses(increasing the parameter "iter").
```R
ParameterEstimation(original_methyl, terminal_methyl, iter = 50, cell_cycle = 2)
```

The output results would include the `estimated_parameters` and the `predicted_matrix` which is the calculated transition matrix using the estimated parameters.

### ```MethylCalculation()```

You can use this function by running
```R
MethylCalculation(original_classes, u, d, p, cell_cycle = 1)
```
- **original_classes** - The original methylation classes
- **u** - The paramater that describing the methylation probablity on CpG site
- **d** - The paramater that describing the de-methylation probablity on 5mCpG site
- **p** - The paramater that describing the methylation probablity on semi-CpG site
- **cell_cycle** - The cell cycle times

Next is a simple example.
We first simulated a proportion of five cases of DNA methylation states at each CpG pair: Unmethylated (C1=0.1), Quarter-methylated (C2=0.2), Half-methylated (C3=0.3), Three-quarter-methylated (C4=0.1) and Fully-methylated (C5=0.1). Thus the `original_classes` would be `c(0.1,0.2,0.3,0.1,0.1)`. With a given group of parameters(*u*=0.01, *d*=0.2, *p*=0.8), the termial proportion of DNA methylation states can calculated by
```R
MethylCalculation(c(0.1,0.2,0.3,0.1,0.1), u = 0.01, d = 0.2, p = 0.8, cell_cycle = 1)
```

The output results would include the `terminal_classes` and the `average_methylation_level` of the terminal DNA methylation states.
