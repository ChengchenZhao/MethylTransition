# MethylTransition overview

MethylTransition is an R package (R version > 3.4) for characterizing the methylation change during one or a few cell cycles at single cell resolution based on a DNA methylation state transition model. MethylTransition relies on the assumption that the transition of DNA methylation state at a CpG site across single cell cycle occurs in three steps: passive demethylation during DNA replication, active DNA methylation transition affected by DNA methylation-modifying enzymes and DNA methylation states combination during the combination of non-sister chromatids. Therefore MethylTransition introduced a methylation state ratio vector, with 5 discrete states, to describe the overall pattern of DNA methylation states for a given cell. To link the two methylation state ratio vectors before and after a cell cycle, MethylTransition conducts a transition matrix comprised by 3 parameters, representing the possibility of DNA methylation maintenance (*u*), active demethylation (*d*) and *de novo* methylation separately (*p*), and it estimates the parameters by a matrix approximation strategy with the Newton-Raphson method. 
<!-- # [Citation: If you use this code, please cite Zhao, C. et.al.(2020). *A DNA methylation state transition model reveals the programmed epigenetic heterogeneity in pre-implantation embryos.* Under revision. ](https://zhanglab.tongji.edu.cn).-->

# Installation

```R
# install.packages("devtools") # run this command if "devtools" is not installed
devtools::install_github("ChengchenZhao/MethylTransition")
library("MethylTransition")
```

# Guidance

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

If this function was not successful to estimated the function, you may try more iterations with different initial guesses (by increasing the parameter "iter").
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
- **u** - The parameter that describing the methylation probability on CpG site
- **d** - The parameter that describing the de-methylation probability on 5mCpG site
- **p** - The parameter that describing the methylation probability on semi-CpG site
- **cell_cycle** - The cell cycle times

Next is a simple example.
We first simulated a proportion of five cases of DNA methylation states at each CpG pair: Unmethylated (Prop<sub>S1</sub> = 0.1), Quarter-methylated (Prop<sub>S2</sub> = 0.2), Half-methylated (Prop<sub>S3</sub> = 0.3), Three-quarter-methylated (Prop<sub>S4</sub> = 0.1) and Fully-methylated (Prop<sub>S5</sub> = 0.1). Thus the `original_classes` would be `c(0.1,0.2,0.3,0.1,0.1)`. With a given group of parameters (*u*=0.01, *d*=0.2, *p*=0.8), the terminal proportion of DNA methylation states can be calculated by
```R
MethylCalculation(c(0.1,0.2,0.3,0.1,0.1), u = 0.01, d = 0.2, p = 0.8, cell_cycle = 1)
```

The output results would include the `terminal_classes` and the `average_methylation_level` of the terminal DNA methylation states.

# Citation
Chengchen Zhao, Naiqian Zhang, Yalin Zhang, Nuermaimaiti Tuersunjiang, Shaorong Gao, Wenqiang Liu, and Yong Zhang. *A DNA methylation state transition model reveals the programmed epigenetic heterogeneity in human pre-implantation embryos*. Genome Biology 2020;21(1):277.