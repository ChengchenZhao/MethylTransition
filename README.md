# MethylTransition overview

This is a R package for the parameter estimation of DNA methylation transition in cell mitosis. The transition matrix of this model describes the changes of DNA methylation during one cell cycle in three steps: passive demethylation by DNA replication, active DNA methylation changes affected by DNA methylation-modifying enzymes and DNA methylation combinations during homologous recombination. This model includes three parameters, *u*, *d* and *p*, to represent the probabilities of three active DNA methylation change types, and those parameters are estimated by minimizing a cost function using the Newton-Raphson method. For now, if you use this code, please cite [the manuscript](https://zhanglab.tongji.edu.cn).

# Install
Install from GitHub : 
```
# install.packages("devtools") # if devtools not be installed
devtools::install_github("ChengchenZhao/MethylTransition")
```

# Step by step guide
