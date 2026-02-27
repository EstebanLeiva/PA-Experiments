# An Exact Method for Reliable Shortest Path Problems with Correlation

This repository contains the experiments for the paper **"An Exact Method for Reliable Shortest Path Problems with Correlation"**. It uses the implementation
of the pulse algorithm in [Julia](https://github.com/EstebanLeiva/PulseAlgorithm.jl.git).

## Authors
* Esteban Leiva (Centro para la Optimización y Probabilidad Aplicada (COPA), Departamento de Ingeniería Industrial, Universidad de los Andes, Bogotá, Colombia; Departamento de Matemáticas, Universidad de los Andes, Bogotá, Colombia)
* Santiago Morales (Department of Mathematics, University of California, Davis, California, USA)
* Daniel Yamín (Tepper School of Business, Carnegie Mellon University, Pittsburgh, Pennsylvania, USA)
* Andrés L. Medaglia (Centro para la Optimización y Probabilidad Aplicada (COPA), Departamento de Ingeniería Industrial, Universidad de los Andes, Bogotá, Colombia)

## Abstract
Shortest path problems often arise in contexts where travel times are uncertain. In these settings, reliable paths are often valued more than paths with lower expected travel times. This has led to several variants of reliable shortest path problems (RSPP) that handle travel time reliability differently. We propose an algorithmic framework for solving RSPPs with non-negatively correlated travel times and resource constraints. By building upon the flexibility of the pulse algorithm, our unified and exact algorithmic framework solves multiple variants of the RSPP: the $\alpha$-reliable shortest path ($\alpha$-RSP), the maximum probability of on-time arrival (MPOAP) problem, and the shortest $\alpha$-reliable path (S-$\alpha$RP). We derive a bound on the reliability of path travel times and incorporate three pruning strategies: bounds, infeasibility, and dominance, leveraging properties of the normal distribution and non-negative correlation structures. Computational experiments on large-scale transportation networks (with up to 33 113 nodes and 75 379 arcs) demonstrate that the framework achieves a ten-fold speed improvement over state-of-the-art methods, highlighting its potential real-world applications and extensions to related problems.

## Citation

If you use this repository or our framework, please cite:

```bibtex
@article{Leiva2026Exact,
author = {Leiva, Esteban and Morales, Santiago and Yamín, Daniel and Medaglia, Andrés L.},
title = {An Exact Method for Reliable Shortest Path Problems With Correlation},
journal = {Networks},
volume = {n/a},
number = {n/a},
pages = {},
keywords = {correlation, optimal paths, pulse algorithm, reliable paths, shortest paths, transportation networks},
doi = {https://doi.org/10.1002/net.70026},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/net.70026},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/net.70026},
abstract = {ABSTRACT Shortest path problems often arise in contexts where travel times are uncertain. In these settings, reliable paths are often valued more than paths with lower expected travel times. This has led to several variants of reliable shortest path problems (RSPP) that handle travel time reliability differently. We propose an algorithmic framework for solving RSPPs with non-negatively correlated travel times and resource constraints. By building upon the flexibility of the pulse algorithm, our unified and exact algorithmic framework solves multiple variants of the RSPP: the \$\$ \alpha \$\$-reliable shortest path (\$\$ \alpha \$\$-RSP), the maximum probability of on-time arrival (MPOAP) problem, and the shortest \$\$ \alpha \$\$-reliable path (S-\$\$ \alpha \mathrm{RP} \$\$). We derive a bound on the reliability of path travel times and incorporate three pruning strategies: bounds, infeasibility, and dominance, leveraging properties of the normal distribution and non-negative correlation structures. Computational experiments on large-scale transportation networks (with up to 33 113 nodes and 75 379 arcs) demonstrate that the framework achieves a ten-fold speed improvement over state-of-the-art methods, highlighting its potential real-world applications and extensions to related problems.}
}
```
