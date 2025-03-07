# Approximation of aggregate curves in the electricity market relying on node selection

**Purpose** - Propose an approximation procedure to efficiently represent aggregate supply and demand curves in the electricity market. 

**Design/methodology/approach** - Two approximation procedures based on one-step functions are designed using the $\mathcal{L}_1$ and $\mathcal{L}_2$ metrics as error measures. For the metric $\mathcal{L}_2$ a closed-form solution is obtained to reduce the global error over the entire domain. For the metric $\mathcal{L}_1$, which lacks a closed-form solution, a linear programming problem is solved to improve the local approximation within intervals defined by the nodes. The dependence on the node locations is addressed by different node selection strategies. In particular, a heuristic strategy is proposed that combines descriptive information from the offers with a dyadic search procedure to minimize the approximation error.

**Findings** - The performance of our proposals is evaluated and compared using curves from the day-ahead Spanish electricity market. Our procedure achieves a promising approximation performance compared with existing approaches.

**Originality** - Existing procedures for representing step supply and demand curves often involve high computational costs or, alternatively, make assumptions of smoothness and differentiability of the curves that contradict the nature of these step curves. Our proposals address these challenges, obtaining approximations that resemble real curves in a parsimonious and practical way.

**Practical implications** The proposed procedures will allow the development of efficient methods for forecasting supply and demand curves and, consequently, market clearing prices. Having these forecasts is of interest to both producers and consumers.

**keywords:** Functional approximation; aggregate curve; electricity market

## About the curve data
The dataset is published on the open science platform, Zenodo, at https://zenodo.org/records/14558660. Overcoming the limitation of Github on the size of files falls out of our research scope:)
