This repository contains simulated data and real world data used in the manuscript
entitled **Zero Inﬂated Poisson Model with Clustered Regression Coefﬁcients: an
Application to Heterogeneity Learning of Field Goal Attempts of Professional
Basketball Players** by Guanyu Hu, Hou-Cheng Yang, Yishu Xue and Yuqing Pan.

* `nmf_basis.rds`: a 1175x5 matrix with each column being one basis function, visualized in Figure 2 of the manuscript.
* `real_data.rds`: a 191x1175 matrix with each row being the shots made by one player on the 1175 blocks on the court; selected players' data visualized in Figure 3 of the manuscript.
* `simulated_balance.rds`: simulated data under the balanced design.
* `simulated_imbalance.rds`: simulated data under the imbalanced design.
* `run_simu_balanced.R`: code that runs simulation for the balanced design.
* `docs/index.html`: a rendered version with explanation for code, view at [https://ys-xue.github.io/MFM-ZIP-Basketball-Supplemental/](https://ys-xue.github.io/MFM-ZIP-Basketball-Supplemental/).