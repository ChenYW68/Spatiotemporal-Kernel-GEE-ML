# Estimating semiparametric models with complex spatiotemporal dependence

This Github page provides code and data for reproducing the results in the manuscript: ``Estimating semiparametric models with complex spatiotemporal dependence'' by Y. Chen, T. Chu, C. Zhou, Y. Shen, and H. Huang.

## Summary
For spatiotemporal data, existing semiparametric models often do not fully account for its complex dependence structure, leading to potential efficiency losses and methodological challenges in inference and computation. To address these issues, we propose a novel framework that incorporates spatiotemporal dependence into the estimation procedures for both parametric and nonparametric components via a K-nearest-neighbor weighted local linear regression. In addition, the proposed approach does not require assuming non-stationarity of the space-time covariance, which makes the modeling more flexible. Under an increasing-domain asymptotic framework, we show that the bias and efficiency of the nonparametric estimator can be significantly improved, while the parametric counterpart still achieves consistency and attains optimal efficiency if the covariance is correctly specified. Simulation studies further demonstrate the finite-sample performance and robustness of the proposed methods, even when the covariance or data-generating process is misspecified. An application to air pollution data illustrates its practical effectiveness.

## Software package
We developed an R package to implement the proposed K-nearest-neighbor weighted local linear regression (KNN-WLLR) in `knnWLLR`; see [`knnWLLR_1.0.zip`](knnWLLR_1.0.zip). All simulation code can be found in the code folder to facilitate reproducibility.


## Simulation results
Figure 1 provides a thorough insight into the performance of each method by showing trajectory plots of the absolute bias, standard deviation, and MSE at each sampling point when estimating the nonparametric functions. Our proposed KNN-WLLR consistently outperforms the competing methods across all three metrics at nearly all sampling points.
<figure id="Figure1">
  <table align="center">
    <tr>
      <td><img src="./figure/Fig2_trajectory.jpg" width="800px"></td>
    </tr>
  </table>
  <figcaption align="center">
    <strong>Figure 1:</strong>  Pointwise absolute bias, standard deviation and mean square error for each method in estimating nonparametric functions.
  </figcaption>
</figure>


## Real data analysis
We analyze daily concentration data of PM2.5 in China's BTH region from November 1 to November 30, 2015. The PM2.5 concentrations are from two sources: (i) readings at 68 spatially sparse monitoring sites and (ii) 
outputs of the Community Multiscale Air Quality (CMAQ) model, a widely used numerical modeling system. The detailed description of the data can be found in the [paper](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-18/issue-2/Efficient-and-effective-calibration-of-numerical-model-outputs-using-hierarchical/10.1214/23-AOAS1823.short) and the [paper](https://www.sciencedirect.com/science/article/abs/pii/S016794732300110X?casa_token=k4sdX2uq82AAAAAA:d7fIbx359tQG0p0V10O3OhEeT19oNDgnq0erS1fBSm97WkAY_o0Tc7Sqy53fhH_HYZOypBWfezw). These data are publicly available on [GitHub](https://github.com/ChenYW68/HDCM). Nonstationary patterns in space and time are observed based on the estimated covariance via the machine learning method.
<figure id="Figure2">
  <table align="center">
    <tr>
      <td><img src="./figure/BTH_Covariance.jpg" width="800px"></td>
    </tr>
  </table>
  <figcaption align="center">
    <strong>Figure 1:</strong> Spatiotemporal covariance estimated via the deeper learning method.
  </figcaption>
</figure>
