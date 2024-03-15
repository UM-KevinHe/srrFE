# srrFE

The `srrFE` software package is designed to implement the SerBIN method proposed by [Wu et al. (2022)](doi:10.1002/sim.9387), specifically tailored to address the computational challenges inherent in large-scale provider data. Additionally, the package is equipped with functions designed to estimate standardized ratios and rates and conduct hypothesis testing.

## Introduction

Provider profiling is acknowledged for its utility in monitoring healthcare quality, facilitating inter-provider care coordination, and enhancing medical cost-effectiveness. Existing methods often use generalized linear models with fixed provider effects, especially when profiling dialysis facilities. However, as the number of providers increases, the computational burden becomes challenging. To address this issue, we introduce a serial blockwise inversion Newton (SerBIN) algorithm that capitalizes on the block structure of the information matrix. Additionally, we propose a shared-memory divide-and-conquer algorithm to further enhance computational efficiency.

Beyond the computational challenge, the current literature lacks a suitable inferential approach for identifying providers with outlier performance, especially in cases involving small providers with extreme outcomes. Traditional score and Wald tests, relying on large-sample distributions of test statistics, inaccurately approximate small-sample properties in this context. To address this inferential gap, we present an exact test of provider effects using exact finite-sample distributions, with the Poisson-binomial distribution serving as a special case when the outcome is binary. This methodological advancement aims to provide a more accurate and reliable assessment of provider performance, considering both computational and inferential challenges.

## Installation

* To install the latest development version from GitHub:
```
require("devtools")
require("remotes")
remotes::install_github("UM-KevinHe/srrFE", ref = "main")
```

## Getting started

See ["Getting Started with srrFE"](https://um-kevinhe.github.io/srrFE/articles/srrFE.html)

## Details of the algorithms

See ["Models"](https://um-kevinhe.github.io/srrFE/articles/Models.html) or [original articles](https://um-kevinhe.github.io/srrFE/articles/Articles.html)

## Getting Help

If you encounter any problems or bugs, please contact us at: [ybshao\@umich.edu](mailto:ybshao@umich.edu){.email}, [kevinhe\@umich.edu](mailto:kevinhe@umich.edu){.email}, [Wenbo.Wu\@nyulangone.org](mailto:Wenbo.Wu@nyulangone.org){.email},
[taoxu\@med.umich.edu](mailto:taoxu@med.umich.edu){.email}

## References

* [Wu, W., Yang, Y., Kang, J., & He, K. (2022). Improving large‐scale estimation and inference for profiling health care providers. Statistics in Medicine, 41(15), 2840-2853.](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9387)
* [He, K., Kalbfleisch, J. D., Li, Y., & Li, Y. (2013). Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects. Lifetime Data Analysis, 19, 490-512.](https://link.springer.com/article/10.1007/s10985-013-9264-6)
