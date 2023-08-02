
# Repo for the paper “A note on joint calibration estimators for totals and quantiles”

## Preliminaries

Install `jointCalib` package

``` r
install.packages("remotes")
remotes::install_github("ncn-foreigners/jointCalib")
remotes::install_github("ncn-foreigners/nonprobsvy") ## for IPW with calibration constraints
```

Then, restart R and load the package

``` r
library(jointCalib)
library(nonprobsvy)
```

## Notebooks

1.  [Main simulation from the paper]() based on Chen, Y., Li, P., and
    Wu, C. (2020). *Doubly robust inference with nonprobability survey
    samples*. Journal of the American Statistical Association,
    115(532):2011–2021
2.  [Tutorial for the `jointCalib`
    package](https://github.com/ncn-foreigners/jointCalib)
3.  [A minimal example of joint calibration in R, Python and
    stata](https://htmlpreview.github.io/?https://raw.githubusercontent.com/ncn-foreigners/paper-note-joint-calibration/main/notebooks/qcalib-r-python-stata.html)

## Structure of the repo

- `codes/` – code for simulation
- `figs/` – figures that correspondst to table in the paper
- `notebooks/` – notebooks with processing results
- `results/` – files with tables based on simulation

## Funding

Work on this package is supported by the the National Science Center,
OPUS 22 grant no. 2020/39/B/HS4/00941.
