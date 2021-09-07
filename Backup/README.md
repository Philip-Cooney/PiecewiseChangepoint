---
output:
  html_document: default
  pdf_document: default
---

# PiecewiseChangepoint

<!-- badges: start -->
<!-- badges: end -->

The goal of PiecewiseChangepoint is to estimate the number and locations of change-points in pieceise exponential models. 

## Installation

You can install the released version of PiecewiseChangepoint from [GitHub](https://github.com/Philip-Cooney/PiecewiseChangepoint) with:

``` r
devtools::install_github("Philip-Cooney/PiecewiseChangepoint")
```

## Example

First we load the package and simulate some piecewise exponential data. 

``` {r}
library(PiecewiseChangepoint)
## basic example code

set.seed(123)
n_obs =300
n_events_req=300
max_time =  2

rate = c(0.75,0.25)
t_change =1

df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)
                   
head(df)
```

``` r

Collapsing.Model <- collapsing.model(df,
                                     n.iter = 5000,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)


```

