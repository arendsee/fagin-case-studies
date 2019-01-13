[![unstable](http://badges.github.io/stability-badges/dist/unstable.svg)](http://github.com/badges/stability-badges)

# A case study for `fagin`

This repo contains the code needed to generate the tables and figures used in
the synder paper and its supplementary files. And more.

See [main fagin page here](https://github.com/arendsee/fagin)

# Dependencies

CRAN R packages:
 * knitr
 * magrittr
 * rmonad

GitHub R Packages:
 * arendsee/fagin
 * arendsee/synder

The CRAN packages can be installed from an R shell with as so

``` R
install.packages("knitr")
```

and the GitHub packages can be installed with devtools:

``` R
library(devtools)
install_github("arendsee/fagin")
```

# Data Retrieval

All data needed to run this program is available on DataHub. To retrieve this
data, first install the data CLI tool by following the instructions available
[here](https://datahub.io/download). Then run the command:

``` sh
make init
```

This will download all the input dataset from DataHub and then organize it the
way `fagin` wants.
