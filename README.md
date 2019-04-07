# twoe

`twoe` (2e) aims at estimating the demographic parameters of tropical tree
species from permanent forest plot data. It includes functions to estimate growth,
mortality and recruitment of tree species. Species are treated as random effects
which allow estimating parameters for rare species.

## Compilation

This package uses the [Scythe Statistical Library](http://scythe.lsa.umich.edu/), an open source C++ library for statistical computation. We suggest using the GCC compiler 4.0 or greater. The current package has been tested using GCC 4.0 on Linux. 

## Installation

`library(devtools)`    
`install_github("ghislainv/twoe")`
