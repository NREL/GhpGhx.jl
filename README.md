# GhpGhx.jl

## Description
This package is used for sizing a ground heat exchanger to serve the specified heating and cooling loads. This package is called by [REopt.jl](https://github.com/NREL/REopt.jl) when GHP is considered in a REopt evaluation.

## License: restricted use
This julia package, repository, and in particular the executable (shared object, compiled fortran code) [tess.so](https://github.com/NREL/GhpGhx.jl/blob/main/ghxmodel/tess.so) file that contains ground heat exchanger (GHX) calculations are governed by a more restricive [license](https://github.com/NREL/GhpGhx.jl/blob/main/LICENSE.md) compared to REopt.jl, and it dictates that this GHX model should only be used in the context of REopt analysis, and it cannot be used for commercial purposes.
