# functional data analysis for matlab

scripts for functional data analysis in matlab
Uses the scripts from the project above to compute spline basis, and then uses custom matlab scripts to conduct functional data analysis on a one-dimensional timeseries
Includes fast vectorised computations of T and F statistic functions, and permutations are computed in parallel.
For a short review on what this is about, see presentation in `docs` folder.

The demo script needs to be updated with some simulated data.

# version history

v0.3 15.03.2023 removed external depedencies for spline computation.

v0.2 14.02.2022 added cluster option for multiple comparisons. Expanded null distribution and data plotting options.

v0.1 included vectorised functional T-test (one, paired or indepedent samples) and functional 2x2 mixed ANOVA. Example script and plotting options.

----------------------------------------------------------------------------------

Functional data analysis code v1.1
Copyright (C) 2022  Dimitris Askitis
dimitrios.askitis@psy.ku.dk

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program. If not, see <http://www.gnu.org/licenses/>.


----------------------------------------------------------------------------------