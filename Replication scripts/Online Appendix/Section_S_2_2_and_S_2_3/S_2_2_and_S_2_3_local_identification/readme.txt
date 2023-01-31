List of files in the folder:



Matlab scripts:

Sections_S2_2_and_S2_3_results.m - a master script for replicating results in the Sections 2.2  and 2.3 of the appendix. For completeness, all the points discussed in the appendix and the main paper are considered at both full spectrum and business cycle frequencies.

local_analysis.m - script that conducts local identification checks for all considered cases in double, quadruple and 50-digit precision

local_analysis_parallel.m - same as local_analysis.m above, but with the parallelized computation of the G-matrix used

local_analysis_indeterminacy.m  - local identification checks in double, quadruple and 50-digit precision for 1000 simulated points from the indeterminacy region.

Simulate_ind.m - simulates points from the indeterminacy region used in local_analysis_indeterminacy.m, saved in local_ind.mat.

Matlab data files:

local_ind.mat - 1000 simulated points for checking local identification in the indeterminacy region. Used in local_analysis_indeterminacy.m.


