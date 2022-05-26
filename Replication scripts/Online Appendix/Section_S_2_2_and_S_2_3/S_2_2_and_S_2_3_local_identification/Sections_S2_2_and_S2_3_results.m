%This script replicates the main results in Subsections S.2.2 and S.2.3 of the online appendix
% (local identification of the small scale model)

clear




%Local analysis for all the main parameter values

%For serial computation:
tic;
local_analysis
toc;

clear


%Same as above, but the G-matrix code has some parallelized steps:
tic;
local_analysis_parallel
toc;

clear

%Simulate and save points from the indeterminacy region
tic;
Simulate_ind
toc;

clear

%Check local identification at the simulated indeterminate parameter values
tic;
local_analysis_indeterminacy
toc;