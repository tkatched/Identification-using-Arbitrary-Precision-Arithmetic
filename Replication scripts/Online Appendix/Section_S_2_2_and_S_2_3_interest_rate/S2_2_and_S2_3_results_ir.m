%This script replicates the local identification  results in Section 5.2 and Subsections
%2.2-2.3 adding interest rate as observable.

clear




%Local analysis for all the main parameter values

%For serial computation:
tic;
local_analysis_ir
toc;

clear


%Check local identification at the simulated indeterminate parameter values
tic;
local_analysis_indeterminacy_ir
toc;