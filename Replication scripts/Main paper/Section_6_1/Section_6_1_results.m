%This script replicates the main results in Subsection 6.1 (local identification
%of the medium scale model)

clear

%% Local identification using observables in levels

sgu_local_levels %utilizes serial G-matrix computation

clear

sgu_local_levels_par %utilizes parallel G-matrix computation

clear
%% Local identification using the original SGU(2012) observables

sgu_local_original_obs %utilizes serial G-matrix computation

clear

sgu_local_original_obs_par %utilizes parallel G-matrix computation