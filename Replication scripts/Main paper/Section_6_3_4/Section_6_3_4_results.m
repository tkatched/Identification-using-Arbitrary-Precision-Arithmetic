%% Replication of results in  Section 6.3.4: 
% searching for the equivalent news shock model of SGU(2012) with a noise
% model of Chahrour and Jurado (2018) as benchmark.

%The search is performed over 21 shock parameters since noise and news are
%equivalent shock representations under conditions in Chahrour and Jurado (2018).

% Search for equivalent news representation in double precision:
runs_sgul_noise_orig_shk

clear

%Verify results in arbitrary precision:

run_sgul_noise_orig_shk_ga_mp