%This script replicates global identification results for the Leeper et al.
%(2017) medium scale model using business cycle frequencies discussed in 
%Online Appendix Section S.2.5


%Search the AMPF region with PMAF as the benchmark model
run_leeper17_ampf_pmaf_bc



clear

%Search the PMAF region with AMPF as the benchmark model:
run_leeper17_pmaf_ampf_bc


clear

%Search the indeterminacy region with PMAF as the benchmark model:
run_leeper17_pmaf_pmpf_bc


clear

%Search the indeterminacy region with AMPF as the benchmark model:
run_leeper17_ampf_pmpf_bc

