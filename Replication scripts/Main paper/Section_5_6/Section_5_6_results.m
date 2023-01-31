%This script replicates global identification results for the Leeper et al.
%(2017) medium scale model using full spectrum discussed in Section 5.6 and
%Online Appendix Section S.2.5


%Search the AMPF region with PMAF as the benchmark model
run_leeper17_ampf_pmaf

clear

%Verify results above using arbitrary precision
run_leeper17_ampf_pmaf_mp

clear

%Search the PMAF region with AMPF as the benchmark model:
run_leeper17_pmaf_ampf

clear

%Verify results above using arbitrary precision
run_leeper17_pmaf_ampf_mp

clear

%Search the indeterminacy region with PMAF as the benchmark model:
run_leeper17_pmaf_pmpf

clear

%Verify results above using arbitrary precision
run_leeper17_pmaf_pmpf_mp

clear

%Search the indeterminacy region with AMPF as the benchmark model:
run_leeper17_ampf_pmpf

clear

%Verify results above using arbitrary precision
run_leeper17_ampf_pmpf_mp