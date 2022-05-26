%This script replicates the main results in Subsection 5.3  and 5.4 (global identification
%of the small scale model and the corresponding nonidentification curves, empirical distances)

clear

%% Observational equivalence between AMPF and PMPF regimes

%Search the indeterminacy region with AMPF1 as the benchmark model
run_ampf1_pmpf

clear

%Verify results above using arbitrary precision
run_ampf1_pmpf_mp

clear

%Search the indeterminacy region with AMPF1 as the benchmark model imposing
%zero restrictions on the last two parameters
run_ampf1_pmpf_restr

clear

%Verify results above using arbitrary precision
run_ampf1_pmpf_restr_mp

clear

%Extend the nonidentification curve from Table 1

curve_ampf_dir1 %direction 1

clear

curve_ampf_dir2 %direction 2

clear
%Verify the values on the nonidentification curve in Table 1 result in
%observational equivalence

ampf_pmpf_endcurve_dir1 %direction1

clear

ampf_pmpf_endcurve_dir2 %direction 2

clear
%Extend  the nonidentification curve from Table 2

curve_ampf2_dir1 %direction 1

clear

curve_ampf2_dir2 %direction 2

clear



%Verify the values on the nonidentification curve in Table 2 result in
%observational equivalence

ampf_pmpf_endcurve2_dir1 %direction1

clear

ampf_pmpf_endcurve2_dir2 %direction 2

clear
%% Observational equivalence between PMAF and PMPF regimes

%Search the indeterminacy region with PMAF1 as the benchmark model
run_pmaf1_pmpf

clear

%Search the indeterminacy region with AMPF1 as the benchmark model imposing
%zero restrictions on M_r and sunspot variance.
run_pmaf1_pmpf_restr_mp

clear

%% Observational equivalence between AMPF and PMAF regimes

%Search the PMAF region with AMPF1 as the benchmark model, report empirical
%distances
run_ampf1_pmaf

clear

%Verify results above using arbitrary precision
run_ampf1_pmaf_mp

clear

%Search the AMPF region with PMAF1 as the benchmark model, report empirical
%distances
run_pmaf1_ampf

clear

%Verify results above using arbitrary precision
run_pmaf1_ampf_mp