%This scripts replicates results from the Online Appendix Section S.2.4

%% The nonidentification curves from Tables S2-S5

%Extend the nonidentification curve from Table S2
curve_pmaf1_dir1 %direction 1

clear

curve_pmaf1_dir2 %direction 2

clear

%Verify the values on the nonidentification curve in Table S2 result in
%observational equivalence

pmaf1_pmpf_endcurve_dir1 %direction1

clear

pmaf1_pmpf_endcurve_dir2 %direction 2

clear



%Extend the nonidentification curve from Table S3
curve_pmaf2_dir1 %direction 1

clear

curve_pmaf2_dir2 %direction 2

clear

%Verify the values on the nonidentification curve in Table S3 result in
%observational equivalence

pmaf2_pmpf_endcurve_dir1 %direction1

clear

pmaf2_pmpf_endcurve_dir2 %direction 2

clear

%Extend the nonidentification curve from Table S4
curve_pmaf3_dir1 %direction 1

clear

curve_pmaf3_dir2 %direction 2

clear

%Verify the values on the nonidentification curve in Table S4 result in
%observational equivalence

pmaf3_pmpf_endcurve_dir1 %direction1

clear

pmaf3_pmpf_endcurve_dir2 %direction 2

clear

%% Global identification analysis considering the noninvertible MA process for er_t

%Search the indeterminacy region with AMPF1 as the benchmark model imposing
%zero restrictions on the last two parameters and noninvertibility of er_t
run_ampf1_pmpf_noninv_restr


%Verify results above using arbitrary precision
run_ampf1_pmpf_noninv_restr_mp

clear


%Verify results above using arbitrary precision
run_ampf1_pmpf_restr_mp


clear

%Search the indeterminacy region with PMAF1 as the benchmark model imposing
%zero restrictions on the last two parameters and noninvertibility of er_t
run_pmaf1_pmpf_noninv_restr

clear

%Verify results above using arbitrary precision
run_pmaf1_pmpf_noninv_restr_mp

