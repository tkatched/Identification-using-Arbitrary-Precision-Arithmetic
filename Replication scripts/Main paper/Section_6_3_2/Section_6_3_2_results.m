%% Identification of anticipated shocks (Subsection 6.3.2, Appendix Tables S12-S13)


run_sgul_nonews %find the closest model with all the news shocks shut down

clear

run_sgul_nonews_mp %verify the results above using quadruple precision

clear

run_sgul_nonews4 %find the closest model with all the 4-period news shocks shut down

clear

run_sgul_nonews4_mp %verify the results above using quadruple precision

clear

run_sgul_nonews8 %find the closest model with all the 8-period news shocks shut down

clear

run_sgul_nonews8_mp %verify the results above using quadruple precision

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shutting down news shock-by-shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_sgul_nonews_g %find the closest model with the gov. spending shock shut down

clear

run_sgul_nonews_g_mp %verify the results above using quadruple precision

clear

run_sgul_nonews_mua %find the closest model with the nonstationary investment-specific productivity shock shut down

clear

run_sgul_nonews_mua_mp %verify the results above using quadruple precision

clear

run_sgul_nonews_mur %find the closest model with the wage markup shock shut down

clear

run_sgul_nonews_mur_mp %verify the results above using quadruple precision

clear

run_sgul_nonews_mux %find the closest model with the nonstationary neutral productivity shock shut down

clear

run_sgul_nonews_mux_mp %verify the results above using quadruple precision

clear

run_sgul_nonews_z %find the closest model with the stationary neutral productivity shock shut down

clear

run_sgul_nonews_z_mp %verify the results above using quadruple precision

clear

run_sgul_nonews_zet %find the closest model with the preference shock shut down

clear

run_sgul_nonews_zet_mp %verify the results above using quadruple precision

clear

run_sgul_nonews_zi %find the closest model with the stationary investment-specific productivity shock shut down

clear

run_sgul_nonews_zi_mp %verify the results above using quadruple precision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shutting down 8-period news shock-by-shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_sgul_non8_1 % %find the closest model with 8-period news in the nonstationary investment-specific productivity shock shut down

clear

run_sgul_non8_1_mp %verify the results above using quadruple precision

clear

run_sgul_non8_2 %find the closest model with 8-period news in the nonstationary neutral productivity shock shut down

clear

run_sgul_non8_2_mp %verify the results above using quadruple precision

clear

run_sgul_non8_3 %find the closest model with 8-period news in the stationary investment-specific productivity shock shut down

clear

run_sgul_non8_3_mp %verify the results above using quadruple precision

clear

run_sgul_non8_4 %find the closest model with 8-period news in the stationary neutral productivity shock shut down

clear

run_sgul_non8_4_mp %verify the results above using quadruple precision

clear

run_sgul_non8_5 %find the closest model with 8-period news in the wage markup shock shut down

clear

run_sgul_non8_5_mp %verify the results above using quadruple precision

clear

run_sgul_non8_6 %find the closest model with 8-period news in the gov. spending shock shut down

clear

run_sgul_non8_6_mp %verify the results above using quadruple precision

clear

run_sgul_non8_7 %find the closest model with 8-period news in the preference shock shut down

clear

run_sgul_non8_7_mp %verify the results above using quadruple precision