%% Script reproducing global identification analysis at the ML estimate
%from SGU (2012)Table II

%Script naming conventions follow Section 6.2 

%% Neighborhood search with c=0.1, 0.5,  1
run_sgul_nsearch_ml

clear

run_sgul_nsearch_1f_ml % neighborhood search with weakly identified parameter sig_zet_4 fixed

clear

%% Sensitivuty to shutting down news shocks

run_sgul_nonews_ml

clear

run_sgul_nonews4_ml

clear

run_sgul_nonews8_ml

clear

%Shutting down 8-period news 1-by-1:

run_sgul_non8_1_ml

clear

run_sgul_non8_2_ml

clear

run_sgul_non8_3_ml

clear

run_sgul_non8_4_ml

clear

run_sgul_non8_5_ml

clear

run_sgul_non8_6_ml

clear

run_sgul_non8_7_ml

clear

%% Shutting down news shocks 1-by-1 

run_sgul_nonews_g_ml

clear

run_sgul_nonews_mua_ml 

clear

run_sgul_nonews_mur_ml

clear

run_sgul_nonews_mux_ml

clear

run_sgul_nonews_z_ml 

clear

run_sgul_nonews_zet_ml

clear

run_sgul_nonews_zi_ml

%% Replacing news shocks by ARMA shocks

run_sgul_ma2_inv2pns_ml
clear

run_sgul_ma2_msel2_ml 

clear

run_sgul_ma2_msel3_ml

clear

run_sgul_ma2_msel4_ml

clear

run_sgul_ma2_msel5_ml

clear

run_sgul_ma2_msel6_ml

clear

run_sgul_ma2_msel7_ml

clear

run_sgul_ma2_msel8_ml

