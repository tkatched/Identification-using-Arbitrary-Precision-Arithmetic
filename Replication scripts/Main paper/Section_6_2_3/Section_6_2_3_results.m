%% Distinguishability of models with ARMA shocks (Subsection 6.2.3)

run_sgul_ma2_inv2pns %find the closest model with 4th and 8th MA lags shocks

clear

run_sgul_ma2_inv2pns_mp %verify the results above using quadruple precision

clear

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replacing news shocks with MA lags one by one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_sgul_ma2_msel2 %find the closest model with the MA nonstationary investment-specific productivity shock

clear

run_sgul_ma2_msel3 %find the closest model with the MA nonstationary neutral productivity shock

clear

run_sgul_ma2_msel4 %find the closest model with the MA stationary investment-specific productivity shock

clear

run_sgul_ma2_msel5 %find the closest model with the MA stationary neutral productivity shock

clear

run_sgul_ma2_msel6 %find the closest model with the MA wage markup shock

clear

run_sgul_ma2_msel7 %find the closest model with the MA gov. spending shock

clear


run_sgul_ma2_msel8 %find the closest model with the MA preference shock

clear

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance decompositions of 4-period shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SGUvar  %compute variance decomposition for the wage markup 
% 4-period news shock of SGU (2012) at the posterior median and the 
%decomposition of the MA wage markup shock

clear
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse responses to anticipated and unanticipated markup shocks 
% in the baseline SGU (2012) vs the 4-period MA specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SGUimp %compute impulse responses for the wage markup 
% unanticipated, 4-period and 8-period news shocks of SGU (2012) at the 
% posterior median and the impulse response of the MA wage markup shock