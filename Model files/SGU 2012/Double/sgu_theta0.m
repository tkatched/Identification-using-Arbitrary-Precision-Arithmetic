%% Parameter values from SGU(2012) Table II posterior median

thet      = 4.74; %Frisch elasticity of labor supply (when gam=b=0)
gam	      = 0.0019; %governs magnitude of wealth elasticity of labor supply
%NB: small values of gamma crucial for importance of news shocks
kap	      = 9.11; %investment adj. cost parameter
del2_del1 = 0.34; %ratio of sensitivity of capacity utilization to variations to rental rate of capital (del2)
%to the steady state capacity utilization
b	      = 0.91; %habit formation parameter
rhoxg     = 0.72; %smoothness in trend of gov. spending

rhoa      = 0.48; %nonstationary investment-specific productivity shock AR
rhox      = 0.38; %neutral productivity shock
rhozi     = 0.47; %investment-specific productivity shock AR
rhoz      = 0.92; %productivity shock AR
rhom      = 0.98; %wage markup shock
rhog      = 0.96; %gov. spending shock
rhozet    = 0.17; %preference shock




%standard deviations of the 21 shocks
sig_mua0 = 0.21; %
sig_mua4 = 0.16;
sig_mua8 = 0.16;
sig_mux0 = 0.38;
sig_mux4 = 0.08;
sig_mux8 = 0.1;
sig_zi0  = 11.72;
sig_zi4  = 1.93;
sig_zi8  = 5.50;
sig_z0   = 0.65;
sig_z4   = 0.11;
sig_z8   = 0.09;
sig_mu0  = 0.5;
sig_mu4  = 4.79;
sig_mu8  = 0.51;
sig_g0   = 0.62;
sig_g4   = 0.57;
sig_g8   = 0.37;
sig_zet0 = 4.03;
sig_zet4 = 1.89;
sig_zet8 = 2.21;

%standard deviation of the measurement error
sig_me   = 0.30;

theta0=[thet;gam;kap;del2_del1;b;rhoxg;rhoa;rhox;rhozi;rhoz;rhom;rhog;rhozet;sig_mua0;sig_mua4;
    sig_mua8;sig_mux0;sig_mux4;sig_mux8;sig_zi0;sig_zi4;sig_zi8;sig_z0;sig_z4;sig_z8;sig_mu0;sig_mu4;sig_mu8;
    sig_g0; sig_g4;sig_g8;sig_zet0;sig_zet4;sig_zet8;sig_me];

clearvars -except theta0
save theta0_sgu
