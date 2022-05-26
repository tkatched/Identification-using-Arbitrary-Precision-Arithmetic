%% Parameter values from SGU(2012) Table II posterior median

thet      = mp('4.74'); %Frisch elasticity of labor supply (when gam=b=0)
gam	      = mp('0.0019'); %governs magnitude of wealth elasticity of labor supply
%NB: small values of gamma induce RC(1,1)=0 in gensys, which violates
%existence condition, although results are still produced. The threshold
%for this at the posterior median of SGU(2012) seems to be around
%gam=0.0066. Smaller values produce RC=[0;1]. 
%HS(2014) seemingly set gam=0.1 in their code when using the SGU(2012) posterior
%median parameter value

kap	      = mp('9.11'); %investment adj. cost parameter
del2_del1 = mp('0.34'); %ratio of sensitivity of capacity utilization to variations to rental rate of capital (del2)
%to the steady state capacity utilization
b	      = mp('0.91'); %habit formation parameter
rhoxg     = mp('0.72'); %smoothness in trend of gov. spending

rhoa      = mp('0.48'); %nonstationary investment-specific productivity shock AR
rhox      = mp('0.38'); %neutral productivity shock
rhozi     = mp('0.47'); %investment-specific productivity shock AR
rhoz      = mp('0.92'); %productivity shock AR
rhom      = mp('0.98'); %wage markup shock
rhog      = mp('0.96'); %gov. spending shock
rhozet    = mp('0.17'); %preference shock




%standard deviations of the 21 shocks
sig_mua0 = mp('0.21'); %
sig_mua4 = mp('0.16');
sig_mua8 = mp('0.16');
sig_mux0 = mp('0.38');
sig_mux4 = mp('0.08');
sig_mux8 = mp('0.1');
sig_zi0  = mp('11.72');
sig_zi4  = mp('1.93');
sig_zi8  = mp('5.50');
sig_z0   = mp('0.65');
sig_z4   = mp('0.11');
sig_z8   = mp('0.09');
sig_mu0  = mp('0.5');
sig_mu4  = mp('4.79');
sig_mu8  = mp('0.51');
sig_g0   = mp('0.62');
sig_g4   = mp('0.57');
sig_g8   = mp('0.37');
sig_zet0 = mp('4.03');
sig_zet4 = mp('1.89');
sig_zet8 = mp('2.21');

%standard deviation of the measurement error
sig_me   = mp('0.30');

theta0mp=[thet;gam;kap;del2_del1;b;rhoxg;rhoa;rhox;rhozi;rhoz;rhom;rhog;rhozet;sig_mua0;sig_mua4;
    sig_mua8;sig_mux0;sig_mux4;sig_mux8;sig_zi0;sig_zi4;sig_zi8;sig_z0;sig_z4;sig_z8;sig_mu0;sig_mu4;sig_mu8;
    sig_g0; sig_g4;sig_g8;sig_zet0;sig_zet4;sig_zet8;sig_me];

clearvars -except theta0mp
save theta0mp_sgu
