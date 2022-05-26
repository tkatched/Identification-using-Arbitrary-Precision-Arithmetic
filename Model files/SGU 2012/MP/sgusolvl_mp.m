%function [T1,TC,T0,RC]=sgusolv(para)   
function [T1,TC,T0,TETA,RC] =sgusolvl_mp(para)

%% parameter vector (must be an MP object)
thet      = para(1);
gam	      = para(2);
kap	      = para(3);
del2_del1 = para(4);
b	      = para(5);
rhoxg     = para(6);
rhoa      = para(7);
rhox      = para(8);
rhozi     = para(9);
rhoz      = para(10);
rhom      = para(11);
rhog      = para(12);
rhozet    = para(13);


%%
neq=91; %no of equations
neta=14; %no of expectational errors
neps=21; %no of shocks
%% Compute SS parameters
sgusspar_mp
%% Label variables, equations, shocks etc.
sgueqvarshocksl
%% Put the model into gensys form
sgusetsysmatl_mp
%% Solve model using gensys

[T1,TC,T0,TY,M,TZ,TETA,GEV,RC]=gensys_mp2(GAM0, GAM1, C, PSI, PPI, onep); %use Sims-based code, also uses ordqz



end