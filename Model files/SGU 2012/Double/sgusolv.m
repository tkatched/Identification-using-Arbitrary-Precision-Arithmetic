%function [T1,TC,T0,RC]=sgusolv(para)   
function [T1,TC,T0,TETA,RC] =sgusolv(para)

%% parameter vector
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
neq=98; %no of equations
neta=14; %no of expectational errors
neps=21; %no of shocks
%% Compute SS parameters
sgusspar
%% Label variables, equations, shocks etc.
sgueqvarshocks
%% Put the model into gensys form
sgusetsysmat
%% Solve model using gensys

[T1,TC,T0,TY,M,TZ,TETA,GEV,RC]=gensys_mod2(GAM0, GAM1, C, PSI, PPI, 1); %use Sims-based code, also uses ordqz

%[T1,TC,T0,M,TZ,TY,GEV,RC] = gensys(GAM0,GAM1,C,PSI,PPI);

end