
function [T1,TC,T0,TETA,RC] =sgusolvl3(para,msel)


%Allows solving for MA(1) shock specification. 
%MA(1) shocks reparameterized as Y_t=theta1*e_t+theta2*e_t-1
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

if msel==2 %ma*4 corrseponds to t-4 lag, ma*8 - to t-8 lag
    mamua0=para(14);
    mamua4=para(15);
    mamua8=para(16);
    mamux0=para(17);
    mamux4=para(18);
    mamux8=para(19);
    mazi0=para(20);
    mazi4=para(21);
    mazi8=para(22);
    maz0=para(23);
    maz4=para(24);
    maz8=para(25);
    mamu0=para(26);
    mamu4=para(27);
    mamu8=para(28);
    mag0=para(29);
    mag4=para(30);
    mag8=para(31);
    mazet0=para(32);
    mazet4=para(33);
    mazet8=para(34);
end
%%
neq=91; %no of equations
neta=14; %no of expectational errors
if msel==1
neps=21; %no of shocks
elseif msel==2
    neps=7;
end
%% Compute SS parameters
sgusspar
%% Label variables, equations, shocks etc.
sgueqvarshocksl2
%% Put the model into gensys form
sgusetsysmatl3
%% Solve model using gensys

[T1,TC,T0,TY,M,TZ,TETA,GEV,RC]=gensys_mod2(GAM0, GAM1, C, PSI, PPI, 1); %use Sims-based code, also uses ordqz

%[T1,TC,T0,M,TZ,TY,GEV,RC] = gensys(GAM0,GAM1,C,PSI,PPI);

end