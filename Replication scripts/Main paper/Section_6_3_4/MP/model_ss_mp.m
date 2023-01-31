% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)


function [ss,param,set] = model_ss_mp(param,set)

%Upack parameters object
param_unpack

onep=mp('1');
twop=mp('2');
zp=mp('0');

who
%BEGIN_EXTRACT_HERE
zY     = onep;
zI     = onep;
zeta   = onep;
u      = onep;
q      = onep;
y_k    = onep/(alphak*muy)*(muy^sig/bet - mua*(onep-delta0));
delta1 = alphak*y_k*mui;
delta2 = rdelta*delta1;
k      = (y_k*mui^alphak/h^alphah*L^(alphak+alphah-onep))^(onep/(alphak-onep));
yd     = y_k*k;%0.395274610790875
iv     = (onep - (onep-delta0)/mui)*k;% 0.068352497249448
g   = G_Y*yd/xg;%0.080394529501488
c      = yd - iv - g*xg; % 0.247867191383252
s      = c*(onep - b/muy)*muy^((gam-onep)/gam);% 0.002427174753172

aux1   = h^theta/(onep-bet*(onep-gam)*muy^(onep-sig));
aux2   = alphah*yd/muw/(theta*h^theta*s);
aux3   = aux2*(onep-bet*b*muy^(-sig));
aux4   = gam*s/(c*(onep-b/muy));
psii    = aux3/(onep+aux1*aux3*aux4);
%term1  = gam*h^theta*s/(c*(1-b/muy))/(1-bet*(1-gam)*muy^(1-sig));
%term2  = theta*h^theta*s*muw/(alphah*yd*(1-bet*b/muy^sig));
%psii    = 1/(term1 + term2); %8.569365092071612e+003
v      = c*(onep-b/muy) - psii*h^theta*s; %0.015306293463611
pai    = psii*h^theta/(onep - bet*(onep-gam)*muy^(onep-sig))/v^sig; %7.849043973959898e+003
lambda = psii*theta*h^theta*s*muw/(alphah*yd*v^sig);%5.384564765310337

% Implied steady-state values
%mua    = muabar;
muy    = muybar; %1.004500000000000
%muw    = muwbar;
mux    = muy*mua^(alphak/(onep-alphak));
%muxbar = mux;
mui    = mux*mua^(onep/(alphak-onep));
h      = hbar;
xg     = muy^(onep/(rhoxg-onep));
delta  = delta0 + delta1*(u-onep) + delta2/twop*(u-onep)^twop;
phi    = zp;
dphi   = zp;

% Observables
gy   = muy;
gc   = muy;
gai  = muy;
gh   = onep;
gg   = muy;
gtfp = mux^(onep-alphak);
ga   = mua;
gbar = g;


E1PT = -(- bet*b*(zeta*(v*muy)^(-sig) - gam*pai/muy^sig*(s/(c*muy - b*c))^(onep-gam)));
E2PT = (bet*(onep-gam)*pai*muy^(-sig)*(s/(c*muy - b*c))^(-gam));
E3PT = (bet*muy^(-sig)*(alphak*yd/k*muy + q*mua*(onep-delta)));
E4PT = (bet*muy^(-sig)*q*mua*zI*dphi*(iv*mui/iv)^twop);

Xss  = [yd,c,iv,h,g,zY,q,s,xg,k,onep,muxbar,muabar,gbar,zI,zeta,muwbar, mp(ones(1,200))];
Yss  = [gy,gc,gai,gh,gg,gtfp,ga,yd,c,iv,h,v,q,lambda,pai,u,E1PT,E2PT,E3PT,E4PT];

ss = [Yss Xss];

%END_EXTRACT_HERE
set.mcbar = mcbar;
set.ikbar = ikbar;
set.irbar = irbar;