% MODEL_SS - Return the steady state of the model computed analytically
%
% usage:
% 
% [ss, parameters] =model_ss(param)


function [ss,param,set] = model_ss(param,set)

%Upack parameters object
param_unpack

who
%BEGIN_EXTRACT_HERE
zY     = 1;
zI     = 1;
zeta   = 1;
u      = 1;
q      = 1;
y_k    = 1/(alphak*muy)*(muy^sig/bet - mua*(1-delta0));
delta1 = alphak*y_k*mui;
delta2 = rdelta*delta1;
k      = (y_k*mui^alphak/h^alphah*L^(alphak+alphah-1))^(1/(alphak-1));
yd     = y_k*k;%0.395274610790875
iv     = (1 - (1-delta0)/mui)*k;% 0.068352497249448
g   = G_Y*yd/xg;%0.080394529501488
c      = yd - iv - g*xg; % 0.247867191383252
s      = c*(1 - b/muy)*muy^((gam-1)/gam);% 0.002427174753172

aux1   = h^theta/(1-bet*(1-gam)*muy^(1-sig));
aux2   = alphah*yd/muw/(theta*h^theta*s);
aux3   = aux2*(1-bet*b*muy^(-sig));
aux4   = gam*s/(c*(1-b/muy));
psii    = aux3/(1+aux1*aux3*aux4);
%term1  = gam*h^theta*s/(c*(1-b/muy))/(1-bet*(1-gam)*muy^(1-sig));
%term2  = theta*h^theta*s*muw/(alphah*yd*(1-bet*b/muy^sig));
%psii    = 1/(term1 + term2); %8.569365092071612e+003
v      = c*(1-b/muy) - psii*h^theta*s; %0.015306293463611
pai    = psii*h^theta/(1 - bet*(1-gam)*muy^(1-sig))/v^sig; %7.849043973959898e+003
lambda = psii*theta*h^theta*s*muw/(alphah*yd*v^sig);%5.384564765310337

% Implied steady-state values
%mua    = muabar;
muy    = muybar; %1.004500000000000
%muw    = muwbar;
mux    = muy*mua^(alphak/(1-alphak));
%muxbar = mux;
mui    = mux*mua^(1/(alphak-1));
h      = hbar;
xg     = muy^(1/(rhoxg-1));
delta  = delta0 + delta1*(u-1) + delta2/2*(u-1)^2;
phi    = 0;
dphi   = 0;

% Observables
gy   = muy;
gc   = muy;
gai  = muy;
gh   = 1;
gg   = muy;
gtfp = mux^(1-alphak);
ga   = mua;
gbar = g;


E1PT = -(- bet*b*(zeta*(v*muy)^(-sig) - gam*pai/muy^sig*(s/(c*muy - b*c))^(1-gam)));
E2PT = (bet*(1-gam)*pai*muy^(-sig)*(s/(c*muy - b*c))^(-gam));
E3PT = (bet*muy^(-sig)*(alphak*yd/k*muy + q*mua*(1-delta)));
E4PT = (bet*muy^(-sig)*q*mua*zI*dphi*(iv*mui/iv)^2);

Xss  = [yd,c,iv,h,g,zY,q,s,xg,k,1,muxbar,muabar,gbar,zI,zeta,muwbar, ones(1,200)];
Yss  = [gy,gc,gai,gh,gg,gtfp,ga,yd,c,iv,h,v,q,lambda,pai,u,E1PT,E2PT,E3PT,E4PT];

ss = [Yss Xss];

%END_EXTRACT_HERE
set.mcbar = mcbar;
set.ikbar = ikbar;
set.irbar = irbar;