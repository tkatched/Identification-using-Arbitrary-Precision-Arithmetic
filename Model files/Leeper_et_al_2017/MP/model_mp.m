function [g0,g1,CC,Psi,Pie,H,C,E,ss,nvar,nexog] = model(param,obs,RegimeM_I)
% function [g0,g1,CC,Psi,Pie,H,C,E,ss,nvar,nexog] = model(param,obs,RegimeM_I)
% For given parameter values, this function puts the model in Chris Sims' 
% Gensys' canonical form and sets up the observation equation.
% The solution will take the form:  x(t) = G * x(t-1) + M * epsilon(t)
%                                   y(t) = H * x(t) + C
%                                   epsilon(t) ~ N(0,E*E')
% INPUT:    
%       param: parameter values
%       obs: vector defining observable variables
%       RegimeM_I: indicator for Regime M versus Regime F estimation

% -------------------------------------------------------------------------
% Initializations 
% ------------------------------------------------------------------------- 
nvar = 58;          % number of variables
nexog = 8;          % number of shocks
nforerrors = 9;     % number of forecast errors

% -------------------------------------------------------------------------
% Call index for endogenous variables
% -------------------------------------------------------------------------
variables


%% ------------------------------------------------------------------------
% Set Parameters
% -------------------------------------------------------------------------

% Calibrated parameters
AD = mp('20');                % average duration of government debt
bet = mp('0.99');           	% discount factor
delt = mp('0.025');           % private capital depreciation rate
alph = mp('0.33');            % share of capital in the prod. function
etaw = mp('0.14');            % elasticity of substitution b/w labor talents
etap = mp('0.14');            % elasticity of substitution b/w intermediate goods
sgc = mp('0.11');             % steady state government consumption to model Y ratio; fiscal variables are US data averages over 1954Q3 to 2014Q2
sb = mp('1.47');              % steady state government debt to model Y ratio
taul = mp('0.186');           % steady state labor tax rate
tauk = mp('0.218');           % steady state capital tax rate
tauc = mp('0.023');           % steady state consumption tax rate

% Estimated parameters 
param_index = 1;
%% steady state growth rate of technology, multiplied by 100
gamm100 = param(param_index);

%% inverse Frisch elasticity
param_index = param_index + 1;
xi =  param(param_index);

%% fraction of non-savers in population
muHH = mp('0');

%% prob. of unions not resetting price
param_index = param_index + 1;
omegaw =  param(param_index);

%% prob. of firms not resetting price
param_index = param_index + 1;
omegap =  param(param_index);

%% K utilization cost parameter
param_index = param_index + 1;
gpsi =  param(param_index); 

%% s" in the investment adjustment cost function
param_index = param_index + 1;
s =  param(param_index);

%% percentage of unions' wage adjustment when they cannot optimize
param_index = param_index + 1;
chiw =  param(param_index);

%% percentage of intermediate firms' price adjustment when they cannot optimize
param_index = param_index + 1;
chip =  param(param_index);

%% Taylor rule inflation coefficient
param_index = param_index + 1;
phipi =  param(param_index);

%% Taylor rule output coefficient
param_index = param_index + 1;
phiy =  param(param_index);

%% response of government consumption to debt
param_index = param_index + 1;
gammgc =  param(param_index);

%% response of capital tax to debt
gammtk = mp('0');

%% response of labor tax to debt
gammtl = mp('0');

%% response of saver transfers to debt
param_index = param_index + 1;
gammz =  param(param_index);

%% serial correlation on technological productiviety growth
param_index = param_index + 1;
rhoa =  param(param_index);

%% serial correlation on preference shock
param_index = param_index + 1;
rhob =  param(param_index);

%% serial correlation on interest rate in Taylor rule
param_index = param_index + 1;
rhor =  param(param_index);

%% serial correlation on investment preference shock
param_index = param_index + 1;
rhoi =  param(param_index);

%% serial correlation on wage markup shock
param_index = param_index + 1;
rhow =  param(param_index);

%% serial correlation on price markup shock
param_index = param_index + 1;
rhop =  param(param_index);

%% serial correlation on gov. consumption shock
param_index = param_index + 1;
rhogc =  param(param_index);

%% serial correlation on capital tax shock
rhotk = mp('0');

%% serial correlation on labor tax shock
rhotl = mp('0');

%% serial correlation on consumption tax shock
rhotc = mp('0');

%% serial correlation on savers transfers shock
if RegimeM_I == 1
    rhoz = mp('0.98'); 
else
    param_index = param_index + 1;
    rhoz = param(param_index);
end

%% technology shock standard deviation
param_index = param_index + 1;
siga =  param(param_index);

%% preference shock standard deviation
param_index = param_index + 1;
sigb =  param(param_index);

%% monetary policy shock standard deviation
param_index = param_index + 1;
sigm =  param(param_index);

%% investment preference shock standard deviation
param_index = param_index + 1;
sigi =  param(param_index);

%% wage markup shock standard deviation
param_index = param_index + 1;
sigw =  param(param_index);

%% price markup shock standard deviation
param_index = param_index + 1;
sigp =  param(param_index);

%% gov. consumption shock standard deviation
param_index = param_index + 1;
siggc =  param(param_index);

%% saver transfer shock standard deviation
param_index = param_index + 1;
sigz =  param(param_index);

%% habit formation
param_index = param_index + 1;
thet =  param(param_index);

%% substitutability (>0) or complementarity (< 0) of private and public consumption              
param_index = param_index + 1;
alphag = param(param_index);    

%% AR coefficients on policy shocks
param_index = param_index + 1;
rhoem = param(param_index);
   
param_index = param_index + 1;
rhoeg = param(param_index); 

if RegimeM_I == 1
    rhoez = mp('0.8'); 
else
    param_index = param_index + 1;
    rhoez = param(param_index); 
end

%% Mean values of observables
param_index = param_index + 1;
Lbar = param(param_index);
param_index = param_index + 1;
Pibar = param(param_index);

                  
%% ------------------------------------------------------------------------
% Standard deviation matrix 
% -------------------------------------------------------------------------

mphd=mp('100');

onep=mp('1');

twop=mp('2');

E = diag([siggc/mphd sigz/mphd siga/mphd sigb/mphd sigm/mphd sigi/mphd sigw/mphd sigp/mphd]);


%% -------------------------------------------------------------------------
% Computation of the steady state
% -------------------------------------------------------------------------
Rho = (onep - (onep/AD))*(onep/bet);
gamm = gamm100/mphd;
expg = exp(gamm);
Rss = expg/bet; % real or nominal interest rate since pi = 1 in steady state
Rbar = (Rss - onep)*mphd;
Pb = onep/(Rss - Rho);
Rk = (exp(gamm)/bet - onep + delt)/(onep-tauk);
psi1 = Rk*(onep-tauk);
mc = onep/(onep + etap);
w = (mc*((onep-alph)^(onep-alph))*(alph^alph)*(Rk^(-alph)))^(onep/(onep-alph));
KL = (w/Rk)*alph/(onep-alph);
OmegL = (KL^alph) - Rk*(KL) - w; 
YL = (KL^alph) - OmegL;
IL = (onep-(onep-delt)*exp(-gamm))*expg*KL;
CL = YL*(onep - sgc) - IL;
rk = Rk;
ZL = ((onep - Rss*exp(-gamm))*sb - sgc)*YL + (tauc)*CL + taul*w + tauk*rk*KL;
ZnL = ZL;
CnL = ((onep-taul)*w + ZnL)/(onep+tauc);
CsL = (CL - muHH*CnL)/(onep-muHH);
CstarL = CsL + alphag*sgc*YL;
l = ((w*(onep-taul)/((onep+tauc)*(onep+etaw)))*(onep/((onep-thet*exp(-gamm))*CstarL)))^(onep/(xi + onep));
cs = CsL*l;
cn = CnL*l;
y = YL*l;
k = KL*l;
Omeg = OmegL*l;
c = CL*l;
inv = IL*l;
z = ZL*l;
b = sb*y;
gc = sgc*y;
ky = k/y;
cy = c/y;
ly = l/y;
TK = tauk*rk*k;
TL = taul*l*w;
TC = tauc*c;
S = tauk*rk*k + taul*l*w + tauc*c - gc - z;
ss = [Rss gc c y inv TK TL TC z bet expg Rho rhogc];

%% -------------------------------------------------------------------------
% System Matrices 
% -------------------------------------------------------------------------
g0 = mp(zeros(nvar,nvar));
g1 = g0;
CC = mp(zeros(nvar,1));
Psi = mp(zeros(nvar,nvar));
Pie = mp(zeros(nvar,nforerrors));


%--------------------------------------------
%(1)	Production Function
%--------------------------------------------
g0(Nl,Ny) = onep;
g0(Nl,Nk) = -((y+Omeg)/y)*alph;
g0(Nl,Nl) = -((y+Omeg)/y)*(onep-alph);


%--------------------------------------------
%(2)	Production Factors
%--------------------------------------------
g0(Nrk,Nrk) = onep;
g0(Nrk,Nw) = -onep;
g0(Nrk,Nk) = onep;
g0(Nrk,Nl) = -onep;


%--------------------------------------------
%(3)	Marginal Cost
%--------------------------------------------
g0(Nmc,Nmc) = onep;
g0(Nmc,Nrk) = -alph;
g0(Nmc,Nw) = alph-onep;

%--------------------------------------------
%(4)	Phillips Equation
%--------------------------------------------
lamprice = ((onep+bet*chip)*omegap)/((onep-bet*omegap)*(onep-omegap));
g0(Npi,Npi) = lamprice;
g0(Npi,Nxpi) = -lamprice*bet/(onep+bet*chip);
g1(Npi,Npi) = lamprice*chip/(onep+bet*chip);
g0(Npi,Nmc) = -onep;
g0(Npi,Nup) = -lamprice;        
                          
%--------------------------------------------
%(5)	Savers' Lagrange Multiplier
%--------------------------------------------
g0(Nlambda,Nlambda) = onep;
g0(Nlambda,Nua) = thet/(expg-thet);
g0(Nlambda,Ncstar) = expg/(expg-thet);
g1(Nlambda,Ncstar) = thet/(expg-thet);
g0(Nlambda,Nub) = -onep;
g0(Nlambda,Ntauc) = (tauc/(onep+tauc));

%--------------------------------------------
%(6)	Long-Run Real Interest Rate
%--------------------------------------------
g0(NrL,NrL) = onep;
g0(NrL,NPb) = onep;
g0(NrL,NxrL) = -bet*Rho/expg;
g0(NrL,NxPb) = -bet*Rho/expg;
g0(NrL,Nxpi) = onep;

%--------------------------------------------
%(7)	Long-Run Inflation Rate
%--------------------------------------------
g0(NpiL,NpiL) = onep;
g0(NpiL,NPb) = onep;
g0(NpiL,NrL) = onep;

%--------------------------------------------
%(8)	Consumption in Utility
%--------------------------------------------
g0(Ncstar,Ncstar) = onep;
g0(Ncstar,Ncs) = -cs/(cs + alphag*gc);
g0(Ncstar,Ngc) = -alphag*gc/(cs + alphag*gc);

%--------------------------------------------
%(9)	Euler Equation
%--------------------------------------------
g0(Ncs,Nlambda) = onep; 
g0(Ncs,NR) = -onep;
g0(Ncs,Nxpi) = onep;
g0(Ncs,Nxlambda) = -onep;
g0(Ncs,Nua) = rhoa;

%--------------------------------------------
%(10)	Capacity Utilization
%--------------------------------------------
g0(Nv,Nrk) = ((onep-gpsi)/gpsi);
g0(Nv,Nv) = -onep;
g0(Nv,Ntauk) = -(((onep-gpsi)/gpsi)*tauk/(onep-tauk));

%--------------------------------------------
%(11)	Capital FOC
%--------------------------------------------
g0(Nq,Nq) = onep;
g0(Nq,NR) = onep;
g0(Nq,Nxpi) = -onep;
g0(Nq,Nxq) = -bet*exp(-gamm)*(onep-delt);
g0(Nq,Nxrk) = -bet*exp(-gamm)*Rk*(onep-tauk);
g0(Nq,Nxtauk) = tauk*exp(-gamm)*bet*Rk;

%--------------------------------------------
%(12)	Investment FOC
%--------------------------------------------
g0(Ni,Nq) = -onep/((onep+bet)*s*(exp(twop*gamm)));
g0(Ni,Ni) = onep;
g0(Ni,Nxi) = -bet/(onep+bet);
g0(Ni,Nua)=(onep-bet*rhoa)/(onep+bet);
g1(Ni,Ni) = onep/(onep+bet);
g0(Ni,Nui) = -onep; % Note: shock normalized by 1/((1+bet)*s*(expg^2))

%--------------------------------------------
%(13)	Effective Capital
%--------------------------------------------
g0(Nk,Nk) = onep;
g0(Nk,Nv) = -onep;
g0(Nk,Nua) = onep;
g1(Nk,Nkbar) = onep;

%--------------------------------------------
%(14)	Law of Motion for Capital
%--------------------------------------------
g0(Nkbar,Nkbar) = onep;
g0(Nkbar,Nui) = -(onep - (onep-delt)*exp(-gamm))*((onep+bet)*s*(exp(twop*gamm))); % Note: extra terms due to normalization
g0(Nkbar,Ni) = -(onep - (onep-delt)*exp(-gamm));
g0(Nkbar,Nua) = (onep-delt)*exp(-gamm);
g1(Nkbar,Nkbar) = (onep-delt)*exp(-gamm);

%--------------------------------------------
%(15)	Wage Equation
%--------------------------------------------
lamwage = (omegaw*(onep+bet)*(onep+xi*(onep+(onep/etaw))))/((onep-omegaw*bet)*(onep-omegaw));
g0(Nw,Nw) = onep + lamwage;
g0(Nw,Nxw) = -lamwage*(bet/(onep+bet));
g0(Nw,Npi) = (lamwage*(onep + bet*chiw)/(onep + bet));
g0(Nw,Nxpi) = -(lamwage*bet/(onep+bet));
g0(Nw,Nl) = -xi;
g0(Nw,Nlambda) = onep;
g0(Nw,Nua) = lamwage*((onep + bet*chiw - rhoa*bet)/(onep+bet));
g0(Nw,Ntaul) = -(taul/(onep-taul));
g0(Nw,Nuw) = -lamwage;
g0(Nw,Nub) = -onep;
g1(Nw,Nw) = (lamwage/(onep+bet));
g1(Nw,Npi) = (lamwage*chiw/(onep + bet));
g1(Nw,Nua) = (lamwage*chiw/(onep+bet));

%--------------------------------------------
%(16)	Monetary Policy Rule
%--------------------------------------------
g0(NR,NR) = onep;
g1(NR,NR) = rhor;
g0(NR,Npi) = -(onep-rhor)*phipi;
g0(NR,Ny) = -(onep-rhor)*phiy;
g0(NR,Num) = -onep;

%--------------------------------------------
%(17)	Aggregate Resource Constraint
%--------------------------------------------
g0(Ny,Nc) = c;
g0(Ny,Ni) = inv;
g0(Ny,Ny) = -y;
g0(Ny,Ngc) = sgc*y;
g0(Ny,Nv) = psi1*k;

%--------------------------------------------
%(18)	Non-Savers Household's Budget
%--------------------------------------------
g0(Ncn,Ncn) = cn*(onep+tauc);
g0(Ncn,Ntauc) = tauc*cn;
g0(Ncn,Nw) = -w*l*(onep-taul);
g0(Ncn,Nl) = -w*l*(onep-taul);
g0(Ncn,Ntaul) = w*l*taul; 
g0(Ncn,Nz) = -z;

%--------------------------------------------
%(19)	Consumption Aggregation
%--------------------------------------------
g0(Nc,Nc) = c;
g0(Nc,Ncs) = -(onep-muHH)*cs;
g0(Nc,Ncn) = -muHH*cn;

%--------------------------------------------
%(20)	Maturity Structure
%--------------------------------------------
g0(NPb,NR) = onep;
g0(NPb,NxPb) = -Rho*Pb/(onep+Rho*Pb);
g0(NPb,NPb) = onep;

%--------------------------------------------
%(21)	Government Budget Constraint
%--------------------------------------------
g0(Nb,Nb) = sb;
g0(Nb,Ngc) = -sgc;
g0(Nb,Nz) = -(z/y); 
g0(Nb,Ntauk) = tauk*rk*ky;
g0(Nb,Nrk) = tauk*rk*ky;
g0(Nb,Nk) = tauk*rk*ky;
g0(Nb,Nua) = sb/bet;
g0(Nb,Ntaul) = taul*w*ly;
g0(Nb,Nw) = taul*w*ly;
g0(Nb,Nl) = taul*w*ly;
g0(Nb,Nc) = (tauc*cy);
g0(Nb,Ntauc) = (tauc*cy);
g1(Nb,Nb) = sb/bet;
g0(Nb,NPb) = -sb*Rho*exp(-gamm);
g1(Nb,NPb) = -sb/bet; 
g0(Nb,Npi) = sb/bet; 

%--------------------------------------------
%(22)    gc Rule
%--------------------------------------------
g0(Ngc,Ngc) = onep;
g1(Ngc,Nsb) = -(onep-rhogc)*gammgc;
g1(Ngc,Ngc) = rhogc;
g0(Ngc,Nugc) = -onep;

%--------------------------------------------
%(23)    Capital Tax Rate Rule
%--------------------------------------------
g0(Ntauk,Ntauk) = onep;
g1(Ntauk,Nsb) = (onep-rhotk)*gammtk;
g1(Ntauk,Ntauk) = rhotk;

%--------------------------------------------
%(24)    Labor Tax Rate Rule
%--------------------------------------------
g0(Ntaul,Ntaul) = onep;
g1(Ntaul,Nsb) = (onep-rhotl)*gammtl;
g1(Ntaul,Ntaul) = rhotl;

%--------------------------------------------
%(25)    Consumption Tax Rate Rule
%--------------------------------------------
g0(Ntauc,Ntauc) = onep;
g1(Ntauc,Ntauc) = rhotc;

%--------------------------------------------
%(26)    Z Rule
%--------------------------------------------
g0(Nz,Nz) = onep;
g1(Nz,Nsb) = -(onep-rhoz)*gammz; 
g1(Nz,Nz) = rhoz;
g0(Nz,Nuz) = -onep;

%--------------------------------------------
%(27)	Fisher Equation
%--------------------------------------------
g0(Nr,Nr) = onep;
g0(Nr,NR) = -onep;
g0(Nr,Nxpi) = onep;

%--------------------------------------------
%(28)	sb Defined
%--------------------------------------------
g0(Nsb,Nsb) = onep;
g0(Nsb,Ny) = onep;
g0(Nsb,Nb) = -onep;

%--------------------------------------------
%(29)	Consumption Tax Revenue
%--------------------------------------------
g0(NTc,NTc) = onep;
g0(NTc,Ntauc) = -onep;
g0(NTc,Nc) = -onep;

%--------------------------------------------
%(30)	Capital Tax Revenue
%--------------------------------------------
g0(NTk,NTk) = onep;
g0(NTk,Ntauk) = -onep;
g0(NTk,Nrk) = -onep;
g0(NTk,Nk) = -onep;

%--------------------------------------------
%(31)	rb defined
%--------------------------------------------
g0(Nrb,Nrb) = onep;
g0(Nrb,NPb) = -Rho*bet/expg;
g1(Nrb,NPb) = -onep;
g0(Nrb,Npi) = onep;

%--------------------------------------------
%(32)	S defined
%--------------------------------------------
g0(NS,NS) = onep;
g0(NS,Ntauk) = -tauk*rk*k/S;
g0(NS,Nrk) = -tauk*rk*k/S;
g0(NS,Nk) = -tauk*rk*k/S;
g0(NS,Ntaul) = -taul*w*l/S;
g0(NS,Nw) = -taul*w*l/S;
g0(NS,Nl) = -taul*w*l/S;
g0(NS,Ntauc) = -tauc*c/S;
g0(NS,Nc) = -tauc*c/S;
g0(NS,Nz) = z/S;
g0(NS,Ngc) = gc/S;

%--------------------------------------------
%(33)	Labor Tax Revenue
%--------------------------------------------
g0(NTl,NTl) = onep;
g0(NTl,Ntaul) = -onep;
g0(NTl,Nw) = -onep;
g0(NTl,Nl) = -onep;

%--------------------------------------------
%(34)	Define Xpi
%--------------------------------------------
g0(Nxpi,Npi)=onep; 
g1(Nxpi,Nxpi)=onep; 
Pie(Nxpi,1)=onep;

%--------------------------------------------
%(35)	Define Xq
%--------------------------------------------
g0(Nxq,Nq)=onep; 
g1(Nxq,Nxq)=onep; 
Pie(Nxq,2)=onep;

%--------------------------------------------
%(36)	Define Xrk
%--------------------------------------------
g0(Nxrk,Nrk)=onep; 
g1(Nxrk,Nxrk)=onep; 
Pie(Nxrk,3)=onep;

%--------------------------------------------
%(37)	Define XI
%--------------------------------------------
g0(Nxi,Ni)=onep; 
g1(Nxi,Nxi)=onep; 
Pie(Nxi,4)=onep;

%--------------------------------------------
%(38)	Define Xtauk
%--------------------------------------------
g0(Nxtauk,Ntauk)=onep; 
g1(Nxtauk,Nxtauk)=onep; 
Pie(Nxtauk,5)=onep;

%--------------------------------------------
%(39)	Define Xw
%--------------------------------------------
g0(Nxw,Nw)=onep; 
g1(Nxw,Nxw)=onep; 
Pie(Nxw,6)=onep;

%--------------------------------------------
%(40)	Define Xlambda
%--------------------------------------------
g0(Nxlambda,Nlambda)=onep; 
g1(Nxlambda,Nxlambda)=onep; 
Pie(Nxlambda,7)=onep;

%--------------------------------------------
%(41)	Define XPb
%--------------------------------------------
g0(NxPb,NPb)=onep; 
g1(NxPb,NxPb)=onep; 
Pie(NxPb,8)=onep;

%--------------------------------------------
%(42)	Define XrL
%--------------------------------------------
g0(NxrL,NrL)=onep; 
g1(NxrL,NxrL)=onep; 
Pie(NxrL,9)=onep;

%--------------------------------------------
%(43)	C Obs
%--------------------------------------------
g0(Ncobs,Ncobs) = onep;
g0(Ncobs,Nc) = -mphd;
g1(Ncobs,Nc) = -mphd;
g0(Ncobs,Nua) = -mphd;

%--------------------------------------------
%(44)	I Obs
%--------------------------------------------
g0(Niobs,Niobs) = onep;
g0(Niobs,Ni) = -mphd;
g1(Niobs,Ni) = -mphd;
g0(Niobs,Nua) = -mphd;

%--------------------------------------------
%(45)	gc Obs
%--------------------------------------------
g0(Ngcobs,Ngcobs) = onep;
g0(Ngcobs,Ngc) = -mphd;
g1(Ngcobs,Ngc) = -mphd;
g0(Ngcobs,Nua) = -mphd;

%--------------------------------------------
%(46)	w Obs
%--------------------------------------------
g0(Nwobs,Nwobs) = onep;
g0(Nwobs,Nw) = -mphd;
g1(Nwobs,Nw) = -mphd;
g0(Nwobs,Nua) = -mphd;

%--------------------------------------------
%(47)	b Obs
%--------------------------------------------
g0(Nbobs,Nbobs) = onep;
g0(Nbobs,Nb) = -mphd;
g1(Nbobs,Nb) = -mphd;
g0(Nbobs,Nua) = -mphd;

%--------------------------------------------
%(48)	R Obs
%--------------------------------------------
g0(NRobs,NRobs) = onep;
g0(NRobs,NR) = -mphd;

%--------------------------------------------
%(49)	Pi Obs
%--------------------------------------------
g0(NPiobs,NPiobs) = onep;
g0(NPiobs,Npi) = -mphd;

%--------------------------------------------
%(50)	L Obs
%--------------------------------------------
g0(NLobs,NLobs) = onep;
g0(NLobs,Nl) = -mphd;

%--------------------------------------------
%(51)	gc Spending Shock
%--------------------------------------------
g0(Nugc,Nugc) = onep;
g1(Nugc,Nugc) = rhoeg; 
Psi(Nugc,Nugc) = onep;

%--------------------------------------------
%(52)	Zs Shock
%--------------------------------------------
g0(Nuz,Nuz)=onep; 
g1(Nuz,Nuz) = rhoez; 
Psi(Nuz,Nuz) = onep;

%--------------------------------------------
%(53)	Growth Rate of Technology shock
%--------------------------------------------
g0(Nua,Nua) = onep; 
g1(Nua,Nua) = rhoa; 
Psi(Nua,Nua) = onep;

%--------------------------------------------
%(54)	Preference Shock
%--------------------------------------------
g0(Nub,Nub) = onep; 
g1(Nub,Nub) = rhob; 
Psi(Nub,Nub) = onep;

%--------------------------------------------
%(55)	Monetary Policy Shock
%--------------------------------------------
g0(Num,Num) = onep; 
g1(Num,Num) = rhoem; 
Psi(Num,Num) = onep;

%--------------------------------------------
%(56)	Investment Shock
%--------------------------------------------
g0(Nui,Nui) = onep; 
g1(Nui,Nui) = rhoi; 
Psi(Nui,Nui) = onep;

%--------------------------------------------
%(57)	wage markup Shock
%--------------------------------------------
g0(Nuw,Nuw) = onep; 
g1(Nuw,Nuw) = rhow; 
Psi(Nuw,Nuw) = onep;

%--------------------------------------------
%(58)	price markup Shock
%--------------------------------------------
g0(Nup,Nup) = onep; 
g1(Nup,Nup) = rhop; 
Psi(Nup,Nup) = onep;


%--------------------------------------------------------------------------
% Observation Equation
%--------------------------------------------------------------------------
if isempty(obs) == 1
    H = [];
    C = [];
else
nobs = size(obs,onep);

% Relating observables to model variables
H = mp(zeros(nobs,nvar));
H(1,obs(1,1)) = onep;
H(2,obs(2,1)) = onep;
H(3,obs(3,1)) = onep;
H(4,obs(4,1)) = onep;
H(5,obs(5,1)) = onep;
H(6,obs(6,1)) = onep;
H(7,obs(7,1)) = onep;
H(8,obs(8,1)) = onep;

% Constants in the observation equation
C = mp(zeros(nobs,1));
C(1) = gamm100;
C(2) = gamm100;
C(3) = gamm100;
C(4) = gamm100;
C(5) = gamm100;
C(6) = Lbar;
C(7) = Pibar;
C(8) = Pibar + Rbar;
end

