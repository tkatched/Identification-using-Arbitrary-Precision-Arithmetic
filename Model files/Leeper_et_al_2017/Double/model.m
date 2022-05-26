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
AD = 20;                % average duration of government debt
bet = 0.99;           	% discount factor
delt = 0.025;           % private capital depreciation rate
alph = 0.33;            % share of capital in the prod. function
etaw = 0.14;            % elasticity of substitution b/w labor talents
etap = 0.14;            % elasticity of substitution b/w intermediate goods
sgc = 0.11;             % steady state government consumption to model Y ratio; fiscal variables are US data averages over 1954Q3 to 2014Q2
sb = 1.47;              % steady state government debt to model Y ratio
taul = 0.186;           % steady state labor tax rate
tauk = 0.218;           % steady state capital tax rate
tauc = 0.023;           % steady state consumption tax rate

% Estimated parameters 
param_index = 1;
%% steady state growth rate of technology, multiplied by 100
gamm100 = param(param_index);

%% inverse Frisch elasticity
param_index = param_index + 1;
xi =  param(param_index);

%% fraction of non-savers in population
muHH = 0;

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
gammtk = 0;

%% response of labor tax to debt
gammtl = 0;

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
rhotk = 0;

%% serial correlation on labor tax shock
rhotl = 0;

%% serial correlation on consumption tax shock
rhotc = 0;

%% serial correlation on savers transfers shock
if RegimeM_I == 1
    rhoz = 0.98; 
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
    rhoez = 0.8; 
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
E = diag([siggc/100 sigz/100 siga/100 sigb/100 sigm/100 sigi/100 sigw/100 sigp/100]);


%% -------------------------------------------------------------------------
% Computation of the steady state
% -------------------------------------------------------------------------
Rho = (1 - (1/AD))*(1/bet);
gamm = gamm100/100;
expg = exp(gamm);
Rss = expg/bet; % real or nominal interest rate since pi = 1 in steady state
Rbar = (Rss - 1)*100;
Pb = 1/(Rss - Rho);
Rk = (exp(gamm)/bet - 1 + delt)/(1-tauk);
psi1 = Rk*(1-tauk);
mc = 1/(1 + etap);
w = (mc*((1-alph)^(1-alph))*(alph^alph)*(Rk^(-alph)))^(1/(1-alph));
KL = (w/Rk)*alph/(1-alph);
OmegL = (KL^alph) - Rk*(KL) - w; 
YL = (KL^alph) - OmegL;
IL = (1-(1-delt)*exp(-gamm))*expg*KL;
CL = YL*(1 - sgc) - IL;
rk = Rk;
ZL = ((1 - Rss*exp(-gamm))*sb - sgc)*YL + (tauc)*CL + taul*w + tauk*rk*KL;
ZnL = ZL;
CnL = ((1-taul)*w + ZnL)/(1+tauc);
CsL = (CL - muHH*CnL)/(1-muHH);
CstarL = CsL + alphag*sgc*YL;
l = ((w*(1-taul)/((1+tauc)*(1+etaw)))*(1/((1-thet*exp(-gamm))*CstarL)))^(1/(xi + 1));
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
g0 = zeros(nvar,nvar);
g1 = g0;
CC = zeros(nvar,1);
Psi = zeros(nvar,nvar);
Pie = zeros(nvar,nforerrors);


%--------------------------------------------
%(1)	Production Function
%--------------------------------------------
g0(Nl,Ny) = 1;
g0(Nl,Nk) = -((y+Omeg)/y)*alph;
g0(Nl,Nl) = -((y+Omeg)/y)*(1-alph);


%--------------------------------------------
%(2)	Production Factors
%--------------------------------------------
g0(Nrk,Nrk) = 1;
g0(Nrk,Nw) = -1;
g0(Nrk,Nk) = 1;
g0(Nrk,Nl) = -1;


%--------------------------------------------
%(3)	Marginal Cost
%--------------------------------------------
g0(Nmc,Nmc) = 1;
g0(Nmc,Nrk) = -alph;
g0(Nmc,Nw) = alph-1;

%--------------------------------------------
%(4)	Phillips Equation
%--------------------------------------------
lamprice = ((1+bet*chip)*omegap)/((1-bet*omegap)*(1-omegap));
g0(Npi,Npi) = lamprice;
g0(Npi,Nxpi) = -lamprice*bet/(1+bet*chip);
g1(Npi,Npi) = lamprice*chip/(1+bet*chip);
g0(Npi,Nmc) = -1;
g0(Npi,Nup) = -lamprice;        
                          
%--------------------------------------------
%(5)	Savers' Lagrange Multiplier
%--------------------------------------------
g0(Nlambda,Nlambda) = 1;
g0(Nlambda,Nua) = thet/(expg-thet);
g0(Nlambda,Ncstar) = expg/(expg-thet);
g1(Nlambda,Ncstar) = thet/(expg-thet);
g0(Nlambda,Nub) = -1;
g0(Nlambda,Ntauc) = (tauc/(1+tauc));

%--------------------------------------------
%(6)	Long-Run Real Interest Rate
%--------------------------------------------
g0(NrL,NrL) = 1;
g0(NrL,NPb) = 1;
g0(NrL,NxrL) = -bet*Rho/expg;
g0(NrL,NxPb) = -bet*Rho/expg;
g0(NrL,Nxpi) = 1;

%--------------------------------------------
%(7)	Long-Run Inflation Rate
%--------------------------------------------
g0(NpiL,NpiL) = 1;
g0(NpiL,NPb) = 1;
g0(NpiL,NrL) = 1;

%--------------------------------------------
%(8)	Consumption in Utility
%--------------------------------------------
g0(Ncstar,Ncstar) = 1;
g0(Ncstar,Ncs) = -cs/(cs + alphag*gc);
g0(Ncstar,Ngc) = -alphag*gc/(cs + alphag*gc);

%--------------------------------------------
%(9)	Euler Equation
%--------------------------------------------
g0(Ncs,Nlambda) = 1; 
g0(Ncs,NR) = -1;
g0(Ncs,Nxpi) = 1;
g0(Ncs,Nxlambda) = -1;
g0(Ncs,Nua) = rhoa;

%--------------------------------------------
%(10)	Capacity Utilization
%--------------------------------------------
g0(Nv,Nrk) = ((1-gpsi)/gpsi);
g0(Nv,Nv) = -1;
g0(Nv,Ntauk) = -(((1-gpsi)/gpsi)*tauk/(1-tauk));

%--------------------------------------------
%(11)	Capital FOC
%--------------------------------------------
g0(Nq,Nq) = 1;
g0(Nq,NR) = 1;
g0(Nq,Nxpi) = -1;
g0(Nq,Nxq) = -bet*exp(-gamm)*(1-delt);
g0(Nq,Nxrk) = -bet*exp(-gamm)*Rk*(1-tauk);
g0(Nq,Nxtauk) = tauk*exp(-gamm)*bet*Rk;

%--------------------------------------------
%(12)	Investment FOC
%--------------------------------------------
g0(Ni,Nq) = -1/((1+bet)*s*(exp(2*gamm)));
g0(Ni,Ni) = 1;
g0(Ni,Nxi) = -bet/(1+bet);
g0(Ni,Nua)=(1-bet*rhoa)/(1+bet);
g1(Ni,Ni) = 1/(1+bet);
g0(Ni,Nui) = -1; % Note: shock normalized by 1/((1+bet)*s*(expg^2))

%--------------------------------------------
%(13)	Effective Capital
%--------------------------------------------
g0(Nk,Nk) = 1;
g0(Nk,Nv) = -1;
g0(Nk,Nua) = 1;
g1(Nk,Nkbar) = 1;

%--------------------------------------------
%(14)	Law of Motion for Capital
%--------------------------------------------
g0(Nkbar,Nkbar) = 1;
g0(Nkbar,Nui) = -(1 - (1-delt)*exp(-gamm))*((1+bet)*s*(exp(2*gamm))); % Note: extra terms due to normalization
g0(Nkbar,Ni) = -(1 - (1-delt)*exp(-gamm));
g0(Nkbar,Nua) = (1-delt)*exp(-gamm);
g1(Nkbar,Nkbar) = (1-delt)*exp(-gamm);

%--------------------------------------------
%(15)	Wage Equation
%--------------------------------------------
lamwage = (omegaw*(1+bet)*(1+xi*(1+(1/etaw))))/((1-omegaw*bet)*(1-omegaw));
g0(Nw,Nw) = 1 + lamwage;
g0(Nw,Nxw) = -lamwage*(bet/(1+bet));
g0(Nw,Npi) = (lamwage*(1 + bet*chiw)/(1 + bet));
g0(Nw,Nxpi) = -(lamwage*bet/(1+bet));
g0(Nw,Nl) = -xi;
g0(Nw,Nlambda) = 1;
g0(Nw,Nua) = lamwage*((1 + bet*chiw - rhoa*bet)/(1+bet));
g0(Nw,Ntaul) = -(taul/(1-taul));
g0(Nw,Nuw) = -lamwage;
g0(Nw,Nub) = -1;
g1(Nw,Nw) = (lamwage/(1+bet));
g1(Nw,Npi) = (lamwage*chiw/(1 + bet));
g1(Nw,Nua) = (lamwage*chiw/(1+bet));

%--------------------------------------------
%(16)	Monetary Policy Rule
%--------------------------------------------
g0(NR,NR) = 1;
g1(NR,NR) = rhor;
g0(NR,Npi) = -(1-rhor)*phipi;
g0(NR,Ny) = -(1-rhor)*phiy;
g0(NR,Num) = -1;

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
g0(Ncn,Ncn) = cn*(1+tauc);
g0(Ncn,Ntauc) = tauc*cn;
g0(Ncn,Nw) = -w*l*(1-taul);
g0(Ncn,Nl) = -w*l*(1-taul);
g0(Ncn,Ntaul) = w*l*taul; 
g0(Ncn,Nz) = -z;

%--------------------------------------------
%(19)	Consumption Aggregation
%--------------------------------------------
g0(Nc,Nc) = c;
g0(Nc,Ncs) = -(1-muHH)*cs;
g0(Nc,Ncn) = -muHH*cn;

%--------------------------------------------
%(20)	Maturity Structure
%--------------------------------------------
g0(NPb,NR) = 1;
g0(NPb,NxPb) = -Rho*Pb/(1+Rho*Pb);
g0(NPb,NPb) = 1;

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
g0(Ngc,Ngc) = 1;
g1(Ngc,Nsb) = -(1-rhogc)*gammgc;
g1(Ngc,Ngc) = rhogc;
g0(Ngc,Nugc) = -1;

%--------------------------------------------
%(23)    Capital Tax Rate Rule
%--------------------------------------------
g0(Ntauk,Ntauk) = 1;
g1(Ntauk,Nsb) = (1-rhotk)*gammtk;
g1(Ntauk,Ntauk) = rhotk;

%--------------------------------------------
%(24)    Labor Tax Rate Rule
%--------------------------------------------
g0(Ntaul,Ntaul) = 1;
g1(Ntaul,Nsb) = (1-rhotl)*gammtl;
g1(Ntaul,Ntaul) = rhotl;

%--------------------------------------------
%(25)    Consumption Tax Rate Rule
%--------------------------------------------
g0(Ntauc,Ntauc) = 1;
g1(Ntauc,Ntauc) = rhotc;

%--------------------------------------------
%(26)    Z Rule
%--------------------------------------------
g0(Nz,Nz) = 1;
g1(Nz,Nsb) = -(1-rhoz)*gammz; 
g1(Nz,Nz) = rhoz;
g0(Nz,Nuz) = -1;

%--------------------------------------------
%(27)	Fisher Equation
%--------------------------------------------
g0(Nr,Nr) = 1;
g0(Nr,NR) = -1;
g0(Nr,Nxpi) = 1;

%--------------------------------------------
%(28)	sb Defined
%--------------------------------------------
g0(Nsb,Nsb) = 1;
g0(Nsb,Ny) = 1;
g0(Nsb,Nb) = -1;

%--------------------------------------------
%(29)	Consumption Tax Revenue
%--------------------------------------------
g0(NTc,NTc) = 1;
g0(NTc,Ntauc) = -1;
g0(NTc,Nc) = -1;

%--------------------------------------------
%(30)	Capital Tax Revenue
%--------------------------------------------
g0(NTk,NTk) = 1;
g0(NTk,Ntauk) = -1;
g0(NTk,Nrk) = -1;
g0(NTk,Nk) = -1;

%--------------------------------------------
%(31)	rb defined
%--------------------------------------------
g0(Nrb,Nrb) = 1;
g0(Nrb,NPb) = -Rho*bet/expg;
g1(Nrb,NPb) = -1;
g0(Nrb,Npi) = 1;

%--------------------------------------------
%(32)	S defined
%--------------------------------------------
g0(NS,NS) = 1;
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
g0(NTl,NTl) = 1;
g0(NTl,Ntaul) = -1;
g0(NTl,Nw) = -1;
g0(NTl,Nl) = -1;

%--------------------------------------------
%(34)	Define Xpi
%--------------------------------------------
g0(Nxpi,Npi)=1; 
g1(Nxpi,Nxpi)=1; 
Pie(Nxpi,1)=1;

%--------------------------------------------
%(35)	Define Xq
%--------------------------------------------
g0(Nxq,Nq)=1; 
g1(Nxq,Nxq)=1; 
Pie(Nxq,2)=1;

%--------------------------------------------
%(36)	Define Xrk
%--------------------------------------------
g0(Nxrk,Nrk)=1; 
g1(Nxrk,Nxrk)=1; 
Pie(Nxrk,3)=1;

%--------------------------------------------
%(37)	Define XI
%--------------------------------------------
g0(Nxi,Ni)=1; 
g1(Nxi,Nxi)=1; 
Pie(Nxi,4)=1;

%--------------------------------------------
%(38)	Define Xtauk
%--------------------------------------------
g0(Nxtauk,Ntauk)=1; 
g1(Nxtauk,Nxtauk)=1; 
Pie(Nxtauk,5)=1;

%--------------------------------------------
%(39)	Define Xw
%--------------------------------------------
g0(Nxw,Nw)=1; 
g1(Nxw,Nxw)=1; 
Pie(Nxw,6)=1;

%--------------------------------------------
%(40)	Define Xlambda
%--------------------------------------------
g0(Nxlambda,Nlambda)=1; 
g1(Nxlambda,Nxlambda)=1; 
Pie(Nxlambda,7)=1;

%--------------------------------------------
%(41)	Define XPb
%--------------------------------------------
g0(NxPb,NPb)=1; 
g1(NxPb,NxPb)=1; 
Pie(NxPb,8)=1;

%--------------------------------------------
%(42)	Define XrL
%--------------------------------------------
g0(NxrL,NrL)=1; 
g1(NxrL,NxrL)=1; 
Pie(NxrL,9)=1;

%--------------------------------------------
%(43)	C Obs
%--------------------------------------------
g0(Ncobs,Ncobs) = 1;
g0(Ncobs,Nc) = -100;
g1(Ncobs,Nc) = -100;
g0(Ncobs,Nua) = -100;

%--------------------------------------------
%(44)	I Obs
%--------------------------------------------
g0(Niobs,Niobs) = 1;
g0(Niobs,Ni) = -100;
g1(Niobs,Ni) = -100;
g0(Niobs,Nua) = -100;

%--------------------------------------------
%(45)	gc Obs
%--------------------------------------------
g0(Ngcobs,Ngcobs) = 1;
g0(Ngcobs,Ngc) = -100;
g1(Ngcobs,Ngc) = -100;
g0(Ngcobs,Nua) = -100;

%--------------------------------------------
%(46)	w Obs
%--------------------------------------------
g0(Nwobs,Nwobs) = 1;
g0(Nwobs,Nw) = -100;
g1(Nwobs,Nw) = -100;
g0(Nwobs,Nua) = -100;

%--------------------------------------------
%(47)	b Obs
%--------------------------------------------
g0(Nbobs,Nbobs) = 1;
g0(Nbobs,Nb) = -100;
g1(Nbobs,Nb) = -100;
g0(Nbobs,Nua) = -100;

%--------------------------------------------
%(48)	R Obs
%--------------------------------------------
g0(NRobs,NRobs) = 1;
g0(NRobs,NR) = -100;

%--------------------------------------------
%(49)	Pi Obs
%--------------------------------------------
g0(NPiobs,NPiobs) = 1;
g0(NPiobs,Npi) = -100;

%--------------------------------------------
%(50)	L Obs
%--------------------------------------------
g0(NLobs,NLobs) = 1;
g0(NLobs,Nl) = -100;

%--------------------------------------------
%(51)	gc Spending Shock
%--------------------------------------------
g0(Nugc,Nugc) = 1;
g1(Nugc,Nugc) = rhoeg; 
Psi(Nugc,Nugc) = 1;

%--------------------------------------------
%(52)	Zs Shock
%--------------------------------------------
g0(Nuz,Nuz)=1; 
g1(Nuz,Nuz) = rhoez; 
Psi(Nuz,Nuz) = 1;

%--------------------------------------------
%(53)	Growth Rate of Technology shock
%--------------------------------------------
g0(Nua,Nua) = 1; 
g1(Nua,Nua) = rhoa; 
Psi(Nua,Nua) = 1;

%--------------------------------------------
%(54)	Preference Shock
%--------------------------------------------
g0(Nub,Nub) = 1; 
g1(Nub,Nub) = rhob; 
Psi(Nub,Nub) = 1;

%--------------------------------------------
%(55)	Monetary Policy Shock
%--------------------------------------------
g0(Num,Num) = 1; 
g1(Num,Num) = rhoem; 
Psi(Num,Num) = 1;

%--------------------------------------------
%(56)	Investment Shock
%--------------------------------------------
g0(Nui,Nui) = 1; 
g1(Nui,Nui) = rhoi; 
Psi(Nui,Nui) = 1;

%--------------------------------------------
%(57)	wage markup Shock
%--------------------------------------------
g0(Nuw,Nuw) = 1; 
g1(Nuw,Nuw) = rhow; 
Psi(Nuw,Nuw) = 1;

%--------------------------------------------
%(58)	price markup Shock
%--------------------------------------------
g0(Nup,Nup) = 1; 
g1(Nup,Nup) = rhop; 
Psi(Nup,Nup) = 1;


%--------------------------------------------------------------------------
% Observation Equation
%--------------------------------------------------------------------------
if isempty(obs) == 1
    H = [];
    C = [];
else
nobs = size(obs,1);

% Relating observables to model variables
H = zeros(nobs,nvar);
H(1,obs(1,1)) = 1;
H(2,obs(2,1)) = 1;
H(3,obs(3,1)) = 1;
H(4,obs(4,1)) = 1;
H(5,obs(5,1)) = 1;
H(6,obs(6,1)) = 1;
H(7,obs(7,1)) = 1;
H(8,obs(8,1)) = 1;

% Constants in the observation equation
C = zeros(nobs,1);
C(1) = gamm100;
C(2) = gamm100;
C(3) = gamm100;
C(4) = gamm100;
C(5) = gamm100;
C(6) = Lbar;
C(7) = Pibar;
C(8) = Pibar + Rbar;
end

