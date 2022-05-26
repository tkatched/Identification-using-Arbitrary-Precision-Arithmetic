%% Compute steady states from the structural parameters 
    
    
%% Calibrated parameter values 

bet	  = 0.99; %disc factor
sig	  = 1; %intert. elast. sub.
alpk  = 0.225; %capital share
alph  = 0.675; %labor share
del0  = 0.025; %SS depreciation rate when capacity utilization is 1
muxss = 1.0032; %SS growth rate of the permanent productivity shock
muass = 0.9957; %SS gross growth rate of price of investment


g_y	  = 0.2; %SS share of gov consumption in gdp
hss	  = 0.2; %SS hours
muss  = 1.15; %SS gross wage markup

%     Alternatively calibrate muyss which implies muxss:
%     muyss=1.0045; %SS gross growth rate of GDP/capita (as in SGU)
%     muxss=muyss/(muxss^(alpk/(alpk-1)));


%% Implied and auxiliary parameters
mukss = muass^(1/(alpk-1))*muxss;% steady state growth for captial
muyss = muass^(alpk/(alpk-1))*muxss;	% steady state output growth rate %ok
xgss  = ((1/muyss)^(1/(1 - rhoxg))); % SS Xg process %ok

%Golden ratios:
g_y   = g_y/xgss; %detrended SS G/Y ratio %ok
y_k   = (1/(bet*muass*muyss^(-sig)) - (1-del0))/(alpk*mukss); %SS Y/K ratio %ok

%other SS ratios
i_k   = 1 - (1-del0)/mukss; %ok
i_y   = i_k/y_k; %ok
c_y   = 1 - g_y*xgss - i_y; %ok
y_c   = 1/c_y; %ok

betbar= bet*b*muyss^(-sig); %ok auxiliary beta-bar
badj  = b/muyss; %auxiliary b_adj

ppsi   = ((1-betbar)*alph/(muss)*y_c)/(hss^thet*(1/muyss)^((1-gam)/gam)*(thet*(1-badj)+gam/(1-bet*(1-gam)*muyss^(1-sig))*alph/(muss)*y_c*(1-betbar))); %scaling parameter psi

del1  = alpk*y_k*mukss; %determination of delta_1
del2  = del2_del1*del1; %impute delta 2
xi    = 1-bet*(1-gam)*muyss^(1-sig); %auxiliary xi
eta   = ppsi*hss^thet*gam*(1/muyss)^((1-gam)/gam)/xi; %auxiliary eta