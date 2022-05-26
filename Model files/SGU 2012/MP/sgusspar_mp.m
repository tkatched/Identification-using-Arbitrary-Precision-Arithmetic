%% Compute steady states from the structural parameters 
    
    
%% Calibrated parameter values 

onep=mp('1'); %1 in mp

bet	  = mp('0.99'); %disc factor
sig	  = mp('1'); %intert. elast. sub.
alpk  = mp('0.225'); %capital share
alph  = mp('0.675'); %labor share
del0  = mp('0.025'); %SS depreciation rate when capacity utilization is 1
muxss = mp('1.0032'); %SS growth rate of the permanent productivity shock
muass = mp('0.9957'); %SS gross growth rate of price of investment


g_y	  = mp('0.2'); %SS share of gov consumption in gdp
hss	  = mp('0.2'); %SS hours
muss  = mp('1.15'); %SS gross wage markup

%     Alternatively calibrate muyss which implies muxss:
%     muyss=mp('1.0045'); %SS gross growth rate of GDP/capita (as in SGU)
%     muxss=muyss/(muxss^(alpk/(alpk-1)));


%% Implied and auxiliary parameters
mukss = muass^(onep/(alpk-onep))*muxss;% steady state growth for captial
muyss = muass^(alpk/(alpk-onep))*muxss;	% steady state output growth rate %ok
xgss  = ((onep/muyss)^(onep/(onep - rhoxg))); % SS Xg process %ok

%Golden ratios:
g_y   = g_y/xgss; %detrended SS G/Y ratio %ok
y_k   = (onep/(bet*muass*muyss^(-sig)) - (onep-del0))/(alpk*mukss); %SS Y/K ratio %ok

%other SS ratios
i_k   = onep - (onep-del0)/mukss; %ok
i_y   = i_k/y_k; %ok
c_y   = onep - g_y*xgss - i_y; %ok
y_c   = onep/c_y; %ok

betbar= bet*b*muyss^(-sig); %ok auxiliary beta-bar
badj  = b/muyss; %auxiliary b_adj

ppsi   = ((onep-betbar)*alph/(muss)*y_c)/(hss^thet*(onep/muyss)^((onep-gam)/gam)*(thet*(onep-badj)+gam/(onep-bet*(onep-gam)*muyss^(onep-sig))*alph/(muss)*y_c*(onep-betbar))); %scaling parameter psi
del1  = alpk*y_k*mukss; %determination of delta_1
del2  = del2_del1*del1; %impute delta 2
xi    = onep-bet*(onep-gam)*muyss^(onep-sig); %auxiliary xi
eta   = ppsi*hss^thet*gam*(onep/muyss)^((onep-gam)/gam)/xi; %auxiliary eta