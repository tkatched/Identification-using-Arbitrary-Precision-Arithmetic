%******************************************************************************
%*      Solves the cashless version of the Leeper(1991) model considered, e.g., in Tan and Walker (2015),
%          using multiple precision (MP) computation.
%**************************************************************************/

function [T1,TC,T0,TETA,RC] = lpsolv_mp2(para,msel)

%Input para must be an MP-object.

% The model features MA(2) monetary and fiscal shocks

%Standard deviations are included into PSI - no need to create covariance
%matrix when computing spectrum.

%% Setup
onep=mp('1');
zp=mp('0');

%%
%********************************************************************
%**      assign names to parameters
%********************************************************************
alpha      = para(1); %monetary rule parameter
beta       = para(2); %discount factor
gam        = para(3); %fiscal rule parameter
%sigt       = para(4); %std. dev. of monetary shock
%sigpsi     = para(5); %std. dev. of fiscal shock

if msel==1 %MA(1) shocks
    
    mat1=para(6); %MA for monetary shock no 1
    mapsi1=para(7); %MA for fiscal shock no 1
    rhot=zp;
    rhopsi=zp;
    mat2=zp;
    mapsi2=zp;
elseif msel==2 %MA(2) shocks
        mat1=para(6);
        mat2=para(7); %MA for monetary shock no 2
        mapsi1=para(8);
        mapsi2=para(9);
        rhot=zp;
        rhopsi=zp;
elseif msel==3 %ARMA(1,1) shocks
    rhot=para(6);
    mat1=para(7);
    rhopsi=para(8);
    mapsi1=para(9);
    mat2=zp;
    mapsi2=zp;
end


%********************************************************************
%**      matrices of the canonical system
%********************************************************************

%** equation indices **
eq_inf   = 1;    % Inflation **
eq_b     = 2;  % Real debt **
eq_theta = 3;  % Monetary AR shock **
eq_psi   = 4;  % Fiscal AR shock **
eq_Epi   = 5;  % Expected inflation **
eq_tma1   = 6; %MA state for theta 1
eq_tma2  = 7; %MA state for theta 2
eq_psima1   = 8; % MA state for psi 1
eq_psima2   = 9; % MA state for psi 2

%** variable indices **
v_pi  = 1; %inflation
v_b   = 2; %real debt
v_theta   = 3; %monetary  shock
v_psi   = 4; %fiscal shock
v_Epi = 5; %inflation expectation
v_tma1 = 6; %MA state 1 for monetary
v_tma2 = 7; %MA state 2 for monetary
v_psima1 = 8; %MA state 1 for fiscal
v_psima2 = 9; %MA state 2 for fiscal

%shock indices
e_tma = 1; % Monetary shock
e_psima = 2; %Fiscal shock




%** expectation error indices (neta) **
n_pi = 1;


%** summary **
neq  = 9; % no of equations
neps = 2; %no of shocks
neta = 1; %no of expectation errors

%** initialize matrices **
GAM0 = mp(zeros(neq,neq));
GAM1 = mp(zeros(neq,neq));
C = mp(zeros(neq,1));
PSI = mp(zeros(neq,neps));
PPI = mp(zeros(neq,neta));


%**********************************************************
%**      1. Inflation
%**********************************************************
GAM0(eq_inf,v_pi)   = -alpha;
GAM0(eq_inf,v_Epi) =  onep;
GAM0(eq_inf,v_theta)   =  -onep;


%**********************************************************
%**      2. Real Debt
%**********************************************************

GAM0(eq_b,v_pi)  = onep/beta;
GAM0(eq_b,v_b) =  onep;
GAM0(eq_b,v_psi)   =  onep/beta - onep;

GAM1(eq_b,v_pi)  = alpha/beta;
GAM1(eq_b,v_b) =  onep/beta - gam*(onep/beta-onep);
GAM1(eq_b,v_theta)   =  onep/beta;



%**********************************************************
%**      Shock processes
%**********************************************************
%** theta **
GAM0(eq_theta,v_theta) = onep;
GAM0(eq_theta,v_tma1) = -onep;
GAM0(eq_theta,v_tma2) = -mat1;
GAM1(eq_theta,v_theta) = rhot;
GAM1(eq_theta,v_tma2) = mat2;


%PSI(eq_theta,e_theta) = 1; %if no MA

%** psi **
GAM0(eq_psi,v_psi) = onep;
GAM0(eq_psi,v_psima1) = -onep;
GAM0(eq_psi,v_psima2) = -mapsi1;
GAM1(eq_psi,v_psima2) = mapsi2;
GAM1(eq_psi,v_psi) = rhopsi;
%PSI(eq_psi,e_psi) = 1; %if no MA

%MA part 
GAM0(eq_tma1,v_tma1)   = onep;
GAM0(eq_tma2,v_tma2)   = onep;
GAM1(eq_tma2,v_tma1)   = onep; 
GAM0(eq_psima1,v_psima1)   = onep;  
GAM0(eq_psima2,v_psima2)   = onep;
GAM1(eq_psima2,v_psima1)   = onep;
%NB: standard deviations not included:
PSI(eq_tma1,e_tma)=onep;
PSI(eq_psima1,e_psima)=onep;





%**********************************************************
%**      Expectation error
%**********************************************************
%** E[pi] **
GAM0(eq_Epi,v_pi)  = onep;
GAM1(eq_Epi,v_Epi) = onep;
PPI(eq_Epi,n_pi)  = onep;

%********************************************************************
%**      QZ(generalized Schur) decomposition by GENSYS
%********************************************************************
 [T1,TC,T0,TY,M,TZ,TETA,GEV,RC]=gensys_mp2(GAM0, GAM1, C, PSI, PPI, onep); %use Sims-based code, also uses ordqz

end