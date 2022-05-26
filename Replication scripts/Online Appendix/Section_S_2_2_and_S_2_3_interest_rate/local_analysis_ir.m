% Run the local identification analysis of the Leeper (1991) model with 
% MA(1) shocks across parameter values in the three regimes.


%% Gaussian quadrature nodes and weights

n=500; %number of nodes
nmp=mp('500');

[x,wg] = quadcomp(n,0);
[xbc,wgbc] = quadcomp(n,1);

[xmp,wgmp] = quadcomp(nmp,0);
[xmpbc,wgmpbc] = quadcomp(nmp,1);


%% Setup model and numerical parameters

msel=1; %model selection: 1=MA(1) shocks, 2=MA(2) shocks, 3=ARMA(1,1) shocks

%Numerical differentiation step sizes:
hd=1e-07; %double precision
hmp=mp('1e-10'); %quadruple precision

%% Base parameter values in the three regimes:


%keep only one side of the quadrature abscissae and weights due to symmetry
%about zero
x=x(1:n/2);
xbc=xbc(1:n/2);
xmp=xmp(1:n/2);
xmpbc=xmpbc(1:n/2);
wg=wg(1:n/2);
wgbc=wgbc(1:n/2);
wgmp=wgmp(1:n/2);
wgmpbc=wgmpbc(1:n/2);

%Parameters not dependent on regimes:
beta=0.9804; %discount factor
sigr=1; %std. dev. of the monetary shock
sigt=1; %std. dev. of the fiscal shock
phir=0.5; %MA coefficient of the monetary shock
phit=0.5; %MA coefficient of the fiscal shock


%Parameter values for AMPF regime, generic value:
alpha1=1.5;
gamma1=1.2;

%Parameter values for PMAF regime, special value:
alpha2=0.3;
gamma2=0.1;


% Special parameter values for AMPF, generic value:
z1m=1/alpha1;
z2m=1/(1/beta-gamma1*(1/beta-1));
a1=-1/z1m; %MA coefficient for monetary shock satisfying equivalence restriction 
c1=-1/z2m; %MA coefficient for fiscal shock satisfying equivalence restriction

% Special parameter values for PMAF, special value:
z1f=1/alpha2;
z2f=1/(1/beta-gamma2*(1/beta-1));
b1=-1/z1f; %MA coefficient for monetary shock satisfying equivalence restriction 
d1=-1/z2f; %MA coefficient for fiscal shock satisfying equivalence restriction

%Indeterminacy regime parameters:
gamma3=1.5;
%sunspot parameters:
m_eta_r=0.3;
m_eta_t=0.3;
sigeta=1; %sunspot shock std. dev.

%% Form parameter vectors:


theta_ampf1=[alpha1,beta,gamma1,sigr,sigt,phir,phit]; % base parameter vector, AMPF regime, generic

theta_ampf2=[alpha1,beta,gamma1,sigr,sigt,a1,c1]; % special parameter vector, AMPF regime, special

theta_pmaf1=[alpha2,beta,gamma2,sigr,sigt,phir,phit]; %base parameter vector, PMAF regime, generic

theta_pmaf2=[alpha2,beta,gamma2,sigr,sigt,b1,d1]; %base parameter vector, PMAF regime, special

thetaind=[alpha2,beta,gamma3,sigr,sigt,phir,phit,m_eta_r,m_eta_t,sigeta]; %base parameter vector, indeterminacy 

%Convert parameter values to MP:
theta_ampf1_mp=mp(theta_ampf1);
theta_ampf2_mp=mp(theta_ampf2);
theta_pmaf1_mp=mp(theta_pmaf1);
theta_pmaf2_mp=mp(theta_pmaf2);
thetaind_mp=mp(thetaind);


ntheta=length(theta_ampf1); %no of parameters under determinacy
nthetaind=length(thetaind); %no pf parameters under indeterminacy

%Store parameter names as strings:
parvecd=['alpha';'beta ';'gamma';'sigr ';'sigt ';'phir ';'phit '];
parvecind=['alpha  ';'beta   ';'gamma  ';'sigr   ';'sigt   ';'phir   ';'phit   ';'m_eta_r';'m_eta_t';'sigeta '];

format longE %set display format to long scientific

%% Run local analysis for the 3 regimes
for kk=1:2
%% Choose the numerical differentiation rule

dif=kk; %numerical derivative selection: 1=symmetric quotient, 2=4-point method

%% Double precision: AMPF, generic value

%Full spectrum analysis:

tic;
Gampf1=gmatlpv_ir(theta_ampf1,x,wg,hd,dif,msel); %compute criterion matrix
tampf1=toc; %save computational time

valam1=eig(Gampf1); %eigenvalues
s = svd(Gampf1);
tolam1 = max(size(Gampf1))*eps(max(s)); %default tolerance
rankam1=rank(Gampf1); %rank based on default tolerance
[subindam1 subeigam1 subnameam1]=nsubs(Gampf1,ntheta,tolam1,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_ampf1):'
disp(rankam1)

disp 'Smallest 2 eigenvalues:'  
disp(valam1(1:2))

disp 'Default tolerance level:'
disp(tolam1)


%Business cycle frequency only analysis:

tic;
Gampf1bc=gmatlpv_ir(theta_ampf1,xbc,wgbc,hd,dif,msel);
tampf1bc=toc;

valam1bc=eig(Gampf1bc);
s = svd(Gampf1bc);
tolam1bc = max(size(Gampf1bc))*eps(max(s));
rankam1bc=rank(Gampf1bc);
[subindam1bc subeigam1bc subnameam1bc]=nsubs(Gampf1bc,ntheta,tolam1bc,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_ampf1) at BC frequencies:'
disp(rankam1bc)

disp 'Smallest eigenvalue:'  
disp(valam1bc(1))

disp 'Default tolerance level:'
disp(tolam1bc)

%% Arbitrary precision: AMPF, generic value

%Full spectrum analysis:

tic;
Gampf1_mp=gmatlp_mpv_ir(theta_ampf1_mp,xmp,wgmp,hmp,dif,msel);
tampf1_mp=toc;
valam1_mp=eig(Gampf1_mp);

disp 'Two smallest eigenvalues in quadruple precision:'
disp(valam1_mp(1:2))


%Business cycle frequency only analysis:

tic;
Gampf1bc_mp=gmatlp_mpv_ir(theta_ampf1_mp,xmpbc,wgmpbc,hmp,dif,msel);
tampf1bc_mp=toc;
valam1bc_mp=eig(Gampf1bc_mp);

disp 'Smallest eigenvalue in quadruple precision, BC frequencies:'
disp(valam1bc_mp(1))



% Check using 50-digit precision:

mp.Digits(50); %increase precision
tic;
Gampf1_mp50=gmatlp_mpv_ir(theta_ampf1_mp,xmp,wgmp,hmp,dif,msel);
tampf1_mp50=toc;
valam1_mp50=eig(Gampf1_mp50);

disp 'Two smallest eigenvalues in 50-digit precision:'
disp(valam1_mp50(1:2))


tic;
Gampf1bc_mp50=gmatlp_mpv_ir(theta_ampf1_mp,xmpbc,wgmpbc,hmp,dif,msel);
tampf1bc_mp50=toc;
valam1bc_mp50=eig(Gampf1bc_mp50);

disp 'Smallest eigenvalue in 50-digit precision, BC frequencies:'
disp(valam1bc_mp50(1))


mp.Digits(34); %return back to quadruple precision

%% Double precision: AMPF, special value

%Full spectrum analysis:

tic;
Gampf2=gmatlpv_ir(theta_ampf2,x,wg,hd,dif,msel); %compute criterion matrix
tampf2=toc; %save computational time

valam2=eig(Gampf2); %eigenvalues
s = svd(Gampf2);
tolam2 = max(size(Gampf2))*eps(max(s)); %default tolerance
rankam2=rank(Gampf2); %rank based on default tolerance
[subindam2 subeigam2 subnameam2]=nsubs(Gampf2,ntheta,tolam2,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_ampf2):'
disp(rankam2)

disp 'Smallest 3 eigenvalues:'  
disp(valam2(1:3))

disp 'Default tolerance level:'
disp(tolam2)


%Business cycle frequency only analysis:

tic;
Gampf2bc=gmatlpv_ir(theta_ampf2,xbc,wgbc,hd,dif,msel);
tampf2bc=toc;

valam2bc=eig(Gampf2bc);
s = svd(Gampf2bc);
tolam2bc = max(size(Gampf2bc))*eps(max(s));
rankam2bc=rank(Gampf2bc);
[subindam2bc subeigam2bc subnameam2bc]=nsubs(Gampf2bc,ntheta,tolam2bc,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_ampf2) at BC frequencies:'
disp(rankam2bc)

disp 'Smallest eigenvalue:'  
disp(valam2bc(1:3))

disp 'Default tolerance level:'
disp(tolam2bc)

%% Arbitrary precision: AMPF, special value

%Full spectrum analysis:

tic;
Gampf2_mp=gmatlp_mpv_ir(theta_ampf2_mp,xmp,wgmp,hmp,dif,msel);
tampf2_mp=toc;
valam2_mp=eig(Gampf2_mp);

disp 'Three smallest eigenvalues in quadruple precision:'
disp(valam2_mp(1:3))


%Business cycle frequency only analysis:

tic;
Gampf2bc_mp=gmatlp_mpv_ir(theta_ampf2_mp,xmpbc,wgmpbc,hmp,dif,msel);
tampf2bc_mp=toc;
valam2bc_mp=eig(Gampf2bc_mp);

disp 'Smallest eigenvalues in quadruple precision, BC frequencies:'
disp(valam2bc_mp(1:3))

% Check using 50-digit precision:

mp.Digits(50); %increase precision

tic;
Gampf2_mp50=gmatlp_mpv_ir(theta_ampf2_mp,xmp,wgmp,hmp,dif,msel);
tampf2_mp50=toc;
valam2_mp50=eig(Gampf2_mp50);

disp 'Three smallest eigenvalues in quadruple precision:'
disp(valam2_mp50(1:3))

tic;
Gampf2bc_mp50=gmatlp_mpv_ir(theta_ampf2_mp,xmpbc,wgmpbc,hmp,dif,msel);
tampf2bc_mp50=toc;
valam2bc_mp50=eig(Gampf2bc_mp50);

disp 'Smallest 3 eigenvalues in 50-digit precision, BC frequencies:'
disp(valam2bc_mp50(1:3))
 

mp.Digits(34); %return back to quadruple precision

%% Double precision: PMAF, generic value

%Full spectrum analysis:

tic;
Gpmaf1=gmatlpv_ir(theta_pmaf1,x,wg,hd,dif,msel);
tpmaf1=toc;

valpm1=eig(Gpmaf1);
s = svd(Gpmaf1);
tolpm1 = max(size(Gpmaf1))*eps(max(s));
rankpm1=rank(Gpmaf1);
[subindpm1 subeigpm1 subnamepm1]=nsubs(Gpmaf1,ntheta,tolpm1,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_pmaf1):'
disp(rankpm1)

disp 'Smallest 2 eigenvalues:'  
disp(valpm1(1:2))

disp 'Default tolerance level:'
disp(tolpm1)


%Business cycle frequency only analysis:

tic;
Gpmaf1bc=gmatlpv_ir(theta_pmaf1,xbc,wgbc,hd,dif,msel);
tpmaf1bc=toc;

valpm1bc=eig(Gpmaf1bc);
s = svd(Gpmaf1bc);
tolpm1bc = max(size(Gpmaf1bc))*eps(max(s));
rankpm1bc=rank(Gpmaf1bc);
[subindpm1bc subeigpm1bc subnamepm1bc]=nsubs(Gpmaf1bc,ntheta,tolpm1bc,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_pmaf1) at BC frequencies:'
disp(rankpm1bc)

disp 'Two smallest eigenvalues, BC frequencies:'  
disp(valpm1bc(1:2))

disp 'Default tolerance level:'
disp(tolpm1bc)

%% Arbitrary precision: PMAF, generic value

%Full spectrum analysis:

tic;
Gpmaf1_mp=gmatlp_mpv_ir(theta_pmaf1_mp,xmp,wgmp,hmp,dif,msel);
tpmaf1_mp=toc;
valpm1_mp=eig(Gpmaf1_mp);

disp 'Two smallest eigenvalues in quadruple precision:'
disp(valpm1_mp(1:2))


%Business cycle frequency only analysis:

tic;
Gpmaf1bc_mp=gmatlp_mpv_ir(theta_pmaf1_mp,xmpbc,wgmpbc,hmp,dif,msel);
tpmaf1bc_mp=toc;
valpm1bc_mp=eig(Gpmaf1bc_mp);

disp 'Two smallest eigenvalues in quadruple precision, BC frequencies:'
disp(valpm1bc_mp(1:2))


% Check using 50-digit precision:

mp.Digits(50); %increase precision

tic;
Gpmaf1_mp50=gmatlp_mpv_ir(theta_pmaf1_mp,xmp,wgmp,hmp,dif,msel);
tpmaf1_mp50=toc;
valpm1_mp50=eig(Gpmaf1_mp50);

disp 'Two smallest eigenvalues in 50-digit precision:'
disp(valpm1_mp50(1:2))


tic;
Gpmaf1bc_mp50=gmatlp_mpv_ir(theta_pmaf1_mp,xmpbc,wgmpbc,hmp,dif,msel);
tpmaf1bc_mp50=toc;
valpm1bc_mp50=eig(Gpmaf1bc_mp50);

disp 'Two smallest eigenvalues in 50-digit precision, BC frequencies:'
disp(valpm1bc_mp50(1:2))


mp.Digits(34); %return back to quadruple precision

%% Double precision: PMAF, special value

%Full spectrum analysis:

tic;
Gpmaf2=gmatlpv_ir(theta_pmaf2,x,wg,hd,dif,msel);
tpmaf2=toc;

valpm2=eig(Gpmaf2);
s = svd(Gpmaf2);
tolpm2 = max(size(Gpmaf2))*eps(max(s));
rankpm2=rank(Gpmaf2);
[subindpm2 subeigpm2 subnamepm2]=nsubs(Gpmaf2,ntheta,tolpm2,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_pmaf2):'
disp(rankpm2)

disp 'Smallest 3 eigenvalues:'  
disp(valpm2(1:3))

disp 'Default tolerance level:'
disp(tolpm2)


%Business cycle frequency only analysis:

tic;
Gpmaf2bc=gmatlpv_ir(theta_pmaf2,xbc,wgbc,hd,dif,msel);
tpmaf2bc=toc;

valpm2bc=eig(Gpmaf2bc);
s = svd(Gpmaf2bc);
tolpm2bc = max(size(Gpmaf2bc))*eps(max(s));
rankpm2bc=rank(Gpmaf2bc);
[subindpm2bc subeigpm2bc subnamepm2bc]=nsubs(Gpmaf2bc,ntheta,tolpm2bc,parvecd); %find nonidentified subsets

disp 'Rank of G(theta_pmaf2) at BC frequencies:'
disp(rankpm2bc)

disp 'Smallest eigenvalue:'  
disp(valpm2bc(1:3))

disp 'Default tolerance level:'
disp(tolpm2bc)

%% Arbitrary precision: PMAF, special value

%Full spectrum analysis:

tic;
Gpmaf2_mp=gmatlp_mpv_ir(theta_pmaf2_mp,xmp,wgmp,hmp,dif,msel);
tpmaf2_mp=toc;
valpm2_mp=eig(Gpmaf2_mp);

disp 'Three smallest eigenvalues in quadruple precision:'
disp(valpm2_mp(1:3))


%Business cycle frequency only analysis:

tic;
Gpmaf2bc_mp=gmatlp_mpv_ir(theta_pmaf2_mp,xmpbc,wgmpbc,hmp,dif,msel);
tpmaf2bc_mp=toc;
valpm2bc_mp=eig(Gpmaf2bc_mp);

disp 'Three smallest eigenvalues in quadruple precision, BC frequencies:'
disp(valpm2bc_mp(1:3))


% Check using 50-digit precision:

mp.Digits(50); %increase precision

tic;
Gpmaf2_mp50=gmatlp_mpv_ir(theta_pmaf2_mp,xmp,wgmp,hmp,dif,msel);
tpmaf2_mp50=toc;
valpm2_mp50=eig(Gpmaf2_mp50);

disp 'Three smallest eigenvalues in 50-digit precision:'
disp(valpm2_mp50(1:3))

tic;
Gpmaf2bc_mp50=gmatlp_mpv_ir(theta_pmaf2_mp,xmpbc,wgmpbc,hmp,dif,msel);
tpmaf2bc_mp50=toc;
valpm2bc_mp50=eig(Gpmaf2bc_mp50);

disp 'Three smallest eigenvalues in 50-digit precision, BC frequencies:'
disp(valpm2bc_mp50(1:3))


mp.Digits(34); %return back to quadruple precision

%% Double precision: PMPF

%Full spectrum analysis:

tic;
Gind=gmatlpv_ir(thetaind,x,wg,hd,dif,msel);
tind=toc;

valind=eig(Gind);
s = svd(Gind);
tolind = max(size(Gind))*eps(max(s));
rankind=rank(Gind);
[subindind subeigind subnameind]=nsubs(Gind,nthetaind,tolind,parvecind); %find nonidentified subsets

disp 'Rank of G(theta_ind):'
disp(rankind)

disp 'Smallest eigenvalue:'  
disp(valind(1))

disp 'Default tolerance level:'
disp(tolind)


%Business cycle frequency only analysis:

tic;
Gindbc=gmatlpv_ir(thetaind,xbc,wgbc,hd,dif,msel);
tindbc=toc;

valindbc=eig(Gindbc);
s = svd(Gindbc);
tolindbc= max(size(Gindbc))*eps(max(s));
rankindbc=rank(Gindbc);

disp 'Rank of G(theta_ind), BC frequencies:'
disp(rankindbc)

disp 'Smallest eigenvalue, BC frequencies:'  
disp(valindbc(1))

disp 'Default tolerance level, BC frequencies:'
disp(tolindbc)
%% Arbitrary precision: indeterminacy

%Full spectrum analysis:

tic;
Gind_mp=gmatlp_mpv_ir(thetaind_mp,xmp,wgmp,hmp,dif,msel);
tind_mp=toc;
valind_mp=eig(Gind_mp);


disp 'Smallest eigenvalue in quadruple precision:'  
disp(valind_mp(1))

%The smallest eigenvalue remains stable and does not
%shrink to zero, confirming that the result in double precision was not
%correct.

%Business cycle frequency only analysis:

tic;
Gindbc_mp=gmatlp_mpv_ir(thetaind_mp,xmpbc,wgmpbc,hd,dif,msel);
tindbc_mp=toc;
valindbc_mp=eig(Gindbc_mp);

disp 'Smallest eigenvalue in quadruple precision, BC frequencies:'  
disp(valindbc_mp(1))

% Check using 50-digit precision:

mp.Digits(50); %increase precision

tic;
Gind_mp50=gmatlp_mpv_ir(thetaind_mp,xmp,wgmp,hmp,dif,msel);
tind_mp50=toc;
valind_mp50=eig(Gind_mp50);

disp 'Smallest eigenvalue in 50-digit precision:'  
disp(valind_mp50(1))

tic;
Gindbc_mp50=gmatlp_mpv_ir(thetaind_mp,xmpbc,wgmpbc,hd,dif,msel);
tindbc_mp50=toc;
valindbc_mp50=eig(Gindbc_mp50);

disp 'Smallest eigenvalue in 50-digit precision, BC frequencies:'  
disp(valindbc_mp50(1))
mp.Digits(34); %return back to quadruple precision
%% save results
filename=['local_point_analysis_ir_dif',num2str(dif)];
save(filename)
end