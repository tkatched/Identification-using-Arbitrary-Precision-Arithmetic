% Run the local identification analysis of the Leeper (1991) model with 
% MA(1) shocks across parameter values in the AMPF and PMPF regimes to
% obtain Result 1 and Result 2 in Section 5.2


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



%Indeterminacy regime parameters:
alpha2=0.3;
gamma3=1.5;
%sunspot parameters:
m_eta_r=0.3;
m_eta_t=0.3;
sigeta=1; %sunspot shock std. dev.

%% Form parameter vectors:


theta_ampf1=[alpha1,beta,gamma1,sigr,sigt,phir,phit]; % base parameter vector, AMPF regime, generic


thetaind=[alpha2,beta,gamma3,sigr,sigt,phir,phit,m_eta_r,m_eta_t,sigeta]; %base parameter vector, indeterminacy 

%Convert parameter values to MP:
theta_ampf1_mp=mp(theta_ampf1);

thetaind_mp=mp(thetaind);


ntheta=length(theta_ampf1); %no of parameters under determinacy
nthetaind=length(thetaind); %no pf parameters under indeterminacy

%Store parameter names as strings:
parvecd=['alpha';'beta ';'gamma';'sigr ';'sigt ';'phir ';'phit '];
parvecind=['alpha  ';'beta   ';'gamma  ';'sigr   ';'sigt   ';'phir   ';'phit   ';'m_eta_r';'m_eta_t';'sigeta '];

format longE %set display format to long scientific

%% Run local analysis for the AMPF and PMPF regimes
for kk=1:2
%% Choose the numerical differentiation rule

dif=kk; %numerical derivative selection: 1=symmetric quotient, 2=4-point method

%% Double precision: AMPF, generic value

%Full spectrum analysis:

tic;
Gampf1=gmatlpv(theta_ampf1,x,wg,hd,dif,msel); %compute criterion matrix
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

disp 'Nonidentified subsets detected in double precision:'
disp '2-element subsets:'
disp(subnameam1.nsub1_1) %this subset is identified spuriously
disp '3-element subsets:'
disp(subnameam1.nsub2_1)

disp 'Smallest eigenvalues corresponding to nonidentified subsets'
disp(subeigam1)


%% Arbitrary precision: AMPF, generic value

%Full spectrum analysis:

tic;
Gampf1_mp=gmatlp_mpv(theta_ampf1_mp,xmp,wgmp,hmp,dif,msel);
tampf1_mp=toc;
valam1_mp=eig(Gampf1_mp);

disp 'Two smallest eigenvalues in quadruple precision:'
disp(valam1_mp(1:2))

disp 'Eigenvalues of nonidentified sets in quadruple precision'
subeigam1_mp.nev1=min(eig(Gampf1_mp(subindam1.nsub1,subindam1.nsub1)));
subeigam1_mp.nev2=min(eig(Gampf1_mp(subindam1.nsub2,subindam1.nsub2)));
disp(subeigam1_mp.nev1) %we see this 2-element subset is identified, eigenvalue remains the same order E-09
disp(subeigam1_mp.nev2)



% Check using 50-digit precision:

mp.Digits(50); %increase precision
tic;
Gampf1_mp50=gmatlp_mpv(theta_ampf1_mp,xmp,wgmp,hmp,dif,msel);
tampf1_mp50=toc;
valam1_mp50=eig(Gampf1_mp50);

disp 'Two smallest eigenvalues in 50-digit precision:'
disp(valam1_mp50(1:2))

disp 'Eigenvalues of nonidentified sets in quadruple precision'
subeigam1_mp50.nev1=min(eig(Gampf1_mp50(subindam1.nsub1,subindam1.nsub1)));
subeigam1_mp50.nev2=min(eig(Gampf1_mp50(subindam1.nsub2,subindam1.nsub2)));
disp(subeigam1_mp50.nev1)
disp(subeigam1_mp50.nev2)


mp.Digits(34); %return back to quadruple precision


%% Double precision: PMPF

%Full spectrum analysis:

tic;
Gind=gmatlpv(thetaind,x,wg,hd,dif,msel);
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

disp 'Nonidentified subsets detected in double precision:'
disp '5-element subsets:'
disp(subnameind.nsub1_1)
disp '--------------------'
disp(subnameind.nsub1_2)
disp '--------------------'
disp(subnameind.nsub1_3)
disp '--------------------'
disp(subnameind.nsub1_4)
disp '--------------------'
disp '6-element subsets:'
disp(subnameind.nsub2_1)
disp '--------------------'
disp '7-element subsets:'
disp(subnameind.nsub3_1)
disp '--------------------'
disp(subnameind.nsub3_2)
disp '--------------------'
disp(subnameind.nsub3_3)
disp '--------------------'
disp(subnameind.nsub3_4)

disp 'Smallest eigenvalues corresponding to nonidentified subsets'
disp(subeigind)

% It can be seen all eigenvalues are not as small as in previous cases of
% identification failure, and are close to the tolerance level. Here
% failure to detect that the G() matrix has full rank leads to seemingly
% many nonidentified subsets.


%% Arbitrary precision: PMPF

%Full spectrum analysis:

tic;
Gind_mp=gmatlp_mpv(thetaind_mp,xmp,wgmp,hmp,dif,msel);
tind_mp=toc;
valind_mp=eig(Gind_mp);


disp 'Smallest eigenvalue in quadruple precision:'  
disp(valind_mp(1))

%The smallest eigenvalue remains stable at about 1.6E-08 and does not
%shrink to zero, confirming that the result in double precision was not
%correct.


% Check using 50-digit precision:

mp.Digits(50); %increase precision

tic;
Gind_mp50=gmatlp_mpv(thetaind_mp,xmp,wgmp,hmp,dif,msel);
tind_mp50=toc;
valind_mp50=eig(Gind_mp50);

disp 'Smallest eigenvalue in 50-digit precision:'  
disp(valind_mp50(1))


mp.Digits(34); %return back to quadruple precision
%% save results
filename=['Section_5_2_results_',num2str(dif)];
save(filename)
end