% Run the local identification analysis of the SGU (2012) model using
% observables in levels

%% Setup model and numerical differentiation parameters


dif=1; %numerical derivative selection: 1=symmetric quotient, 2=4-step method

%Numerical differentiation step sizes:
hd=1e-07; %double precision
hmp=mp('1e-10'); %quadruple precision


%% Gaussian quadrature

n=500; %no of quadrature points
nmp=mp('500');

[x,wg] = quadcomp(n,0);
[xbc,wgbc] = quadcomp(n,1);

[xmp,wgmp] = quadcomp(nmp,0);
[xmpbc,wgmpbc] = quadcomp(nmp,1);

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

%% Load parameters
load theta0_sgu
load theta0mp_sgu

theta0=theta0(1:34); %drop measurement error
theta0mp=theta0mp(1:34);

ntheta=length(theta0);
format longE %set display format to long scientific
%% Double precision,  business cycle frequency only analysis:

tic;
Gbc=gmatsgulv(theta0,xbc,wgbc,hd,dif);
t_dbc=toc;

valgbc=eig(Gbc);
s = svd(Gbc);
tolgbc = max(size(Gbc))*eps(max(s));
rankgbc=rank(Gbc);

disp 'Rank of G(theta0) at BC frequencies:'
disp(rankgbc)

disp 'Smallest eigenvalues:'  
disp(valgbc(1:5))

disp 'Default tolerance level:'
disp(tolgbc)


%% Arbitrary precision, business cycle frequency only analysis:

tic;
Gbc_mp=gmatsgul_mpv(theta0mp,xmpbc,wgmpbc,hmp,dif);
t_mpbc=toc;
valgbc_mp=eig(Gbc_mp);

disp 'Smallest eigenvalue in quadruple precision, BC frequencies:'
disp(valgbc_mp(1:5))

mp.Digits(50); %increase precision
hmp=mp('1e-20'); %reduce differentiation step size

tic;
Gbc_mp50=gmatsgul_mpv(theta0mp,xmpbc,wgmpbc,hmp,dif);
t_mp50bc=toc;
valgbc_mp50=eig(Gbc_mp50);

disp 'Five smallest eigenvalues in 50-digit precision, BC frequencies:'
disp(valgbc_mp50(1:5))

mp.Digits(34); %return back to quadruple precision

save sgu_local_results_levels_bc