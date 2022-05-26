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

%% Double precision, full spectrum

%Full spectrum analysis:

tic;
G=gmatsgulv(theta0,x,wg,hd,dif); %compute criterion matrix
t_d=toc; %save computational time

valg=eig(G); %eigenvalues
s = svd(G);
tolg = max(size(G))*eps(max(s)); %default tolerance
rankg=rank(G); %rank based on default tolerance
%[subind subeig subname]=nsubs(G,ntheta,tolg); %find nonidentified subsets

disp 'Rank of G(theta0):'
disp(rankg)

disp 'Smallest 5 eigenvalues:'  
disp(valg(1:5))

disp 'Default tolerance level:'
disp(tolg)

%% Arbitrary precision

%Full spectrum analysis:

tic;
G_mp=gmatsgul_mpv(theta0mp,xmp,wgmp,hmp,dif);
t_mp=toc;
valg_mp=eig(G_mp);

disp 'Five smallest eigenvalues in quadruple precision:'
disp(valg_mp(1:5))

mp.Digits(50); %increase precision
hmp=mp('1e-20'); %reduce differentiation step size

tic;
G_mp50=gmatsgul_mpv(theta0mp,xmp,wgmp,hmp,dif);
t_mp50=toc;
valg_mp50=eig(G_mp50);

disp 'Five smallest eigenvalues in 50-digit precision:'
disp(valg_mp50(1:5))

mp.Digits(34); %return back to quadruple precision

save sgu_local_results_levels