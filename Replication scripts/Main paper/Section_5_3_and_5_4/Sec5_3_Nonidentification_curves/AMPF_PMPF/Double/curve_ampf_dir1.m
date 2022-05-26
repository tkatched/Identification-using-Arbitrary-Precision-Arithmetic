filename='Lpr_kl_br1_ar3_sb1_sa1_inv1_par1  2  3  4  5  6  7  8_mp34';

load(filename,'thetaindmp')

%Set base parameters
theta0mp=[thetaindmp,mp(zeros(1,2))];
theta0=[double(theta0mp)];

%% Setup model and numerical differentiation parameters

msel=1; %model selection: 1=MA(1) shocks, 2=MA(2) shocks, 3=ARMA(1,1) shocks

dif=1; %numerical derivative selection: 1=symmetric quotient, 2=4-step method

%Numerical differentiation step sizes:
hd=1e-07; %double precision
hmp=mp('1e-10'); %quadruple precision

ntheta=4;
parvecind=['alpha  ';'beta   ';'gamma  ';'sigr   ';'sigt   ';'phir   ';'phit   ';'m_eta_r';'m_eta_t';'sigeta '];
%% Compute quadrature weights and nodes
n=500;
nmp=mp('500');

[x,wg] = quadcomp(n,0);
[xmp,wgmp] = quadcomp(nmp,0);

%% Compute G-matrix and nonidentified subsets in double precision
Gd=gmatlpv(theta0,x,wg,hd,dif,msel);
[subind subeig subname]=nsubs(Gd,ntheta,[],parvecind); %find nonidentified subsets

%% Validate in quadruple precision
% Gmp = gmatlp_mpv(theta0mp,xmp,wgmp,hmp,dif,msel);
% subeig_mp.nev1=min(eig(Gmp(subind.nsub2,subind.nsub2)));
% subeig_mp.nev2=min(eig(Gmp(subind.nsub3,subind.nsub3)));
% disp(subeig_mp.nev1)
% disp(subeig_mp.nev2)

%% Nonidentification curve direction 1

nit=100000; %predefine necessary iterations

nt=length(theta0);
he=1e-04; %step size in Euler approximation 

%blanks
G=zeros(ntheta,ntheta);
theta2=zeros(nit+1,nt); %matrix collecting parameters on the curve
theta2(1,:)=theta0;
rankvec2=zeros(nit+1,1);
ctheta2=zeros(ntheta,nit+1);
eigcheck2=zeros(ntheta,nit+1);

filenamedir1='curve_ampf_dir1';
alpha=theta0(1);
k=2;
tic;
while alpha<1 
if mod(k,1000)==0
        k %display counter at each 1000th iteration
        save(filenamedir1)
       
    end
    %G=gmatlpi(theta2(k-1,:),x,wg,hd,dif,msel,subind.nsub3);
    G=infmatlpi(theta2(k-1,:),x,wg,hd,dif,msel,subind.nsub3);

rankvec2(k,1)= rank(G); %keep track of rank
[vect, val] = eig(G);
val=abs(val); %take absolute values
b=sort(diag(val));
a=find(diag(val)==b(1)); %choose smallest eigenvalue
c0=real(vect(:,a));
if c0(1,1)<0
    c0=c0*(-1);
end
ctheta2(:,k)=c0;
eigcheck2(:,k)=diag(val);
clear c0 a
theta2(k,:)=theta2(k-1,:)+ [he*ctheta2(1,k)',zeros(1,2),he*ctheta2(2,k),0,he*ctheta2(3,k),0,he*ctheta2(4,k),zeros(1,2)]; %update parm. vector
alpha=theta2(k,1);
k=k+1;
end
toc;
save(filenamedir1)
%%
% check KL and Empirical distances for equally spaced points - direction 1
ind=[2583:2583:25830]';
thetacheck=theta2(ind,:);
klind=zeros(size(thetacheck,1),1);
normind=klind;
edind=zeros(size(thetacheck,1),2);
for i=1:size(thetacheck,1)
    klind(i,1)= kllp(theta2(1,:),thetacheck(i,:),1,1,500,0); %check KL
    edind(i,:)=[pfhlp(theta2(1,:),thetacheck(i,:),1,1,0.05,100,500,0);pfhlp(theta2(1,:),thetacheck(i,:),1,1,0.05,1000,500,0)];
    normind(i,1)=norm(theta2(1,:)-thetacheck(i,:)); 
    
end
save(filenamedir1)