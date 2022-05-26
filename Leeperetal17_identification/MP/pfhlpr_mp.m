function [P1 KL1 P2 KL2] = pfhlpr_mp(theta1,theta2,breg,areg,a,T,n,bc)

%Function to compute the empirical distance measure (pfh) for significance
%level a, for sample size T, Leeper et al. (2017) model, multiple precision

%Inputs:
% theta1 - parameter vector of the null model
% msel1 - selects shock specification for the null model

% theta2 - parameter vector of the alternative model
% msel2 - selects shock specification for the alternative model 

% a - significance level. E.g., a=0.05 for 5% level.

% T - sample size

% n - number of points for Gaussian quadrature

%bc - business cycle frequency index: =0 for full spectrum, =1 for BC frequencies only 

onep=mp('1');
twop=mp('2');
mpz=mp('0');


z=norminv(onep-a,mpz,onep); %Normal critical value
%% Some setup
ny=mp('8'); %NB: number of observables has to be the same for both models

% Define variables from DSGE model
variables;

% Set observable positions
obs = zeros(8,1);
obs(1,1) = Ncobs; % 34, consumption, 100 times log difference
obs(2,1) = Niobs; % 35, investment,100 times log difference
obs(3,1) = Nwobs; % 37, wage, 100 times log difference
obs(4,1) = Ngcobs;% 36, gov spending,100 times log difference
obs(5,1) = Nbobs; %38, debt,100 times log difference
obs(6,1) = NLobs; %41, hours worked, level
obs(7,1) = NPiobs; %40, inflation, level
obs(8,1) = NRobs; % 39, interest,level

%% Baseline model solution
[g0,g1,CC,Psi,Pie,H,C,E,ss,nvar,nexog] = model_mp(theta1,obs,breg);

[G,TC,M,TY,fmat,TZ,TETA,GEV,eu]=gensys_mp2(g0, g1, CC, Psi, Pie, 1); %use Sims-based code, also uses ordqz

% Check whether model 1 has solutions
if isempty(G)
    disp('empty solution: nonexistence/numerical issue')
    KL1=[];
    KL2=[];
    P1=[];
    P2=[];
    return
end

id1=mp(eye(nvar)); %state vector dimension

M=M(:,nvar-nexog+1:end); % need to select the part that multiply the structural shocks

if eu(1,1)==1 && eu(2,1)==1 % determinacy
    
    Q = M*E*E'*M';% this is the error covariance matrix
    
elseif eu(1,1)==1 && eu(2,1)==0 %indeterminacy
    %shocks written as sunshock_t=proj_coef*et+sigma*eta_t; eta_t N(0,1)
    TETA=rref(TETA').'; %reduced column echelon form for TETA
    MM=[M,TETA];
    rotmat1 = [mp(eye(8)),mp(zeros(8,1))];
    rotmat = [rotmat1; theta1(end-8:end)]; %last nine elements are projection coefs followed by standard deviation
    Q1 = [E*E',mp(zeros(8,1));mp(zeros(1,8)),onep]; %structural shock variance and 1
    Q = MM*rotmat*Q1*rotmat'*MM';
else%no equilibrium exists/numerical problems.
    KL1=[];
    KL2=[]
    P1=[];
    P2=[];
    return
end

%% Alternative model solution
[g0_a,g1_a,CC_a,Psi_a,Pie_a,H_a,C_a,E_a,ss_a,nvar_a,nexog_a] = model_mp(theta2,obs,areg);

[G_a,TC_a,M_a,TY_a,fmat_a,TZ_a,TETA_a,GEV_a,eu_a]=gensys_mp2(g0_a, g1_a, CC_a, Psi_a, Pie_a, 1); %use Sims-based code, also uses ordqz

% Check whether model 2 has solutions
if isempty(G_a)
    disp('empty solution: nonexistence/numerical issue')
    KL1=[];
    KL2=[];
    P1=[];
    P2=[];
    return
end

id1_a=mp(eye(nvar_a)); %state vector dimension



M_a=M_a(:,nvar_a-nexog_a+1:end); % need to select the part that multiply the structural shocks

if eu_a(1,1)==1 && eu_a(2,1)==1 % determinacy
    
    Q_a = M_a*E_a*E_a'*M_a';% this is the error covariance matrix
    
elseif eu_a(1,1)==1 && eu_a(2,1)==0 %indeterminacy
    %shocks written as sunshock_t=proj_coef*et+sigma*eta_t; eta_t N(0,1)
    TETA_a=rref(TETA_a').'; %reduced column echelon form for TETA
    MM_a=[M_a,TETA_a];
    rotmat1_a = [mp(eye(8)),mp(zeros(8,1))];
    rotmat_a = [rotmat1_a; theta2(end-8:end)]; %last nine elements are projection coefs followed by standard deviation
    Q1_a = [E_a*E_a',mp(zeros(8,1));mp(zeros(1,8)),onep]; %structural shock variance and 1
    Q_a = MM_a*rotmat_a*Q1_a*rotmat_a'*MM_a';
else%no equilibrium exists/numerical problems.
    KL1=[];
    KL2=[];
    P1=[];
    P2=[];
    return
end
%% Compute Gaussian quadrature abscissae and nodes

[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well 
    x=x(1:floor(n/2)+1);
    wg=wg(1:floor(n/2)+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end



%% Compute KL divergences and asymptotic variances

%initialize quantities
kl12=mpz;
kl21=mpz;
vfh12=mpz;
vfh21=mpz;

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end


tempeye1=mp(eye(nvar));
tempeye2=mp(eye(nvar_a));


for i=1:nx %compute spectrum
    exe=exp(mp('-1i')*x(i));
    
    mat11=(tempeye1-G*exe)\tempeye1;
    mat21=mat11';  %note that ' gives the conjugate transpose
    spec1=(H*mat11*Q*mat21*H')/(mp('2*pi')); %compute spectrum using representation in the paper
    
    mat12=(tempeye2-G_a*exe)\tempeye2;
    mat22=mat12';  %note that ' gives the conjugate transpose
    spec2=(H_a*mat12*Q_a*mat22*H_a')/(mp('2*pi')); %compute spectrum using representation in the paper
    
    
    temp1=spec2\spec1;
    temp2=spec1\spec2;
    temp3 = (mp(eye(ny)) - (temp1))^twop; %(I-f_theta*inv(h_phi))
    temp4 = (mp(eye(ny)) - (temp2))^twop; %(I-f_theta*inv(h_phi))
    
    kl12=kl12+2*(trace(temp1) - log(det(temp1))-ny)*wg(i);
    kl21=kl21+2*(trace(temp2) - log(det(temp2))-ny)*wg(i);
    vfh12=vfh12+2*trace(temp3)*wg(i);
    vfh21=vfh21+2*trace(temp4)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    exe=exp(mp('-1i')*x(end));
    
    mat11=(tempeye1-G*exe)\tempeye1;
    mat21=mat11';  %note that ' gives the conjugate transpose
    spec1=(H*mat11*Q*mat21*H')/(mp('2*pi')); %compute spectrum using representation in the paper
    
    mat12=(tempeye2-G_a*exe)\tempeye2;
    mat22=mat12';  %note that ' gives the conjugate transpose
    spec2=(H_a*mat12*Q_a*mat22*H_a')/(mp('2*pi')); %compute spectrum using representation in the paper
    
    
    temp1=spec2\spec1;
    temp2=spec1\spec2;
    temp3 = (mp(eye(ny)) - (temp1))^twop; %(I-f_theta*inv(h_phi))
    temp4 = (mp(eye(ny)) - (temp2))^twop; %(I-f_theta*inv(h_phi))
    
    kl12=kl12+2*(trace(temp1) - log(det(temp1))-ny)*wg(end);
    kl21=kl21+2*(trace(temp2) - log(det(temp2))-ny)*wg(end);
    vfh12=vfh12+2*trace(temp3)*wg(end);
    vfh21=vfh21+2*trace(temp4)*wg(end);
end

%final answers.  real() necessary since small complex residual of order e-15i sometimes remains
kl12=real(kl12)/mp('4*pi'); 
kl21=real(kl21)/mp('4*pi'); 
vfh12=real(vfh12)/mp('4*pi');
vfh21=real(vfh21)/mp('4*pi');
%% Compute q-alpha
q=(-sqrt(T)*kl12) + sqrt(vfh12)*z; %q-alpha
%%
temped1=(q-sqrt(T)*kl21)/(sqrt(vfh21)); %term inside brackets
temped2=(q-sqrt(T)*kl12)/(sqrt(vfh12)); %term inside brackets
P1=onep-normcdf(temped1); %final probability
KL1=kl12;
P2=onep-normcdf(temped2); %final probability
KL2=kl21;
end