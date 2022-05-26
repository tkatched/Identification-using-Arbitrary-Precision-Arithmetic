function P = pfhsgulma3ns(theta1,theta2,msel1,msel2,a,T,n,bc)

%Function to compute the empirical distance measure (pfh) for significance
%level a, for sample size T, SGU(2012) model, including specifications with
%MA unanticipated shocks.

%Inputs:
% theta1 - parameter vector of the null model
% msel1 - selects shock specification for the null model

% theta2 - parameter vector of the alternative model
% msel2 - selects shock specification for the alternative model 

% a - significance level. E.g., a=0.05 for 5% level.

% T - sample size

% n - number of points for Gaussian quadrature

%bc - business cycle frequency index: =0 for full spectrum, =1 for BC frequencies only 

z=norminv(1-a,0,1); %Normal critical value
%% Solution
ny=7; %NB: number of observables has to be the same for both models

%solve both models
[TT1,TC1,TEPS1,TETA1,RC1] =  sgusolvl4(theta1,msel1); %solve model using Sims algorithm
[TT2,TC2,TEPS2,TETA2,RC2] =  sgusolvl4(theta2,msel2);

%check whether both models have solutions
if isempty(TT1) || isempty(TT2)
    disp('empty solution: nonexistence/numerical issue')
    KL=[];
    return
end

neq1=size(TT1,2); %number of equations - model 1
neq2=size(TT2,2); %number of equations - model 2

intsel=[[14,15,16];[17,18,19];[20,21,22];[23,24,25];[26,27,28];[29,30,31];[32,33,34]]; %ma parameter indices for various selections
%% Compute Gaussian quadrature abscissae and nodes

[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well 
    x=x(1:floor(n/2)+1);
    wg=wg(1:floor(n/2)+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end

%% Selection and covariance matrices
A1=zeros(ny,neq1); %selection matrix A1

%vector of observables: [y; cons;iinv;hrs;g;tfp;100*a]

eq_y = 1;
eq_c = 2;
eq_i = 3;
eq_h = 4;
eq_g = 5;
eq_z = 6;
eq_agr = 7;

%variable indices for selection
v_i	     = 69;
v_c	     = 70;
v_y	     = 71;
v_h	     = 73;
v_g	     = 49;
v_z	     = 31;
v_mua    = 4;


%fill selection matrix
A1(eq_y,v_y)=1; 

A1(eq_c,v_c)=1; 

A1(eq_i,v_i)=1; 

A1(eq_h,v_h)=1; 

A1(eq_g,v_g)=1; 

A1(eq_z,v_z)=1; 

A1(eq_agr,v_mua)=100;
A2=A1; %selection matrix A2


QQ1=zeros(21,21);
for i=14:34
    QQ1(i-13,i-13)=(theta1(i))^2; %put variances on the diagonal 
end

if msel1>1
    intsel1=intsel(msel1-1,:);
    QQ1(intsel1(1)-13,intsel1(1)-13)=1;
    QQ1(intsel1(2:3)-13,:)=[];
    QQ1(:,intsel1(2:3)-13)=[];
end



QQ2=zeros(21,21);
for i=14:34
    QQ2(i-13,i-13)=(theta2(i))^2; %put variances on the diagonal 
end

if msel2>1
    intsel2=intsel(msel2-1,:);
    QQ2(intsel2(1)-13,intsel2(1)-13)=1;
    QQ2(intsel2(2:3)-13,:)=[];
    QQ2(:,intsel2(2:3)-13)=[];
end

%% Adjustment for indeterminacy

if RC1(1,1)==1 && RC1(2,1)==1 %determinacy
    
    RR1 = TEPS1;
    
    
else
    disp('indeterminacy/no equilibrium exists/numerical problems in model 1')
    KL=[];
    return
end

if RC2(1,1)==1 && RC2(2,1)==1 %determinacy
    
    RR2 = TEPS2;
    
    
else
    disp('indeterminacy/no equilibrium exists/numerical problems in model 2')
    KL=[];
    return
end

%% Compute KL divergences and asymptotic variances

%initialize quantities
kl12=0;
kl21=0;
vfh12=0;
vfh21=0;

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end

for i=1:nx %compute spectrum
    exe=exp(-1i*x(i));
    
    mat11=(eye(neq1)-TT1*exe)\eye(neq1);
    mat21=mat11';  %note that ' gives the conjugate transpose
    spec1=(A1*mat11*RR1*QQ1*RR1'*mat21*A1')/(2*pi); %compute spectrum using representation in the paper
    
    mat12=(eye(neq2)-TT2*exe)\eye(neq2);
    mat22=mat12';  %note that ' gives the conjugate transpose
    spec2=(A2*mat12*RR2*QQ2*RR2'*mat22*A2')/(2*pi); %compute spectrum using representation in the paper
    
    
    temp1=spec2\spec1;
    temp2=spec1\spec2;
    temp3 = (eye(ny) - (temp1))^2; %(I-f_theta*inv(h_phi))
    temp4 = (eye(ny) - (temp2))^2; %(I-f_theta*inv(h_phi))
    
    kl12=kl12+2*(trace(temp1) - log(det(temp1))-ny)*wg(i);
    kl21=kl21+2*(trace(temp2) - log(det(temp2))-ny)*wg(i);
    vfh12=vfh12+2*trace(temp3)*wg(i);
    vfh21=vfh21+2*trace(temp4)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    exe=exp(-1i*x(end));
    
    mat11=(eye(neq1)-TT1*exe)\eye(neq1);
    mat21=mat11';  %note that ' gives the conjugate transpose
    spec1=(A1*mat11*RR1*QQ1*RR1'*mat21*A1')/(2*pi); %compute spectrum using representation in the paper
    
    mat12=(eye(neq2)-TT2*exe)\eye(neq2);
    mat22=mat12';  %note that ' gives the conjugate transpose
    spec2=(A2*mat12*RR2*QQ2*RR2'*mat22*A2')/(2*pi); %compute spectrum using representation in the paper
    
    temp1=spec2\spec1;
    temp2=spec1\spec2;
    temp3 = (eye(ny) - (temp1))^2; %(I-f_theta*inv(h_phi))
    temp4 = (eye(ny) - (temp2))^2; %(I-f_theta*inv(h_phi))
    
    kl12=kl12+2*(trace(temp1) - log(det(temp1))-ny)*wg(end);
    kl21=kl21+2*(trace(temp2) - log(det(temp2))-ny)*wg(end);
    vfh12=vfh12+2*trace(temp3)*wg(end);
    vfh21=vfh21+2*trace(temp4)*wg(end);
end

%final answers.  real() necessary since small complex residual of order e-15i sometimes remains
kl12=real(kl12)/(4*pi); 
kl21=real(kl21)/(4*pi); 
vfh12=real(vfh12)/(4*pi);
vfh21=real(vfh21)/(4*pi);
%% Compute q-alpha
q=(-sqrt(T)*kl12) + sqrt(vfh12)*z; %q-alpha
%%
temp=(q-sqrt(T)*kl21)/(sqrt(vfh21)); %term inside brackets

P=1-normcdf(temp); %final probability

end