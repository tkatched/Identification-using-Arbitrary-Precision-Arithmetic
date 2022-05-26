function V = vfhlp(theta1,theta2,msel1,msel2,n,bc)
%computes asymptotic variance for LR test between parameterizations (theta1,theta2) of the 
%Leeper(1991) type models

%Inputs:

% theta1 - benchmark parameter vector
% theta2 -alternative parameter vector

%msel1 - select model 1 shock specification
%msel2 - select model 2 shock specification

%n - number of points for Gaussian quadrature 

%bc - business cycle frequency index: =0 for full spectrum, =1 for BC frequencies only 
ny=2; %NB: number of observables has to be the same for both models

%solve both models
[TT1,TC1,TEPS1,TETA1,RC1] =  lpsolv2(theta1,msel1); %solve model using Sims algorithm
[TT2,TC2,TEPS2,TETA2,RC2] =  lpsolv2(theta2,msel2);

%check whether both models have solutions
if isempty(TT1) || isempty(TT2)
    disp('empty solution: nonexistence/numerical issue')
    KL=[];
    return
end

neq1=size(TT1,2); %number of equations - model 1
neq2=size(TT2,2); %number of equations - model 2
%% Compute Gaussian quadrature abscissae and nodes

[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well 
    x=x(1:n/2+1);
    wg=wg(1:n/2+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end
%% Selection and covariance matrices
A1=[eye(ny),zeros(ny,neq1-ny)]; %selection matrix A1
A2=[eye(ny),zeros(ny,neq2-ny)]; %selection matrix A2

QQ1=createcov_lp(theta1,msel1,RC1);
QQ2=createcov_lp(theta2,msel2,RC2);
%% Adjustment for indeterminacy

if RC1(1,1)==1 && RC1(2,1)==1 %determinacy
    
    RR1 = [TEPS1,zeros(neq1,1)];
    
elseif RC1(1,1)==1 && RC1(2,1)==0  %indeterminacy
    TETA1=rref(TETA1').'; %reduced column echelon form for TETA
    RR1 = [TEPS1,TETA1];
    
else
    disp('no equilibrium exists/numerical problems in model 1')
    KL=[];
    return
end

if RC2(1,1)==1 && RC2(2,1)==1 %determinacy
    
    RR2 = [TEPS2,zeros(neq2,1)];
    
elseif RC2(1,1)==1 && RC2(2,1)==0  %indeterminacy
    TETA2=rref(TETA2').'; %reduced column echelon form for TETA
    RR2 = [TEPS2,TETA2];
    
else
    disp('no equilibrium exists/numerical problems in model 2')
    KL=[];
    return
end

%% Compute asymptotic variance

V=0;


if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end

for i=1:nx %compute spectrum
    exe=exp(-1i*x(i));
    
    mat11=(eye(neq1)-TT1*exe)\eye(neq1);
    mat21=mat11';  %note that ' gives the conjugate transpose
    spec1=A1*mat11*RR1*QQ1*RR1'*mat21*A1'/(2*pi); %compute spectrum using representation in the paper
    
    mat12=(eye(neq2)-TT2*exe)\eye(neq2);
    mat22=mat12';  %note that ' gives the conjugate transpose
    spec2=A2*mat12*RR2*QQ2*RR2'*mat22*A2'/(2*pi); %compute spectrum using representation in the paper
    
    
    temp = (eye(ny) - (spec1/spec2))^2; %(I-f_theta*inv(h_phi))
    V=V+2*trace(temp)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    exe=exp(-1i*x(end));
    
    mat11=(eye(neq1)-TT1*exe)\eye(neq1);
    mat21=mat11';  %note that ' gives the conjugate transpose
    spec1=A1*mat11*RR1*QQ1*RR1'*mat21*A1'/(2*pi); %compute spectrum using representation in the paper
    
    mat12=(eye(neq2)-TT2*exe)\eye(neq2);
    mat22=mat12';  %note that ' gives the conjugate transpose
    spec2=A2*mat12*RR2*QQ2*RR2'*mat22*A2'/(2*pi); %compute spectrum using representation in the paper
    
    
    temp = (eye(ny) - (spec1/spec2))^2; %(I-f_theta*inv(h_phi))
    V=V+trace(temp)*wg(end);
end


V=real(V)*(1/4/pi); %final answer: KL. real() necessary since
%small complex residual of order e-15i sometimes remains

end