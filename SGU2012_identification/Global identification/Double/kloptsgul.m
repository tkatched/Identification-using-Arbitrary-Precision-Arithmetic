function X = kloptsgul(thetainput,thetab,thetaa,truespec,x,wg,bc,parsel,con,wgt,nrm,indp)

%Function to compute the KL distance between SGU (2012) models (observables in levels)

% Inputs:

%thetainput - candidate solution vector

%w - vector of frequencies between [-pi,pi]

%bc=0 corresponds to the full spectrum case; when bc=1 the function uses only the business cycle frequencies

%parsel - parameter selection index.

% con - constraint evaluation index. If con = 1, the function will evaluate
% the constraint given by constraintlp.m first.

%% Misc. setup
warning('off','all');

n=length(x); %number of abscissae supplied
%% Evaluate constraint and apply flat penalty if violated
if con(1)==1
    conval=constraintsgu(thetainput,thetab,thetaa,parsel,con(2),wgt,nrm,indp);
    if max(conval)>1e-10 % tolerance level for constraint violation
        X=1e10;
        return
    end
end

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=thetainput; %replace parameters that are varying in optimization


%% Solve the model


ny=7; %number of observables

[TT,TC,TEPS,TETA,RC] = sgusolvl(theta0); %solve model using Sims algorithm

if isempty(TT)
    disp('empty solution: nonexistence/numerical issue')
    X=1e12;
    return
end
neq=size(TT,2); %number of equations
dck=rcond(eye(neq)-TT);

if dck<1e-10
    X=1e12;
    return
end

A=zeros(ny,neq); %selection matrix A

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
A(eq_y,v_y)=1; 

A(eq_c,v_c)=1; 

A(eq_i,v_i)=1; 

A(eq_h,v_h)=1; 

A(eq_g,v_g)=1; 

A(eq_z,v_z)=1; 

A(eq_agr,v_mua)=100;

QQ=zeros(21,21);
for i=14:34
    QQ(i-13,i-13)=(theta0(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
    %QQ(i-13,i-13)=theta0(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
end

if RC(1,1)==1  %determinacy
    
    RR = TEPS;
    
else%no equilibrium exists/numerical problems.
    X=1e12;
    return
end

%% Compute KL divergence from the benchmark spectrum


%preparations to compute spectrum
X=0;
if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end

for i=1:nx %compute spectrum
    exe=exp(-1i*x(i));
    
    mat1=(eye(neq)-TT*exe)\eye(neq);
    mat2=mat1';  %note that ' gives the conjugate transpose
    spec=(A*mat1*RR*QQ*RR'*mat2*A')/(2*pi); %compute spectrum using representation in the paper
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+2*(trace(temp) - log(det(temp))-ny)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    i=i+1;
    exe=exp(-1i*x(i));
    
    mat1=(eye(neq)-TT*exe)\eye(neq);
    mat2=mat1';  %note that ' gives the conjugate transpose
    spec=(A*mat1*RR*QQ*RR'*mat2*A')/(2*pi); %compute spectrum using representation in the paper
    
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+(trace(temp) - log(det(temp))-ny)*wg(i);
end

X=real(X)*(1/4/pi); %final answer: KL. real() necessary since
%small complex residual of order e-15i sometimes remains

X=X*10000;

% if X > -1e-06 %deal with tiny negative values due to approximation when KL close to zero
%     X=abs(X);
% end

if isnan(X) | isinf(X) | isempty(X) | X<0
    
    X=1e12;
    return
end


end