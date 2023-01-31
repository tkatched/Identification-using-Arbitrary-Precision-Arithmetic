function X = kloptsguln(thetainput,thetaa,truespec,msel,x,wg,bc,parsel,con,wgt,nrm,indp)

%Function to compute the KL distance between Leeper(1991) models

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

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=thetainput; %replace parameters that are varying in optimization


%% Solve the model


[gx,hx,eta,exitflag,npara] = sgusolvn(theta0,msel);

if exitflag ~=1
    X=1e12;
    return
end



neq=size(hx,2); %number of equations

ny=7; %no of observables


dck=rcond(eye(neq)-hx);

if dck<1e-10
    X=1e12;
    return
end


A=zeros(ny,22); %selection matrix A

%vector of observables: [y; cons;inv;hrs;g;tfp;a_growth]

eq_y = 1;
eq_c = 2;
eq_i = 3;
eq_h = 4;
eq_g = 5;
eq_z = 6;
eq_agr = 7;

%variable indices for selection
v_i	     = 10;
v_c	     = 9;
v_y	     = 8;
v_h	     = 11;
v_g	     = 21;
v_z	     = 22;
v_mua    = 7;
        
        %fill selection matrix
A(eq_y,v_y)=1; 

A(eq_c,v_c)=1; 

A(eq_i,v_i)=1; 

A(eq_h,v_h)=1; 

A(eq_g,v_g)=1; 

A(eq_z,v_z)=1; 

A(eq_agr,v_mua)=100;

gx2=gx;zeros(2,neq); %add selection of levels of g and z
gx2(21,5)=1;
gx2(22,6)=1;

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
    
    mat1=(eye(neq)-hx*exe)\eye(neq);
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    spec=(A*gx2*mat1*eta*eta'*mat2*gx2'*A')/(2*pi); %compute spectrum using representation in the paper
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+2*(trace(temp) - log(det(temp))-ny)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    i=i+1;
    exe=exp(-1i*x(i));
    
    mat1=(eye(neq)-hx*exe)\eye(neq);
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    spec=(A*gx2*mat1*eta*eta'*mat2*gx2'*A')/(2*pi); %compute spectrum using representation in the paper
    
    
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