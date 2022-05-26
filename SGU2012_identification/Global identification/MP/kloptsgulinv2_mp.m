function X = kloptsgulinv2_mp(thetainput,thetab,thetaa,truespec,x,wg,bc,parsel,con,wgt,nrm,indp,lb,ub)

%Function to compute the KL distance between SGU (2012) models (observables
%in levels) featuring all news or all MA unanticipated shocks in arbitrary
%precision.

% Inputs:

%thetainput - candidate solution vector

%w - vector of frequencies between [-pi,pi]

%bc=0 corresponds to the full spectrum case; when bc=1 the function uses only the business cycle frequencies

%parsel - parameter selection index.

% con - constraint evaluation index. If con = 1, the function will evaluate
% the constraint given by constraintlp.m first.

%% Misc. setup
warning('off','all');

thetat=transfmpsgu(thetainput,2,lb,ub); %transform parameters from unconstrained parameterization

n=length(x); %number of abscissae supplied

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=thetat; %replace parameters that are varying in optimization


%% Solve the model

ny=mp('7'); %number of observables

[TT,TC,TEPS,TETA,RC] = sgusolvl3_mp(theta0,2); %solve model using Sims algorithm

if isempty(TT)
    disp('empty solution: nonexistence/numerical issue')
    X=mp('1e12');
    return
end
neq=size(TT,2); %number of equations
dck=rcond(mp(eye(neq))-TT);

if dck<mp('1e-10')
    X=mp('1e12');
    return
end

A=mp(zeros(ny,neq)); %selection matrix A

%vector of observables: [y; cons;iinv;hrs;g;tfp;100*a]
% Some auxiliary quantities defined in mp
ci=mp('-1i'); %imaginary -i
pi2=mp('2*pi'); %2*pi in mp
onep=mp('1');
twop=mp('2');
hhmp=mp('100');

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
A(eq_y,v_y)=onep; 

A(eq_c,v_c)=onep; 

A(eq_i,v_i)=onep; 

A(eq_h,v_h)=onep; 

A(eq_g,v_g)=onep; 

A(eq_z,v_z)=onep; 

A(eq_agr,v_mua)=hhmp;

tempeye=mp(eye(neq)); %auxiliary identity matrix

%QQ omitted since equals identity matrix

if RC(1,1)==1  %determinacy
    
    RR = TEPS;
    
else%no equilibrium exists/numerical problems.
    X=mp('1e12');
    return
end

%% Compute KL divergence from the benchmark spectrum


%preparations to compute spectrum
X=mp('0');
if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end

for i=1:nx %compute spectrum
    exe=exp(ci*x(i));
    
    mat1=(tempeye-TT*exe)\tempeye;
    mat2=mat1';  %note that ' gives the conjugate transpose
    spec=(A*mat1*RR*RR'*mat2*A')/(pi2); %compute spectrum using representation in the paper
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+twop*(trace(temp) - log(det(temp))-ny)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
  exe=exp(ci*x(i));
    
    mat1=(tempeye-TT*exe)\tempeye;
    mat2=mat1';  %note that ' gives the conjugate transpose
    spec=(A*mat1*RR*RR'*mat2*A')/(pi2); %compute spectrum using representation in the paper
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+(trace(temp) - log(det(temp))-ny)*wg(i);
end

X=real(X)*mp('1/(4*pi)'); %final answer: KL. real() necessary since
%small complex residual of order e-15i sometimes remains


if isnan(X) | isinf(X) | isempty(X) | X<0
    
    X=mp('1e12');
    return
end


end