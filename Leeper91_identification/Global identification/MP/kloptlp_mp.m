function X = kloptlp_mp(thetainput,thetab,thetaa,truespec,x,wg,bc,msela,areg,parsel,invsel,con,wgt,nrm,indp,lb,ub)

%Function to compute the KL distance between Leeper(1991) models

% Inputs:

%thetainput - candidate solution vector

%w - vector of frequencies between [-pi,pi]

%bc=0 corresponds to the full spectrum case; when bc=1 the function uses only the business cycle frequencies

%parsel - parameter selection index. If parsel = 1, the search is over the
%shock parameters only. If parsel = 2, the search is over all structural
%parameters.

% regime - denotes the benchmark regime. If regime = 1, base regime is
% AMPF, so search is over PMAF. If regime = 2, base regime is PMAF, so
% search AMPF.

% con - constraint evaluation index. If con = 1, the function will evaluate
% the constraint given by constraintlp.m first.

%% Setup
warning('off','all');

thetat=transfmplp(thetainput,2,lb,ub); %transform parameters from unconstrained parameterization

n=mp(length(x));
%% Evaluate constraint and apply flat penalty if violated
if con(1)==1
    conval=constraintlp_mp(thetat,thetab,thetaa,msela,parsel,invsel,con(2),wgt,nrm,indp);
    if max(conval)>mp('1e-10') % tolerance level for constraint violation
        X=mp('1e10');
        return
    end
end

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=thetat; %replace parameters that are varying in optimization

%% Solve the model
%Specify the monetary policy rule using msel

ny=2; %number of observables

[TT,TC,TEPS,TETA,RC] = lpsolv_mp2(theta0,msela); %solve model using Sims algorithm

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

A=[mp(eye(ny)),mp(zeros(ny,neq-ny))]; %selection matrix A

QQ=createcov_lpmp(theta0,msela,RC);

if RC(1,1)==1 && RC(2,1)==1 && areg<3 %determinacy
    
    RR = [TEPS,mp(zeros(neq,1))];
    
elseif RC(1,1)==1 && RC(2,1)==0 && areg<3 %indeterminacy when searching determinacy region
    X=mp('1e10');
    return
elseif RC(1,1)==1 && RC(2,1)==0 && areg==3 %indeterminacy
    TETA=rref(TETA').'; %reduced column echelon form for TETA
    RR = [TEPS,TETA];
elseif RC(1,1)==1 && RC(2,1)==1 && areg==3 %%determinacy when searching indeterminacy region
    X=mp('1e10');
    return
else%no equilibrium exists/numerical problems.
    X=mp('1e12');
    return
end

%% Compute KL divergence from the benchmark spectrum


X=mp('0');
sqi=mp('-1i');
id=mp(eye(neq));
cc=mp('2*pi');

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end

for i=1:nx
    exe=exp(sqi*x(i));
    mat1=(id-TT*exe)\id; %inv(1-T1L)
    mat2=mat1';  %note that ' gives the conjugate transpose
    spec=A*mat1*RR*QQ*RR'*mat2*A'/cc; %spectral density matrix
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+mp('2')*(trace(temp) - log(det(temp))-mp('2'))*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    i=i+1;
    
    exe=exp(sqi*x(i));
    mat1=(id-TT*exe)\id; %inv(1-T1L)
    mat2=mat1';  %note that ' gives the conjugate transpose
    spec=A*mat1*RR*QQ*RR'*mat2*A'/cc; %spectral density matrix
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+(trace(temp) - log(det(temp))-mp('2'))*wg(i);
end

X=real(X)*mp('1/4/pi'); %final answer: KL. real() necessary since
%small complex residual of order e-15i sometimes remains

if isnan(X) | isinf(X) | isempty(X) | X<0
    
    X=mp('1e12');
    return
end


end