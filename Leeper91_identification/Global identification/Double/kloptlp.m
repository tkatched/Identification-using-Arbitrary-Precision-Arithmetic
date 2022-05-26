function X = kloptlp(thetainput,thetab,thetaa,truespec,x,wg,bc,msela,areg,parsel,invsel,con,wgt,nrm,indp)

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

%% Misc. setup
warning('off','all');

n=length(x); %number of abscissae supplied
%% Evaluate constraint and apply flat penalty if violated
if con(1)==1
    conval=constraintlp(thetainput,thetab,thetaa,msela,parsel,invsel,con(2),wgt,nrm,indp);
    if max(conval)>1e-10 % tolerance level for constraint violation
        X=1e10;
        return
    end
end

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=thetainput; %replace parameters that are varying in optimization


%% Solve the model
%Specify the monetary policy rule using msel

ny=2; %number of observables

[TT,TC,TEPS,TETA,RC] = lpsolv2(theta0,msela); %solve model using Sims algorithm

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

A=[eye(ny),zeros(ny,neq-ny)]; %selection matrix A

QQ=createcov_lp(theta0,msela,RC);

if RC(1,1)==1 && RC(2,1)==1 && areg<3 %determinacy
    
    RR = [TEPS,zeros(neq,1)];
    
elseif RC(1,1)==1 && RC(2,1)==0 && areg<3 %indeterminacy when searching determinacy region
    X=1e10;
    return
elseif RC(1,1)==1 && RC(2,1)==0 && areg==3 %indeterminacy
    TETA=rref(TETA').'; %reduced column echelon form for TETA
    RR = [TEPS,TETA];
elseif RC(1,1)==1 && RC(2,1)==1 && areg==3 %%determinacy when searching indeterminacy region
    X=1e10;
    return
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
        spec=A*mat1*RR*QQ*RR'*mat2*A'/(2*pi); %compute spectrum using representation in the paper
        
        temp=spec\truespec(((i-1)*ny+1):i*ny,:);
        X=X+2*(trace(temp) - log(det(temp))-ny)*wg(i);
        
    end
    
    if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
        i=i+1;
        exe=exp(-1i*x(i));
        
        mat1=(eye(neq)-TT*exe)\eye(neq);
        mat2=mat1';  %note that ' gives the conjugate transpose
        spec=A*mat1*RR*QQ*RR'*mat2*A'/(2*pi); %compute spectrum using representation in the paper
        
        
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