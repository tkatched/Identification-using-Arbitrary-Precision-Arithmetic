function [c, ceq] = constraintsgu(xest,thetab,thetaa,parsel,ns,wgt,nrm,indp)

%Constraint function that ensures the neighborhood constraint around theta0

%Inputs:

% xest - candidate solution point

%parsel - parameter selection index.

%ns - neighborhood size (denoted "c" in the paper)
%wgt - user-specified vector of weights (dimension=dim(indp))
%nrm - specified norm, takes inputs 1 (for L-1 norm), 2(for L-2 norm) and Inf (for infinity norm)
%indp - vector of indices for parameters to be constrained, in increasing
%order

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=xest; %replace parameters that are varying in optimization

%% This part checks the neighborhood constraint if argument ns is different from 0.

if ns>0 %if neighborhood size specified as nonzero, enforce neighborhood constraint
    xcon=theta0(indp); %select the constrained parameters
    
    c2 = ns-norm((xcon-thetab(indp))./wgt,nrm);
else
    c2=[];
end
%c=[c1;c2]; %collect the constraint evaluations into 1 vector
c=c2;
ceq = [];
end