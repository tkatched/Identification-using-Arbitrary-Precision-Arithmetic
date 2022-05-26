function [c, ceq] = constraintlp(xest,thetab,thetaa,msela,parsel,invsel,ns,wgt,nrm,indp)

%Constraint function that ensures 2 out of 3 invertibility conditions for
%MA(2) shock processes (the third is handled directly by bounds).

%Inputs:

% xest - candidate solution point

%parsel - parameter selection index. If parsel = 1, the search is over the
%shock parameters only. If parsel = 2, the search is over all structural
%parameters.

%ns - neighborhood size (denoted "c" in the paper)
%wgt - user-specified vector of weights (dimension=dim(indp))
%nrm - specified norm, takes inputs 1 (for L-1 norm), 2(for L-2 norm) and Inf (for infinity norm)
%indp - vector of indices for parameters to be constrained, in increasing
%order

%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=xest; %replace parameters that are varying in optimization

%% This part enforces invertibility conditions specified by invsel if it is
% supplied as a nonempty value.
if ~isempty(invsel)

    if invsel==1
        c1=[];
    end
    
    if msela==1
        
            ma11=theta0(6);
            ma21=theta0(7);
        
    elseif msela==2
        
            ma11=theta0(6);
            ma12=theta0(7);
            ma21=theta0(8);
            ma22=theta0(9);  
    end
    % Compute invertibility constraints
    if msela==1 %MA(1) case
        if invsel==2
            c1=1.001-abs(ma21);
        elseif invsel==3
            c1=1.001-abs(ma11);
        elseif invsel==4
            c1=[1.001-abs(ma11);1.001-abs(ma21)];
        end
    elseif msela==2
        
        %Evaluate the absolute values of the MA polynomial roots for each shock process
        cond11=abs((-ma11+sqrt(ma11^2-4*ma12))/(2*ma12));
        cond12=abs((-ma11-sqrt(ma11^2-4*ma12))/(2*ma12));
        
        cond21=abs((-ma21+sqrt(ma21^2-4*ma22))/(2*ma22));
        cond22=abs((-ma21-sqrt(ma21^2-4*ma22))/(2*ma22));
        
        if invsel==1
            c1=[1.001-cond11;1.001-cond12;1.001-cond21;1.001-cond22]; %all roots outside unit circle
        elseif invsel==2
            c1=[1.001-cond11;1.001-cond12;min([cond21-0.999;cond22-0.999])]; %monetary ma roots outside the unit circle, at least 1 fiscal inside
        elseif invsel==3
            c1=[min([cond11-0.999;cond12-0.999]);1.001-cond21;1.001-cond22]; % at least 1 monetary ma root inside the unit circle, fiscal outside
        elseif invsel==4
            c1=[min([cond11-0.999;cond12-0.999]);min([cond21-0.999;cond22-0.999])]; %at least 1 root inside the unit circle for each MA
        end
    end
else
    c1=[];
end

%% This part checks the neighborhood constraint if argument ns is different from 0.

if ns>0 %if neighborhood size specified as nonzero, enforce neighborhood constraint
    xcon=theta0(indp); %select the constrained parameters
    
    c2 = ns-norm((xcon-thetab(indp))./wgt,nrm);
else
    c2=[];
end
c=[c1;c2]; %collect the constraint evaluations into 1 vector
ceq = [];
end