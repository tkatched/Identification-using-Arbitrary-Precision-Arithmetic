function X = transfmpsgu(theta,dir,lb,ub)
%Transformation of parameters for SGU(2012) model from the unconstrained specification to the
%initial specification and back - (adapted from fminsearchbnd code by John D'Errico)


%Inputs:

%theta - parameter vector, an mp-object

%dir - direction of the transformation: 1 - from initial specification to
%unconstrained; 2 - the other way around


ntheta=length(theta); %no of parameters

%%
onep=mp('1');
twop=mp('2');
pim=mp('pi');

if dir==1 %transforming into unconstrained form (adapted from fminsearchbnd code by John D'Errico)
    
    for i=1:ntheta

        
        
        tmp = twop*(theta(i) - lb(i))/(ub(i)-lb(i)) - onep;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        theta(i) = twop*pim+asin(max(-onep,min(onep,tmp)));
       
    end
    
    
elseif dir==2
    
    for i=1:ntheta


        
      theta(i) = (sin(theta(i))+onep)/twop;
      theta(i) = theta(i)*(ub(i) - lb(i)) + lb(i);
      % just in case of any floating point problems
      theta(i) = max(lb(i),min(ub(i),theta(i)));
       
    end
    
    
end

X=theta;
end