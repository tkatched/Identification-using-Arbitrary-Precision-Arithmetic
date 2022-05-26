function [x,wg] = quadcomp(n,bc)

%Auxiliary function to construct Gauss-Legendre abscissae and weights.

%Inputs: 

%n - an even number of points to approximate the integral (order of the
%quadrature). Should be an mp-object for arbitrary precision computation, 
%double otherwise.

%bc -  0 if the full spectrum is used,  1 if only the business cycle 
%frequencies are considered.
try
mpflag=ismp(n); %detect whether input is double or mp
catch
    disp(' Multiprecision toolbox not on the Matlab path')
    mpflag=0;
end
if mpflag==1 
    if bc==0
        [x,wg] = mp.GaussLegendre(n,mp('-pi'),mp('pi'));
    elseif bc==1
        [x1,wg1] = mp.GaussLegendre(floor(n/2),mp('-pi/3'),mp('-pi/16'));
        [x2,wg2] = mp.GaussLegendre(floor(n/2),mp('pi/16'),mp('pi/3')); %construct abscissae and weights for bc frequencies     
        x=[x1;x2];
        wg=[wg1;wg2];
    end
    
else
    if bc==0
        [x,wg] = legendre_rulem(n,-pi,pi);
    elseif bc==1
        [x1,wg1] = legendre_rulem(floor(n/2),-pi/3,-pi/16);
        [x2,wg2] = legendre_rulem(floor(n/2),pi/16,pi/3); %construct abscissae and weights for bc frequencies     
        x=[x1;x2];
        wg=[wg1;wg2];
    end
end

end
