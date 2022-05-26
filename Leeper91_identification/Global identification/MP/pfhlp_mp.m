function P = pfhlp_mp(theta1,theta2,msel1,msel2,a,T,n,bc)

%Function to compute the empirical distance measure (pfh) for significance
%level a, for sample size T, Leeper (1991) model

%Inputs:
% theta1 - parameter vector of the null model
% msel1 - selects shock specification for the null model

% theta2 - parameter vector of the alternative model
% msel2 - selects shock specification for the alternative model 

% a - significance level. E.g., a=0.05 for 5% level.

% T - sample size

% n - number of points for Gaussian quadrature (must be an MP object)

%bc - business cycle frequency index: =0 for full spectrum, =1 for BC frequencies only 

z=norminv(1-a,0,1); %Normal critical value
z=mp(z);
%% Compute q-alpha
q=(-sqrt(T)*kllp_mp(theta1,theta2,msel1,msel2,n,bc)) + sqrt(vfhlp_mp(theta1,theta2,msel1,msel2,n,bc))*z; %q-alpha
%%
temp=(q-sqrt(T)*kllp_mp(theta2,theta1,msel2,msel1,n,bc))/(sqrt(vfhlp_mp(theta2,theta1,msel2,msel1,n,bc))); %term inside brackets

P=mp('1')-normcdf(double(temp)); %final probability

end