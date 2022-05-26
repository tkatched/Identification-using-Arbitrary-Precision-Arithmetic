%INVERT_REP - get noise variance parameters from a vector of news variances
%for a single process of of the SGU form
function [sigu,sigv] = invert_rep(sige2)

N = length(sige2)-1;
sigu = sum(sige2);

%Sum of inverse of remaining variances
S = (1/sige2(1) - 1/(sigu));

%Recursive formula for solving
sigv = zeros(1,N);
for jj = 1:N
    sigv(jj) = (1/sige2(jj+1)*(1/sigu+S)^-1 + 1)/(1/sigu + S);
    S = S-1/sigv(jj);
end
