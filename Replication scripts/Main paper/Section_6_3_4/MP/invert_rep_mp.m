%INVERT_REP - get noise variance parameters from a vector of news variances
%for a single process of of the SGU form
function [sigu,sigv] = invert_rep_mp(sige2)

onep=mp('1');

N = length(sige2)-1;
sigu = sum(sige2);

%Sum of inverse of remaining variances
S = (onep/sige2(1) - onep/(sigu));

%Recursive formula for solving
sigv = mp(zeros(1,N));
for jj = 1:N
    sigv(jj) = (onep/sige2(jj+1)*(onep/sigu+S)^-onep + onep)/(onep/sigu + S);
    S = S-onep/sigv(jj);
end