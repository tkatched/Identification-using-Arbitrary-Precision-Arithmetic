% NOISE_PROCESS - Generate H noisy signal kalman represetnation for an iid 
% news process.
%
% Usage:
%
% [M,NU,EC,ED,H,G] = noise_process(sigu,sigv)
%

function [out,F,ET,H,GG0,GG1,GG2,K] = noise_process(sigu,sigv)

%Largest horizon of news.
N = length(sigv);

%Only active signals
sig_idx = find(sigv);
sigv = sigv(sigv~=0);
nsig = length(sigv);

%State process
F = [zeros(N+1+nsig,1),[eye(N);zeros(1+nsig,N)],zeros(N+1+nsig,nsig)]; %state evol, last entries are the shocks
L = [zeros(nsig+1,N),eye(nsig+1)]';
H2 = diag([sigu,sigv]);

ET = L*sqrt(H2);
Q = L*H2*L';

%Observable process
H = zeros(1+nsig,length(F));
H(1,1) = 1; %See current prod shock
H(2:nsig+1,sig_idx+1) = eye(nsig); %See signal about shock in the pipeline
H(2:nsig+1,end-nsig+1:end) = eye(nsig);%Add noise to signals
H = H';

R = zeros(size(H,2));

%Compute the Ricatti solution and SS Kalman Gain via iteration
P = .001*eye(length(F));
crit = 1;
while crit > 1e-14
    P_new = F*P*F' + Q - F*P*H*((H'*P*H + R)\(H'*P*F'));
    crit = max(max(abs(P_new-P)));
    P = P_new;
end

%Depdendence on xi(t-1)
GG0 = P*H*((H'*P*H+R)\(H'*F));

%Depedence on xi(t-1,t-1)
GG1 = F + P*H*((H'*P*H+R)\(-H'*F));

%Depepdence on v(t), state shocks
GG2 = P*H*((H'*P*H+R)\(H'))*L*sqrt(H2);

%Depepdence on w(t), obs shocks
GG3 = P*H/(H'*P*H+R)*sqrt(R);

%K = (F*P*H)/(H'*P*H+R);
K = (P*H)/(H'*P*H+R);
sqQ = sqrt(Q);
sqR = sqrt(R);


out = [GG0(:);GG1(:);GG2(:);GG3(:);F(:);sqQ(:)];


