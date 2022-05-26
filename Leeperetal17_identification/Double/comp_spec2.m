
function truespec = comp_spec2(param,x,RegimeM_I)

% -------------------------------------------------------------------------
% Define Observables
% -------------------------------------------------------------------------
% Define variables from DSGE model
variables;

% Set observable positions
obs = zeros(8,1);
obs(1,1) = Ncobs; % 34, consumption, 100 times log difference
obs(2,1) = Niobs; % 35, investment,100 times log difference
obs(3,1) = Nwobs; % 37, wage, 100 times log difference
obs(4,1) = Ngcobs;% 36, gov spending,100 times log difference
obs(5,1) = Nbobs; %38, debt,100 times log difference
obs(6,1) = NLobs; %41, hours worked, level
obs(7,1) = NPiobs; %40, inflation, level
obs(8,1) = NRobs; % 39, interest,level

% consumption, investment, wage, gov spending, debt, hours worked,
% inflation, interest

% -------------------------------------------------------------------------
% Load Data
% -------------------------------------------------------------------------
%[obsdata] = get_data(d1,d2);
%cd ./mat_files; save logpostinfo.mat obs obsdata; cd ..; 


[g0,g1,CC,Psi,Pie,H,C,E,ss,nvar,nexog] = model(param,obs,RegimeM_I);
% function [g0,g1,CC,Psi,Pie,H,C,E,ss,nvar,nexog] = model(param,obs,RegimeM_I)
% For given parameter values, this function puts the model in Chris Sims' 
% Gensys' canonical form and sets up the observation equation.
% The solution will take the form:  x(t) = G * x(t-1) + M * epsilon(t)
%                                   y(t) = H * x(t) + C
%                                   epsilon(t) ~ N(0,E*E')
% INPUT:    
%       param: parameter values
%       obs: vector defining observable variables
%       RegimeM_I: indicator for Regime M versus Regime F estimation
[G,TC,MM,TY,fmat,TZ,TETA,GEV,RC]=gensys_mod2(g0, g1, CC, Psi, Pie, 1); %use Sims-based code, also uses ordqz

M=MM(:,nvar-nexog+1:end); % need to select the part that multiply the structural shocks
Q = M*E*E'*M';% this is the error covariance matrix
%initX = zeros(nvar,1);
%initP = lyapunov_symm(G,Q,1); % lyapunov_symm is a dynare code that solves a Lyapunov equation symmetrically
%lastlike = kalman(obsdata,G,Q,initP,initX,T,nobs,obs,C);

% Y=H*[(I-G(L))^(-1)]M*epsilon(t)
% spectral density H*[(I-G(exp(-iw))^(-1)]*Q*{[(I-G(exp(-iw))^(-1)]*}*H';
sqi=-1i;
        
id1=eye(nvar); %state vector dimension
cc=2*pi;
        
ny=8; % obs dimension
        
truespec=zeros(ny*length(x),ny);
% spectral density H*[(I-G(exp(-iw))^(-1)]*Q*{[(I-G(exp(-iw))^(-1)]*}*H';
for i=1:length(x)
     exe=exp(sqi*x(i));
     mat1=(id1-G*exe)\id1; %inv(1-T1L)
     mat2=mat1';  %note that ' gives the conjugate transpose
     truespec(((i-1)*ny+1):i*ny,:)=H*mat1*Q*mat2*H'/cc;
end


end





