%**************************************************************
% SOLVE_LINEAR: Solves the full info and incomplete info versions of the
% Schmitt-Grohe & Uribe (2012) model.
%**************************************************************
addpath('../DSGE_tools');

nper =20;
load model_object
[param,set] = parameters;
paramv = struct2array(param);
param_unpack;

%Get indexes to track
load v_idx
vvidx     = [gy_idx, gc_idx, gai_idx, gh_idx];

%*************************************************************************
% FULL INFO MODEL = SGU 2012 News Representation (no inference to solve)
%*************************************************************************
set.fullm = 1; set.sguinfo = 1;
[f, fx, fy, fxp, fyp, eta]=model_prog(paramv,struct2array(set));
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp);
[sigY1,sigX1] = mom(gx(vvidx,:),hx,eta*eta');

%FEVD for figure
surp_only = 1:3:21;
news_only  =[2:3:21,3:3:21];
[Vy_new2, ~, Vy_tot] = fevd(gx(vvidx,:),hx,eta,[1:19,500]); 
fevd_lines = zeros(4,20);
for jj = 1:20
   fevd_lines(:,jj) = sum(Vy_new2(news_only,:,jj)); 
end
save ../info_figure/sgu_fevd fevd_lines

%**************************************************************
% NOISE VERSION: UNLIKE BS AND BLL CODE, AUGMENTING MODEL OUTSIDE
% OF MODEL.M...saves time in symbolic toolbox.
%**************************************************************

%Using results to get noise representation for each of the 7 news shocks
Mi    = cell(1,7);
ETi   = cell(1,7);

zlist  = {'zY','mux','mua','g','zI','zeta','muw'};
sigu_mat = zeros(7,1);
sigv_mat = zeros(7,8);
plm_idx = [];  %will store state indexes that are X(t,t)
alm_idx = [];  %will store state indexes that are X(t)
for jj = 1:7 

    %Get the right noise-rep parameters
    [sigu_mat(jj), sigv_mat(jj,[4,8])]  = invert_rep(eval(['[sig0',zlist{jj}, ' sig4',zlist{jj}, ' sig8', zlist{jj} '].^2']));
    
    %Compute the kalam filter
    [out,F,E,H,GG0,GG1,GG2] = quick_kalman_noise(sigu_mat(jj),sigv_mat(jj,:));
    
    % putting each process in as [XNi(t), XNi(t,t)] = M[(XNi(t-1),XNi(t-1,t-1)] + ETi*enoise(t)
    Mi{jj}  = [F,zeros(length(F));GG0,GG1];
    ETi{jj} = [E;GG2];
    
    %tracking PLM and ALM indexes
    plm_idx = [plm_idx, length(hx)+ (jj-1)*22 + (12:22)];
    alm_idx = [alm_idx, length(hx) + (jj-1)*22 + (1:11)];
end


%Stocking processes for expectations and shocks
M   = blkdiag(Mi{:});
ET  = blkdiag(ETi{:});

%Augmenting state space of full info model
set.fullm = 1; set.sguinfo = 0;
[f, fx, fy, fxp, fyp, eta]=model_prog(paramv,struct2array(set));

%Adding in noise states to model equations 
fx2  = [fx,zeros(size(fx,1),size(M,2)); zeros(size(M,1),size(fx,2)),M];
fy2  = [fy;zeros(size(M,1),size(fy,2))];
fxp2 = [fxp,zeros(size(fx,1),size(M,2)); zeros(size(M,1),size(fx,2)),eye(size(M))];
fyp2 = [fyp;zeros(size(M,1),size(fy,2))];

%Link shocks to economic states
shcks_idx = [length(hx)+(2:22:length(M))];
fx2(87,shcks_idx(1)) = -1;
fx2(88,shcks_idx(2)) = -1;
fx2(89,shcks_idx(3)) = -1;
fx2(90,shcks_idx(4)) = -1;
fx2(91,shcks_idx(5)) = -1;
fx2(92,shcks_idx(6)) = -1;
fx2(93,shcks_idx(7)) = -1;

%Solving full-info case
[gx2,hx2]=gx_hx_alt(fy2,fx2,fyp2,fxp2);

%Adding in noise shocks to eta
eta2 = [0*eta;ET];

%Pull out conjectured LOM
Gxmatn = gx2(e1pt_idx:e4pt_idx,:);
Gxmatn(:,plm_idx) = Gxmatn(:,plm_idx)+Gxmatn(:,alm_idx);
Gxmatn(:,alm_idx) = 0;

%Re-solve the model
set.fullm = 0; set.sguinfo = 0;
[f, fx, fy, fxp, fyp, eta]=model_prog(paramv,struct2array(set));

%Adding in noise states to model equations 
fx3  = [fx,zeros(size(fx,1),size(M,2)); zeros(size(M,1),size(fx,2)),M];
fy3  = [fy;zeros(size(M,1),size(fy,2))];
fxp3 = [fxp,zeros(size(fx,1),size(M,2)); zeros(size(M,1),size(fx,2)),eye(size(M))];
fyp3 = [fyp;zeros(size(M,1),size(fy,2))];

%Link shocks to economic states
fx3(87,shcks_idx(1)) = -1;
fx3(88,shcks_idx(2)) = -1;
fx3(89,shcks_idx(3)) = -1;
fx3(90,shcks_idx(4)) = -1;
fx3(91,shcks_idx(5)) = -1;
fx3(92,shcks_idx(6)) = -1;
fx3(93,shcks_idx(7)) = -1;

%Adding in noise shocks to eta
eta3 = [0*eta;ET];

%Adding PLM to expectation equations: eqs 69-72 contain defs of expectation
%vars.
fy3(69:72,e1pt_idx:e4pt_idx) =  eye(4);
fx3(69:72,:)                 = -Gxmatn;

%Solve impose PLM=ALM
[gx3,hx3]=gx_hx_alt(fy3,fx3,fyp3,fxp3);

%Save result maib tables
save('SGU_solution', 'gx','hx','eta', 'gx3','hx3','eta3', '*idx');

return

%% Checking that model give same variances
[sigY3,sigX3] = mom(gx3(vvidx,:),hx3,eta3*eta3');
diag(sigY1)./diag(sigY3)

diag(sigX1(11:17,11:17))
diag(sigX3(11:17,11:17))

diag(sigX3(plm_idx,plm_idx))./diag(sigX3(alm_idx,alm_idx))

Vy = variance_decomposition(gx2(vvidx,:),hx2,eta2)
Vy = variance_decomposition(gx3(vvidx,:),hx3,eta3)


