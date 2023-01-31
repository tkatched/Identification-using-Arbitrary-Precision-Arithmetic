function [gxm,hxm,etam,exitflag,npara] = sgusolvn_mp(theta0,msel)

onep=mp('1');
twop=mp('2');
zp=mp('0');

exitflag=[];
npara=mp(zeros(1,21));
npind=[[23:25];[17:19];[14:16];[29:31];[20:22];[32:34];[26:28]];

nper =20;
load model_object
[param,set] = parameters_mp(theta0);
paramv = struct2array(param);
param_unpack;

%Get indexes to track
load v_idx
vvidx     = [gy_idx, gc_idx, gai_idx, gh_idx];


%*************************************************************************
% FULL INFO MODEL = SGU 2012 News Representation (no inference to solve)
%*************************************************************************
set.fullm = onep; set.sguinfo = onep;
[f, fx, fy, fxp, fyp, eta]=model_prog_mp(paramv,struct2array(set));
[gx,hx,exitflag1]=gx_hx_alt_mp(fy,fx,fyp,fxp);


if msel==1
    gxm=gx; % done: save solution matrices under news shocks
    hxm=hx;
    etam=eta;
    exitflag=exitflag1;
end


if msel==2
    %**************************************************************
% NOISE VERSION: UNLIKE BS AND BLL CODE, AUGMENTING MODEL OUTSIDE
% OF MODEL.M...saves time in symbolic toolbox.
%**************************************************************

%Using results to get noise representation for each of the 7 news shocks
Mi    = cell(1,7);
ETi   = cell(1,7);

zlist  = {'zY','mux','mua','g','zI','zeta','muw'};
sigu_mat = mp(zeros(7,1));
sigv_mat = mp(zeros(7,8));
plm_idx = [];  %will store state indexes that are X(t,t)
alm_idx = [];  %will store state indexes that are X(t)

for jj = 1:7 

    %Get the right noise-rep parameters
    [sigu_mat(jj), sigv_mat(jj,[4,8])]  = invert_rep_mp(eval(['[sig0',zlist{jj}, ' sig4',zlist{jj}, ' sig8', zlist{jj} '].^twop']));
    
   
    npara(npind(jj,:))=[sigu_mat(jj), sigv_mat(jj,[4,8])]; %fill in equivalent noise parameters for later check
    
    %Compute the kalman filter
    [out,F,E,H,GG0,GG1,GG2] = quick_kalman_noise_mp(sigu_mat(jj),sigv_mat(jj,:));
    
    % putting each process in as [XNi(t), XNi(t,t)] = M[(XNi(t-1),XNi(t-1,t-1)] + ETi*enoise(t)
    Mi{jj}  = [F,mp(zeros(length(F)));GG0,GG1];
    ETi{jj} = [E;GG2];
    
    %tracking PLM and ALM indexes
    plm_idx = [plm_idx, length(hx)+ (jj-1)*22 + (12:22)];
    alm_idx = [alm_idx, length(hx) + (jj-1)*22 + (1:11)];
end


%Stocking processes for expectations and shocks
M   = blkdiag(Mi{:});
ET  = blkdiag(ETi{:});

%Augmenting state space of full info model
set.fullm = onep; set.sguinfo = zp;
[f, fx, fy, fxp, fyp, eta]=model_prog_mp(paramv,struct2array(set));

%Adding in noise states to model equations 
fx2  = [fx,mp(zeros(size(fx,1),size(M,2))); mp(zeros(size(M,1),size(fx,2))),M];
fy2  = [fy;mp(zeros(size(M,1),size(fy,2)))];
fxp2 = [fxp,mp(zeros(size(fx,1),size(M,2))); mp(zeros(size(M,1),size(fx,2))),mp(eye(size(M)))];
fyp2 = [fyp;mp(zeros(size(M,1),size(fy,2)))];

%Link shocks to economic states
shcks_idx = [length(hx)+(2:22:length(M))];
fx2(87,shcks_idx(1)) = -onep;
fx2(88,shcks_idx(2)) = -onep;
fx2(89,shcks_idx(3)) = -onep;
fx2(90,shcks_idx(4)) = -onep;
fx2(91,shcks_idx(5)) = -onep;
fx2(92,shcks_idx(6)) = -onep;
fx2(93,shcks_idx(7)) = -onep;

%Solving full-info case
[gx2,hx2]=gx_hx_alt_mp(fy2,fx2,fyp2,fxp2);

%Adding in noise shocks to eta
eta2 = [zp*eta;ET];

%Pull out conjectured LOM
Gxmatn = gx2(e1pt_idx:e4pt_idx,:);
Gxmatn(:,plm_idx) = Gxmatn(:,plm_idx)+Gxmatn(:,alm_idx);
Gxmatn(:,alm_idx) = zp;

%Re-solve the model
set.fullm = zp; set.sguinfo = zp;
[f, fx, fy, fxp, fyp, eta]=model_prog_mp(paramv,struct2array(set));

%Adding in noise states to model equations 
fx3  = [fx,mp(zeros(size(fx,1),size(M,2))); mp(zeros(size(M,1),size(fx,2))),M];
fy3  = [fy;mp(zeros(size(M,1),size(fy,2)))];
fxp3 = [fxp,mp(zeros(size(fx,1),size(M,2))); mp(zeros(size(M,1),size(fx,2))),mp(eye(size(M)))];
fyp3 = [fyp;mp(zeros(size(M,1),size(fy,2)))];

%Link shocks to economic states
fx3(87,shcks_idx(1)) = -onep;
fx3(88,shcks_idx(2)) = -onep;
fx3(89,shcks_idx(3)) = -onep;
fx3(90,shcks_idx(4)) = -onep;
fx3(91,shcks_idx(5)) = -onep;
fx3(92,shcks_idx(6)) = -onep;
fx3(93,shcks_idx(7)) = -onep;

%Adding in noise shocks to eta
eta3 = [zp*eta;ET];

%Adding PLM to expectation equations: eqs 69-72 contain defs of expectation
%vars.
fy3(69:72,e1pt_idx:e4pt_idx) =  mp(eye(4));
fx3(69:72,:)                 = -Gxmatn;

%Solve impose PLM=ALM
[gx3,hx3]=gx_hx_alt_mp(fy3,fx3,fyp3,fxp3);
gxm=gx3;
hxm=hx3;
etam=eta3;
end

end