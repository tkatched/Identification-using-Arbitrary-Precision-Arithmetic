% MODEL - Generate linearized version of the Schmitt-Grohe and Uribe (2012) model

function [mod,param,set] = model(param,set)

%Name of text files
mod.fname   = 'model_prog.m';
mod.ss_call = 'model_ss.m';

%Declare parameters symbols: parameters are values to be estimated, symbols are values that are "fixed"
PARAM = struct2sym(param);
SET   = struct2sym(set);

syms lambda zeta v pai s c h yd muw q k u zI zY g xg mux mua gy gc gg gai gh gtfp clag ivlag hlag ydlag glag zYlag qlag
syms E1PT E2PT E3PT E4PT

%Declare Needed Symbols 
iv      = sym('iv'); 
ga      = sym('ga');

%The exogenous states and news processes
z  = [zY,mux,mua,g,zI,zeta,muw];
X  = z;
f = sym([]);
etaM = cell(1,7);
for jj = 1:length(z)
    [M,etaM{jj},xs,xsp] = nper_news(['e_' char(z(jj))],[4,8],eval(['[sig4',char(z(jj)), ' sig8' char(z(jj)) ']']));
    f(end+1:end+length(M)) = transpose(log(xsp)) - M*transpose(log(xs));
    X = [X,xs];
    eval(['syms ' char(xs(1))]);
end

% States and controls
X  = [ydlag,clag,ivlag,hlag,glag,zYlag,qlag,s,xg,k,X];
Y  = [gy,gc,gai,gh,gg,gtfp,ga,yd,c,iv,h,v,q,lambda,pai,u,E1PT,E2PT,E3PT,E4PT];
XP = make_prime(X);
YP = make_prime(Y);
make_index([Y,X]);

% Definitions
muy     = mux*mua^(alphak/(alphak-1));
muy_p   = mux_p*mua_p^(alphak/(alphak-1)); 
mui     = muy/mua;
mui_p   = muy_p/mua_p;
muibar  = muxbar*muabar^(1/(alphak-1));
phi     = kappa/2*(iv*mui/ivlag - muibar)^2;
phi_p   = kappa/2*(iv_p*mui_p/ivlag_p - muibar)^2;
dphi    = kappa*(iv*mui/ivlag - muibar);
dphi_p  = kappa*(iv_p*mui_p/ivlag_p - muibar);
delta2  = rdelta*delta1;
delta   = delta0 + delta1*(u-1) + delta2/2*(u-1)^2;
delta_p = delta0 + delta1*(u_p-1) + delta2/2*(u_p-1)^2;
ddelta  = delta1 + delta2*(u-1);

% Eq. conditions with expectations
f(end+1) = -lambda + zeta*v^(-sig) - gam*pai*s_p/(c - b*clag/muy) - E1PT;% - bet*b*(zeta_p*(v_p*muy_p)^(-sig) - gam*pai_p/muy_p^sig*(s_p/(c_p*muy_p - b*c))^(1-gam));
f(end+1) = -pai + psii*zeta*v^(-sig)*h^theta + E2PT;%+ bet*(1-gam)*pai_p*muy_p^(-sig)*(s_p/(c_p*muy_p - b*c))^(-gam);
f(end+1) = -q + E3PT;%+bet*lambda_p/lambda*muy_p^(-sig)*(alphak*yd_p/k_p*muy_p + q_p*mua_p*(1-delta_p));
f(end+1) = -1 +q*zI*(1 - phi - dphi*(iv*mui/ivlag))  + E4PT;%+ bet*lambda_p/lambda*muy_p^(-sig)*q_p*mua_p*zI_p*dphi_p*(iv_p*mui_p/iv)^2;

% Eq. conditions without expectations (other than endog states)
f(end+1) = -v + c - b*clag/muy - psii*h^theta*s_p;
f(end+1) = (c - b*clag/muy)^gam*(s/muy)^(1-gam) -s_p ;
f(end+1) =  alphah*yd/muw*lambda -theta*psii*zeta*v^(-sig)*h^theta*s_p;
f(end+1) = -alphak*yd/k*mui + q*u*ddelta;
f(end+1) = -yd + zY*(u*k/mui)^alphak*h^alphah*L^(1-alphak-alphah);
f(end+1) = -k_p + (1-delta)*k/mui + (1-phi)*zI*iv;
f(end+1) = -yd + c + iv + g*xg_p;
f(end+1) = -xg_p + xg^rhoxg/muy;

%Define expectations terms
f(end+1) = E1PT + fullm*(- bet*b*(zeta_p*(v_p*muy_p)^(-sig) - gam*pai_p/muy_p^sig*(s_p/(c_p*muy_p - b*c))^(1-gam)));
f(end+1) = E2PT - fullm*(bet*(1-gam)*pai_p*muy_p^(-sig)*(s_p/(c_p*muy_p - b*c))^(-gam));
f(end+1) = E3PT - fullm*(bet*lambda_p/lambda*muy_p^(-sig)*(alphak*yd_p/k_p*muy_p + q_p*mua_p*(1-delta_p)));
f(end+1) = E4PT - fullm*(bet*lambda_p/lambda*muy_p^(-sig)*q_p*mua_p*zI_p*dphi_p*(iv_p*mui_p/iv)^2);

% Lagged variables
f(end+1) = -ydlag_p + yd;
f(end+1) = -clag_p  + c;
f(end+1) = -ivlag_p + iv;
f(end+1) = -hlag_p  + h;
f(end+1) = -glag_p  + g;
f(end+1) = -zYlag_p + zY;
f(end+1) = -qlag_p  + q;

% Observables
f(end+1) = -gy   + yd/ydlag*muy;
f(end+1) = -gc   + c/clag*muy;
f(end+1) = -gai  + iv/ivlag*mua*mui;
f(end+1) = -gh   + h/hlag;
f(end+1) = -gg   + g/glag*xg^(rhoxg-1);
f(end+1) = -gtfp + zY/zYlag*mux^(1-alphak);
f(end+1) = -ga   + mua;

% Exogenous processes
f(end+1) = log(zY_p)  - rhozY  *log(zY)                + sguinfo*log(e_zY1);
f(end+1) = log(mux_p/muxbar) - rhomux *log(mux/muxbar) + sguinfo*log(e_mux1);
f(end+1) = log(mua_p/muabar) - rhomua *log(mua/muabar) + sguinfo*log(e_mua1);
f(end+1) = log(g_p/gbar)   - rhog   *log(g/gbar)       + sguinfo*log(e_g1);
f(end+1) = log(zI_p)  - rhozI  *log(zI)                + sguinfo*log(e_zI1);
f(end+1) = log(zeta_p)- rhozeta*log(zeta)              + sguinfo*log(e_zeta1);
f(end+1) = log(muw_p/muwbar) - rhomuw *log(muw/muwbar) + sguinfo*log(e_muw1);

% disp(['Nx:   ', num2str(length(X))]);
% disp(['Neq:  ', num2str(length(f))]);
% disp(['Nvar: ', num2str(length(X)+length(Y))]);

%Log-linear approx (Pure linear if log_var = [])
xlog = true(1,length(X));
ylog = true(1,length(Y)); ylog(1,e4pt_idx) = false;
log_var = [X(xlog) Y(ylog) XP(xlog) YP(ylog)];


mod.f = subs(f, log_var, exp(log_var));
mod.X = X;
mod.XP = XP;
mod.Y = Y;
mod.YP = YP;
mod.PARAM = PARAM;
mod.param = param;
mod.SET = SET;
mod.set = set;
mod.adiff = set.adiff; %Include anaylytical derivatives?
mod.xlog = xlog;
mod.ylog = ylog;


%Standard Errors
nx = length(X);
ny = length(Y);
mod.shck = sym(zeros(nx,7*3));

%Shock order: [zY,mux,mua,g,zI,zeta,muw], then news;
for jj = 1:7
    eval(['mod.shck('   lower(char(z(jj))) '_idx-ny ,3*(jj-1)+1) = sig0' char(z(jj)) ';'])
    eval(['mod.shck(e_' lower(char(z(jj))) '4_idx-ny,3*(jj-1)+2) = sig4' char(z(jj)) ';'])
    eval(['mod.shck(e_' lower(char(z(jj))) '8_idx-ny,3*(jj-1)+3) = sig8' char(z(jj)) ';'])
end

%Measurement Error (parameters are std-devs in param.m file)
mod.me = [];

%Derivatives using numerical toolbox
mod = anal_deriv(mod);

%Save index variables for use in main_prog
!rm -f v_idx.mat
save v_idx *_idx



