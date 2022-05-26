% PARAMETERS - This function returns a parameter structure to use in the model solution.


function [param,set,GG,GG2,Gxmat] = parameters(theta0)
                                                       

set.adiff      = 0;
set.approx_deg = 1;


%Full or partial info
set.fullm = 1;

% Calibrated parameters
set.bet   = 0.99;
set.alphak = 0.225;
set.alphah = 0.675;
set.delta0 = 0.025;
set.muybar = 1.0045;
set.muabar = 0.9957;
set.G_Y    = 0.2;
set.hbar   = 0.2;
set.L      = 1;
set.sig  = 1;
set.muwbar = 1.15;


%Estimated parameters
% load parameter;
% param_in       = median_bayesian;
% param.theta    = param_in(1) + 1; 
% param.gam      = param_in(2);
% param.kappa    = param_in(3);
% param.rdelta   = param_in(4);
% param.b        = param_in(5)*0.99; 
% param.rhoxg    = param_in(6)*0.99; 
% param.rhozY    = param_in(7)*0.99; 
% param.rhomua   = param_in(8)*0.99; 
% param.rhog     = param_in(9)*0.99; 
% param.rhomux   = param_in(10)- 0.5; 
% param.rhomuw   = param_in(11)*0.99; 
% param.rhozeta  = param_in(12)*0.99; 
% param.rhozI    = param_in(13)*0.99; 
% param.sig0zY   = param_in(14);
% param.sig4zY   = param_in(15);
% param.sig8zY   = param_in(16);
% param.sig0mua  = param_in(17);
% param.sig4mua  = param_in(18);
% param.sig8mua  = param_in(19);
% param.sig0g    = param_in(20);
% param.sig4g    = param_in(21);
% param.sig8g    = param_in(22);
% param.sig0mux  = param_in(23);
% param.sig4mux  = param_in(24);
% param.sig8mux  = param_in(25);
% param.sig0muw  = param_in(26);
% param.sig4muw  = param_in(27);
% param.sig8muw  = param_in(28);
% param.sig0zeta = param_in(29);
% param.sig4zeta = param_in(30);
% param.sig8zeta = param_in(31);
% param.sig0zI   = param_in(32);
% param.sig4zI   = param_in(33);
% param.sig8zI   = param_in(34);
% param.sigyme   = param_in(35);

%% Baseline parameterization consistent with our code
param.theta    = theta0(1); 
param.gam      = theta0(2);
param.kappa    = theta0(3);
param.rdelta   = theta0(4);
param.b        = theta0(5); 
param.rhoxg    = theta0(6); 
param.rhozY    = theta0(10); 
param.rhomua   = theta0(7); 
param.rhog     = theta0(12); 
param.rhomux   = theta0(8); 
param.rhomuw   = theta0(11); 
param.rhozeta  = theta0(13); 
param.rhozI    = theta0(9); 
param.sig0zY   = theta0(23);
param.sig4zY   = theta0(24);
param.sig8zY   = theta0(25);
param.sig0mua  = theta0(14);
param.sig4mua  = theta0(15);
param.sig8mua  = theta0(16);
param.sig0g    = theta0(29);
param.sig4g    = theta0(30);
param.sig8g    = theta0(31);
param.sig0mux  = theta0(17);
param.sig4mux  = theta0(18);
param.sig8mux  = theta0(19);
param.sig0muw  = theta0(26);
param.sig4muw  = theta0(27);
param.sig8muw  = theta0(28);
param.sig0zeta = theta0(32);
param.sig4zeta = theta0(33);
param.sig8zeta = theta0(34);
param.sig0zI   = theta0(20);
param.sig4zI   = theta0(21);
param.sig8zI   = theta0(22);
param.sigyme   = theta0(35);

%% Alternative parameterization to check noise shocks (news shock std. dev. changed - noise2)
% param.theta    = 4.74; 
% param.gam      = 0.0019;
% param.kappa    = 9.11;
% param.rdelta   = 0.34;
% param.b        = 0.91; 
% param.rhoxg    = 0.72; 
% param.rhozY    = 0.92; 
% param.rhomua   = 0.48; 
% param.rhog     = 0.96; 
% param.rhomux   = 0.38; 
% param.rhomuw   = 0.98; 
% param.rhozeta  = 0.17; 
% param.rhozI    = 0.47; 
% param.sig0zY   = 0.55;
% param.sig4zY   = 0.2;
% param.sig8zY   = 0.05;
% param.sig0mua  = 0.25;
% param.sig4mua  = 0.10;
% param.sig8mua  = 0.10;
% param.sig0g    = 0.7;
% param.sig4g    = 0.4;
% param.sig8g    = 0.6;
% param.sig0mux  = 0.25;
% param.sig4mux  = 0.2;
% param.sig8mux  = 0.2;
% param.sig0muw  = 0.8;
% param.sig4muw  = 3;
% param.sig8muw  = 0.3;
% param.sig0zeta = 5;
% param.sig4zeta = 1;
% param.sig8zeta = 3;
% param.sig0zI   = 9;
% param.sig4zI   = 3;
% param.sig8zI   = 4;
% param.sigyme   = 0.30;
%%

% Implied steady-state values
set.mua    = set.muabar;
set.muy    = set.muybar; %1.004500000000000
set.muw    = set.muwbar;
set.mux    = set.muy*set.mua^(set.alphak/(1-set.alphak));
set.muxbar = set.mux;
set.mui    = set.mux*set.mua^(1/(set.alphak-1));
set.h      = set.hbar;
set.xg     = set.muy^(1/(param.rhoxg-1));

%Determined in equilbrium
set.ybar  = NaN;      
set.mcbar = NaN;
set.ikbar = NaN;
set.irbar = NaN;
set.delta1 = NaN;
set.delta2 = NaN;
set.psii   = NaN;
set.gbar = NaN;
%**************************************************************************
% CHOOSE INFO STRUCTURE TO SOLVE WITH
%**************************************************************************
set.sguinfo = true;   %true mean use SGU rep; false means use noise rep


