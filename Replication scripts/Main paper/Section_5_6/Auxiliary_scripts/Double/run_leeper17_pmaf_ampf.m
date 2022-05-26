%Searches the alternative region for Leeper et al (2017) model

%Specify parameters for search:
%parstruc(1,1).ind=[1:35]; %all parameters in the region F

%% Setup

% Choose baseline and alternative regimes: regime =1 for AMPF, regime=0 for
% PMAF

breg=2; %benchmark model regime (only supports AMPF=1, PMAF=2)

areg=1; %objective function will search this regime (supports AMPF=1, PMAF=2, PMPF=3)

%parsel=parstruc(1,1).ind; %select parameters to optimize over

bc=0; %0 for full spectrum; 1 for business cycle frequencies only
                
%% Base parameter values and spectrum computation

%Load baseline regime

if breg==1
    load posterior_regimeM.mat; %AMPF
    
else
    
    load posterior_regimeF.mat; %PMAF
end


thetab=meandraw; %set baseline value to posterior mean

clear meandraw postdraw postlike prcdraw_90 % clean up a bit

n=100; %number of points to evaluate the integral
[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    x=x(1:n/2+1);
    wg=wg(1:n/2+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end


truespec = comp_spec2(thetab,x,breg); %baseline spectrum

%% Algorithm selection and Matlab versions
%Select one procedure:
runga=0; %1 if using Genetic algorithm + Multistart combination, 0 otherwise
runpso=1; %1 if using Particle Swarm + Multistart combination, 0 otherwise
%-----Note PSO is available only in R2014b and later----
NonlinCon='penalty'; %this is the preferred option, available in the matlab version 2014b and later; for earlier versions, set to 'auglag'.

%% Set up parallel computation (uncomment and modify as appropriate)

%  numcore=10; %specify the number of cores
% %
% % %for versions before R2014a
% % matlabpool open local numcore
% %
% % %for versions R2014a and above
% ppool=parpool('local',numcore); %create parallel pool object


%% Parameter bounds



%Determine policy parameter bounds according to regimes

if areg==1 %searching AMPF
    load posterior_regimeM.mat;
    
    %Set bounds to 90% interval
    lb=prcdraw_90(1,:);
    ub=prcdraw_90(2,:);
    
    lgn = 2*(ub-lb);
    lb = meandraw-lgn;
    ub = meandraw+lgn;
    
    lb(1)=0.01;
    lb(2)=0.1;
    lb(5)=0.1;
    lb(6)=0.1;
    lb(7)=0.01;
    lb(8)=0.01;
    lb(9)=1.01;
    lb(15)=0.01;
    %ub(19)=0.99;
    %ub(20)=0.99;
    lb(21)=1;
    lb(27)=0.1;
    ub(30)=0.9;
    lb(33)=0;
    clear meandraw postdraw postlike prcdraw_90 % clean up a bit
    
elseif areg==2 %searching PMAF
    load posterior_regimeF.mat;
    
    %Set bounds to 90% interval
    lb=prcdraw_90(1,:);
    ub=prcdraw_90(2,:);
    
    lgn = 2*(ub-lb);
    lb = meandraw-lgn;
    ub = meandraw+lgn;
    
    lb(2)=0.1;
    lb(5)=0.1;
    lb(6)=0.1;
    lb(7)=0.01;
    lb(8)=0.01;
    lb(9)=0.1;
    lb(15)=0.01;
    ub(19)=0.99;
    ub(20)=0.99;
    lb(22)=1;
    lb(28)=0.1;
    lb(35)=0;
    clear meandraw postdraw postlike prcdraw_90 % clean up a bit
    
elseif areg==3
    load posterior_regimeF.mat;
    
    %Set bounds to 90% interval based on the PMAF regime, then enlarged to
    %accomodate greater variation.
    
    lb=prcdraw_90(1,:);
    slb = min(2*lb,lb/2);
    lb = min (slb,lb);
    ub=prcdraw_90(2,:);
    sub = max(2*ub,ub/2);
    ub = max (sub,ub);
    
    %introduce nine additional parameters
    %first eight are projections of sunspot on structural shocks,
    %the last one is the standard deviation of the sunspot shock
    
    %parasun = zeros(9,1);
    %parasun(9) = 1;
    
    %update the lower and upper bounds
    lb = [lb,-2*ones(1,8),0.01];
    ub = [ub,1*ones(1,8),10];
    numpar = length(lb);
    lb(1)=0.001;%make some bounds wider
    ub(3)=0.99;
    ub(4)=0.99; %stickiness
    lb(9)=0.0;
    ub(9)=0.9;
    lb(11)=-0.5;
    ub(11)=0.5;
    lb(12)=-0.5;
    ub(12)=0.9;
    ub(15)=0.9;
    ub(16)=0.99;
    ub(18)=0.99;
    ub(19)=0.99;
    ub(20)=0.99;
    ub(28)=4;
    ub(30)=0.5;
    lb(31)=0.0;
    ub(31)=0.99;
    ub(33)=0.99;
    ub(34)=500;
    lb(35)=0;
    ub(35)=1.5;
    lb(38)=-4;
    %lb(44)=-4;
    lb(43)=-8;
    lb(40)=-4;
    lb(41)=-4;
    lb(42)=-4;
    
    lb(44)=1e-10;
    
    
    clear meandraw postdraw postlike prcdraw_90 % clean up a bit
end


numpar=length(lb);%number of parameters in the objective function


%% Detailed specifications; you should not have to modify them

%Specify neighborhood and parameter constraint characteristics

nrm=[]; %set norm for the constraint function (1,2, or Inf)
indp=[]; %vector of parameter indices to be constrained
wgt=[]; %vector of weights for the constraint function

con1=[0;0]; %constraint indicator for objective function used by GA and fmincon

if bc==0
    resfilename=[['Lpr17_kl_br'],num2str(breg),'_ar',num2str(areg)];
elseif bc==1
    resfilename=[['Lpr17_kl_bc_br'],num2str(breg),'_ar',num2str(areg)];
end

ConstraintFunction=[];

ObjectiveFunction = @(theta0)kloptlpr2(theta0,truespec,x,wg,bc,areg); %set objective function for GA/multistart
ObjectiveFunctionP = @(theta0)kloptlpr2(theta0,truespec,x,wg,bc,areg); %set objective function for PSO


dispalg='iter'; %set whether algorithm iterations are displayed.
dispint=20; %interval between displayed iterations (for PSO)

%% GA algorithm settings
if runga==1
    gen=1000; %max number of generation for GA
    stgenlim=50; %max number of stall generations (i.e., no improvement found)
    
    initpop=[]; %set initial population (if smaller than popsize, MATLAB will
    %randomly draw the rest. If [], the whole population is randomly drawn.
    %Can be a row vector of dimension numpar or a matrix. Each row of the
    %matrix is then a candidate initial value.
    
    popsize=100; %population size
    elcnt=3; %elite count - number of elite individuals retained in population
    
    tolfunga=1e-10; %tolerance level for improvement in the objective for GA
    tolconga=1e-10; %tolerance level for constraint for GA
    
    usepga=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []
    %In later versions of Matlab, 1 and 0 can also be used respectively.
end

%% PSO algorithm settings
if runpso==1
    
    swarmsize=800; %swarm size (similar concept to population size for GA)
    maxitpso=1000; %maximum number of iterations (similar concept to generations for GA)
    stiterlim=500; %max number of stall PSO iterations (i.e., no improvement found)
    initswarm=[]; %set initial population (if smaller than swarmsize, MATLAB will
    %randomly draw the rest. If [], the whole population is randomly drawn.
    %Can be a row vector of dimension numpar or a matrix. Each row of the
    %matrix is then a candidate initial value.
    
    minfn=0.1; %smallest fraction of neighbors for PSO (smallest size of the
    %adaptive neighborhood)
    
    tolfunpso=1e-06; %tolerance level for improvement in the objective for PSO
    
    psoname=[['psolpr_bch_br'],num2str(breg),'_ar',num2str(areg)]; %set name for a temp output file that stores the swarms (problem-based)
    OutFun=@(optimValues,state)psout(optimValues,state,psoname); %output
    %function for extracting swarms from PSO runs for further local optimization.
    
    useppso=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []
    
    
end
%% Multistart algorithm settings

numrpoints=50; %number of random starting points for Multistart

usepms=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []

% settings for fmincon
maxit=10000; % set max number of iterations
maxfev=50000; % set max number of function evaluations
tolfunfmc=1e-10; %tolerance level for improvement in the objective for fmincon
tolconfmc=1e-10; %tolerance level for constraint for fmincon
tolx=1e-10; %tolerance on solution value

localg='active-set'; %set main local algorithm to be used for multistart

%% Run optimization

if runga==1
    timega=tic;
    GA_optim %run GA+Multistart
    timelga=toc(timega); %time taken by GA/Multistart
    save(resfilename) %save intermediate results
end
save(resfilename)
if runpso==1
    timepso=tic;
    PSO_optim %run PSO+Multistart
    timelpso=toc(timepso); %time taken by PSO/Multistart
    save(resfilename) %save intermediate results
end

%% Empirical distance 
T=200; %sample size
a=0.05; %significance level
ed200 = pfhlpr(thetab,xestpso2,breg,areg,a,T,n,bc); %empirical distance at T=200, alpha=5%
save(resfilename)