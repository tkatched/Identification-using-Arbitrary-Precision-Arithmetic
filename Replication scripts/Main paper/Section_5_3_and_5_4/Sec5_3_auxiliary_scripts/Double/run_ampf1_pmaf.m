%Searches the PMAF region with the benchmark model given by AMPF1

%Specify parameters for search:
parstruc(1,1).ind=[1:7];


for jj=1:1 %consider invertible case
clearvars -except jj kk parstruc ppool
%This script performs search for observationally equivalent points for the
%cashless Leeper (1991) model under determinacy
%% Setup

% Choose baseline and alternative regimes: regime =1 for AMPF, regime=2 for
% PMAF, regime =3 for indeterminacy (so far supported only as alternative regime)

breg=1; %benchmark model regime   

areg=2; %objective function will search this regime

% Choose shock specification in base and alternative models (1-MA(1), 2-MA(2), 3-ARMA(1,1))

mselb=1;

msela=1;

%Choose whether to optimize over shock parameters only or all parameters
%(1= all parameters, 2=shock parameters only):

parsel=parstruc(1,1).ind; %select parameters to optimize over

invsel=jj; %Choose whether to search invertible/noninvertible MA in alternative regimes (1=(Inv,Inv); 2=(Inv,Ninv);3=(Ninv,Inv);4=(Ninv,Ninv))

bc=0; %0 for full spectrum; 1 for business cycle frequencies only
%% Base parameter values in the two regimes:
alpha1=1.5;
gamma1=1.2;
alpha2=0.3;
gamma2=0.1;
beta=0.9804;
sigt=1;
sigpsi=1;

 %% Get shock parameters that observe equivalence restriction derived using the method in Tan and Walker (2015)
        
        
        
        %Compute MA parameters and form baseline parameter vectors for each regime:
        
        if mselb==1
            a1=0.5;
            b1=0.5;
            c1=0.5;
            d1=0.5;
            theta1=[alpha1,beta,gamma1,sigt,sigpsi,a1,c1]; % base parameter vector, AMPF regime
            theta2=[alpha2,beta,gamma2,sigt,sigpsi,b1,d1]; %base parameter vector, PMAF regime
        elseif mselb==2
            a1=0.5;
            a2=0.5;
            b1=0.5;
            b2=0.5;
            c1=0.5;
            c2=0.5;
            d1=0.5;
            d2=0.5;
            theta1=[alpha1,beta,gamma1,sigt,sigpsi,a1,a2,c1,c2]; % base parameter vector, AMPF regime
            theta2=[alpha2,beta,gamma2,sigt,sigpsi,b1,b2,d1,d2]; %base parameter vector, PMAF regime
        end
%% Create and save the benchmark model spectrum for use by objective function

%determine baseline parameter vector

eval(['thetab=theta',num2str(breg),';']);
eval(['thetaa=theta',num2str(areg),';']);

n=100; %number of points to evaluate the integral
[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well 
    x=x(1:n/2+1);
    wg=wg(1:n/2+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end

[T1,TC,T0,TETA,RC] = lpsolv2(thetab,mselb);
sqi=-1i;
neq=size(T1,2); %number of equations
id1=eye(neq);
cc=2*pi;

ny=2;

truespec=zeros(ny*length(x),ny);

A=[eye(ny),zeros(ny,neq-ny)]; %selection matrix A

QQ=createcov_lp(thetab,mselb,RC);

RR = [T0,zeros(neq,1)];
for i=1:length(x)
    exe=exp(sqi*x(i));
    mat1=(id1-T1*exe)\id1; %inv(1-T1L)
    mat2=mat1';  %note that ' gives the conjugate transpose
    truespec(((i-1)*ny+1):i*ny,:)=A*mat1*RR*QQ*RR'*mat2*A'/cc;
end

%% Clean up before running the KL minimization

clearvars -except breg areg mselb msela thetab thetaa parsel invsel jj kk parstruc bc n x wg truespec ppool
%% Algorithm selection and Matlab versions
%Select one procedure:
runga=0; %1 if using Genetic algorithm + Multistart combination, 0 otherwise
runpso=1; %1 if using Particle Swarm + Multistart combination, 0 otherwise
          %-----Note PSO is available only in R2014b and later----
NonlinCon='penalty'; %this is the preferred option, available in the matlab version 2014b and later; for earlier versions, set to 'auglag'.
%% Set up parallel computation

%  numcore=14; %specify the number of cores
% % 
% % %for versions before R2014a
% % matlabpool open local numcore
% % 
% % %for versions R2014a and above
% ppool=parpool('local',numcore); %create parallel pool object
%% Set parameter bounds


%Parameter order: 
%[alpha beta gamma sigt sigpsi mat1 mat2  mapsi1 mapsi2 ]


 %Determine policy parameter bounds according to regimes
     
     %Determine policy parameter bounds according to regimes
     if areg==1
         %search active monetary
         
         lb1=[1.01 0.9 1.01 ];
         ub1=[3   0.999  3  ];
         
     elseif areg==2
         lb1=[0.01 0.9 0.01 ];
         ub1=[0.99  0.999 0.99];
     elseif areg==3
         lb1=[0.01 0.9 1.01 ];
         ub1=[0.99  0.999 3];
     end
     




if msela==1
    if invsel==1 %all invertible
        lb2=[0.001 0.001 -0.999 -0.999];
        
        ub2=[10   10      0.999  0.999];
    elseif invsel==2
        lb2=[0.001 0.001 -0.999 -5];
        
        ub2=[10   10      0.999  5];
    elseif invsel==3
        lb2=[0.001 0.001 -5 -0.999];
        
        ub2=[10   10      5  0.999];
    elseif invsel==4
        lb2=[0.001 0.001 -5 -5];
        
        ub2=[10   10      5  5];
    end
    
elseif msela==2
    lb2=[0.001 0.001 -5 -5 -5 -5]; % (non)invertibility will be enforced in constraint
    
    ub2=[10   10      5  5  5  5];
elseif msela==3 %option for ARMA shocks (for now set stationary and invertible)
    lb2=[0.001 0.001 -0.99 -0.999 -0.99 -0.999];
    
    ub2=[10   10      0.99  0.999 0.99  0.999];
end

%append indeterminacy parameters if necessary
if areg==3
    lb3=[-5,-5,0];
    ub3=[ 5, 5, 10];
else
    lb3=[];
    ub3=[];
end


lb=[lb1,lb2,lb3];
ub=[ub1,ub2,ub3];

lb=lb(parsel);
ub=ub(parsel);

numpar=length(lb);%number of parameters in the objective function

%% Detailed specifications; you should not have to modify them

con1=[0;0]; %constraint indicator for objective function used by GA and fmincon

%Specify neighborhood and parameter constraint characteristics
wgt=[];
nrm=[];
indp=[];

%Decide whether to call constraint in the objective for PSO
if invsel>1 || msela==2
    con2=[1;0]; %constraint indicator for objective function used by PSO
    ConstraintFunction=@(xest)constraintlp(xest,thetab,thetaa,msela,parsel,invsel,con2(2),wgt,nrm,indp); %set constraint for the problem

else
    con2=[0;0];
    ConstraintFunction=[];
end


if bc==0
    resfilename=[['Lpr_klbch_br'],num2str(breg),'_ar',num2str(areg),'_sb',num2str(mselb),'_sa',num2str(msela),'_inv',num2str(invsel),'_par',num2str(parsel)];
elseif bc==1
    resfilename=[['Lpr_klbch_bc_br'],num2str(breg),'_ar',num2str(areg),'_sb',num2str(mselb),'_sa',num2str(msela),'_inv',num2str(invsel),'_par',num2str(parsel)];
end

ObjectiveFunction = @(theta0)kloptlp(theta0,thetab,thetaa,truespec,x,wg,bc,msela,areg,parsel,invsel,con1,wgt,nrm,indp); %set objective function for GA/multistart
ObjectiveFunctionP = @(theta0)kloptlp(theta0,thetab,thetaa,truespec,x,wg,bc,msela,areg,parsel,invsel,con2,wgt,nrm,indp); %set objective function for PSO


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
    
swarmsize=300; %swarm size (similar concept to population size for GA)
maxitpso=1000; %maximum number of iterations (similar concept to generations for GA)
stiterlim=100; %max number of stall PSO iterations (i.e., no improvement found)
initswarm=[]; %set initial population (if smaller than swarmsize, MATLAB will
%randomly draw the rest. If [], the whole population is randomly drawn.
%Can be a row vector of dimension numpar or a matrix. Each row of the
%matrix is then a candidate initial value.

minfn=0.1; %smallest fraction of neighbors for PSO (smallest size of the 
%adaptive neighborhood)

tolfunpso=1e-06; %tolerance level for improvement in the objective for PSO

psoname=[['psolp_int_br'],num2str(breg),'_ar',num2str(areg),'_sb',num2str(mselb),'_sa',num2str(msela),'_inv',num2str(invsel),'_par',num2str(parsel)]; %set name for a temp output file that stores the swarms (problem-based)
OutFun=@(optimValues,state)psout(optimValues,state,psoname); %output 
%function for extracting swarms from PSO runs for further local optimization.

useppso=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []


end
%% Multistart algorithm settings

numrpoints=50; %number of random starting points for Multistart

usepms=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []

% settings for fmincon
maxit=1000; % set max number of iterations
maxfev=10000; % set max number of function evaluations
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
%% Empirical distance between regimes
alph=0.05; %significance level

P80 = pfhlp(thetab,xestpso2,mselb,msela,alph,80,100,0);
P150 = pfhlp(thetab,xestpso2,mselb,msela,alph,150,100,0);
P200 = pfhlp(thetab,xestpso2,mselb,msela,alph,200,100,0);
P1K = pfhlp(thetab,xestpso2,mselb,msela,alph,1000,100,0);

disp 'Empirical distance for T = 80, 150, 200, 1000'
[P80;P150;P200;P1K]

disp 'Deviations of empirical distance from 1 for T = 80, 150, 200, 1000'
[1-P80;1-P150;1-P200;1-P1K]
save(resfilename)
    end
