%Searches theta0 vs. MA(1) shocks

        %This script performs search for observationally equivalent points for the
        %cashless Leeper (1991) model under determinacy
        %% Setup
        
        
        parsel=[1:34]; %all news shock variances set to zero
        
        
        bc=1; %0 for full spectrum; 1 for business cycle frequencies only
        
        msel1=1;
        msel2=8;
        %% Base parameter values:
        thet      = 4.74; %Frisch elasticity of labor supply (when gam=b=0)
        gam	      = 0.0019; %governs magnitude of wealth elasticity of labor supply
        %NB: small values of gamma crucial for importance of news shocks
        kap	      = 9.11; %investment adj. cost parameter
        del2_del1 = 0.34; %ratio of sensitivity of capacity utilization to variations to rental rate of capital (del2)
        %to the steady state capacity utilization
        b	      = 0.91; %habit formation parameter
        rhoxg     = 0.72; %smoothness in trend of gov. spending
        
        rhoa      = 0.48; %nonstationary investment-specific productivity shock AR
        rhox      = 0.38; %neutral productivity shock
        rhozi     = 0.47; %investment-specific productivity shock AR
        rhoz      = 0.92; %productivity shock AR
        rhom      = 0.98; %wage markup shock
        rhog      = 0.96; %gov. spending shock
        rhozet    = 0.17; %preference shock
        
        %standard deviations of the 21 shocks
        sig_mua0 = 0.21; %
        sig_mua4 = 0.16;
        sig_mua8 = 0.16;
        sig_mux0 = 0.38;
        sig_mux4 = 0.08;
        sig_mux8 = 0.1;
        sig_zi0  = 11.72;
        sig_zi4  = 1.93;
        sig_zi8  = 5.50;
        sig_z0   = 0.65;
        sig_z4   = 0.11;
        sig_z8   = 0.09;
        sig_mu0  = 0.5;
        sig_mu4  = 4.79;
        sig_mu8  = 0.51;
        sig_g0   = 0.62;
        sig_g4   = 0.57;
        sig_g8   = 0.37;
        sig_zet0 = 4.03;
        sig_zet4 = 1.89;
        sig_zet8 = 2.21;
        
        thetab=[thet;gam;kap;del2_del1;b;rhoxg;rhoa;rhox;rhozi;rhoz;rhom;rhog;rhozet;sig_mua0;sig_mua4;
            sig_mua8;sig_mux0;sig_mux4;sig_mux8;sig_zi0;sig_zi4;sig_zi8;sig_z0;sig_z4;sig_z8;sig_mu0;sig_mu4;sig_mu8;
            sig_g0; sig_g4;sig_g8;sig_zet0;sig_zet4;sig_zet8]';
        
        thetaa=thetab;
        intsel=[32,33,34];
        thetaa(intsel)=0.1;
      
        %% Create and save the benchmark model spectrum for use by objective function
        
        n=400; %number of points to evaluate the integral
        [x,wg] = quadcomp(n,bc);
        
        if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
            x=x(1:n/2+1);
            wg=wg(1:n/2+1);
        else
            x=x(1:n/2);
            wg=wg(1:n/2);
        end
        
        [TT,TC,TEPS,TETA,RC] = sgusolvl4(thetab,msel1);
        neq=size(TT,2); %number of equations

        ny=7;
        
        truespec=zeros(ny*length(x),ny);
        
        A=zeros(ny,neq); %selection matrix A
        
        %vector of observables: [y_growth; cons_growth;iinv_growth;hrs_growth;g_growth;tfp_growth;a_growth]

eq_y = 1;
eq_c = 2;
eq_i = 3;
eq_h = 4;
eq_g = 5;
eq_z = 6;
eq_agr = 7;

%variable indices for selection
v_i	     = 69;
v_c	     = 70;
v_y	     = 71;
v_h	     = 73;
v_g	     = 49;
v_z	     = 31;
v_mua    = 4;
        
        %fill selection matrix
A(eq_y,v_y)=1; 

A(eq_c,v_c)=1; 

A(eq_i,v_i)=1; 

A(eq_h,v_h)=1; 

A(eq_g,v_g)=1; 

A(eq_z,v_z)=1; 

A(eq_agr,v_mua)=100;

QQ=zeros(21,21);
for i=14:34
    QQ(i-13,i-13)=(thetab(i))^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
    %QQ(i-13,i-13)=theta0(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
end

        
        RR = TEPS;
for i=1:length(x) %compute spectrum
    exe=exp(-1i*x(i));
    
    mat1=(eye(neq)-TT*exe)\eye(neq);
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    truespec(((i-1)*ny+1):i*ny,:)=(A*mat1*RR*QQ*RR'*mat2*A')/(2*pi); %compute spectrum using representation in the paper
     
end

        
        %% Clean up before running the KL minimization
        
        %clearvars -except breg areg mselb msela thetab thetaa parsel invsel jj kk ii parstruc ns bc n x wg truespec
        %% Algorithm selection and Matlab versions
        %Select one procedure:
        runga=0; %1 if using Genetic algorithm + Multistart combination, 0 otherwise
        runpso=1; %1 if using Particle Swarm + Multistart combination, 0 otherwise
        %-----Note PSO is available only in R2014b and later----
        NonlinCon='penalty'; %this is the preferred option, available in the matlab version 2014b and later; for earlier versions, set to 'auglag'.
        %% Set up parallel computation
        
%          numcore=8; %specify the number of cores
%         %
%         % %for versions before R2014a
%         % matlabpool open local numcore
%         %
%         % %for versions R2014a and above
%         ppool=parpool('local',numcore); %create parallel pool object
        %% Set parameter bounds
        
        
        %Parameter order:
        
        
        
        lb=[0.001 ,7E-06,0.01,0.01,0.00001,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,zeros(1,21)];
        ub=[15,0.999,35,15,0.999,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,ones(1,21)*20];
        
        lb(intsel)=[0.01,-7,-5];
        ub(intsel)=[20,7,5];
        lb=lb(parsel);
        ub=ub(parsel);
        
        numpar=length(lb);%number of parameters in the objective function
        
        %% Detailed specifications; you should not have to modify them
        
        con1=[0;0]; %constraint indicator for objective function used by GA and fmincon
        con2=[0;0]; %constraint indicator for objective function used by PSO
        
        %Specify neighborhood and parameter constraint characteristics
        
        nrm=[]; %set norm for the constraint function (1,2, or Inf)
        indp=[]; %vector of parameter indices to be constrained
        wgt=[]; %vector of weights for the constraint function
        
        %Specify constraint function:
        
        ConstraintFunction=[]; %set constraint for the problem
        
        
        if bc==0
            resfilename=['Sgul_ma2_msel_',num2str(msel2)];
        elseif bc==1
            resfilename=['Sgul_ma2_msel_',num2str(msel2),'_bc'];
        end
        
        ObjectiveFunction = @(theta0)kloptsgulinv3(theta0,thetab,thetaa,msel2,intsel,truespec,x,wg,bc,parsel,con1,wgt,nrm,indp); %set objective function for GA/multistart
        ObjectiveFunctionP = @(theta0)kloptsgulinv3(theta0,thetab,thetaa,msel2,intsel,truespec,x,wg,bc,parsel,con2,wgt,nrm,indp); %set objective function for PSO
        
        
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
            
            popsize=300; %population size
            elcnt=3; %elite count - number of elite individuals retained in population
            
            tolfunga=1e-10; %tolerance level for improvement in the objective for GA
            tolconga=1e-10; %tolerance level for constraint for GA
            
            usepga=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []
            %In later versions of Matlab, 1 and 0 can also be used respectively.
        end
        
        %% PSO algorithm settings
        if runpso==1
            
            swarmsize=1000; %swarm size (similar concept to population size for GA)
            maxitpso=1000; %maximum number of iterations (similar concept to generations for GA)
            stiterlim=100; %max number of stall PSO iterations (i.e., no improvement found)
            initswarm=[]; %set initial population (if smaller than swarmsize, MATLAB will
            %randomly draw the rest. If [], the whole population is randomly drawn.
            %Can be a row vector of dimension numpar or a matrix. Each row of the
            %matrix is then a candidate initial value.
            
            minfn=0.1; %smallest fraction of neighbors for PSO (smallest size of the
            %adaptive neighborhood)
            
            tolfunpso=1e-06; %tolerance level for improvement in the objective for PSO
            
            psoname=['psosgul_ma2_msel_',num2str(msel2),'_',num2str(bc)]; %set name for a temp output file that stores the swarms (problem-based)
            OutFun=@(optimValues,state)psout(optimValues,state,psoname); %output
            %function for extracting swarms from PSO runs for further local optimization.
            
            useppso=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []
            
            
        end
        %% Multistart algorithm settings
        
        numrpoints=50; %number of random starting points for Multistart
        
        usepms=['Always']; %Set to 'Always' to use parallel computation, otherwise to 'Never' or []
        
        % settings for fmincon
        maxit=1000; % set max number of iterations
        maxfev=20000; % set max number of function evaluations
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
