%This script performs Particle Swarm algorithm optimization with additional local
%optimization by Multistart/fmincon.

%The script uses options and problem descriptions (objective and constraint
%function if applicable)specified in the parent script of the type
%QT_TablesX_part_Y.m

%% Set up and run PSO
rng(0) %set default seed 

%Set PSO options
psopt=optimoptions('particleswarm','UseParallel',useppso,'Display',dispalg,'DisplayInterval',dispint,'SwarmSize',swarmsize,'TolFun',tolfunpso,'MaxIter',maxitpso,'StallIterLimit',stiterlim,'MinFractionNeighbors',minfn,'OutputFcns',OutFun,'InitialSwarm',initswarm);

%Run PSO
tic;
[xestpso,fvalpso,exitflagpso,outputpso]=particleswarm(ObjectiveFunctionP,numpar,lb,ub,psopt); %note penalized objective used
load(psoname) %import data from the output function
delete([psoname,'.mat']) %clean up the output file
timelpso1=toc;
%% Set up and run Multistart

%Set fmincon options

optfmc=optimset('fmincon'); %create default option structure for fmincon
optfmc=optimset('MaxIter',maxit,'MaxFunEvals',maxfev,'TolFun',tolfunfmc,'TolCon',tolconfmc,'TolX',tolx,'Algorithm',localg);

mspso = MultiStart('UseParallel',usepms,'Display',dispalg); %set multistart structure

psopoints1=swarmfinal(1:49,:); %take first 50 points of the final population
psopoints2=swarmfinal([50:((swarmsize-50)/10):swarmsize]',:); %add 10 equally spaced points from the rest
psopoints=CustomStartPointSet([psopoints1;psopoints2]); %stack together starting points for Multistart
rpointspso=RandomStartPointSet('NumStartPoints',numrpoints); %preset random points
spointspso={psopoints,rpointspso};


problempso = createOptimProblem('fmincon','x0',xestpso,...
'objective',ObjectiveFunction,'lb',lb,'ub',ub,'nonlcon',ConstraintFunction,...
'options',optfmc); %create problem structure
tic;
[xestpso2,fvalpso2,exitflagpso2,outputpso2,minimapso] = run(mspso,problempso,spointspso); %run multistart with 50 extra random points
save(resfilename) %save intermediate result before running multistart
timelfmc=toc;