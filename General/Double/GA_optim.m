%This script performs Genetic algorithm optimization with additional local
%optimization by Multistart/fmincon.

%The script uses options and problem descriptions (objective and constraint
%function if applicable)specified in the parent script of the type
%QT_TablesX_part_Y.m

%% Set up and run GA
rng(0) %set default seed 

gaoptions=gaoptimset; %create options structure with default values
gaoptions = gaoptimset(gaoptions,'Generation',gen,'TolFun',tolfunga,'PopulationSize',popsize,'InitialPopulation', initpop,'UseParallel',usepga);
gaoptions = gaoptimset(gaoptions,'Display',dispalg,'StallGenLimit',stgenlim,'TolCon',tolconga,'PopInitRange',[lb;ub],'EliteCount',elcnt,'NonlinConAlgorithm',NonlinCon);

%Run GA
[xestga,fvalga,exitflagga,outputga,populationga,scoresga]=ga(ObjectiveFunction,numpar,[],[],[],[],lb,ub,ConstraintFunction,gaoptions);

%% Set up and run Multistart

%Set fmincon options
optfmc=optimset('fmincon'); %create default option structure for fmincon
optfmc=optimset('MaxIter',maxit,'MaxFunEvals',maxfev,'TolFun',tolfunfmc,'TolCon',tolconfmc,'TolX',tolx,'Algorithm',localg);

msga = MultiStart('UseParallel',usepms,'Display',dispalg); %set multistart structure

gapoints1=populationga(1:49,:); %take first 50 points of the final population
gapoints2=populationga([50:((popsize-50)/10):popsize]',:); %add 10 equally spaced points from the rest
gapoints=CustomStartPointSet([gapoints1;gapoints2]); %stack together starting points for Multistart
rpointsga=RandomStartPointSet('NumStartPoints',numrpoints); %preset random points
spointsga={gapoints,rpointsga};

problemga = createOptimProblem('fmincon','x0',xestga,...
'objective',ObjectiveFunction,'lb',lb,'ub',ub,'nonlcon',ConstraintFunction,...
'options',optfmc); %create problem structure

[xestga2,fvalga2,exitflagga2,outputga2,minimaga] = run(msga,problemga,spointsga); %run multistart with 50 extra random points
save(resfilename) %save intermediate result before running multistart