%Verifies results in arbitrary precision for Leeper et al (2017) model


 %% Specify case and retrieve respective results

%Select benchmark and alternative regimes

breg=2; %benchmark model regime (only supports AMPF=1, PMAF=2)

areg=1; %objective function will search this regime (supports AMPF=1, PMAF=2, PMPF=3)

bc=0; %0 for full spectrum; 1 for business cycle frequencies only


if bc==0
    resfilename=[['Lpr17_kl_br'],num2str(breg),'_ar',num2str(areg)];
elseif bc==1
    resfilename=[['Lpr17_kl_bc_br'],num2str(breg),'_ar',num2str(areg)];
end

resfilename2=[resfilename,'_mp_mins']; %new result file name


%% Create and save the benchmark model spectrum

load(resfilename,'thetab') %load the benchmark parameter

thetab=mp(thetab);

n=mp('100'); %number of points to evaluate the integral
[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    x=x(1:n/2+1);
    wg=wg(1:n/2+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end


truespec = comp_spec2_mp(thetab,x,breg); %baseline spectrum

%% Clean up before running the KL minimization

clearvars -except breg areg thetab resfilename resfilename2 bc n x wg truespec

 %% Recreate parameter bounds in MP

%Determine policy parameter bounds according to regimes

if areg==1 %searching AMPF
    load posterior_regimeM.mat;
    
    %Set bounds to 90% interval
    lb=mp(prcdraw_90(1,:));
    ub=mp(prcdraw_90(2,:));
    
    lgn = mp('2')*(ub-lb);
    lb = mp(meandraw)-lgn;
    ub = mp(meandraw)+lgn;
    
    lb(1)=mp('0.01');
    lb(2)=mp('0.1');
    lb(5)=mp('0.1');
    lb(6)=mp('0.1');
    lb(7)=mp('0.01');
    lb(8)=mp('0.01');
    lb(9)=mp('1.01');
    lb(15)=mp('0.01');
    %ub(19)=mp('0.99');
    %ub(20)=mp('0.99');
    lb(21)=mp('1');
    lb(27)=mp('0.1');
    ub(30)=mp('0.9');
    lb(33)=mp('0');
    clear meandraw postdraw postlike prcdraw_90 % clean up a bit
    
elseif areg==2 %searching PMAF
    load posterior_regimeF.mat;
    
    %Set bounds to 90% interval
    lb=mp(prcdraw_90(1,:));
    ub=mp(prcdraw_90(2,:));
    
    lgn = mp('2')*(ub-lb);
    lb = mp(meandraw)-lgn;
    ub = mp(meandraw)+lgn;
    
    lb(2)=mp('0.1');
    lb(5)=mp('0.1');
    lb(6)=mp('0.1');
    lb(7)=mp('0.01');
    lb(8)=mp('0.01');
    lb(9)=mp('0.1');
    lb(15)=mp('0.01');
    ub(19)=mp('0.99');
    ub(20)=mp('0.99');
    lb(22)=mp('1');
    lb(28)=mp('0.1');
    lb(35)=mp('0');
    clear meandraw postdraw postlike prcdraw_90 % clean up a bit
    
elseif areg==3
    load posterior_regimeF.mat;
    
    %Set bounds to 90% interval based on the PMAF regime, then enlarged to
    %accomodate greater variation.
    
    lb=mp(prcdraw_90(1,:));
    slb = min(mp('2')*lb,lb/mp('2'));
    lb = min (slb,lb);
    ub=mp(prcdraw_90(2,:));
    sub = max(mp('2')*ub,ub/mp('2'));
    ub = max (sub,ub);
    
    %introduce nine additional parameters
    %first eight are projections of sunspot on structural shocks,
    %the last one is the standard deviation of the sunspot shock
    
    %parasun = zeros(9,1);
    %parasun(9) = 1;
    
    %update the lower and upper bounds
    lb = [lb,-mp('2')*mp(ones(1,8)),mp('0.01')];
    ub = [ub,mp('1')*mp(ones(1,8)),mp('10')];
    numpar = length(lb);
    lb(1)=mp('0.001');%make some bounds wider
    ub(3)=mp('0.99');
    ub(4)=mp('0.99'); %stickiness
    lb(9)=mp('0');
    ub(9)=mp('0.9');
    lb(11)=mp('-0.5');
    ub(11)=mp('0.5');
    lb(12)=mp('-0.5');
    ub(12)=mp('0.9');
    ub(15)=mp('0.9');
    ub(16)=mp('0.99');
    ub(18)=mp('0.99');
    ub(19)=mp('0.99');
    ub(20)=mp('0.99');
    ub(28)=mp('4');
    ub(30)=mp('0.5');
    lb(31)=mp('0');
    ub(31)=mp('0.99');
    ub(33)=mp('0.99');
    ub(34)=mp('500');
    lb(35)=mp('0');
    ub(35)=mp('1.5');
    lb(38)=mp('-4');
    %lb(44)=-4;
    lb(43)=mp('-8');
    lb(40)=mp('-4');
    lb(41)=mp('-4');
    lb(42)=mp('-4');
    
    lb(44)=mp('1e-10');
    
    
    clear meandraw postdraw postlike prcdraw_90 % clean up a bit
end


numpar=length(lb);%number of parameters in the objective function

%% Additional local optimization using MP version of fminsearch

load(resfilename,'fvalpso2','xestpso2')

thetamp=mp(xestpso2);

%set fminsearch options
optfms=optimset('fminsearch');
optfms=optimset(optfms,'Display','iter','MaxFunEvals',50000,'MaxIter',5000,'TolFun',1e-30,'TolX',1e-12);

ObjFunmp=@(theta0)kloptlpr2_mp(theta0,truespec,x,wg,bc,areg,lb,ub);

%perform optimization:

tic;
[xestmp,fvalmp,eflagmp,outmp] = fminsearch(ObjFunmp,transfmplpr(thetamp,1,lb,ub),optfms);
timelfms=toc;
        
%% Arrange and save results

thetaindmp=transfmplpr(xestmp,2,lb,ub); %recover the new minimizer
klmp=fvalmp; %record the Kl distance

%% Empirical distance
T=mp('200'); %sample size
a=mp('0.05'); %significance level
ed200 = pfhlpr_mp(thetab,thetaindmp,breg,areg,a,T,n,bc); %empirical distance at T=200, alpha=5%
%% Display results

klmp_o=ObjFunmp(transfmplpr(thetamp,1,lb,ub)); %original KL in mp

disp('-----------------------------------------------')

disp('                 Original minimizer                              New minimizer')
disp([thetamp',thetaindmp']);

display ('                     Original KL                                     New KL');
display ([klmp_o,klmp]);

disp(' Empirical distance at T=200');
disp(ed200);
save(resfilename2) %save results into a new file

