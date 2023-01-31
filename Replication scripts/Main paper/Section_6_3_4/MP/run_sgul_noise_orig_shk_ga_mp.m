%% Script to examine robustness of global identification results in multiple precision

%% Setup


parsel=[14:34]; %search all parameters except measurement error (not present in observables in levels)


bc=0; %0 for full spectrum; 1 for business cycle frequencies only

%Define filenames:
if bc==0
    resfilename=['Sgul_noise_orig_shk'];
elseif bc==1
    resfilename=['Sgul_noise_orig_shk_bc'];
end
resfilename2=[resfilename,'_mp_mins']; %new result file name
%% Create and save the benchmark model spectrum for use by objective function
load(resfilename,'thetab','msel1','msel2','thetabnoise') %load the benchmark parameter

thetab=mp(thetab);
n=mp('400'); %number of points to evaluate the integral
[x,wg] = quadcomp(n,bc);

if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    x=x(1:n/2+1);
    wg=wg(1:n/2+1);
else
    x=x(1:n/2);
    wg=wg(1:n/2);
end


[gx,hx,eta,npara] = sgusolvn_mp(thetab,msel1);
neq=size(hx,2); %number of equations

ny=7; %no of observables

truespec=mp(zeros(ny*length(x),ny));

A=mp(zeros(ny,22)); %selection matrix A

%vector of observables: [y; cons;inv;hrs;g;tfp;a_growth]
% Some auxiliary quantities defined in mp
ci=mp('-1i'); %imaginary -i
pi2=mp('2*pi'); %2*pi in mp
onep=mp('1');
twop=mp('2');
hhmp=mp('100');

eq_y = 1;
eq_c = 2;
eq_i = 3;
eq_h = 4;
eq_g = 5;
eq_z = 6;
eq_agr = 7;

%variable indices for selection
v_i	     = 10;
v_c	     = 9;
v_y	     = 8;
v_h	     = 11;
v_g	     = 21;
v_z	     = 22;
v_mua    = 7;

%fill selection matrix
A(eq_y,v_y)=onep;

A(eq_c,v_c)=onep;

A(eq_i,v_i)=onep;

A(eq_h,v_h)=onep;

A(eq_g,v_g)=onep;

A(eq_z,v_z)=onep;

A(eq_agr,v_mua)=hhmp;

gx2=[gx;mp(zeros(2,neq))]; %add selection of levels of g and z
gx2(21,5)=onep;
gx2(22,6)=onep;

tempeye=mp(eye(neq));

for i=1:length(x) %compute spectrum under noise
    exe=exp(ci*x(i));
    
    mat1=(tempeye-hx*exe)\tempeye;
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    truespec(((i-1)*ny+1):i*ny,:)=(A*gx2*mat1*eta*eta'*mat2*gx2'*A')/(pi2); %compute spectrum using representation in the paper
    
end

%% Recreate parameter bounds in MP
lb=[mp('[0.001 ,7E-06,0.01,0.01,0.00001,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99]'),mp(zeros(1,22))];
ub=[mp('[15,0.999,35,15,0.999,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99]'),mp(ones(1,22))*mp('20')];

lb=lb(parsel);
ub=ub(parsel);
%% Additional local optimization using MP version of fminsearch

%load double precision optimization results (minimizer and obj function value in double precision)

load(resfilename,'thetaa','fvalga2','xestga2','con2','wgt','nrm','indp')
thetaa=mp(thetaa); %mp value of the benchmark parameters of the alternative regime

thetamp=mp(xestga2);

%set fminsearch options
optfms=optimset('fminsearch');
optfms=optimset(optfms,'Display','iter','MaxFunEvals',5000000,'MaxIter',50000,'TolFun',1e-30,'TolX',1e-16);

ObjFunmp=@(theta0)kloptsguln_mp(theta0,thetaa,truespec,msel2,x,wg,bc,parsel,con2,wgt,nrm,indp,lb,ub);

%perform optimization:

tic;
[xestmp,fvalmp,eflagmp,outmp] = fminsearch(ObjFunmp,transfmpsgu(thetamp,1,lb,ub),optfms);
timelfms=toc;

%% Arrange and save results
thetaindmp=transfmpsgu(xestmp,2,lb,ub); %recover the new minimizer
klmp=fvalmp; %record the Kl distance

%% Display results

klmp_o=ObjFunmp(transfmpsgu(thetamp,1,lb,ub)); %original KL in mp

disp('-----------------------------------------------')

disp('                 Original minimizer                              New minimizer')
disp([thetamp',thetaindmp']);

display ('                     Original KL                                     New KL');
display ([klmp_o,klmp]);
save(resfilename2) %save results into a new file