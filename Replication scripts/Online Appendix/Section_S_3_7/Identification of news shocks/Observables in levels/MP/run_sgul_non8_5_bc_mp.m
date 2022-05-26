%% Script to examine robustness of global identification results in multiple precision

                %% Setup
        
        
        parsel=[1:34]; %mua news shock variances set to zero
        parsel(28)=[];
        
        bc=1; %0 for full spectrum; 1 for business cycle frequencies only
        
                %Define filenames:
        if bc==0
            resfilename=['Sgul_non8_5'];
        elseif bc==1
            resfilename=['Sgul_non8_5_bc'];
        end
        resfilename2=[resfilename,'_mp_mins']; %new result file name
                %% Create and save the benchmark model spectrum for use by objective function
                load(resfilename,'thetab') %load the benchmark parameter
        
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
        
        [TT,TC,TEPS,TETA,RC] = sgusolvl_mp(thetab);
        neq=size(TT,2); %number of equations

        ny=7;
        
        truespec=mp(zeros(ny*length(x),ny));
        
        A=mp(zeros(ny,neq)); %selection matrix A
        
        %vector of observables: [y_growth; cons_growth;iinv_growth;hrs_growth;g_growth;tfp_growth;a_growth]
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
v_i	     = 69;
v_c	     = 70;
v_y	     = 71;
v_h	     = 73;
v_g	     = 49;
v_z	     = 31;
v_mua    = 4;
        
        %fill selection matrix
A(eq_y,v_y)=onep; 

A(eq_c,v_c)=onep; 

A(eq_i,v_i)=onep; 

A(eq_h,v_h)=onep; 

A(eq_g,v_g)=onep; 

A(eq_z,v_z)=onep; 

A(eq_agr,v_mua)=hhmp;


tempeye=mp(eye(neq)); %auxiliary identity matrix

QQ=mp(zeros(21,21));
for i=14:34
    QQ(i-13,i-13)=(thetab(i)/hhmp)^twop; %put variances on the diagonal (dividing by 100 in HS(2014) code)
    %QQ(i-13,i-13)=theta0(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
end

        
        RR = TEPS;
for i=1:length(x) %compute spectrum
    exe=exp(ci*x(i));
    
    mat1=(tempeye-TT*exe)\tempeye;
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    truespec(((i-1)*ny+1):i*ny,:)=(A*mat1*RR*QQ*RR'*mat2*A')/(pi2); %compute spectrum using representation in the paper
     
end

%% Recreate parameter bounds in MP
        lb=[mp('[0.001 ,7E-06,0.01,0.01,0.00001,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99,-0.99]'),mp(zeros(1,21))];
        ub=[mp('[15,0.999,35,15,0.999,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99]'),mp(ones(1,21))*mp('20')];
        
                lb=lb(parsel);
        ub=ub(parsel);
                %% Additional local optimization using MP version of fminsearch
        
        %load double precision optimization results (minimizer and obj function value in double precision)
     
        load(resfilename,'thetaa','fvalpso2','xestpso2','con2','wgt','nrm','indp')
         thetaa=mp(thetaa); %mp value of the benchmark parameters of the alternative regime
         
         thetamp=mp(xestpso2);
         
                 %set fminsearch options
        optfms=optimset('fminsearch');
        optfms=optimset(optfms,'Display','iter','MaxFunEvals',50000,'MaxIter',5000,'TolFun',1e-30,'TolX',1e-12);
        
        ObjFunmp=@(theta0)kloptsgul_mp(theta0,thetab,thetaa,truespec,x,wg,bc,parsel,con2,wgt,nrm,indp,lb,ub);
        
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
        save(resfilename2) %save results