%This script verifies that values on the nonidentification curve in Table 2 result
%in observational equivalence (Direction 1)

%Specify parameters for search:
parstruc(1,1).ind=[1:8]; %modify nonzero parameters only


    for jj=1:10
        clearvars -except kk jj parstruc
        %% Specify case and retrieve respective results
        
        %description of search problem
        
        
        breg=1; %benchmark model regime
        areg=3; %objective function will search this regime
        
        % Choose shock specification in base and alternative models (1-MA(1), 2-MA(2), 3-ARMA(1,1))
        mselb=1;
        
        msela=1;
        
        %Choose whether to optimize over shock parameters only or all parameters
        %(1= all parameters, 2=shock parameters only):
        parsel=parstruc(1,1).ind; %select parameters to optimize over
        
        invsel=1; %Choose whether to search invertible/noninvertible MA in alternative regimes (1=(Inv,Inv); 2=(Inv,Ninv);3=(Ninv,Inv);4=(Ninv,Ninv))
        
        
        bc=0; %choose frequencies
        
%         pl=50; %set precision level in digits
%         
%         mp.Digits(pl); %impose precision level
        
        %Define filenames:
        if bc==0
            resfilename=[['Lpr_kl_br'],num2str(breg),'_ar',num2str(areg),'_sb',num2str(mselb),'_sa',num2str(msela),'_inv',num2str(invsel),'_par',num2str(parsel)];
        elseif bc==1
            resfilename=[['Lpr_kl_bc_br'],num2str(breg),'_ar',num2str(areg),'_sb',num2str(mselb),'_sa',num2str(msela),'_inv',num2str(invsel),'_par',num2str(parsel)];
        end
        %resfilename2=[resfilename,'_mp',num2str(pl)]; %new result file name
        resfilename2=[resfilename,'_mp_ampf2_cdir1_',num2str(jj)];
        
        
        
        
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
        
        [T1,TC,T0,TETA,RC] = lpsolv_mp2(thetab,mselb);
        
        neq=size(T1,2); %number of equations
        
        sqi=mp('-1i');
        id1=mp(eye(neq));
        cc=mp('2*pi');
        
        ny=2;
        
        truespec=mp(zeros(ny*length(x),ny));
        
        A=[mp(eye(ny)),mp(zeros(ny,neq-ny))]; %selection matrix A
        
        QQ=createcov_lpmp(thetab,mselb,RC);
        
        RR = [T0,mp(zeros(neq,1))];
        for i=1:length(x)
            exe=exp(sqi*x(i));
            mat1=(id1-T1*exe)\id1; %inv(1-T1L)
            mat2=mat1';  %note that ' gives the conjugate transpose
            truespec(((i-1)*ny+1):i*ny,:)=A*mat1*RR*QQ*RR'*mat2*A'/cc;
        end
        
        %% Clean up before running the KL minimization
        
        clearvars -except breg areg mselb msela thetab parsel invsel resfilename resfilename2 bc kk jj n x wg truespec parstruc
        
        %% Recreate parameter bounds in MP
        %Parameter order:
        %[alpha beta gamma sigt sigpsi mat1 mat2  mapsi1 mapsi2 ]
        
        
        %Determine policy parameter bounds according to regimes
        
        %Determine policy parameter bounds according to regimes
        if areg==1
            %search active monetary
            
            lb1=mp('[1.01 0.9 1.01 ]');
            ub1=mp('[3   0.999  3  ]');
            
        elseif areg==2
            lb1=mp('[0.01 0.9 0.01 ]');
            ub1=mp('[0.99  0.999 0.99]');
        elseif areg==3
            lb1=mp('[0 0.9 1.01 ]');
            ub1=mp('[0.99  0.999 3]');
        end
        
        
        
        
        
        if msela==1
            if invsel==1 %all invertible
                lb2=mp('[0.001 0.001 -1 -0.999]');
                
                ub2=mp('[10   10      0.999  0.999]');
            elseif invsel==2
                lb2=mp('[0.001 0.001 -0.999 -5]');
                
                ub2=mp('[10   10      0.999  5]');
            elseif invsel==3
                lb2=mp('[0.001 0.001 -5 -0.999]');
                
                ub2=mp('[10   10      5  0.999]');
            elseif invsel==4
                lb2=mp('[0.001 0.001 -5 -5]');
                
                ub2=mp('[10   10      5  5]');
            end
            
        elseif msela==2
            lb2=mp('[0.001 0.001 -5 -5 -5 -5]'); % (non)invertibility will be enforced in constraint
            
            ub2=mp('[10   10      5  5  5  5]');
        elseif msela==3 %option for ARMA shocks (for now set stationary and invertible)
            lb2=mp('[0.001 0.001 -0.99 -0.999 -0.99 -0.999]');
            
            ub2=mp('[10   10      0.99  0.999 0.99  0.999]');
        end
        
        %append indeterminacy parameters if necessary
        if areg==3
            lb3=mp('[-5,-5,0]');
            ub3=mp('[ 10, 5, 10]');
        else
            lb3=mp('[]');
            ub3=mp('[]');
        end
        
        
        lb=[lb1,lb2,lb3];
        ub=[ub1,ub2,ub3];
        
        lb=lb(parsel);
        ub=ub(parsel);
        
        %% Additional local optimization using MP version of fminsearch
        
        %load double precision optimization results (minimizer and obj function value in double precision)
        
        load(resfilename,'thetaa','con2','wgt','nrm','indp')
        
        thetaa=mp(thetaa); %mp value of the benchmark parameters of the alternative regime
        
       
        load('curve_ampf2_dir1','thetacheck')
        
        thetamp=mp(thetacheck(jj,parstruc(1,1).ind));
      
      
        
        
        %set fminsearch options
        optfms=optimset('fminsearch');
        optfms=optimset(optfms,'Display','iter','MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-30,'TolX',1e-16);
        
        ObjFunmp=@(theta0)kloptlp_mp(theta0,thetab,thetaa,truespec,x,wg,bc,msela,areg,parsel,invsel,con2,wgt,nrm,indp,lb,ub);
        
        %perform optimization:
        
        tic;
        [xestmp,fvalmp,eflagmp,outmp] = fminsearch(ObjFunmp,transfmplp(thetamp,1,lb,ub),optfms);
        timelfms=toc;
        
        %% Arrange and save results
        thetaindmp=transfmplp(xestmp,2,lb,ub); %recover the new minimizer
        klmp=fvalmp; %record the Kl distance
        
        %% Display results
        
        
        disp('-----------------------------------------------')
        
        disp('                  New minimizer')
        disp([thetaindmp']);
        
        display ('                     New KL');
        display ([klmp]);
        save(resfilename2) %save results into a new file
    end
