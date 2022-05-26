function G = gmatsguv2(theta0,x,wg,h,dif)

%Function to compute local identification condition for the SGU (2012) model at
%a point theta0, using partially parallelized computation.

%Inputs: theta0 - vector of model parameters in the order:
%                 [thet;gam;kap;del2_del1;b;rhoxg;rhoa;rhox;rhozi;rhoz;rhom;rhog;rhozet;sig_mua0;sig_mua4;
%             sig_mua8;sig_mux0;sig_mux4;sig_mux8;sig_zi0;sig_zi4;sig_zi8;sig_z0;sig_z4;sig_z8;sig_mu0;sig_mu4;sig_mu8;
%             sig_g0; sig_g4;sig_g8;sig_zet0;sig_zet4;sig_zet8;sig_me]';
%
% x - Gauss-Legendre quadrature abscissae, only on one side about zero as
% the code exploits symmetry of quadrature abscissae/weights and the
% equivalence of spectral density to its complex conjugate at frequencies 
% symmetric about zero.
% w - Gauss-Legendre quadrature weights, corresponding to abscissae x above
% h - step size in numerical differentiation

%Observables order:  [y_growth; cons_growth;iinv_growth;hrs_growth;g_growth;tfp_growth;a_growth]

%uses auxiliary files, gensys_mod2.m, svdrr.m, sgusolv.m,
%sgusspar.m, sgueqvarshocks.m, sgusetsysmat.m

%The numerical derivative is computed using: 1) symmetric difference quotient
%2-point method if dif=1; 2) symmetric quotient difference 4-point method otherwise

%msel - selects specification of shocks (TBA)

%% Model solution

ntheta=length(theta0); %number of parameters

[TT,TC,TEPS,TETA,RC] = sgusolv(theta0);
nstate = size(TT,2); %number of state variables
dimy=7; %no of observables

        
if RC~=[1;1] %if not determinacy
    G=[];
    disp(RC)
    disp('Indeterminacy,nonexistence or numerical problems')
    return
end



%% Spectrum computation
%set up selection matrix

nfreq=length(x); %no of frequencies used

A=zeros(dimy,nstate); %selection matrix A

%vector of observables: [y_growth; cons_growth;iinv_growth;hrs_growth;g_growth;tfp_growth;a_growth]
eq_ygr = 1;
eq_cgr = 2;
eq_igr = 3;
eq_hgr = 4;
eq_ggr = 5;
eq_tfp = 6;
eq_agr = 7;

%variable indices for selection
v_i	     = 72;
v_c	     = 73;
v_y	     = 74;
v_h	     = 76;
v_g	     = 49;
v_z	     = 31;
v_mua    = 4;
v_muy    = 1;
v_xg     = 3;
v_mux    = 13;
v_ylag   = 81; %lags 
v_clag   = 82;
v_ilag   = 83;
v_hlag   = 84;
v_glag   = 67; 
v_xglag  = 68;
v_zlag   = 69;

alpk  = 0.225; %calibrated capital share

%fill selection matrix
A(eq_ygr,v_y)=100; %(100(1-L) )
A(eq_ygr,v_ylag)=-100; %(100(1-L) )
A(eq_ygr,v_muy)=100;

A(eq_cgr,v_c)=100; %(100(1-L))
A(eq_cgr,v_clag)=-100; %(100(1-L))
A(eq_cgr,v_muy)=100;

A(eq_igr,v_i)=100; %(100(1-L))
A(eq_igr,v_ilag)=-100; %(100(1-L))
A(eq_igr,v_muy)=100;

A(eq_hgr,v_h)=100; %(100(1-L))
A(eq_hgr,v_hlag)=-100; %(100(1-L))

A(eq_ggr,v_g)=100; %(100(1-L))
A(eq_ggr,v_glag)=-100; %(100(1-L))
A(eq_ggr,v_muy)=100;
A(eq_ggr,v_xg)=100; %(100(1-L))
A(eq_ggr,v_xglag)=-100; %(100(1-L))

A(eq_tfp,v_z)=100; %(100(1-L))
A(eq_tfp,v_zlag)=-100; %(100(1-L))
A(eq_tfp,v_mux)=100*(1-alpk);

A(eq_agr,v_mua)=100;
%% Numerical derivative
dimy2=dimy^2;

if dif==1 %use symmetric Newton's quotient
    nsp_mat=zeros(dimy2*nfreq,ntheta); %blank for derivatives (a-d)
    psp_mat=nsp_mat; %blank for derivatives (a+d)
    
    
    
    for j=1:ntheta;
        
        %perturb one parameter at a time
        npara=theta0;
        ppara=theta0;
        npara(j)=theta0(j)-h;
        ppara(j)=theta0(j)+h;
        
        
        
        %recompute the spectral density matrices)
        if j<14 %for parameters affecting gensys solution, resolve the model
        [nT1,nTC,nTEPS,nTETA,nRC] = sgusolv(npara);
        [pT1,pTC,pTEPS,pTETA,pRC] = sgusolv(ppara);
        else %keep the solution at theta0 as shock parameters do not affect gensys solution
            nT1=TT;
            nTC=TC;
            nTEPS=TEPS;
            nTETA=TETA;
            nRC=RC;
            
            pT1=TT;
            pTC=TC;
            pTEPS=TEPS;
            pTETA=TETA;
            pRC=RC;
        end
        
        
        %form the covariance matrix
        nQQ=zeros(21,21);
        for i=14:34
            nQQ(i-13,i-13)=(npara(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %nQQ(i-13,i-13)=npara(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        
        nsme=zeros(dimy,dimy);
        nsme(eq_ygr,eq_ygr)=npara(35)^2; %covariance matrix of measurement error
        
        if nRC==[1;1]; % determinacy
            
            nRR=nTEPS;
            
        else
            G=[];
            disp(nRC)
            disp('Indeterminacy, nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        
        pQQ=zeros(21,21);
        for i=14:34
            pQQ(i-13,i-13)=(ppara(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %pQQ(i-13,i-13)=ppara(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        
        psme=zeros(dimy,dimy);
        psme(eq_ygr,eq_ygr)=ppara(35)^2; %covariance matrix of measurement error
        
        if pRC==[1;1]; % determinacy
            
            pRR = pTEPS;
            
        else
            G=[];
            disp(pRC)
            disp('Indeterminacy, nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        
        
        parfor i=1:nfreq;
            exe=exp(-1i*x(i));
            
            nmat1=(eye(nstate)-nT1*exe)\eye(nstate);
            nmat2=nmat1';
            
            pmat1=(eye(nstate)-pT1*exe)\eye(nstate);
            pmat2=pmat1';
            
            tempn=(A*nmat1*nRR*nQQ*nRR'*nmat2*A'+nsme)/(2*pi); %spectra at perturbed theta0
            tempp=(A*pmat1*pRR*pQQ*pRR'*pmat2*A'+psme)/(2*pi); %spectra at perturbed theta0
            
            tempn2(:,:,i)=tempn(:); %spectra at perturbed theta0
            tempp2(:,:,i)=tempp(:); %spectra at perturbed theta0
            
        end
        
        for i=1:nfreq
            nsp_mat(((i-1)*dimy2+1):(i*dimy2),j)=tempn2(:,:,i); %spectra at perturbed theta0
            psp_mat(((i-1)*dimy2+1):(i*dimy2),j)=tempp2(:,:,i); %spectra at perturbed theta0
        end
    end
    
    %now construct the G matrix, note that the jkth component is simply the derivative with respect to j and k
    %a symmetric difference quotient method is used
    
    Gmat=zeros(ntheta,ntheta); %blank
    
    nd1=(psp_mat-nsp_mat)/(2*h);
    nd2=conj(nd1);
    
    for freq=1:nfreq;
        
        
        dv1=nd1(((freq-1)*dimy2+1):(freq*dimy2),:);
        dv2=nd2(((freq-1)*dimy2+1):(freq*dimy2),:);
        
        Gmat=Gmat+wg(freq)*(dv1'*dv1+dv2'*dv2);
        
    end
else %use 4-point derivative approximation method
    
    nsp_mat=zeros(dimy2*nfreq,ntheta); %blank for derivatives (a-h)
    psp_mat=nsp_mat; %blank for derivatives (a+h)
    n2sp_mat=nsp_mat; %blank for derivatives (a-2h)
    p2sp_mat=nsp_mat; %blank for derivatives (a+2h)
    
    for j=1:ntheta;
        
        %perturb one parameter at a time
        npara=theta0;
        ppara=theta0;
        n2para=theta0;
        p2para=theta0;
        npara(j)=theta0(j)-h;
        ppara(j)=theta0(j)+h;
        n2para(j)=theta0(j)-2*h;
        p2para(j)=theta0(j)+2*h;
        
        
        
        %recompute the spectral density matrices)
        if j<14
        [nT1,nTC,nTEPS,nTETA,nRC] = sgusolv(npara);
        [pT1,pTC,pTEPS,pTETA,pRC] = sgusolv(ppara);
        [n2T1,n2TC,n2TEPS,n2TETA,n2RC] = sgusolv(n2para);
        [p2T1,p2TC,p2TEPS,p2TETA,p2RC] = sgusolv(p2para);
        
        else

            nT1=TT;
            nTC=TC;
            nTEPS=TEPS;
            nTETA=TETA;
            nRC=RC;
            
            n2T1=TT;
            n2TC=TC;
            n2TEPS=TEPS;
            n2TETA=TETA;
            n2RC=RC;
            
            pT1=TT;
            pTC=TC;
            pTEPS=TEPS;
            pTETA=TETA;
            pRC=RC;
            
            p2T1=TT;
            p2TC=TC;
            p2TEPS=TEPS;
            p2TETA=TETA;
            p2RC=RC;
        end
        %form the covariance matrix
        nQQ=zeros(21,21);
        for i=14:34
            nQQ(i-13,i-13)=(npara(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %nQQ(i-13,i-13)=npara(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        
        nsme=zeros(dimy,dimy);
        nsme(eq_ygr,eq_ygr)=npara(35)^2; %covariance matrix of measurement error
        
        if nRC==[1;1]; % determinacy
            
            nRR=nTEPS;
            
        else
            G=[];
            disp(nRC)
            disp('Indeterminacy, nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        %form the covariance matrix
        n2QQ=zeros(21,21);
        for i=14:34
            n2QQ(i-13,i-13)=(n2para(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %n2QQ(i-13,i-13)=n2para(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        
        n2sme=zeros(dimy,dimy);
        n2sme(eq_ygr,eq_ygr)=n2para(35)^2; %covariance matrix of measurement error
        
        if n2RC==[1;1]; % determinacy
            
            n2RR=n2TEPS;
            
        else
            G=[];
            disp(n2RC)
            disp('Indeterminacy, nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        
        pQQ=zeros(21,21);
        for i=14:34
            pQQ(i-13,i-13)=(ppara(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %pQQ(i-13,i-13)=ppara(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        
        psme=zeros(dimy,dimy);
        psme(eq_ygr,eq_ygr)=ppara(35)^2; %covariance matrix of measurement error
        
        if pRC==[1;1]; % determinacy
            
            pRR = pTEPS;
            
        else
            G=[];
            disp(pRC)
            disp('Indeterminacy, nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        p2QQ=zeros(21,21);
        for i=14:34
            p2QQ(i-13,i-13)=(p2para0(i)/100)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %p2QQ(i-13,i-13)=p2para(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        
        p2sme=zeros(dimy,dimy);
        p2sme(eq_ygr,eq_ygr)=p2para(35)^2; %covariance matrix of measurement error
        
        if p2RC==[1;1]; % determinacy
            
            p2RR = p2TEPS;
            
        else
            G=[];
            disp(p2RC)
            disp('Indeterminacy, nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        parfor i=1:nfreq;
            exe=exp(-1i*x(i));
            
            
            nmat1=(eye(nstate)-nT1*exe)\eye(nstate);
            nmat2=nmat1';
            
            n2mat1=(eye(nstate)-n2T1*exe)\eye(nstate);
            n2mat2=n2mat1';
            
            pmat1=(eye(nstate)-pT1*exe)\eye(nstate);
            pmat2=pmat1';
            
            p2mat1=(eye(nstate)-p2T1*exe)\eye(nstate);
            p2mat2=p2mat1';
            
            tempn=(A*nmat1*nRR*nQQ*nRR'*nmat2*A'+nsme)/(2*pi); %spectra at perturbed theta0
            tempn2=(A*n2mat1*n2RR*n2QQ*n2RR'*n2mat2*A'+n2sme)/(2*pi); %spectra at perturbed theta0
            tempp=(A*pmat1*pRR*pQQ*pRR'*pmat2*A'+psme)/(2*pi); %spectra at perturbed theta0
            tempp2=(A*p2mat1*p2RR*p2QQ*p2RR'*p2mat2*A'+p2sme)/(2*pi); %spectra at perturbed theta0
            
            junkn(:,:,i)=tempn(:);
            junkn2(:,:,i)=tempn2(:);
            junkp(:,:,i)=tempp(:);
            junkp2(:,:,i)=tempp2(:);
            
        end
        
        for i=1:nfreq
            nsp_mat(((i-1)*dimy2+1):(i*dimy2),j)=junkn(:,:,i); %spectra at perturbed theta0
            n2sp_mat(((i-1)*dimy2+1):(i*dimy2),j)=junkn2(:,:,i); %spectra at perturbed theta0
            psp_mat(((i-1)*dimy2+1):(i*dimy2),j)=junkp(:,:,i); %spectra at perturbed theta0
            p2sp_mat(((i-1)*dimy2+1):(i*dimy2),j)=junkp2(:,:,i); %spectra at perturbed theta0
        end
    end
    
    %now construct the G matrix, note that the jkth component is simply the derivative with respect to j and k
    %a 4-point method is used
    
    Gmat=zeros(ntheta,ntheta); %blank
    
    nd1=(8*psp_mat-8*nsp_mat+n2sp_mat-p2sp_mat)/(12*h); %numerical derivatives stacked
    nd2=conj(nd1);
    
    for freq=1:nfreq;
        
        
        dv1=nd1(((freq-1)*dimy2+1):(freq*dimy2),:);
        dv2=nd2(((freq-1)*dimy2+1):(freq*dimy2),:);
        
        Gmat=Gmat+wg(freq)*(dv1'*dv1+dv2'*dv2);
        
    end
end
G=real(Gmat);



end