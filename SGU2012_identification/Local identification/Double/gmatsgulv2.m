function G = gmatsgulv2(theta0,x,wg,h,dif)

%Function to compute local identification condition for the SGU (2012) model at
%a point theta0, using partially parallelized computation.

%Inputs: theta0 - vector of model parameters in the order:
% [thet;gam;kap;del2_del1;b;rhoxg;rhoa;rhox;rhozi;rhoz;rhom;rhog;rhozet;sig_mua0;sig_mua4;
%             sig_mua8;sig_mux0;sig_mux4;sig_mux8;sig_zi0;sig_zi4;sig_zi8;sig_z0;sig_z4;sig_z8;sig_mu0;sig_mu4;sig_mu8;
%             sig_g0; sig_g4;sig_g8;sig_zet0;sig_zet4;sig_zet8]';
%
% x - Gauss-Legendre quadrature abscissae, only on one side about zero as
% the code exploits symmetry of quadrature abscissae/weights and the
% equivalence of spectral density to its complex conjugate at frequencies 
% symmetric about zero.
% w - Gauss-Legendre quadrature weights, corresponding to abscissae x above
% h - step size in numerical differentiation

%Observables order:  [y; cons;iinv;hrs;g;tfp;100*a]

%uses auxiliary files gensys_mod2.m, svdrr.m, sgusolvl.m,
%sgusspar.m, sgueqvarshocksl.m, sgusetsysmatl.m

%The numerical derivative is computed using: 1) symmetric difference quotient
%2-point method if dif=1; 2) symmetric quotient difference 4-point method otherwise


%% Model solution

ntheta=length(theta0); %number of parameters

[TT,TC,TEPS,TETA,RC] = sgusolvl(theta0);
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
        [nT1,nTC,nTEPS,nTETA,nRC] = sgusolvl(npara);
        [pT1,pTC,pTEPS,pTETA,pRC] = sgusolvl(ppara);
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
            
            tempn=(A*nmat1*nRR*nQQ*nRR'*nmat2*A')/(2*pi); %spectra at perturbed theta0
            tempp=(A*pmat1*pRR*pQQ*pRR'*pmat2*A')/(2*pi); %spectra at perturbed theta0
            
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
        [nT1,nTC,nTEPS,nTETA,nRC] = sgusolvl(npara);
        [pT1,pTC,pTEPS,pTETA,pRC] = sgusolvl(ppara);
        [n2T1,n2TC,n2TEPS,n2TETA,n2RC] = sgusolvl(n2para);
        [p2T1,p2TC,p2TEPS,p2TETA,p2RC] = sgusolvl(p2para);
        
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
            
            tempn=(A*nmat1*nRR*nQQ*nRR'*nmat2*A')/(2*pi); %spectra at perturbed theta0
            tempn2=(A*n2mat1*n2RR*n2QQ*n2RR'*n2mat2*A')/(2*pi); %spectra at perturbed theta0
            tempp=(A*pmat1*pRR*pQQ*pRR'*pmat2*A')/(2*pi); %spectra at perturbed theta0
            tempp2=(A*p2mat1*p2RR*p2QQ*p2RR'*p2mat2*A')/(2*pi); %spectra at perturbed theta0
            
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