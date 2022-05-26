function G = gmatlpv2(theta0,x,wg,h,dif,msel)

%Function to compute local identification condition for the Leeper (1991) model at
%a point theta0 under both determinacy and indeterminacy.

%Parallelization is implemented when computing spectra for numerical
%derivatives.

%Inputs: theta0 - vector of model parameters in the order:
%(msel=1): theta=%[alpha beta gam sigt sigpsi mat1 mapsi1]';

%(msel=2): theta=%[alpha beta gam sigt sigpsi mat1 mat2 mapsi1 mapsi2]';

%(msel=3): theta=%[alpha beta gam sigt sigpsi rhot1 mat1 rhopsi1 mapsi1]';
%
% x - Gauss-Legendre quadrature abscissae, only on one side about zero as
% the code exploits symmetry of quadrature abscissae/weights and the
% equivalence of spectral density to its complex conjugate at frequencies
% symmetric about zero.
% w - Gauss-Legendre quadrature weights, corresponding to abscissae x above
% h - step size in numerical differentiation

%Observables order: 1) inflation, 2) real debt

%uses auxiliary files createcov_lp.m, leepsolv.m, gensys_mod.m, svdrr.m,

%The numerical derivative is computed using: 1) symmetric difference quotient
%2-point method if dif=1; 2) symmetric quotient difference 4-point method otherwise

%msel - selects specification of shocks

%% Model solution

ntheta=length(theta0); %number of parameters

[TT,TC,TEPS,TETA,RC] = lpsolv2(theta0,msel);
nstate = size(TT,2); %number of state variables
dimy=2; %no of observables

QQ=createcov_lp(theta0,msel,RC); %create covariance matrix

if RC==[1;1] %if determinacy
    RR = [TEPS,zeros(nstate,1)];
    
    
elseif RC==[1;0] %if indeterminacy
    TETA=rref(TETA').'; %reduced column echelon form for TETA
    RR = [TEPS,TETA]; %note - to conform with Gauss/LS gensys
    
else
    G=[];
    disp('Nonexistence or numerical problems')
    return
end



%% Spectrum computation
%set up selection matrix

nfreq=length(x); %no of frequencies used

A=[eye(dimy),zeros(dimy,nstate-dimy)]; %selection matrix A


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
        [nT1,nTC,nTEPS,nTETA,nRC] = lpsolv2(npara,msel);
        [pT1,pTC,pTEPS,pTETA,pRC] = lpsolv2(ppara,msel);
        
        nQQ=createcov_lp(npara,msel,nRC); %augmented covariance matrix
        
        if nRC==[1;1]; % determinacy
            
            nRR=[nTEPS,zeros(nstate,1)];
            
        elseif nRC==[1;0]; % indeterminacy
            nTETA=rref(nTETA').'; %reduced column echelon form for nTETA
            nRR = [nTEPS,nTETA]; %note - to conform with Gauss/LS gensys
            
        else
            G=[];
            disp('Nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        
        pQQ=createcov_lp(ppara,msel,pRC); %augmented covariance matrix
        
        if pRC==[1;1]; % determinacy
            
            pRR = [pTEPS,zeros(nstate,1)];
            
        elseif pRC==[1;0]; % indeterminacy
            
            pRR = [pTEPS,pTETA];
            pTETA=rref(pTETA').'; %reduced column echelon form for nTETA
            
        else
            G=[];
            disp('Nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        
        
        parfor i=1:nfreq;
            exe=exp(-1i*x(i));
            
            nmat1=(eye(nstate)-nT1*exe)\eye(nstate);
            nmat2=nmat1';
            
            pmat1=(eye(nstate)-pT1*exe)\eye(nstate);
            pmat2=pmat1';
            
            tempn=A*nmat1*nRR*nQQ*nRR'*nmat2*A'/(2*pi);
            tempp=A*pmat1*pRR*pQQ*pRR'*pmat2*A'/(2*pi);
            
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
        [nT1,nTC,nTEPS,nTETA,nRC] = lpsolv2(npara,msel);
        [pT1,pTC,pTEPS,pTETA,pRC] = lpsolv2(ppara,msel);
        [n2T1,n2TC,n2TEPS,n2TETA,n2RC] = lpsolv2(n2para,msel);
        [p2T1,p2TC,p2TEPS,p2TETA,p2RC] = lpsolv2(p2para,msel);
        
        nQQ=createcov_lp(npara,msel,nRC); %augmented covariance matrix
        
        if nRC==[1;1]; % determinacy
            
            nRR=[nTEPS,zeros(nstate,1)];
            
        elseif nRC==[1;0]; % indeterminacy
            nTETA=rref(nTETA').'; %reduced column echelon form for nTETA
            nRR = [nTEPS,nTETA]; %note - to conform with Gauss/LS gensys
            
        else
            G=[];
            disp('Nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        n2QQ=createcov_lp(n2para,msel,n2RC); %augmented covariance matrix
        
        if n2RC==[1;1]; % determinacy
            
            n2RR=[n2TEPS,zeros(nstate,1)];
            
        elseif n2RC==[1;0]; % indeterminacy
            n2TETA=rref(n2TETA').'; %reduced column echelon form for nTETA
            n2RR = [n2TEPS,n2TETA]; %note - to conform with Gauss/LS gensys
            
        else
            G=[];
            disp('Nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        
        pQQ=createcov_lp(ppara,msel,pRC); %augmented covariance matrix
        
        if pRC==[1;1]; % determinacy
            
            pRR = [pTEPS,zeros(nstate,1)];
            
        elseif pRC==[1;0]; % indeterminacy
            
            pRR = [pTEPS,pTETA];
            pTETA=rref(pTETA').'; %reduced column echelon form for nTETA
            
        else
            G=[];
            disp('Nonexistence or numerical problems - cannot compute derivative')
            return
            
        end
        
        p2QQ=createcov_lp(p2para,msel,p2RC); %augmented covariance matrix
        
        if p2RC==[1;1]; % determinacy
            
            p2RR = [p2TEPS,zeros(nstate,1)];
            
        elseif p2RC==[1;0]; % indeterminacy
            
            p2RR = [p2TEPS,p2TETA];
            p2TETA=rref(p2TETA').'; %reduced column echelon form for nTETA
            
        else
            G=[];
            disp('Nonexistence or numerical problems - cannot compute derivative')
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
            
            tempn=A*nmat1*nRR*nQQ*nRR'*nmat2*A'/(2*pi); %spectra at perturbed theta0
            tempn2=A*n2mat1*n2RR*n2QQ*n2RR'*n2mat2*A'/(2*pi);
            tempp=A*pmat1*pRR*pQQ*pRR'*pmat2*A'/(2*pi);
            tempp2=A*p2mat1*p2RR*p2QQ*p2RR'*p2mat2*A'/(2*pi);
            
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