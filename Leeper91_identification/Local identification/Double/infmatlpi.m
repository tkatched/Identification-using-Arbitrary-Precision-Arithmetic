function G = infmatlpi(theta0,x,wg,h,dif,msel,parind)

%Function to compute local identification condition via the Information matrix for the Leeper (1991) model at
%a point theta0 under both determinacy and indeterminacy.

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

%parind - index of parameters for computing the criterion matrix. If empty, the
%default option of using all parameters is used.
%% Check parameter index
if isempty(parind)
    parind=(1:length(theta0)); %by default, compute G-matrix for the full parameter vector
end

npar=length(parind); %number of perturbed parameters
%% Model solution


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
sp_mat=zeros(dimy*nfreq,dimy); %blank for spectra at all frequencies

for i=1:nfreq
    exe=exp(-1i*x(i));
    
    mat1=(eye(nstate)-TT*exe)\eye(nstate);
    mat2=mat1';
    
    sp_mat(((i-1)*dimy+1):(i*dimy),:)=A*mat1*RR*QQ*RR'*mat2*A'/(2*pi); %spectra at perturbed theta0
    
end

%% Numerical derivative

if dif==1 %use symmetric Newton's quotient
    nsp_mat=zeros(dimy*nfreq,dimy*npar); %blank for derivatives (a-d)
    psp_mat=nsp_mat; %blank for derivatives (a+d)
    
    
    
    for j=1:npar
        
        %perturb one parameter at a time
        npara=theta0;
        ppara=theta0;
        npara(parind(j))=theta0(parind(j))-h;
        ppara(parind(j))=theta0(parind(j))+h;
        
        
        
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
        
        
        
        for i=1:nfreq;
            exe=exp(-1i*x(i));
            
            nmat1=(eye(nstate)-nT1*exe)\eye(nstate);
            nmat2=nmat1';
            
            pmat1=(eye(nstate)-pT1*exe)\eye(nstate);
            pmat2=pmat1';
            
            nsp_mat(((i-1)*dimy+1):(i*dimy),((j-1)*dimy+1):(j*dimy))=A*nmat1*nRR*nQQ*nRR'*nmat2*A'/(2*pi); %spectra at perturbed theta0
            psp_mat(((i-1)*dimy+1):(i*dimy),((j-1)*dimy+1):(j*dimy))=A*pmat1*pRR*pQQ*pRR'*pmat2*A'/(2*pi); %spectra at perturbed theta0
            
        end
    end
    
    %now construct the G matrix, note that the jkth component is simply the derivative with respect to j and k
    %a symmetric difference quotient method is used
    
    Gmat=zeros(npar,npar); %blank
    
    for j=1:npar
        
        for k=j:npar
            aj1=(psp_mat(:,((j-1)*dimy+1):(j*dimy))-nsp_mat(:,((j-1)*dimy+1):(j*dimy)))/(2*h);
            aj2=conj(aj1);
            ak1=(psp_mat(:,((k-1)*dimy+1):(k*dimy))-nsp_mat(:,((k-1)*dimy+1):(k*dimy)))/(2*h);
            ak2=conj(ak1);

            for freq=1:nfreq;
                %Gmat(j,k)=Gmat(j,k)+wg(freq)*(sum(diag(sp_mat(((freq-1)*dimy+1):(freq*dimy),:)\aj1(((freq-1)*dimy+1):(freq*dimy),:)*sp_mat(((freq-1)*dimy+1):(freq*dimy),:)\ak1(((freq-1)*dimy+1):(freq*dimy),:)))+sum(diag(conj(sp_mat(((freq-1)*dimy+1):(freq*dimy),:))\aj2(((freq-1)*dimy+1):(freq*dimy),:)*conj(sp_mat(((freq-1)*dimy+1):(freq*dimy),:))\ak2(((freq-1)*dimy+1):(freq*dimy),:))));
                Gmat(j,k)=Gmat(j,k)+wg(freq)*(sum(diag(pinv(sp_mat(((freq-1)*dimy+1):(freq*dimy),:))*aj1(((freq-1)*dimy+1):(freq*dimy),:)*pinv(sp_mat(((freq-1)*dimy+1):(freq*dimy),:))*ak1(((freq-1)*dimy+1):(freq*dimy),:)))+sum(diag(pinv(conj(sp_mat(((freq-1)*dimy+1):(freq*dimy),:)))*aj2(((freq-1)*dimy+1):(freq*dimy),:)*pinv(conj(sp_mat(((freq-1)*dimy+1):(freq*dimy),:)))*ak2(((freq-1)*dimy+1):(freq*dimy),:))));
            end
            
        end
        
    end
    Gmat=Gmat+(triu(Gmat)-diag(diag(Gmat)))'; %use Hermitian property to construct lower off-diagonal elements as conjugates of their upper counterparts
else %use 4-point derivative approximation method
    
    nsp_mat=zeros(dimy*nfreq,dimy*npar); %blank for derivatives (a-h)
    psp_mat=nsp_mat; %blank for derivatives (a+h)
    n2sp_mat=zeros(dimy*nfreq,dimy*npar); %blank for derivatives (a-2h)
    p2sp_mat=nsp_mat; %blank for derivatives (a+2h)
    
    for j=1:npar;
        
        %perturb one parameter at a time
        npara=theta0;
        ppara=theta0;
        n2para=theta0;
        p2para=theta0;
        npara(parind(j))=theta0(parind(j))-h;
        ppara(parind(j))=theta0(parind(j))+h;
        n2para(parind(j))=theta0(parind(j))-2*h;
        p2para(parind(j))=theta0(parind(j))+2*h;
        
        
        
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
        
        
        
        for i=1:nfreq;
            exe=exp(-1i*x(i));
            
            nmat1=(eye(nstate)-nT1*exe)\eye(nstate);
            nmat2=nmat1';
            
            n2mat1=(eye(nstate)-n2T1*exe)\eye(nstate);
            n2mat2=n2mat1';
            
            pmat1=(eye(nstate)-pT1*exe)\eye(nstate);
            pmat2=pmat1';
            
            p2mat1=(eye(nstate)-p2T1*exe)\eye(nstate);
            p2mat2=p2mat1';
            
            nsp_mat(((i-1)*dimy+1):(i*dimy),((j-1)*dimy+1):(j*dimy))=A*nmat1*nRR*nQQ*nRR'*nmat2*A'/(2*pi); %spectra at perturbed theta0
            n2sp_mat(((i-1)*dimy+1):(i*dimy),((j-1)*dimy+1):(j*dimy))=A*n2mat1*n2RR*n2QQ*n2RR'*n2mat2*A'/(2*pi); %spectra at perturbed theta0
            psp_mat(((i-1)*dimy+1):(i*dimy),((j-1)*dimy+1):(j*dimy))=A*pmat1*pRR*pQQ*pRR'*pmat2*A'/(2*pi); %spectra at perturbed theta0
            p2sp_mat(((i-1)*dimy+1):(i*dimy),((j-1)*dimy+1):(j*dimy))=A*p2mat1*p2RR*p2QQ*p2RR'*p2mat2*A'/(2*pi); %spectra at perturbed theta0
            
        end
    end
    
    %now construct the G matrix, note that the jkth component is simply the derivative with respect to j and k
    %a 4-point method is used
    
    Gmat=zeros(npar,npar); %blank
    
    for j=1:npar
        
        for k=j:npar
            aj1=(8*psp_mat(:,((j-1)*dimy+1):(j*dimy))-8*nsp_mat(:,((j-1)*dimy+1):(j*dimy))+n2sp_mat(:,((j-1)*dimy+1):(j*dimy))-p2sp_mat(:,((j-1)*dimy+1):(j*dimy)))/(12*h);
            aj2=conj(aj1);
            ak1=(8*psp_mat(:,((k-1)*dimy+1):(k*dimy))-8*nsp_mat(:,((k-1)*dimy+1):(k*dimy))+n2sp_mat(:,((k-1)*dimy+1):(k*dimy))-p2sp_mat(:,((k-1)*dimy+1):(k*dimy)))/(12*h);
            ak2=conj(ak1);
            for freq=1:nfreq;
                 Gmat(j,k)=Gmat(j,k)+wg(freq)*(sum(diag(aj1(((freq-1)*dimy+1):(freq*dimy),:)*ak1(((freq-1)*dimy+1):(freq*dimy),:)))+sum(diag(aj2(((freq-1)*dimy+1):(freq*dimy),:)*ak2(((freq-1)*dimy+1):(freq*dimy),:))));
                
            end
            
        end
        
    end
    Gmat=Gmat+(triu(Gmat)-diag(diag(Gmat)))'; %use Hermitian property to construct lower off-diagonal elements as conjugates of their upper counterparts
end
G=real(Gmat)*(1/4/pi);



end