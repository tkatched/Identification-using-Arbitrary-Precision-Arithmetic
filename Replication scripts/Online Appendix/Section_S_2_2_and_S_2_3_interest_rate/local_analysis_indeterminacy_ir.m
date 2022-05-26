%% Run the local analysis on points from the indeterminacy region

% Load simulated values
load local_ind theta_ind

%% Gaussian quadrature nodes and weights

n=500; %number of nodes
nmp=mp('500');

[x,wg] = quadcomp(n,0);
[xbc,wgbc] = quadcomp(n,1);

[xmp,wgmp] = quadcomp(nmp,0);
[xmpbc,wgmpbc] = quadcomp(nmp,1);

%keep only one side of the quadrature abscissae and weights due to symmetry
%about zero
x=x(1:n/2);
xbc=xbc(1:n/2);
xmp=xmp(1:n/2);
xmpbc=xmpbc(1:n/2);
wg=wg(1:n/2);
wgbc=wgbc(1:n/2);
wgmp=wgmp(1:n/2);
wgmpbc=wgmpbc(1:n/2);

nmat=1000; %number of points considered

nthetaind=length(theta_ind(1,:)); %no of parameters
msel=1; %model selection: 1=MA(1) shocks, 2=MA(2) shocks, 3=ARMA(1,1) shocks

dif=1; %numerical derivative rule (1=symmetric quotient,2=4-point rule)

filename=['local_ind_ir_dif',num2str(dif)]; %name for result file

%Numerical differentiation step sizes for double and arbitrary precision:
hd=1e-07;
hmp=mp('1e-10');

%% Preallocate blanks

G_ind=zeros(nthetaind,nthetaind,nmat);
Gmp_ind=mp(zeros(nthetaind,nthetaind,nmat));

G_indbc=zeros(nthetaind,nthetaind,nmat);
Gmp_indbc=mp(zeros(nthetaind,nthetaind,nmat));

valind=zeros(nthetaind,nmat);
seigind=zeros(nmat,1);
tolind=seigind;
rankind=seigind;

valindbc=zeros(nthetaind,nmat);
seigindbc=zeros(nmat,1);
tolindbc=seigind;
rankindbc=seigind;

valind_mp=mp(zeros(nthetaind,nmat));
seigind_mp=mp(zeros(nmat,1));
tolind_mp=seigind_mp;
rankind_mp=seigind_mp;

valindbc_mp=mp(zeros(nthetaind,nmat));
seigindbc_mp=mp(zeros(nmat,1));
tolindbc_mp=seigind_mp;
rankindbc_mp=seigind_mp;


tind=seigind;
tindbc=seigind;


tind_mp=seigind_mp;
tindbc_mp=seigind_mp;


%% Run the local analysis in double and quadruple precision
for k=1:nmat
    if mod(k,10)==0
        save(filename)
        disp(k)
    end
    

    thetaind=theta_ind(k,:);

    thetaind_mp=mp(thetaind);
    
    
    tic;
    Gind=gmatlpv_ir(thetaind,x,wg,hd,dif,msel);
    tind(k,1)=toc;
    
    valind(:,k)=eig(Gind);
    seigind(k,1)=min(abs(valind(:,k)));
    s = svd(Gind);
    tolind(k,1) = max(size(Gind))*eps(max(s));
    rankind(k,1)=rank(Gind);
    
    tic;
    Gindbc=gmatlpv_ir(thetaind,xbc,wgbc,hd,dif,msel);
    tindbc(k,1)=toc;
    
    valindbc(:,k)=eig(Gindbc);
    seigindbc(k,1)=min(abs(valindbc(:,k)));
    s = svd(Gindbc);
    tolindbc(k,1) = max(size(Gindbc))*eps(max(s));
    rankindbc(k,1)=rank(Gindbc);
    
    
    %Multiple precision computations:
    
    tic;
    Gind_mp=gmatlp_mpv_ir(thetaind_mp,xmp,wgmp,hmp,dif,msel);
    tind_mp(k,1)=toc;
    
    valind_mp(:,k)=eig(Gind_mp);
    seigind_mp(k,1)=min(abs(valind_mp(:,k)));
    s_mp = svd(Gind_mp);
    tolind_mp(k,1) = max(size(Gind_mp))*eps(max(s_mp));
    rankind_mp(k,1)=rank(Gind_mp);
    
    tic;
    Gindbc_mp=gmatlp_mpv_ir(thetaind_mp,xmpbc,wgmpbc,hmp,dif,msel);
    tindbc_mp(k,1)=toc;
    
    valindbc_mp(:,k)=eig(Gindbc_mp);
    seigindbc_mp(k,1)=min(abs(valindbc_mp(:,k)));
    s_mp = svd(Gindbc_mp);
    tolindbc_mp(k,1) = max(size(Gindbc_mp))*eps(max(s_mp));
    rankindbc_mp(k,1)=rank(Gindbc_mp);
    
    

    
    
    
    G_ind(:,:,k)=Gind;
    G_indbc1(:,:,k)=Gindbc;
    Gmp_ind(:,:,k)=Gind_mp;
    Gmp_indbc(:,:,k)=Gindbc_mp;
    
end
disp 'Smallest eigenvalue range in double precision:'
disp([min(seigind),max(seigind)])

disp 'Smallest eigenvalue range in double precision, BC frequencies:'
disp([min(seigindbc),max(seigindbc)])

disp 'Smallest eigenvalue range in quadruple precision:'
disp([min(seigind_mp),max(seigind_mp)])

disp 'Smallest eigenvalue range in quadruple precision, BC frequencies:'
disp([min(seigindbc_mp),max(seigindbc_mp)])

save(filename)
