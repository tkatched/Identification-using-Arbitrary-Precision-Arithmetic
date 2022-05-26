function X = kloptlpr2_mp(thetainput,truespec,x,wg,bc,areg,lb,ub)

%Function to compute the KL distance between Leeper et al (2017) model in
%multiple precision


% Inputs:

%thetainput - candidate solution vector

%x - vector of frequencies between [-pi,pi]

%wg - vector of Gaussian quadrature weights corresponding to x

%bc=0 corresponds to the full spectrum case; when bc=1 the function uses only the business cycle frequencies


% areg - denotes the alternative regime. If regime = 1, alternate regime is
% AMPF, if areg = 2, alt regime is PMAF, if areg=3, alt regime is PMPF.


%truespec - base regime spectrum.

%% Misc. setup
warning('off','all');

n=mp(length(x)); %number of abscissae supplied


theta0=transfmplpr(thetainput,2,lb,ub); %transform parameters from the unconstrained parameterization


%% Solve the model
%Specify the monetary policy rule using msel

% Define variables from DSGE model
variables;

% Set observable positions
obs = zeros(8,1);
obs(1,1) = Ncobs; % 34, consumption, 100 times log difference
obs(2,1) = Niobs; % 35, investment,100 times log difference
obs(3,1) = Nwobs; % 37, wage, 100 times log difference
obs(4,1) = Ngcobs;% 36, gov spending,100 times log difference
obs(5,1) = Nbobs; %38, debt,100 times log difference
obs(6,1) = NLobs; %41, hours worked, level
obs(7,1) = NPiobs; %40, inflation, level
obs(8,1) = NRobs; % 39, interest,level

[g0,g1,CC,Psi,Pie,H,C,E,ss,nvar,nexog] = model_mp(theta0,obs,areg);

[G,TC,M,TY,fmat,TZ,TETA,GEV,eu]=gensys_mp2(g0, g1, CC, Psi, Pie, 1); %use Sims-based code, also uses ordqz

if isempty(G)
    disp('empty solution: nonexistence/numerical issue')
    X=mp('1e12');
    return
end

id1=mp(eye(nvar)); %state vector dimension

dck=rcond(id1-G);

if dck<mp('1e-10')
    X=mp('1e12');
    return
end


M=M(:,nvar-nexog+1:end); % need to select the part that multiply the structural shocks


if eu(1,1)==1 && eu(2,1)==1 && areg<3 % searching determinacy
    
    Q = M*E*E'*M';% this is the error covariance matrix
    
elseif eu(1,1)==1 && eu(2,1)==0 && areg<3 %indeterminacy when searching determinacy region
    X=mp('1e12');
    return
elseif eu(1,1)==1 && eu(2,1)==0 && areg==3 %searching indeterminacy
    %shocks written as sunshock_t=proj_coef*et+sigma*eta_t; eta_t N(0,1)
    TETA=rref(TETA').'; %reduced column echelon form for TETA
    MM=[M,TETA];
    rotmat1 = [mp(eye(8)),mp(zeros(8,1))];
    rotmat = [rotmat1; theta0(end-8:end)]; %last nine elements are projection coefs followed by standard deviation
    Q1 = [E*E',mp(zeros(8,1));mp(zeros(1,8)),mp('1')]; %structural shock variance and 1
    Q = MM*rotmat*Q1*rotmat'*MM';
elseif eu(1,1)==1 && eu(2,1)==1 && areg==3 %%determinacy when searching indeterminacy region
    X=mp('1e12');
    return
else%no equilibrium exists/numerical problems.
    X=mp('1e12');
    return
end
  



        

cc=mp('2*pi');
sqi=mp('-1i');

ny=mp('8'); % obs dimension
%% Compute KL divergence from the benchmark spectrum


    %preparations to compute spectrum
    X=mp('0');
    if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
        nx =length(x)-1;
    else
        nx = length(x); %number of frequencies used
    end
    
    for i=1:nx %compute spectrum
        exe=exp(sqi*x(i));
        
        mat1=(id1-G*exe)\id1;
        mat2=mat1';  %note that ' gives the conjugate transpose
        spec=H*mat1*Q*mat2*H'/cc; %compute spectrum using representation in the paper
        
        temp=spec\truespec(((i-1)*ny+1):i*ny,:);
        X=X+mp('2')*(trace(temp) - log(det(temp))-ny)*wg(i);
        
    end
    
    if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
        i=i+1;
        exe=exp(sqi*x(i));
        
        mat1=(id1-G*exe)\id1;
        mat2=mat1';  %note that ' gives the conjugate transpose
        spec=H*mat1*Q*mat2*H'/cc; %compute spectrum using representation in the paper
        
        
        temp=spec\truespec(((i-1)*ny+1):i*ny,:);
        X=X+(trace(temp) - log(det(temp))-ny)*wg(i);
    end
    
    X=real(X)*mp('1/4/pi'); %final answer: KL. real() necessary since
    %small complex residual of order e-15i sometimes remains
    

% if X > -1e-06 %deal with tiny negative values due to approximation when KL close to zero
%     X=abs(X);
% end

if isnan(X) | isinf(X) | isempty(X) | X<0 
    
    X=mp('1e12');
    return
end


end