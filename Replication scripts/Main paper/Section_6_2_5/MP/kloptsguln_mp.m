function X = kloptsguln_mp(thetainput,thetaa,truespec,msel,x,wg,bc,parsel,con,wgt,nrm,indp,lb,ub)

%Function to compute the KL distance between Leeper(1991) models

% Inputs:

%thetainput - candidate solution vector

%w - vector of frequencies between [-pi,pi]

%bc=0 corresponds to the full spectrum case; when bc=1 the function uses only the business cycle frequencies

%parsel - parameter selection index.

% con - constraint evaluation index. If con = 1, the function will evaluate
% the constraint given by constraintlp.m first.

%% Misc. setup
warning('off','all');

thetat=transfmpsgu(thetainput,2,lb,ub); %transform parameters from unconstrained parameterization

n=length(x); %number of abscissae supplied


%% Insert parameters specified by parsel into benchmark value

theta0=thetaa; %benchmark parameter values

theta0(parsel)=thetat; %replace parameters that are varying in optimization




%% Solve the model
%Specify the monetary policy rule using msel

ny=mp('7'); %number of observables

[gx,hx,eta,exitflag,npara] = sgusolvn_mp(theta0,msel);


if exitflag ~=1
    X=mp('1e12');
    return
end


neq=size(hx,2); %number of equations

A=mp(zeros(ny,22)); %selection matrix A

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


tempeye=mp(eye(neq)); %auxiliary identity matrix

gx2=[gx;mp(zeros(2,neq))]; %add selection of levels of g and z
gx2(21,5)=onep;
gx2(22,6)=onep;
%% Compute KL divergence from the benchmark spectrum


%preparations to compute spectrum
X=mp('0');
if bc==0 && mod(n,2)==1 %if not even, keep frequency 0 as well
    nx =length(x)-1;
else
    nx = length(x); %number of frequencies used
end

for i=1:nx %compute spectrum
    exe=exp(ci*x(i));
    
    mat1=(tempeye-hx*exe)\tempeye;
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    spec=(A*gx2*mat1*eta*eta'*mat2*gx2'*A')/(pi2); %compute spectrum using representation in the paper
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+twop*(trace(temp) - log(det(temp))-ny)*wg(i);
    
end

if bc==0 && mod(n,2)==1 %if not even, compute at frequency 0 separately
    i=i+1;
    exe=exp(ci*x(i));
    
    mat1=(tempeye-hx*exe)\tempeye;
    mat2=mat1';  %note that ' gives the conjugate transpose
    
    spec=(A*gx2*mat1*eta*eta'*mat2*gx2'*A')/(pi2); %compute spectrum using representation in the paper
    
    
    temp=spec\truespec(((i-1)*ny+1):i*ny,:);
    X=X+(trace(temp) - log(det(temp))-ny)*wg(i);
end

X=real(X)*mp('1/(4*pi)'); %final answer: KL. real() necessary since
%small complex residual of order e-15i sometimes remains


if isnan(X) | isinf(X) | isempty(X) | X<0
    
    X=mp('1e12');
    return
end


end