function [G1,C,impact,ywt,fmat,fwt,loose,gev,eu]=gensys_mp2(g0,g1,c,psi,pi,div) 

%This is a modification of the matlab program of Sims to use mp 
%computation - input matrices should be mp objects.

%It also allows for indeterminacy following the development 
%in Lubik and Schorfheide (2003,2004).

%(Sims' original program is availabe at http://sims.princeton.edu/yftp/gensys/mfiles/).

%It uses the matlab commands qz and ordqz to obtain the generalized Schur
%decomposition

% System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Christopher A. Sims
% Corrected 10/28/96 by CAS

eu=[0;0];
veta2=[];
realsmall=mp('1e-6');
fixdiv=(nargin==6);
n=size(g0,1);
[a b q z v]=qz(g0,g1); %initial QZ decomposition
if ~fixdiv, div=mp('1.01'); end
nunstab=0;
zxz=0;

for i=1:n
% ------------------div calculation------------
   if ~fixdiv
      if abs(a(i,i)) > 0
         divhat=abs(b(i,i))/abs(a(i,i));
         if 1+realsmall<divhat & divhat<=div
            div=.5*(1+divhat);
         end
      end
   end
% ----------------------------------------
   nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
   if abs(a(i,i))<realsmall & abs(b(i,i))<realsmall
      zxz=1;
   end
end
div ;
nunstab;
if ~zxz
    
  ordering = abs(diag(b)) <= (abs(diag(a))*div); %produce logical vector with 0's where unstable roots are  and
  %1's otherwise
  [a,b,q,z] = ordqz(a,b,q,z,ordering); %use Matlab ordqz routine to move the unstable eigenvalue to bottom right corner.

end
gev=[diag(a) diag(b)]; %for computing generalized eigenvalues
if zxz
   disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
   eu=[-2;-2];
   G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];loose=[]; %added empty loose here
   return
end
q1=q(1:n-nunstab,:); %Q1.
q2=q(n-nunstab+1:n,:); %Q2.
z1=z(:,1:n-nunstab)'; %Z1
z2=z(:,n-nunstab+1:n)'; %Z2
a2=a(n-nunstab+1:n,n-nunstab+1:n); %Lambda22
b2=b(n-nunstab+1:n,n-nunstab+1:n); %Omega22
etawt=q2*pi; %Q2.*Pi

neta = size(pi,2); %no of forecast errors
if nunstab == 0
  etawt == zeros(0,neta);
  ueta = zeros(0,0);
  deta = zeros(0,0);
  veta = zeros(neta,0);
  bigev = 0;
else
    %Sims' code - no V_vull (V.2) computed 
  [ueta,deta,veta]=svd(etawt); %singular value decomposition of Q.*Pi into U,D,V
  md=min(size(deta));
  bigev=find(diag(deta(1:md,1:md))>realsmall);
  ueta=ueta(:,bigev); %U.1
  
  veta=veta(:,bigev); %V.1
  deta=deta(bigev,bigev); %D11

%Schorfheide's code translation to get U.1, D11, V.1, V.2
[ueta,deta,veta,veta2,R]=svdrr_mp(etawt);
end

eu(1) = length(bigev)>=nunstab;


if nunstab == n
  etawt1 = zeros(0,neta);
  bigev =0;
  ueta1 = zeros(0, 0);
  veta1 = zeros(neta,0);
  deta1 = zeros(0,0);

else
  etawt1 = q1 * pi;
  ndeta1 = min(n-nunstab,neta);
  [ueta1,deta1,veta1]=svd(etawt1);
  md=min(size(deta1));
  bigev=find(diag(deta1(1:md,1:md))>realsmall);
  ueta1=ueta1(:,bigev);
  veta1=veta1(:,bigev);
  deta1=deta1(bigev,bigev);
end

if isempty(veta1)
	unique=1;
else
	loose = veta1-veta*veta'*veta1;
	[ul,dl,vl] = svd(loose);
	nloose = sum(abs(diag(dl)) > realsmall*n);
	unique = (nloose == 0);
end
if unique
   %disp('solution unique');
   eu(2)=1;
else
   %fprintf(1,'Indeterminacy.  %d loose endog errors.\n',nloose);

end

tempeye1=mp(eye(n-nunstab));
tempeye2=mp(eye(nunstab));
tempeye3=mp(eye(n));

tempz1=mp(zeros(nunstab,n-nunstab));
tempz2=mp(zeros(nunstab,n));

tmat = [tempeye1 -(ueta*(deta\veta')*veta1*deta1*ueta1')'];
G0= [tmat*a; tempz1 tempeye2];
G1= [tmat*b; tempz2];

G0I=G0\tempeye3;

%check if G0 is singular

[warnmsg, msgid] = lastwarn; %save the last warning (will be 'MATLAB:singularMatrix'
%if there is a singularity problem)

if strcmp(msgid,'MATLAB:singularMatrix')
 
    eu=[-3;-3];
  lastwarn(''); %reset the warning
  disp('G0 singular');
  G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];loose=[]; %added empty loose here
  return
end

G1=G0I*G1;
usix=n-nunstab+1:n;
C=G0I*[tmat*q*c;(a(usix,usix)-b(usix,usix))\q2*c];
impact=G0I*[tmat*q*psi;mp(zeros(nunstab,size(psi,2)))];
fmat=b(usix,usix)\a(usix,usix);
fwt=-b(usix,usix)\q2*psi;
ywt=G0I(:,usix);

if eu==[1;1]
loose = G0I * [etawt1 * (mp(eye(neta)) - veta * veta');mp(zeros(nunstab, neta))]; % TETA under determinacy - Sims
elseif eu(1)==1 && eu(2)==0 && isempty(veta2)==0
loose=G0I * [etawt1 * (mp(eye(neta)) - veta * veta')*veta2;mp(zeros(nunstab, size(veta2,2)))]; %TETA augmented for indeterminacy
elseif eu(1)==1 && eu(2)==0 && isempty(veta2)==1
   %disp('problematic parameter value(s)')
   %G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];loose=[]; %modification to consider only [1;1] and [1;0] return codes
   %return
   loose=G0I * [etawt1 * (mp(eye(neta)) - veta * veta');mp(zeros(nunstab, size(veta2,2)))]; %change to account for indeterminacy with 0 unstable roots
else
disp('existence with only non s.c. z')
   G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];loose=[]; %modification to consider only [1;1] and [1;0] return codes
   return    
end
% -------------------- above are output for system in terms of z'y -------
G1=real(z*G1*z');
C=real(z*C);
impact=real(z*impact);
loose = real(z * loose);


ywt=z*ywt;
end