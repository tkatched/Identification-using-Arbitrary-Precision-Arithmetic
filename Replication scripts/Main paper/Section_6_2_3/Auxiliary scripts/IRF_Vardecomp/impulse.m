
%% This function computes impulse responses on levels for a model 

function irf=impulse(T0,T1,A0,A1,order);
%Note, we only need to compute (A0+A1*L)inv(1-T1*L)*T0
%first create a matrix that contains powers of T1
n=rows(T1);
%ny=cols(T0);
PT1=zeros(n,n*(order+1));
PT1(1:n,1:n)=eye(n);
for i=1:order;
PT1(:,i*n+1:(i+1)*n)=PT1(:,(i-1)*n+1:i*n)*T1;
end

%now multiplying by T0
PT=zeros(n,cols(T0)*(order+1));

for i=1:order+1;
PT(:,(i-1)*cols(T0)+1:i*cols(T0))=PT1(:,(i-1)*n+1:i*n)*T0;
end
%finally multiply by A0
factor1=A0*PT;
%factor2=horzcat(zeros(cols(T0),cols(T0)),(A1*PT(:,1:cols(PT)-cols(T0))));
irf=factor1; %+factor2;
end