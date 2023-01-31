function V=variance_decomp(T0,T1,A0,A1,QQ,order);
%Note, data: (A0+A1*L)inv(1-T1*L)*T0
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
%finally multiply by A0+A1*L
factor1=A0*PT;
%size (zeros(cols(T0),cols(T0)))
factor2=A1*[zeros(n,cols(T0)),PT(:,1:cols(PT)-cols(T0))];


factor =factor1+factor2;

V = 0;
for i=1:order+1;

V = V+factor(:,(i-1)*cols(T0)+1:i*cols(T0))*QQ*factor(:,(i-1)*cols(T0)+1:i*cols(T0))';

end
end