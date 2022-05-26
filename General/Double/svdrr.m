function  [Urr,Drr,Vrr,V_null,R]=svdrr(A)

%Translation of Lubik and Schorfheide (2004) Gauss routine svdrr
realsmall = 1E-6;

[U, D, V] = svd(A);
V = V';
temp=abs(diag(D));
s=size(temp);
for i=1:s(1)
    for j=1:s(2)
        if temp(i,j) > realsmall
            temp(i,j)=1;
        end
    end
end
R = sum(temp);
if R == 0
   Urr = U(:,1);
   Drr = D(1,1);
   Vrr = V(1,:);
else
   Urr = U(:,1:R);
   Drr = D(1:R,1:R);
   Vrr = V(1:R,:);
end

if R == size(V,1);
   V_null = zeros(1,size(V,2));
else
   V_null = V(R+1:size(V,1),:);
end

Vrr=Vrr';
V_null = V_null';

end