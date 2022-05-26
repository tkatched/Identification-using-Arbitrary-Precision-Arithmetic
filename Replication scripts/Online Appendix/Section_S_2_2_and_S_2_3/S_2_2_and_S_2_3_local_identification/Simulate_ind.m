%% This file draws parameters from the indeterminacy region of the Leeper (1991) model

%% Draw parameters from the indeterminacy region for a "generic" local identification checks
rng(0)

nmat=1000; %number of draws

lbind=[0.01 0.9  1.01 0.005 0.005 -0.99 -0.99 -3 -3 0.005];
ubind=[0.99 0.999 2  1   1      0.99  0.99    3   3  1];

nthetaind=length(lbind);

%Preallocate blanks
theta_ind=zeros(nmat,nthetaind);


for k=1:nmat
    
    thetaind=zeros(nthetaind,1);

    for i=1:nthetaind
        thetaind(i)=(ubind(i)-lbind(i))*rand+lbind(i);
    end
    
    theta_ind(k,:)=thetaind;
end

save local_ind