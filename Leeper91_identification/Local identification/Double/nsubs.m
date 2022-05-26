function [SBS SE SBN] = nsubs(G,num,tol,parvec)

%%Function for analysis of nonidentified parameter subsets of the criterion
%%matrix G(theta).

%The code is precision independent - if G, tol provided are mp object, all
%relevant computation automatically proceeds in MP with the Advanpix
%toolbox installed.

%Inputs: G - the criterion matrix
% num - consider subset of up to 'num' elements, must be <= no of
% parameters in the model

% tol (optional) - tolerance level for determining rank of submatrices corresponding
% to subsets. If not provided, the default Matlab rank tolerance for
% G(theta) is used.
% For arbitrary precision computation, it is recommended to supply tolerance
% as the default computation may produce a number that is too small to be reasonable.


% parvec (optional) - a column vector containing parameter names. Used to
% display parameter names of the nonidentified subsets in addition to their
% indices.
;


ntheta=size(G,1); %no of parameters

if num>ntheta %check max. subsert size doesn't exceed parameter dimension
    disp 'Error: maximum subset size exceeds parameter dimension'
    return
end

%Check if the tolerance level supplied, if not calculate default:
if ~exist('tol','var') || isempty(tol)
    s=svd(G);
    tol = ntheta*eps(max(s)); %default tolerance
end


%create structures to store results
index=struct;
eigv=struct;
seigv=struct;
grank=struct;
Gs=struct;
elimpool=struct;
%%

for kk=1:num
    index.(['ind',num2str(kk)])=nchoosek([1:1:ntheta]',kk); %all possible subsets of kk elements
    
    if kk>1 %remove subsets containing those previously detected as nonidentified as proper subsets
        for ii=1:kk-1
            ps=size(elimpool.(['ep',num2str(ii)]),1);
            for m=1:ps
            temp1=ismember(index.(['ind',num2str(kk)]),elimpool.(['ep',num2str(ii)])(m,:));
            temp2=sum(temp1');
            temp3=find(temp2==ii);
            index.(['ind',num2str(kk)])(temp3,:)=[];
            end
        end
    end
    
    ncombs=size(index.(['ind',num2str(kk)]),1); %no of feasible combinations of kk-element subsets
    Gs.(['G',num2str(kk)])=zeros(kk,kk,ncombs); %blanks for kk-element subset G matrices
    grank.(['rk',num2str(kk)])=zeros(ncombs,1); %blanks for storing ranks based on tol
    eigv.(['ev',num2str(kk)])=zeros(kk,ncombs); %blanks for storing eigenvalues
    for k=1:ncombs
        for i=1:kk
            for j=1:kk %assemble G-submatrices for each subset
                Gs.(['G',num2str(kk)])(i,j,k)=G(index.(['ind',num2str(kk)])(k,i),index.(['ind',num2str(kk)])(k,j));
            end
        end
        grank.(['rk',num2str(kk)])(k,1)=rank(Gs.(['G',num2str(kk)])(:,:,k),tol); %check rank using tol
        eigv.(['ev',num2str(kk)])(:,k)=eig(Gs.(['G',num2str(kk)])(:,:,k)); % save eigenvalues
        
    end
    seigv.(['sev',num2str(kk)])=min(abs(eigv.(['ev',num2str(kk)])),[],1); %save smallest eigenvalues
    temp4=find(seigv.(['sev',num2str(kk)])<tol); %identify subsets with one zero eigenvalue according to tol
    elimpool.(['ep',num2str(kk)])=index.(['ind',num2str(kk)])(temp4,:); %save the detected nonidentified subset (as parameter indices)
    neigv.(['nev',num2str(kk)])=seigv.(['sev',num2str(kk)])(temp4); %save the corresponding smallest eigenvalues
end
%%
c=0; %initialize counter for detected subsets
for v=1:num
    if ~isempty(elimpool.(['ep',num2str(v)]))
        c=c+1;
        SBS.(['nsub',num2str(c)])=elimpool.(['ep',num2str(v)]); %save nonidentified subset indices if any
        SE.(['nev',num2str(c)])=neigv.(['nev',num2str(v)]); %save corresponding smallest eigenvalues
    end
end

if exist('parvec','var') %if parvec supplied, match parameter names to detected subsets in
    for vv=1:c
        temp5=size(SBS.(['nsub',num2str(vv)]),1);
        for w=1:temp5
        SBN.(['nsub',num2str(vv),'_',num2str(w)])=parvec(SBS.(['nsub',num2str(vv)])(w,:),:);
        end
    end
else
    SBN=[];
end

%Make sure all outputs are allocated even if there are no nonidentified
%subsets:
if exist('SBN')==0
    SBN=[];
end

if exist('SE')==0
    SE=[];
end

if exist('SBS')==0
    SBS=[];
end
end
