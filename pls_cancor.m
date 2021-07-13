function [pls,cancor]=pls_cancor(x,y,nperm,nboot,CI)
% [pls,cancor]=pls_cancor(x,y,nperm,nboot,CI)
%Input 
% x and y are the data assume x is the higher dimension matrix
% option arguments
% nperm = number of permutation resamplings
% nboot = number of bootstrap resamplings
% confidence interval - default 95
%Output
% pls.u, pls.v, pls.s are from the SVD of Rxy
% pls.perm contains matrix of all singular values for each permutation
% pls.boot.u pls.boot.us bootstrapped U and U*S
% pls.boot.v pls.boot.vs boostrapped V and V*S
% pls.boot.ul_u, pls.boot.ll_u upper and lower percentile of bootstrap dist for U
% pls.boot.ul_v, pls.boot.ll_v upper and lower percentile of bootstrap dist for V
% cancor.u, cancor.v, cancor.s are from the SVD of Rxy
% cancor.perm contains matrix of all singular values for each permutation
%
% cancor.boot.u cancor.boot.us bootstrapped U and U*S
% cancor.boot.v cancor.boot.vs boostrapped V and V*S
% cancor.boot.ul_u, cancor.boot.ll_u upper and lower percentile of bootstrap dist for U
% cancor.boot.ul_v, cancor.boot.ll_v upper and lower percentile of bootstrap dist for V
% cancor.a, cancor.b, cancor.sa, cancor.sb are standardize weights and structure coefficients
% 
% cancor.boot.a bootstrapped standardized weights
% cancor.boot.b bootstrapped standardized weights
% cancor.boot.ul_a, cancor.boot.ll_a upper and lower percentile of bootstrap dist
% cancor.boot.ul_b, cancor.boot.ll_b upper and lower percentile of bootstrap dist
%
% cancor.boot.sa bootstrapped structure weights
% cancor.boot.sb bootstrapped structure weights
% cancor.boot.ul_sa, cancor.boot.ll_sa upper and lower percentile of bootstrap dist
% cancor.boot.ul_sb, cancor.boot.ll_sb upper and lower percentile of bootstrap dist
%
% UPDATE - added checks for badconditioned data for CCA and for bootstrap
%
%  Written by ARMcIntosh, December 2020
%
% Dependencies: requires pls_cmd resampling code: rri_boot_order
% 


cancor.badcondition=0; %counter in case Rxx or Ryy is rank deficient

rxy=corr(x,y);
rxx=corr(x);
ryy=corr(y);
[pls.u,pls.s,pls.v]=svd(rxy,0);
if rank(rxx)==length(rxx) & rank(ryy)==length(ryy)
rx1=chol(rxx);
ry1=chol(ryy);
omega=inv(rx1)'*rxy*inv(ry1);
[cancor.u,cancor.s,cancor.v]=svd(omega,0);
cancor.a=inv(rx1)*cancor.u;  %standardized canonical weights
cancor.sa=rxx*cancor.a; %canonical structure coeffcients
cancor.b=inv(ry1)*cancor.v;  %standardized canonical weights
cancor.sb=ryy*cancor.b; %canonical structure coeffcients
else
    cancor.badcondition=cancor.badcondition+1;
    cancor.warning=['Unable to run cca because either X or Y is rank deficient'];
end

n=length(y);
[r,c]=size(x);
[r,c2]=size(y);




if nargin>2
  %perm_order=rri_perm_order(n,1,nperm);
cancor.badcondition_perm=0;
  disp("Permutation loop");
for i=1:nperm
xperm=x(randperm(n),:); %keep the permutation simple for now
%     for j=1:c
%         xperm(:,j)=x(randperm(n),j); %permute columns independently
%     end
%     for j=1:c2
%                 yperm(:,j)=y(randperm(n),j);
%     end
%     pls.xperm{i}=xperm;
%     pls.yperm{i}=yperm;
    rxPy=corr(xperm,y);
    rxP=corr(xperm);
    ryP=corr(y);
    pls.perms(:,i)=svd(rxPy);
    
    if cancor.badcondition==0 %prevents CCA loop if the X or Y are not full rank
   
    if rank(rxx)==length(rxP) 
        Pomega=inv(chol(rxP))'*rxPy*inv(chol(ryP));
        cancor.perms(:,i)=svd(Pomega);
    else
        cancor.badcondition_perm=cancor.badcondition_perm+1;
    end
    end
    
end
 %use rri_boot_order - part of plscmd release
    [boot_order,~]=rri_boot_order(n,1,nboot);
    disp("Bootstrap loop");
    pls.boot.u=zeros([size(pls.u),nboot]);
    pls.boot.v=zeros([size(pls.v),nboot]);
    pls.boot.us=zeros([size(pls.u),nboot]);
    pls.boot.vs=zeros([size(pls.v),nboot]);
    pls.boot.s=zeros([size(pls.s),nboot]);
    
    if cancor.badcondition==0 %prevents CCA loop if the X or Y are not full rank
    cancor.boot.u=zeros([size(cancor.u),nboot]);
    cancor.boot.v=zeros([size(cancor.v),nboot]);
    cancor.boot.us=zeros([size(cancor.u),nboot]);
    cancor.boot.vs=zeros([size(cancor.v),nboot]);
    cancor.boot.s=zeros([size(cancor.s),nboot]);
    end

    pls.boot.badboot=0;
for i=1:nboot
    %disp(i);
    xboot=x(boot_order(:,i),:);
    yboot=y(boot_order(:,i),:);
    %check for zero variance in bootstrap resample
    test_x=std(xboot);
    test_y=std(yboot);
    test_x=sum(test_x==0);
    test_y=sum(test_y==0);
    
    if test_x==0 & test_y==0
               
    rxBy=corr(xboot,yboot);
    rxB=corr(xboot);
    ryB=corr(yboot);
    
    [pls.boot.u(:,:,i),pls.boot.s(:,:,i),pls.boot.v(:,:,i)]=svd(rxBy,0);
    if cancor.badcondition==0 %prevents CCA loop if the X or Y are not full rank
        
    if rank(rxB)==length(rxB) & rank(ryB)==length(ryB)
    Bomega=inv(chol(rxB))'*rxBy*inv(chol(ryB));
    [cancor.boot.u(:,:,i),cancor.boot.s(:,:,i),cancor.boot.v(:,:,i)]=svd(Bomega,0);
    flip_idx=[];
    flip_idx=find(diag(cancor.boot.u(:,:,i)'*cancor.u)<0);
    if isempty(flip_idx)==0
        cancor.boot.u(:,flip_idx,i)=cancor.boot.u(:,flip_idx,i)*-1;
        cancor.boot.v(:,flip_idx,i)=cancor.boot.v(:,flip_idx,i)*-1;
    end
        cancor.boot.us(:,:,i)=cancor.boot.u(:,:,i)*cancor.boot.s(:,:,i);
        cancor.boot.vs(:,:,i)=cancor.boot.v(:,:,i)*cancor.boot.s(:,:,i);
        cancor.boot.a(:,:,i)=inv(chol(rxB))*cancor.boot.u(:,:,i);
        cancor.boot.sa(:,:,i)=rxB*cancor.boot.a(:,:,i);
        cancor.boot.b(:,:,i)=inv(chol(ryB))*cancor.boot.v(:,:,i);
        cancor.boot.sb(:,:,i)=ryB*cancor.boot.b(:,:,i);
    else
        cancor.badcondition_perm=cancor.badcondition_perm+1;
    end
    
    end %outter if
    
   
    %check for sign flips on a bootstrap iteration 
    flip_idx=[];
    flip_idx=find(diag(pls.boot.u(:,:,i)'*pls.u)<0);
    if isempty(flip_idx)==0
        pls.boot.u(:,flip_idx,i)=pls.boot.u(:,flip_idx,i)*-1;
        pls.boot.v(:,flip_idx,i)=pls.boot.v(:,flip_idx,i)*-1;
    end
    
    %save scaled singular vectors for CI calculations
    pls.boot.us(:,:,i)=pls.boot.u(:,:,i)*pls.boot.s(:,:,i);
    pls.boot.vs(:,:,i)=pls.boot.v(:,:,i)*pls.boot.s(:,:,i);
    
    else
        pls.boot.badboot=pls.boot.badboot+1;
    end

  
end

%compute percentile CI's for U and V
[ru,cu]=size(pls.u);
[rv,cv]=size(pls.v);

if nargin>3
    CI=95;
end

pls.boot.ul_u=zeros(ru,cu);
pls.boot.ll_u=zeros(ru,cu);
pls.boot.ul_v=zeros(rv,cv);
pls.boot.ll_v=zeros(rv,cv);

if cancor.badcondition==0 %prevents CCA loop if the X or Y are not full rank
cancor.boot.ul_u=zeros(ru,cu);
cancor.boot.ll_u=zeros(ru,cu);
cancor.boot.ul_v=zeros(rv,cv);
cancor.boot.ll_v=zeros(rv,cv);
cancor.boot.ul_a=zeros(ru,cu);
cancor.boot.ll_a=zeros(ru,cu);
cancor.boot.ul_b=zeros(rv,cv);
cancor.boot.ll_b=zeros(rv,cv);
cancor.boot.ul_sa=zeros(ru,cu);
cancor.boot.ll_sa=zeros(ru,cu);
cancor.boot.ul_sb=zeros(rv,cv);
cancor.boot.ll_sb=zeros(rv,cv);
end

%Loop to step through and get upper and lower bounds of bootstrap dist
for i=1:cu
    for j=1:ru
        pls.boot.ul_u(j,i)=prctile(pls.boot.us(j,i,:),CI);
        pls.boot.ll_u(j,i)=prctile(pls.boot.us(j,i,:),100-CI);
          if cancor.badcondition==0 %prevents CCA loop if the X or Y are not full rank
        cancor.boot.ul_u(j,i)=prctile(cancor.boot.us(j,i,:),CI);
        cancor.boot.ll_u(j,i)=prctile(cancor.boot.us(j,i,:),100-CI);
        cancor.boot.ul_a(j,i)=prctile(cancor.boot.a(j,i,:),CI);
        cancor.boot.ll_a(j,i)=prctile(cancor.boot.a(j,i,:),100-CI);
        cancor.boot.ul_sa(j,i)=prctile(cancor.boot.sa(j,i,:),CI);
        cancor.boot.ll_sa(j,i)=prctile(cancor.boot.sa(j,i,:),100-CI);
          end
    end
end

for i=1:cv
    for j=1:rv
        pls.boot.ul_v(j,i)=prctile(pls.boot.vs(j,i,:),CI);
        pls.boot.ll_v(j,i)=prctile(pls.boot.vs(j,i,:),100-CI);
          if cancor.badcondition==0 %prevents CCA loop if the X or Y are not full rank
        cancor.boot.ul_v(j,i)=prctile(cancor.boot.vs(j,i,:),CI);
        cancor.boot.ll_v(j,i)=prctile(cancor.boot.vs(j,i,:),100-CI);
        cancor.boot.ul_b(j,i)=prctile(cancor.boot.b(j,i,:),CI);
        cancor.boot.ll_b(j,i)=prctile(cancor.boot.b(j,i,:),100-CI);
        cancor.boot.ul_sb(j,i)=prctile(cancor.boot.sb(j,i,:),CI);
        cancor.boot.ll_sb(j,i)=prctile(cancor.boot.sb(j,i,:),100-CI);
          end
    end
end

end

