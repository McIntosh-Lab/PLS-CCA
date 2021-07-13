function [pls_repro,cca_repro,splitflag]=split_half_PLSCCA_TrainTest(x,y,numsplits,null_flag)
% [pls_repro,cca_repro,splitflag]=split_half_PLSCCA_TrainTest(x,y,numsplits,null_flag)
% numsplits: number of split half samples
% null_flag: save distributions 1=true default=0
%
% Computes "expected" singular values from 1/2 of a sample left out from
% training sample.
%
% Output
% pls_repro.pls_s_train Distribution singular values from training samples
% pls_repro.pls_s_test  Distribution singular values from corresponding
% test samples
%
% cca_repro.cca_s_train Distribution singular values from training sample
% cca_repro.cca_s_test  Distribution singular values from corresponding
% test samples
%
% pls_repro.z Z-values for each test singular value
% (mean(s_test)/std(s_test))
%
% cca_repro.z Z-values for each test singular value
% (mean(s_test)/std(s_test))
%
% pls_repro.pls_s_train_null Distribution singular values from training
% sample where rows for one data block are permuted
%
% pls_repro.pls_s_test_null Distribution singular values from test
% sample where rows for one data block are permuted
%
% cca_repro.cca_s_train_null Distribution singular values from training
% sample where rows for one data block are permuted
%
% cca_repro.cca_s_test_null Distribution singular values from test
% sample where rows for one data block are permuted
%
% Dependencies pls_cancor
%
% Written ARMcIntosh Dec 2020
%

if nargin==3
    null_flag=0;
end

[n,p]=size(x);
nsplit=floor(n/2);
[~,q]=size(y);

d=min(p,q);


splitflag=[];

for i=1:numsplits
    %disp(i);
     idx=randperm(n);
     idx_1=idx(1:nsplit);
     idx_2=idx(ceil(n/2)+1:n);  %controls in case n is an odd number
     x1=x(idx_1,:);
     y1=y(idx_1,:);
     x2=x(idx_2,:);
     y2=y(idx_2,:);
     %test for zero variance in x or y
     test_std=[std(x1) std(x2) std(y1) std(y2)];

    
     if all((test_std)~=0)
         rxy=corr(x2,y2);
         rxx=corr(x2);
         ryy=corr(y2);
         if rank(rxx)==length(rxx) && rank(ryy)==length(ryy)
             rx1=chol(rxx);
             ry1=chol(ryy);
             omega=inv(rx1)'*rxy*inv(ry1);
         end
         [pls1,cca1]=pls_cancor(x1,y1);
         pls_s_train(:,:,i)=pls1.s;
                 
         pls_s_test(:,:,i)=pls1.u'*rxy*pls1.v;
         if isfield(cca1,'u') && rank(rxx)==length(rxx) && rank(ryy)==length(ryy)
             cca_s_train(:,:,i)=cca1.s;
             cca_s_test(:,:,i)=cca1.u'*omega*cca1.v;
         else
             cca_s_train(:,:,i)=NaN(q,q);
             cca_s_test(:,:,i)=NaN(q,q);

         end
    else
        splitflag=[splitflag,i];
         pls_s_train(:,:,i)=NaN(q,q);
         pls_s_test(:,:,i)=NaN(q,q);
         cca_s_train(:,:,i)=NaN(q,q);
         cca_s_test(:,:,i)=NaN(q,q);
    end
     
    
    

end



if null_flag==1
for i=1:numsplits %create a null distribution - scramble X-Y
    %disp(i);
     idx=randperm(n);
     idx_1=idx(1:nsplit);
     idx_2=idx(ceil(n/2)+1:n);
     x1=x(idx_1,:);
     y1=y(idx_2,:);
     x2=x(idx_2,:);
     y2=y(idx_1,:);
     
     %test for zero variance in x or y
     test_std=[std(x1) std(x2) std(y1) std(y2)];

    
    if all((test_std)~=0)
        rxy=corr(x2,y2);
         rxx=corr(x2);
         ryy=corr(y2);
         if rank(rxx)==length(rxx) && rank(ryy)==length(ryy)
             rx1=chol(rxx);
             ry1=chol(ryy);
             omega=inv(rx1)'*rxy*inv(ry1);
         end
         [pls1_null,cca1_null]=pls_cancor(x1,y1);
         pls_s_train_null(:,:,i)=pls1_null.s;
                 
         pls_s_test_null(:,:,i)=pls1_null.u'*rxy*pls1_null.v;
         
         if isfield(cca1_null,'u') 
             cca_s_train_null(:,:,i)=cca1_null.s;
             cca_s_test_null(:,:,i)=cca1_null.u'*omega*cca1_null.v;

         end
    end

end

end

pls_repro.pls_s_train=pls_s_train;
pls_repro.pls_s_test=pls_s_test;
cca_repro.cca_s_train=cca_s_train;
cca_repro.cca_s_test=cca_s_test;

for i=1:d
    pls_repro.z(i)=mean(pls_repro.pls_s_test(i,i,:),'omitnan')/std(pls_repro.pls_s_test(i,i,:),'omitnan');
    cca_repro.z(i)=mean(cca_repro.cca_s_test(i,i,:),'omitnan')/std(cca_repro.cca_s_test(i,i,:),'omitnan');
end

if null_flag==1
pls_repro.pls_s_train_null=pls_s_train_null;
pls_repro.pls_s_test_null=pls_s_test_null;
cca_repro.cca_s_train_null=cca_s_train_null;
cca_repro.cca_s_test_null=cca_s_test_null;
end

