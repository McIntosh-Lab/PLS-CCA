function [pls_repro,cca_repro,splitflag]=split_half_PLSCCA(x,y,numsplits,lv,CI,dist_flag)
% [pls_repro,cca_repro,splitflag]=split_half_PLSCCA(x,y,numsplits,lv,CI,dist_flag)
%
% computes the cosines between singular vectors U and V from numsplit
% random splits of X and Y data
%
% nsplits: number of split half samples
% lv: largest number of LV to be evaluated - e.g. lv=3 means 1,2,3 are
% assesed
% CI: confidence interval percentile
% dist_flag: save distributions 1=true default=0
%
%OUTPUT
% pls_repro.pls_rep_mean_u  average of cosines for u distribution from split-half
% pls_repro.pls_rep_mean_v average of cosines for v distribution from split-half
% pls_repro.pls_rep_z_u  Z-value for v distribition (mean_u/std_u)
% pls_repro.pls_rep_z_v  Z-value for v distribition (mean_v/std_v)
% pls_repro.pls_rep_ul_u=pls_rep_ul_u upper bound of u distribution  
% pls_repro.pls_rep_ll_u=pls_rep_ll_u lower bound of u distribution 
% pls_repro.pls_rep_ul_v=pls_rep_ul_v upper bound of v distribution 
% pls_repro.pls_rep_ll_v=pls_rep_ll_v upper bound of v distribution  
% pls_repro.pls_null_mean_u=pls_null_mean_u average of null u distribution
% created by permutation 
% pls_repro.pls_null_mean_v=pls_null_mean_v average of null v distribution
% pls_repro.pls_null_z_u=pls_null_u_z Z-value for null u distribition 
% pls_repro.pls_null_z_v=pls_null_v_z Z-value for null v distribition 
% pls_repro.pls_null_ul_u=pls_null_ul_u upper bound of null u distribution
% pls_repro.pls_null_ll_u=pls_null_ll_u; lower bound of null u distribution
% pls_repro.pls_null_ul_v=pls_null_ul_v; upper bound of null v distribution
% pls_repro.pls_null_ll_v=pls_null_ll_v lower bound of null u distribution
% pls_repro.pls_dist_u  full distribution of u cosines if dist_flag=1
% pls_repro.pls_dist_v  full distribution of v cosines if dist_flag=1
% pls_repro.pls_dist_null_u full distribution of null_u cosines if dist_flag=1
% pls_repro.pls_dist_null_v full distribution of null_v cosines if dist_flag=1
%
% This is repeated for CCA
% cca_repro.cca_rep_mean_u
% cca_repro.cca_rep_mean_v
% cca_repro.cca_rep_z_u
% cca_repro.cca_rep_z_v
% cca_repro.cca_rep_ul_u
% cca_repro.cca_rep_ll_u
% cca_repro.cca_rep_ul_v 
% cca_repro.cca_rep_ll_v 
% cca_repro.cca_null_mean_u
% cca_repro.cca_null_mean_v
% cca_repro.cca_null_z_u
% cca_repro.cca_null_z_v
% cca_repro.cca_null_ul_u 
% cca_repro.cca_null_ll_u
% cca_repro.cca_null_ul_v 
% cca_repro.cca_null_ll_v
% %structure coefficients
% cca_repro.cca_rep_mean_sa
% cca_repro.cca_rep_mean_sb
% cca_repro.cca_rep_z_sa
% cca_repro.cca_rep_z_sb
% cca_repro.cca_rep_ul_sa 
% cca_repro.cca_rep_ll_sa
% cca_repro.cca_rep_ul_sb 
% cca_repro.cca_rep_ll_sb 
% cca_repro.cca_null_mean_sa
% cca_repro.cca_null_mean_sb
% cca_repro.cca_null_z_sa
% cca_repro.cca_null_z_sb;
% cca_repro.cca_null_ul_sa 
% cca_repro.cca_null_ll_sa
% cca_repro.cca_null_ul_sb 
% cca_repro.cca_null_ll_sb
% cca_repro.cca_dist_sa
% cca_repro.cca_dist_sb
% cca_repro.cca_dist_null_sa
% cca_repro.cca_dist_null_sb
%
% Dependencies: calls pls_cancor
%
% Written ARMcIntosh December 2020


if nargin==5
    dist_flag=0;
end

[n,p]=size(x);
nsplit=floor(n/2);
[~,q]=size(y);

d=min(p,q);

pls_u_repro=zeros(d,d,numsplits);
pls_v_repro=zeros(d,d,numsplits);
cca_u_repro=zeros(d,d,numsplits);
cca_v_repro=zeros(d,d,numsplits);
cca_sa_repro=zeros(d,d,numsplits);
cca_sb_repro=zeros(d,d,numsplits);
pls_u_null=zeros(d,d,numsplits);
pls_v_null=zeros(d,d,numsplits);
cca_u_null=zeros(d,d,numsplits);
cca_v_null=zeros(d,d,numsplits);
cca_sa_null=zeros(d,d,numsplits);
cca_sb_null=zeros(d,d,numsplits);
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
     [pls1,cca1]=pls_cancor(x1,y1);
     [pls2,cca2]=pls_cancor(x2,y2);
     pls_u_repro(:,:,i)=pls1.u'*pls2.u;
     pls_v_repro(:,:,i)=pls1.v'*pls2.v;
     if isfield(cca1,'u') && isfield(cca2,'u')
     cca_u_repro(:,:,i)=cca1.u'*cca2.u;
     cca_v_repro(:,:,i)=cca1.v'*cca2.v;
     cca_sa_repro(:,:,i)=normalize(cca1.sa')*normalize(cca2.sa);
     cca_sb_repro(:,:,i)=normalize(cca1.sa)'*normalize(cca2.sa);
     end
    else
        splitflag=[splitflag,i];
    end
     
    
    

end

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
     [pls1,cca1]=pls_cancor(x1,y1);
     [pls2,cca2]=pls_cancor(x2,y2);
     pls_u_null(:,:,i)=pls1.u'*pls2.u;
     pls_v_null(:,:,i)=pls1.v'*pls2.v;
     
     if isfield(cca1,'u') && isfield(cca2,'u')
     cca_u_null(:,:,i)=cca1.u'*cca2.u;
     cca_v_null(:,:,i)=cca1.v'*cca2.v;
     cca_sa_null(:,:,i)=normalize(cca1.sa')*normalize(cca2.sa);
     cca_sb_null(:,:,i)=normalize(cca1.sa)'*normalize(cca2.sa);
     end
     
    end
    

end

for i=1:lv
    pls_rep_mean_u(i)=mean(abs(pls_u_repro(i,i,:)));
    pls_rep_std_u(i)=std(abs(pls_u_repro(i,i,:)));
    pls_rep_u_z(i)=pls_rep_mean_u(i)/pls_rep_std_u(i);
    pls_rep_ul_u(i)=prctile(abs(pls_u_repro(i,i,:)),CI);
    pls_rep_ll_u(i)=prctile(abs(pls_u_repro(i,i,:)),100-CI);
    pls_rep_mean_v(i)=mean(abs(pls_v_repro(i,i,:)));
    pls_rep_std_v(i)=std(abs(pls_v_repro(i,i,:)));
    pls_rep_v_z(i)=pls_rep_mean_v(i)/pls_rep_std_v(i);
    pls_rep_ul_v(i)=prctile(abs(pls_v_repro(i,i,:)),CI);
    pls_rep_ll_v(i)=prctile(abs(pls_v_repro(i,i,:)),100-CI);
    %if isfield(cca1,'u')==1 & isfield(cca2,'u')==1
    cca_rep_mean_u(i)=mean(abs(cca_u_repro(i,i,:)));
    cca_rep_std_u(i)=std(abs(cca_u_repro(i,i,:)));
    cca_rep_u_z(i)=cca_rep_mean_u(i)/cca_rep_std_u(i);
    cca_rep_ul_u(i)=prctile(abs(cca_u_repro(i,i,:)),CI);
    cca_rep_ll_u(i)=prctile(abs(cca_u_repro(i,i,:)),100-CI);
    cca_rep_mean_v(i)=mean(abs(cca_v_repro(i,i,:)));
    cca_rep_std_v(i)=std(abs(cca_v_repro(i,i,:)));
    cca_rep_v_z(i)=cca_rep_mean_v(i)/cca_rep_std_v(i);
    cca_rep_ul_v(i)=prctile(abs(cca_v_repro(i,i,:)),CI);
    cca_rep_ll_v(i)=prctile(abs(cca_v_repro(i,i,:)),100-CI);
    %for structure coefficients
    cca_rep_mean_sa(i)=mean(abs(cca_sa_repro(i,i,:)));
    cca_rep_std_sa(i)=std(abs(cca_sa_repro(i,i,:)));
    cca_rep_sa_z(i)=cca_rep_mean_sa(i)/cca_rep_std_sa(i);
    cca_rep_ul_sa(i)=prctile(abs(cca_sa_repro(i,i,:)),CI);
    cca_rep_ll_sa(i)=prctile(abs(cca_sa_repro(i,i,:)),100-CI);
    cca_rep_mean_sb(i)=mean(abs(cca_sb_repro(i,i,:)));
    cca_rep_std_sb(i)=std(abs(cca_sb_repro(i,i,:)));
    cca_rep_sb_z(i)=cca_rep_mean_sb(i)/cca_rep_std_sb(i);
    cca_rep_ul_sb(i)=prctile(abs(cca_sb_repro(i,i,:)),CI);
    cca_rep_ll_sb(i)=prctile(abs(cca_sb_repro(i,i,:)),100-CI);
    %end
end
 
for i=1:lv
    pls_null_mean_u(i)=mean(abs(pls_u_null(i,i,:)));
    pls_null_std_u(i)=std(abs(pls_u_null(i,i,:)));
    pls_null_u_z(i)=pls_null_mean_u(i)/pls_null_std_u(i);
    pls_null_ul_u(i)=prctile(abs(pls_u_null(i,i,:)),CI);
    pls_null_ll_u(i)=prctile(abs(pls_u_null(i,i,:)),100-CI);
    pls_null_mean_v(i)=mean(abs(pls_v_null(i,i,:)));
    pls_null_std_v(i)=std(abs(pls_v_null(i,i,:)));
    pls_null_v_z(i)=pls_null_mean_v(i)/pls_null_std_v(i);
    pls_null_ul_v(i)=prctile(abs(pls_v_null(i,i,:)),CI);
    pls_null_ll_v(i)=prctile(abs(pls_v_null(i,i,:)),100-CI);
    %if isfield(cca1,'u')==1 & isfield(cca2,'u')==1
    cca_null_mean_u(i)=mean(abs(cca_u_null(i,i,:)));
    cca_null_std_u(i)=std(abs(cca_u_null(i,i,:)));
    cca_null_u_z(i)=cca_null_mean_u(i)/cca_null_std_u(i);
    cca_null_ul_u(i)=prctile(abs(cca_u_null(i,i,:)),CI);
    cca_null_ll_u(i)=prctile(abs(cca_u_null(i,i,:)),100-CI);
    cca_null_mean_v(i)=mean(abs(cca_v_null(i,i,:)));
    cca_null_std_v(i)=std(abs(cca_v_null(i,i,:)));
    cca_null_v_z(i)=cca_null_mean_v(i)/cca_null_std_v(i);
    cca_null_ul_v(i)=prctile(abs(cca_v_null(i,i,:)),CI);
    cca_null_ll_v(i)=prctile(abs(cca_v_null(i,i,:)),100-CI);
    %for structure coefficients
    cca_null_mean_sa(i)=mean(abs(cca_sa_null(i,i,:)));
    cca_null_std_sa(i)=std(abs(cca_sa_null(i,i,:)));
    cca_null_sa_z(i)=cca_null_mean_sa(i)/cca_null_std_sa(i);
    cca_null_ul_sa(i)=prctile(abs(cca_sa_null(i,i,:)),CI);
    cca_null_ll_sa(i)=prctile(abs(cca_sa_null(i,i,:)),100-CI);
    cca_null_mean_sb(i)=mean(abs(cca_sb_null(i,i,:)));
    cca_null_std_sb(i)=std(abs(cca_sb_null(i,i,:)));
    cca_null_sb_z(i)=cca_null_mean_sb(i)/cca_null_std_sb(i);
    cca_null_ul_sb(i)=prctile(abs(cca_sb_null(i,i,:)),CI);
    cca_null_ll_sb(i)=prctile(abs(cca_sb_null(i,i,:)),100-CI);
    %end
end
 
pls_repro.pls_rep_mean_u=pls_rep_mean_u;
pls_repro.pls_rep_mean_v=pls_rep_mean_v;
pls_repro.pls_rep_z_u=pls_rep_u_z;
pls_repro.pls_rep_z_v=pls_rep_v_z;
pls_repro.pls_rep_ul_u=pls_rep_ul_u; 
pls_repro.pls_rep_ll_u=pls_rep_ll_u;
pls_repro.pls_rep_ul_v=pls_rep_ul_v; 
pls_repro.pls_rep_ll_v=pls_rep_ll_v; 
pls_repro.pls_null_mean_u=pls_null_mean_u;
pls_repro.pls_null_mean_v=pls_null_mean_v;
pls_repro.pls_null_z_u=pls_null_u_z;
pls_repro.pls_null_z_v=pls_null_v_z;
pls_repro.pls_null_ul_u=pls_null_ul_u; 
pls_repro.pls_null_ll_u=pls_null_ll_u;
pls_repro.pls_null_ul_v=pls_null_ul_v; 
pls_repro.pls_null_ll_v=pls_null_ll_v;
if dist_flag==1
    pls_repro.pls_dist_u=pls_u_repro;
    pls_repro.pls_dist_v=pls_v_repro;
    pls_repro.pls_dist_null_u=pls_u_null;
    pls_repro.pls_dist_null_v=pls_v_null;
end

%if isfield(cca1,'u')==1 & isfield(cca2,'u')==1
cca_repro.cca_rep_mean_u=cca_rep_mean_u;
cca_repro.cca_rep_mean_v=cca_rep_mean_v;
cca_repro.cca_rep_z_u=cca_rep_u_z;
cca_repro.cca_rep_z_v=cca_rep_v_z;
cca_repro.cca_rep_ul_u=cca_rep_ul_u; 
cca_repro.cca_rep_ll_u=cca_rep_ll_u;
cca_repro.cca_rep_ul_v=cca_rep_ul_v; 
cca_repro.cca_rep_ll_v=cca_rep_ll_v; 
cca_repro.cca_null_mean_u=cca_null_mean_u;
cca_repro.cca_null_mean_v=cca_null_mean_v;
cca_repro.cca_null_z_u=cca_null_u_z;
cca_repro.cca_null_z_v=cca_null_v_z;
cca_repro.cca_null_ul_u=cca_null_ul_u; 
cca_repro.cca_null_ll_u=cca_null_ll_u;
cca_repro.cca_null_ul_v=cca_null_ul_v; 
cca_repro.cca_null_ll_v=cca_null_ll_v;
%structure coefficients
cca_repro.cca_rep_mean_sa=cca_rep_mean_sa;
cca_repro.cca_rep_mean_sb=cca_rep_mean_sb;
cca_repro.cca_rep_z_sa=cca_rep_sa_z;
cca_repro.cca_rep_z_sb=cca_rep_sb_z;
cca_repro.cca_rep_ul_sa=cca_rep_ul_sa; 
cca_repro.cca_rep_ll_sa=cca_rep_ll_sa;
cca_repro.cca_rep_ul_sb=cca_rep_ul_sb; 
cca_repro.cca_rep_ll_sb=cca_rep_ll_sb; 
cca_repro.cca_null_mean_sa=cca_null_mean_sa;
cca_repro.cca_null_mean_sb=cca_null_mean_sb;
cca_repro.cca_null_z_sa=cca_null_sa_z;
cca_repro.cca_null_z_sb=cca_null_sb_z;
cca_repro.cca_null_ul_sa=cca_null_ul_sa; 
cca_repro.cca_null_ll_sa=cca_null_ll_sa;
cca_repro.cca_null_ul_sb=cca_null_ul_sb; 
cca_repro.cca_null_ll_sb=cca_null_ll_sb;
if dist_flag==1
    cca_repro.cca_dist_sa=cca_sa_repro;
    cca_repro.cca_dist_sb=cca_sb_repro;
    cca_repro.cca_dist_null_sa=cca_sa_null;
    cca_repro.cca_dist_null_sb=cca_sb_null;
end
%else
%    cca_repro.badcondition={"Cannot do CCA because X or Y are ill-conditioned"};
%end