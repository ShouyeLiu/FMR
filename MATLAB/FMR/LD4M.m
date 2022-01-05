function  LD4Mout = LD4M( chisq,l2,l4,annot,varargin)
%LD4M (stratified LD 4th moments regression).
%
%   REQUIRED INPUT ARGUMENTS:
%   chisq:  Mx1 vector of chi^2 statistics;
%   l2: Mtot x P matrix of LD scores (LD second moments), where P is no. annotations
%   and Mtot is the number of reference SNPs;
%   l4: Mtot x P matrix of LD 4th moments, same size as l2;
%   annot: Mtot x P matrix of annotation values;

%   OPTIONAL INPUT ARGUMENTS as name-value pairs:
%   RegressionIndices: Mx1 vector of indices corresponding to the reference SNPs that
%   will be used in the regression (regression SNPs must be a subset of
%   reference SNPs);
%   l2Weights: Mx1 vector of weights for LDSC;
%   l4Weights: Mx1 vector of weights for LD4M;
%   NoJackknifeBlocks: number of jackknife blocks (default 100);
%   OutputSubspace: which annotations to report estimates for. This can be
%   either a list of P' indices <=P, or a MxP' projection matrix whose columns
%   correspond to linear combinations of input annotations.
%   ConditioningSubspace: which annotations to condition on. This can be
%   either a list of P' indices <=P, or a MxP' projection matrix whose columns
%   correspond to linear combinations of input annotations.
%   NoBootstrapIter: optionally, run block bootstrap instead of block
%   jackknife (default: 0 -> run block jackknife)
%   SeedBootstrap: optionally, set random seed for bootstrap resampling
%   MixedModel: 1 if desired to treat jackknife blocks containing
%   large-effect loci as fixed effects. Other SNPs in the same block will
%   be ignored.
%   MixedModelThreshold: chi^2 threshold such that if any SNP in a
%   block is larger than this, it is treated as a fixed effect
%
%   OUTPUT ARGUMENTS: all outputs are matrices of size NoJacknifeBlocks x
%   P', ie each row = a leave-one-block-out estimate and each column = an
%   annotation.
%   Qannot: excess heritability overlap compared with null model including
%   all annotations
%   Qsub: excess heritability overlap compared with null model including
%   subset of annotations, if OutputSubspace is specified
%   Me1: effective number of independent SNPs for trait 1
%   Me2: effective number of independent SNPs for trait 2
%   rg: genetic correlation
%   h21: "heritability" for trait 1, if input summary stats are in
%   standardized units
%   h22: "heritability" for trait 2, if input summary stats are in
%   standardized units

mm_regression=length(chisq);
mean_chisq=mean(chisq);
[mm_tot,pp_l2]=size(annot);

checkvector=@(x)isreal(x) && isvector(x) && length(x)==mm_regression;
check4matrix=@(x)isreal(x) && size(x,1)==mm_regression && size(x,2)==4;
checkmatrix=@(x)isreal(x) && size(x,1)==mm_tot && size(x,2)==pp_l2;
checkmatrix_l=@(x)isreal(x) && size(x,1)==mm_tot && (size(x,2)==pp_l2 || size(x,2)==pp_l2+1);
checkIndexlist=@(x)(isvector(x) && all(mod(x,1)==0) && all(x<=pp_l2) && all(x>0));
checkProjectionmatrix=@(x)isreal(x) && (size(x,1)==pp_l2 && rank(x)==size(x,2));
p=inputParser;

addRequired(p,'chisq',checkvector)
addRequired(p,'l2',checkmatrix_l)
addRequired(p,'l4',checkmatrix_l)
addRequired(p,'annot',checkmatrix)
addParameter(p,'RegressionIndices',1:mm_tot,@(x)checkvector(x) && all(mod(x,1)==0) && all(x<=mm_tot))
addParameter(p,'l2Weights',[],checkvector)
addParameter(p,'l4Weights',[],checkvector)
addParameter(p,'NoJackknifeBlocks',100,@(x)isscalar(x) && all(mod(x,1)==0) && all(x<=mm_regression))
addParameter(p,'JackknifeBlocks',[],@(x)true)
addParameter(p,'NoBootstrapIter',0,@(x)isscalar(x) && all(mod(x,1)==0) && all(x>=0))
addParameter(p,'SeedBootstrap',0,@(x)isscalar(x) && all(mod(x,1)==0) && all(x>=0))
addParameter(p,'OutputSubspace',[],@(x)checkProjectionmatrix(x) || checkIndexlist(x))
addParameter(p,'NoiseVariance',[],@(x)isscalar(x));
addParameter(p,'MixedModel',0,@(x)isscalar(x) && (x==0 || x==1));
addParameter(p,'MixedModelThreshold',mean_chisq*100,@(x)isscalar(x) && x>0 );
addParameter(p,'NoRefSNPs',[],@(x)isscalar(x) && x>=mm_regression );

parse(p,chisq,l2,l4,annot,varargin{:});
if checkIndexlist(p.Results.OutputSubspace)
    output_subspace=sparse(p.Results.OutputSubspace,1:length(p.Results.OutputSubspace),ones(length(p.Results.OutputSubspace),1),pp_l2,length(p.Results.OutputSubspace));
else
    output_subspace=p.Results.OutputSubspace;
end

no_blocks=p.Results.NoJackknifeBlocks;
idx_regression=p.Results.RegressionIndices;

if ~isempty(p.Results.NoRefSNPs)
    mm_ref = p.Results.NoRefSNPs;
else
    mm_ref=mm_tot;
end

if isempty(p.Results.l2Weights)
    l2weights=1./max(1,l2(idx_regression,1));
else
    l2weights=p.Results.l2Weights;% 2nd-moment weights matrix
end
WW2_mat=diag(sparse(l2weights));


if isempty(p.Results.l4Weights)
    l4weights=1./max(1,l4(idx_regression,1));
else
    l4weights=p.Results.l4Weights;% 2nd-moment weights matrix
end
WW4_mat=diag(sparse(l4weights));

NoiseVariance=p.Results.NoiseVariance;

run_bootstrap=p.Results.NoBootstrapIter>0;
if  run_bootstrap
    no_iters=p.Results.NoBootstrapIter;
    if p.Results.SeedBootstrap==0
        rng('shuffle');
    else
        rng(p.Results.SeedBootstrap)
    end
else
    no_iters=no_blocks;
end

% Jackknife blocks (ref SNPs)
if checkIndexlist(p.Results.JackknifeBlocks)
    reference_blocks=p.Results.JackknifeBlocks;
    no_blocks=length(reference_blocks);
else
    blocksize=floor(mm_tot/no_blocks);
    reference_blocks=cell(no_blocks,1);
    for jk=1:no_blocks
        reference_blocks{jk}=(jk-1)*blocksize+1:jk*blocksize;
    end
end

% Jackknife blocks (regression SNPs)
regression_blocks=reference_blocks;
for jk=1:no_blocks
    regression_blocks{jk}=...
        find(idx_regression >= min(reference_blocks{jk}),1,'first') : ...
        find(idx_regression <= max(reference_blocks{jk}),1,'last');
end

% Mixed model setup
fixedeffect_blocks=false(no_blocks,1);
if p.Results.MixedModel
    MixedModelThreshold=p.Results.MixedModelThreshold;
    
    for ii=1:no_blocks
        if max(chisq(regression_blocks{ii}))>MixedModelThreshold
            fixedeffect_blocks(ii)=true;
        end
    end
    fprintf('Treating %d out of %d blocks as fixed effects\n',...
        sum(fixedeffect_blocks),no_blocks);
    if ~run_bootstrap
        no_iters=no_iters-sum(fixedeffect_blocks);
    end
    
end
clear p varargin


if isempty(output_subspace);  output_subspace=eye(pp_l2);end

output_annot=annot*output_subspace;

% Handling LD score intercept options
if size(l2,2)==pp_l2 && isempty(NoiseVariance)
    l2=[ones(mm_tot,1) l2];
end
if size(l2,2)==pp_l2+1 && any(l2(:,1)~=1)
    warning('If size(l2,2)==size(annot,2)+1, then first column of l2 should usually be all ones')
end
if size(l2,2)==pp_l2+1 && ~isempty(NoiseVariance)
    warning('If size(l2,2)==size(annot,2)+1, then NoiseCovariance input is ignored')
end

l2_intercept=size(l2,2)-pp_l2;
l4_intercept=size(l4,2)-pp_l2;

pp_l2=size(l2,2);
pp_l4=size(l4,2);
pp_annot=size(annot,2);
pp_output=size(output_annot,2);

blocksize=floor(mm_regression/no_blocks);

ell_ell=zeros(pp_l2,pp_l2,no_blocks);
ell_ld4m=zeros(pp_l2,pp_l4,no_blocks);
ld4m_ld4m=zeros(pp_l4,pp_l4,no_blocks);
annot_annot=zeros(pp_annot,pp_output,no_blocks);
%ell_annot=zeros(pp_annot,pp_output,no_blocks);
ell_ld4m_aa=zeros(pp_l2,pp_l4,no_blocks);
annot_annot_aa=zeros(pp_annot,pp_output,no_blocks);
ell_aa=zeros(pp_l2,no_blocks);
ld4m_aa=zeros(pp_l4,no_blocks);
ld4m_aaaa=zeros(pp_l4,no_blocks);
ell_mean=zeros(pp_l2,no_blocks);
ld4m_mean=zeros(pp_l4,no_blocks);
aa=zeros(no_blocks,1);
% Compute various moments for each jackknife block
for jk=1:no_blocks
    ind=regression_blocks{jk};
    ind2=idx_regression(regression_blocks{jk});
    ind3=reference_blocks{jk};
    ell_ell(:,:,jk)=l2(ind2,:)'*WW2_mat(ind,ind)*l2(ind2,:)/sum(l2weights(ind));
    
    ld4m_ld4m(:,:,jk)=l4(ind2,:)'*WW4_mat(ind,ind)*l4(ind2,:)/sum(l4weights(ind));
    ell_ld4m_aa(:,:,jk)=((l2(ind2,:).*chisq(ind))'*WW4_mat(ind,ind)*l4(ind2,:))/sum(l4weights(ind));
    ell_aa(:,jk)=l2(ind2,:)'*WW2_mat(ind,ind)*(chisq(ind))/sum(l2weights(ind));%
    ld4m_aa(:,jk)=l4(ind2,:)'*WW4_mat(ind,ind)*(chisq(ind))/sum(l4weights(ind));
    annot_annot_aa(:,:,jk)=(annot(ind2,:).*chisq(ind))'*output_annot(ind2,:)./sum(output_annot(ind2,:));
    
    ell_ld4m(:,:,jk)=(l2(ind2,:))'*WW4_mat(ind,ind)*l4(ind2,:)/sum(l4weights(ind));
    ld4m_aaaa(:,jk)=l4(ind2,:)'*WW4_mat(ind,ind)*chisq(ind).^2/sum(l4weights(ind));%
    
    ld4m_mean(:,jk)=mean(WW4_mat(ind,ind)*l4(ind2,:))/mean(l4weights(ind));
    ell_mean(:,jk)=mean(WW2_mat(ind,ind)*l2(ind2,:))/mean(l2weights(ind));
    
    annot_annot(:,:,jk)=annot(ind3,:)'*output_annot(ind3,:)./sum(output_annot(ind3,:));
    %ell_annot(:,:,jk)=l2(ind3,l2_intercept+1:end)'*output_annot(ind3,:)./sum(output_annot(ind3,:));
    
    aa(jk)=mean(WW2_mat(ind,ind)*chisq(ind))/mean(l2weights(ind));
end

% Combine jackknife estimates to obtain leave-one-out regression
% coefficients
covarcoef=zeros(pp_annot,no_iters+1);
intercept=zeros(1,no_iters+1);
kurtcoef=zeros(pp_annot,no_iters+1);
kurtexcess=zeros(pp_output,no_iters+1);
covar=zeros(pp_output,no_iters+1);

ind=cell(no_iters+1,1);
included_blocks=find(~fixedeffect_blocks);
for jk=1:no_iters
    if run_bootstrap
        ind{jk}=randsample(included_blocks,length(included_blocks),true);
    else
        ind{jk}=included_blocks([1:jk-1,jk+1:no_iters]);
    end
end
ind{no_iters+1}=included_blocks;

% Leave-one-out estimates
for jk=1:no_iters+1
    
    %   S-LDSC
    if l2_intercept==1
        temp=mean(ell_ell(:,:,ind{jk}),3)\mean(ell_aa(:,ind{jk}),2);% l2var^-1 * l2a2
        covarcoef(:,jk)=temp(2:end);
        intercept(jk)=temp(1);
    else
        temp=mean(ell_ell(:,:,ind{jk}),3)\(mean(ell_aa(:,ind{jk})-...
            NoiseVariance*ell_mean(:,ind{jk}),2));% l2var^-1 * l2a2
        covarcoef(:,jk)=temp;
        intercept(jk)=NoiseVariance;
    end
    
    % 4th cumulants
    ld4m_ld4m_leave1out=mean(ld4m_ld4m(l4_intercept+1:end,l4_intercept+1:end,ind{jk}),3);
    
    ld4m_aaaa_leave1out=mean(ld4m_aaaa(l4_intercept+1:end,ind{jk}),2);
    
    
    % Correct for sampling noise 
    ld4m_aaaa_leave1out=ld4m_aaaa_leave1out-...
        6*(intercept(jk)*mean(ld4m_aa(l4_intercept+1:end,ind{jk}),2)-...
        intercept(jk)*intercept(jk)*mean(ld4m_mean(l4_intercept+1:end,ind{jk}),2)/2);
    
    % Correct for expected 4th moment 
    ld4m_aaaa_leave1out=ld4m_aaaa_leave1out-...
        6*((mean(ell_ld4m_aa(l2_intercept+1:end,l4_intercept+1:end,ind{jk}),3)'*covarcoef(:,jk)-...
        intercept(jk)*mean(ell_ld4m(l2_intercept+1:end,l4_intercept+1:end,ind{jk}),3)'*covarcoef(:,jk))/2);

    
    % Fourth cumulants of effect-size distribution
    kurtcoef(:,jk)=ld4m_ld4m_leave1out\ld4m_aaaa_leave1out;
    kurtexcess(:,jk)=kurtcoef(:,jk)'*mean(annot_annot(:,:,ind{jk}),3);
    
    % Second cumulants of effect-size distribution
    covar(:,jk)=covarcoef(:,jk)'*mean(annot_annot(:,:,ind{jk}),3);
    
end

avar=mean(aa)-intercept(end);

% fixed-effect h2 and kurt
fixedeffect_blocks=find(fixedeffect_blocks);
fixedeffect_h2=zeros(length(fixedeffect_blocks),1);
for jk=1:length(fixedeffect_blocks)
    fixedeffect_h2(jk)=max(chisq(regression_blocks{fixedeffect_blocks(jk)}));
end
fixedeffect_kurt=fixedeffect_h2.^2;


base_ind=find(sum(annot==1)==mm_tot,1);
if isempty(base_ind)
    warning('No base annotation detected; using first column as base');
    base_ind=1;
end
tempkurt=sum(fixedeffect_kurt)/mm_ref;
tempcov=sum(fixedeffect_h2)/mm_ref;

covcombined=zeros(1,no_blocks);
kurtexcesscombined=covcombined;
for jk=1:no_iters
    kurtexcesscombined(included_blocks(jk))=tempkurt+kurtexcess(base_ind,jk);
    covcombined(included_blocks(jk))=tempcov+covar(base_ind,jk);
end


for jk=1:length(fixedeffect_kurt)
    kurtexcesscombined(fixedeffect_blocks(jk))=tempkurt+kurtexcess(base_ind,end)-...
        fixedeffect_kurt(jk)/mm_ref;
    covcombined(fixedeffect_blocks(jk))=tempcov+covar(base_ind,end)-...
        fixedeffect_h2(jk,:)/mm_ref;
end

Me=3*sum(annot(:,base_ind))*covcombined.^2./(kurtexcesscombined+3*covcombined.*avar);

LD4Mout(1).kurt=kurtexcess(:,1:end-1);
LD4Mout(1).kurtcoef=kurtcoef(:,1:end-1);
LD4Mout(1).cov=covar(:,1:end-1);
LD4Mout(1).covcoef=covarcoef(:,1:end-1);
LD4Mout(1).fixedeffecth2=fixedeffect_h2;
LD4Mout(1).fixedeffectkurt=fixedeffect_kurt;
LD4Mout(1).combinedkurt=kurtexcesscombined;
LD4Mout(1).combinedcov=covcombined;
LD4Mout(1).intercept=intercept(1:end-1);
LD4Mout(1).fixedeffectblocks=fixedeffect_blocks;
LD4Mout(1).log10Me=log10(Me);
LD4Mout(1).kurt=kurtexcess(:,1:end-1);

end




