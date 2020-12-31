
function [ww,ss,LD4Mout,warningflag] = run_FMR(chisq,lF,l2,l4,ss,tt,varargin)
%Run Fourier Mixture Regression to estimate the effect-size distribution
%from GWAS summary statistics and fourier LD scores.
%   Input arguments:
%   chisq: Mx1 GWAS chi^2 stats or effect-size estimates in per-s.d. units
%   lF: MxKtot matrix of Fourier LD scores. Columns correspond to ss-tt
%   pairs, but there can be more pairs than columns.
%   l2: Mx1 vector of LD scores
%   l4: Mx1 vector of LD 4th moments
%   ss: 1xK1 vector of variance parameters for the Gaussian mixture cpts
%   tt: 1xK2 vector of sampling times
%
%   Optional input arguments:
%   l2Weights: regression weights for LDSC and FMR.
%   l4Weights: regression weights for LD4M.
%   JackknifeBlocks: cell array of user-specified jackknife blocks
%   NoJackknifeBlocks: number of jackknife blocks (default 100)
%   NGWAS: GWAS sample size, if it is desired to use fixed-intercept LD
%   score regression, which is not recommended.
%   RescaleParam: parameter that allows the user to shift the mixture model
%   grid relative to to the summary statistics. When it is equal to 1, the
%   sumstats are scaled such that the mean of the HDM is centered over the
%   value 1 in the array ss. When it is larger/smaller, the mean of the HDM
%   is centered over a value that is smaller/larger respectively. If it is
%   set to 0, the scaling is performed to match the maximum chi^2 stat with
%   the maximum value of ss.
%   WeightParam: determines relative weighting of LDSC/LD4M vs FMR moment
%   eqs. By default, a large value is used so that the LDSC/LD4M eqs are
%   treated as constraints (which constrain but do not determine the
%   regression weights)
%
%   Output arguments:
%   ww: estimated mixture weights
%   ss: variance parameters, scaled to match the observed sumstats
%   LD4Mout: output of LD 4th moments regression
%   warningflag: if =1, the method determined that FMR was underpowered

p=inputParser;
addRequired(p,'chisq',@isvector)
addRequired(p,'lF',@ismatrix)
addRequired(p,'l2',@isvector)
addRequired(p,'l4',@isvector)
addRequired(p,'ss',@isvector)
addRequired(p,'tt',@isvector)
no_tvals=length(tt);
mm_regression=length(chisq);

addParameter(p,'l2Weights',1./l2,@isvector)
addParameter(p,'l4Weights',1./l4,@isvector)
addParameter(p,'JackknifeBlocks',[],@iscell)
addParameter(p,'NoJackknifeBlocks',100,@(x)isscalar(x) && all(mod(x,1)==0) && all(x<=mm_regression))
addParameter(p,'NGWAS',[],@(x)isscalar(x));
addParameter(p,'RescaleParam',1,@(x)isscalar(x) && x>=0 );
addParameter(p,'WeightParam',10^4,@(x)isscalar(x) && x>=0 );

parse(p,chisq,lF,l2,l4,ss,tt,varargin{:});

rescale_param=p.Results.RescaleParam;
wt2=p.Results.l2Weights;
wt4=p.Results.l4Weights;
jk_blocks=p.Results.JackknifeBlocks;
NoJackknifeBlocks=p.Results.NoJackknifeBlocks;
nn_GWAS=p.Results.NGWAS;
weight_param=p.Results.WeightParam;

if isempty(nn_GWAS)
    LD4Mout=LD4M(chisq,l2,l4,ones(mm_regression,1),'l2Weights',wt2,'l4Weights',wt4);
    nn_GWAS=1/mean(LD4Mout.intercept);
else
    LD4Mout=LD4M(chisq,l2,l4,ones(mm_regression,1),'l2Weights',wt2,'l4Weights',wt4,'NoiseVar',1/nn_GWAS);
end
perSNPh2=mean(LD4Mout.cov);
Eh2a2=real(perSNPh2*mm_regression./(10.^LD4Mout.log10Me));

if rescale_param==0
    rescale=max(chisq)/max(ss);
else
    rescale=full(mean(Eh2a2))*rescale_param;%perSNPh2*full(mean(ell2));%*ss(1);
end

if perSNPh2<=0
    error('Negative h2 estimate. Sample size is too small, phenotype is not heritable, or there is some error in the data input.')
end

warningflag=mean(LD4Mout.kurt)/std(LD4Mout.kurt)/sqrt(NoJackknifeBlocks+1)<2;
if mean(LD4Mout.kurt)<0
    warning('LD4M produced very noisy results. Using a heuristic to rescale summary statistics. Results probably unreliable.')
    rescale=perSNPh2*mean(l2)*10*rescale_param;
end

perSNPh2=perSNPh2/rescale;
chisq=chisq/rescale;
nn_GWAS=nn_GWAS*rescale;

hh=10^-8;
yyF=cos(tt.*sqrt(chisq))./exp(-1/2*tt.^2/nn_GWAS);
yy2=cos((hh+tt).*sqrt(chisq))./exp(-1/2*(hh+tt).^2/nn_GWAS);
yyF=-(yy2-yyF)/hh;

yy4=1/3*(chisq.^2-6*chisq/nn_GWAS+3/nn_GWAS^2-3*(chisq-1/nn_GWAS).*(l2*perSNPh2));
yy2=chisq-1/nn_GWAS;

l4=l4.*(ss-mean(chisq-1/nn_GWAS));
l2=repmat(l2,1,length(ss));

index_matrix=(1:no_tvals)+(0:no_tvals-1)';
tt_matrix=reshape(ones(no_tvals).*tt',1,no_tvals.^2);
lF=lF(:,index_matrix).*tt_matrix;
lF=reshape(lF,no_tvals*mm_regression,no_tvals);

weights=min(10^30,real(wt2./[mean(yyF.^2) mean(yy2.^2)/weight_param mean(yy4.^2)/weight_param]));

if isempty(jk_blocks)
    jk_blocks=cell(NoJackknifeBlocks,1);
    blocksize=ceil(mm_regression/NoJackknifeBlocks);
    for jk=1:NoJackknifeBlocks
        block=((jk-1)*blocksize+1:min(mm_regression,jk*blocksize)) + ...
            (0:mm_regression:(size(lF,1)-mm_regression))';
        jk_blocks{jk}=block(:);
    end
end
ww=nonnegative_regression_jk(real([lF;l2;l4]),real([yyF(:);yy2;yy4]),weights(:),'JackknifeBlocks',jk_blocks);

ss=ss*rescale;
ww=ww./sum(ww,2);%*implied_perSNPh2;
%implied_log10Me=log10(mm_reference*implied_perSNPh2./(ss*ww/implied_perSNPh2));
%sprintf('Implied effect-size variance = %f\n Implied log10Me = %.1f',full(implied_perSNPh2*rescale),full(mean(implied_log10Me)))


end




