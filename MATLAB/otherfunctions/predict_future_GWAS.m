function [h2GWAS,numGWAS,NTPR] = predict_future_GWAS(ss,ww,increase_NGWAS,sig_thresh)
%predict_future_GWAS predicts results of a larger GWAS from FMR output
%   
%   Input arguments:
%   ss: scaled variance parameters of Gaussian mixture cpts
%   ww: mixture weights as fractions of h2
%   increase_NGWAS: increase in GWAS power (eg 2 = double effective N)
%   sig_thresh: significance threshold in units that correspond to ss. If
%   ss is estimated from chi^2 statistics, then sig_thresh=30 corresponds
%   to genome-wide significance.
%   
%   Output arguments:
%   h2GWAS: proportion of h2g explained by significant SNPs
%   numGWAS: number of loci that are significant
%   NTPR: not-by-chance true positive rate at this sample size + threshold

if ~exist('sig_thresh')
    sig_thresh=30;
end

powerfn=@(x)normcdf(sqrt(increase_NGWAS)*x-sqrt(sig_thresh))+normcdf(-sqrt(increase_NGWAS)*x-sqrt(sig_thresh));
samesignfn=@(x)normcdf(sqrt(increase_NGWAS)*(x)-sqrt(sig_thresh));

no_samples=1e5;
xx=randn(no_samples,1);
xsum=sum(xx.^2);cpt_power=zeros(1,length(ww));num_sig=cpt_power;
for ii=1:length(ww)
    cpt_power(ii)=dot(powerfn(xx*sqrt(ss(ii))),xx.^2)/xsum;
    if nargout>1
        num_sig(ii)=chi2cdf(sig_thresh/(ss(ii)*increase_NGWAS+1),1,'upper');
    end
    if nargout>2
        ntpr_numer(ii)=mean(samesignfn(abs(xx)*sqrt(ss(ii)))-samesignfn(-abs(xx)*sqrt(ss(ii))));
        ntpr_denom(ii)=mean(samesignfn(abs(xx)*sqrt(ss(ii)))+samesignfn(-abs(xx)*sqrt(ss(ii))));
    end
end
h2GWAS=cpt_power*ww/sum(ww);
if nargout>1
    numGWAS=num_sig*(ww./ss');
end
if nargout>2
    NTPR=ntpr_numer*ww/(ntpr_denom*ww);
end
end

