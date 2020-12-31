function [quantiles,numloci] = estimate_h2_quantiles(ss,ww,qq,Nh2)
%estimate_h2_quantiles estimates the quantiles of fixed-effect heritability
%from FMR output.
%   
%   Input arguments:
%   ss: scaled variance parameters of Gaussian mixture cpts
%   ww: mixture weights as fractions of h2
%   qq: quantiles of the distribution to be estimated
%   Nh2: total genetic variance, on the scale of ss parameters. For
%   example, if ss is estimated from chi^2 statistics, this is equal to N
%   times the observed-scale heritability
%
%   Output arguments:
%   quantiles: effect sizes corresponding to the quantiles of h2
%   numloci: number of loci req'd to explain the quanitles of h2


if ~exist('Nh2')
    Nh2=1;
    warning('Assuming that ss parameter is scaled to be a fraction of h2, which is not the default FMR output')
end

no_samples=1e5;
xx=randn(no_samples,1);
pp=rand(no_samples,1);
ss=reshape(ss,length(ss),1)/Nh2;

for jk=1:size(ww,1)
    wt=cumsum(ww(jk,:))/sum(ww(jk,:));
    ii=discretize(pp,[0 wt]);% Pick a random cpt with probability ww
    [samples,order]=sort(xx.^2.*ss(ii));
    weights=1./ss(ii(order));
    dist=cumsum(samples.*weights)/sum(samples.*weights);
    for kk=1:length(qq)
        tt=find(dist>=qq(kk),1,'first');
        numloci(jk,kk)=sum(weights(tt:end))/sum(weights.*samples);
        quantiles(jk,kk)=samples(tt);
    end
end
end

