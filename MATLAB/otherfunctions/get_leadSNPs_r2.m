function [leadSNPs,h2GWAS] = get_leadSNPs_r2(chisq,RRb,sig_thresh)
%get_leadSNPs_r2 performs greedy pruning + thresholding to select GWAS lead
%SNPs
%   Input: chisq: GWAS sumstats. RRb: boolean sparse LD matrix. SNPs with
%   LD greater than theshold have 'true'. sig_thresh: chisq significance
%   threshold.
%   Output: leadSNPs: boolean vector of whether each SNP was included.
%   h2GWAS: sum of chisq statistics for the included SNPs.

incl=chisq>sig_thresh;
RRb=RRb(incl,incl);
chisq=chisq(incl);

[chisq,ix]=sort(chisq,'descend');
RRb=RRb(ix,ix);

ii=1;
incl=true(size(chisq));
leadSNPs=false(size(chisq));
while any(incl)
    ii=ii+1;
    imax=find(incl,1,'first');
    incl(RRb(:,imax))=false;
    leadSNPs(imax)=true;
    if mod(ii,1000)==0
        disp(ii)
    end
end

h2GWAS=sum(chisq(leadSNPs));
end

