% Run FMR on 32 well-powered traits; also estimate h2GWAS, MGWAS
% Sept 2020

clear all;clc

savepath='../../matfiles/FMRestimates_32traits.mat';

load('../../matfiles/fourierLDscores.base.mat','SNPs','lF','l2','l4')
RefSNPs=vertcat(SNPs{:});
load('../../matfiles/1kg_LD.HM3.window1cm.noblocks.mat','RRb','LDSNPs')
load('../../matfiles/32wellpowered_sumstats.mat','sumstats')

addpath('../FMR')
addpath('../otherfunctions')

mm=length(lF);
no_blocks=100;
t_ratio_step=2;
ss=t_ratio_step.^(-7:5); %sigma^2 values of mixture cpts
tt=sqrt(ss); %sampling times
rel_wt=1;

sig_thresh_array=[30,100,300,1000];



for traitnum=1:32
    disp(sumstats(traitnum).traitname)
    traits{traitnum}=sumstats(traitnum).traitname;
    Neff(traitnum)=sumstats(traitnum).Neff;
    tic
    
    [~,i2,i1]=intersect(RefSNPs,sumstats(traitnum).SNPs_N145k,'stable');
    
    chisq=sumstats(traitnum).chisq_N145k;
    [ww,sigmasq,LD4Mout(traitnum),warningflag(traitnum)] = ...
        run_FMR(chisq(i1),lF(i2,:),l2(i2),l4(i2),...
        ss,tt,'l2Weights',1./l2(i2),'l4Weights',1./l4(i2));
    
    ww_est(traitnum,:)=mean(ww);
    ss_est(traitnum,:)=sigmasq;
    ww_jk(traitnum,:,:)=ww';
    toc
    time(traitnum)=toc;
    
    % Observed h2gwas, MGWAS
    [~,i1,i2]=intersect(sumstats(traitnum).SNPs_N145k,LDSNPs);
    [leadSNPs] = get_leadSNPs_r2(chisq(i1),RRb(i2,i2),min(sig_thresh_array));
    indep_chisq=chisq(i1);
    indep_chisq=indep_chisq(indep_chisq>min(sig_thresh_array));
    indep_chisq=indep_chisq(leadSNPs);
    for ii=1:length(sig_thresh_array)
        h2GWAS=sum(indep_chisq(indep_chisq>sig_thresh_array(ii)));
        h2gs_obs(traitnum,ii)=h2GWAS/(mm*mean(LD4Mout(traitnum).cov));
        numgs_obs(traitnum,ii)=sum(indep_chisq>sig_thresh_array(ii));
    end
end


save(savepath,'*obs','*est','warningflag','LD4Mout','*jk','time','mm','traits','Neff')
