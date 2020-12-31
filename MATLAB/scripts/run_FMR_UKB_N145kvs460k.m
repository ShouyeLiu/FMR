
% Run FMR on 22 UK Biobank traits with N=145k and N=460k sumstats
% September 2020
clear all;clc
special_plot_pheno=[5,8,16:17,19,21];% Must be in order!
special_pheno_names={'Platelet vol','Platelet count','BMI','Height','Diastolic BP','College' };

savepath='../../matfiles/FMRestimates_UKB_N145kvs460k.mat';

load('../../matfiles/fourierLDscores.base.mat','SNPs','lF','l2','l4')
RefSNPs=vertcat(SNPs{:});
load('../../matfiles/1kg_LD.HM3.window1cm.noblocks.mat','RRb','LDSNPs')
load('../../matfiles/UKB_sumstats_N145kvs460k.mat','sumstats')

addpath('../FMR')
addpath('../otherfunctions')

mm=length(lF);
no_blocks=100;
t_ratio_step=2;
ss=t_ratio_step.^(-7:5); %sigma^2 values of mixture cpts
tt=sqrt(ss); %sampling times
rel_wt=1;

sig_thresh_array=[30,100,300,1000];

for traitnum=1:22
    
    disp(sumstats(traitnum).traitname)
    traits{traitnum}=sumstats(traitnum).traitname;
    
    tic
    % Predicted (N=145k)
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
    
    % Observed (N=460k)
    chisq2=sumstats(traitnum).chisq_N460k;

    [~,i2,i1]=intersect(RefSNPs,sumstats(traitnum).SNPs_N460k,'stable');
    LD4Mout2(traitnum)=LD4M(chisq2(i1),l2(i2),l4(i2),ones(length(i2),1),...
        'l2Weights',1./l2(i2),'l4Weights',1./l4(i2));
    
    [~,i1,i2]=intersect(sumstats(traitnum).SNPs_N460k,LDSNPs);
    
    [leadSNPs] = get_leadSNPs_r2(chisq2(i1),RRb(i2,i2),min(sig_thresh_array));
    indep_chisq=chisq2(i1);
    indep_chisq=indep_chisq(indep_chisq>min(sig_thresh_array));
    indep_chisq=indep_chisq(leadSNPs);
    
    % Predicted/observed at different significance thresholds
    for ii=1:length(sig_thresh_array)
        h2GWAS=sum(indep_chisq(indep_chisq>sig_thresh_array(ii)));
        h2gs_obs_thresholds(traitnum,ii)=h2GWAS/(mm*mean(LD4Mout2(traitnum).cov));
        numgs_obs_thresholds(traitnum,ii)=sum(indep_chisq>sig_thresh_array(ii));
    end
    
end

save(savepath,'*est*','*obs*','*jk','LD4Mout*','warningflag','time','traits','mm')
