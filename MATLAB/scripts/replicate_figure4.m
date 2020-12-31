% Replicate Fig 4 results using output from run_FMR_32traits.m
% Sept 2020

clear all;clc

addpath('../FMR')
addpath('../otherfunctions')

load('../../matfiles/FMRestimates_32traits.mat','h2gs_obs','numgs_obs',...
    'ww_est','ss_est','ww_jk','LD4Mout','mm','traits','Neff')


special_plot_pheno=[1 6 14 28 30 31 32];
special_pheno_names={'BMI','Height','SCZ','IBD','BPD','Alz','CAD' };

Mh=4e5;
sig_thresh=[4 10 20 30];

for ii=1:4
    for tt=1:32
        Nh2=mean(LD4Mout(tt).cov)*mm;
        int=mean(LD4Mout(tt).intercept);
        [h2GWAS(ii,tt),MGWAS(ii,tt)] = predict_future_GWAS(ss_est(tt,:),ww_est(tt,:)'*Nh2/int,1,sig_thresh(ii)/int);
        www=ww_est(tt,:)'*mean(LD4Mout(tt).cov*mm/Mh)./ss_est(tt,:)';
        [~,~,NTPR(ii,tt)] = predict_future_GWAS([0 ss_est(tt,:)],[1-sum(www); www],1,sig_thresh(ii));
    end
end

ind=[special_plot_pheno, setdiff(1:32,special_plot_pheno)];

figure;
empty_strings=special_pheno_names;empty_strings(:)={''};
titles={'(a) \chi^2>4','(b) \chi^2>10','(c) \chi^2>20','(d) \chi^2>30'};
for sp=1:4
    subplot(1,4,sp)
    if sp==1
        errorbar_text(h2GWAS(sp,ind)',NTPR(sp,ind)',zeros(length(ind),1),zeros(length(ind),1),special_pheno_names)
    else
        errorbar_text(h2GWAS(sp,ind)',NTPR(sp,ind)',zeros(length(ind),1),zeros(length(ind),1),empty_strings)
    end
    xlim([0 1]);ylim([0 1])
    xlabel('h^2 explained')
    ylabel('NTPR')
    title(titles{sp})
end

