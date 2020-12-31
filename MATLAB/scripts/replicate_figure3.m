% Replicate Fig 3 results using output from run_FMR_32traits.m
% Sept 2020

clear all;clc

addpath('../FMR')
addpath('../otherfunctions')

load('../../matfiles/FMRestimates_32traits2.mat','h2gs_obs','numgs_obs',...
    'ww_est','ss_est','ww_jk','LD4Mout','mm','traits','Neff')


special_plot_pheno=[1 6 14 28 30 31 32];
special_pheno_names={'BMI','Height','SCZ','IBD','BPD','Alz','CAD' };

nn_ratio_array=10.^(-1:.05:5);
sig_thresh=30;

h2GWAS=zeros(length(nn_ratio_array),32);numGWAS=h2GWAS;
for ii=1:length(nn_ratio_array)
    for tt=1:32
        int=mean(LD4Mout(tt).intercept);
        [h2GWAS(ii,tt),numGWAS(ii,tt)] = predict_future_GWAS(ss_est(tt,:),ww_est(tt,:)'*(mm*mean(LD4Mout(tt).cov))/int,nn_ratio_array(ii),sig_thresh/int);
    end
end
for tt=1:32
    ix50=find(h2GWAS(:,tt)>=.5,1,'first');
    ix90=find(h2GWAS(:,tt)>=.9,1,'first');
    NN50(tt)=nn_ratio_array(ix50);
    NN90(tt)=nn_ratio_array(ix90);
    MGWAS_N50(tt)=numGWAS(ix50,tt);
    MGWAS_N90(tt)=numGWAS(ix90,tt);
    Nh2(tt)=mean(LD4Mout(tt).cov)*10^7;
end

ss=2.^(-7:5);
for tt=1:32
    r2fn=@(nn)(ss_est(tt,:)./(1/nn+ss_est(tt,:)))*ww_est(tt,:)';
    r2_50_est(tt)=r2fn(NN50(tt));
end

Nh250=NN50.*Nh2;
Nh290=NN90.*Nh2;

NN50=NN50.*Neff;
NN90=NN90.*Neff;

ind=[special_plot_pheno, setdiff(1:32,special_plot_pheno)];

figure;subplot(1,3,1)
hold on
plot([1e5,2*10^7],[1e6,2*10^8],'color',[.5 .5 .5])
errorbar_text(NN50(ind),NN90(ind),zeros(length(ind),1),zeros(length(ind),1),special_pheno_names)
title('Sample size requirements')
xlabel('N_{eff} for %h^2_{GWAS}=0.5')
ylabel('N_{eff} for %h^2_{GWAS}=0.9')
set(gca,'YScale','log','XScale','log','XTick',10.^(5:7),'YTick',10.^(6:8))
text(2*10^5,2*10^6,'y=10x','color',[.5 .5 .5])

subplot(1,3,2)
hold on
errorbar_text(NN50(ind),MGWAS_N50(ind),zeros(length(ind),1),zeros(length(ind),1),special_pheno_names)
title('Number of discoveries when %h^2_{GWAS}=0.5')
xlabel('N_{eff} for %h^2_{GWAS}=0.5')
ylabel('M_{GWAS}')
set(gca,'XScale','log','XTick',10.^(5:7),'Yscale','log')
%ylim([0 1]);
xlim([1e5 2e7])

subplot(1,3,3)
hold on
errorbar_text(NN50(ind),r2_50_est(ind),zeros(length(ind),1),zeros(length(ind),1),special_pheno_names)
title('PGS r^2 when %h^2_{GWAS}=0.5')
xlabel('N_{eff} for %h^2_{GWAS}=0.5')
ylabel('r^2_{PGS}/h^2_g')
set(gca,'XScale','log','XTick',10.^(5:7),'YTick',[.5:.1:1])
ylim([.5 1]);xlim([1e5 2e7])
