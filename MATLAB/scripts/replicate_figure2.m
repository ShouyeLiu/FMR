% Replicate Fig 2 results using output from run_FMR_UKB_N145kvs460k.m
% Sept 2020

clear all;clc

special_plot_pheno=[5 8 13 14 15 17];
special_pheno_names={'Platelet vol','Platelet count','BMI','Height','Diastolic BP','College' };

addpath('../FMR')
addpath('../otherfunctions')

load('../../matfiles/FMRestimates_UKB_N145kvs460k.mat',...
    'ss_est','ww_est','ww_jk','h2gs_obs_thresholds','numgs_obs_thresholds',...
    'traits','LD4Mout','LD4Mout2','mm')

sig_thresh=30;
nn_ratio=460/145;

for jk=1:100
    for tt=1:22
        est_int=mean(LD4Mout(tt).intercept);
        [h2GWS_jk(jk,tt),numGWS_jk(jk,tt)] = ...
            predict_future_GWAS(ss_est(tt,:),...
            ww_jk(tt,:,jk)'*(mm*LD4Mout(tt).cov(jk))/est_int,...
            nn_ratio,sig_thresh/est_int);
    end
end

h2GWS_est=mean(h2GWS_jk);
h2GWS_err=std(h2GWS_jk)*sqrt(101);
numGWS_est=mean(numGWS_jk);
numGWS_err=std(numGWS_jk)*sqrt(101);

plotpheno=[special_plot_pheno,setdiff(1:22,special_plot_pheno)];

figure;hold on
subplot(1,2,1)
plot([0 1],[0 1],'color',[0 0 0])
errorbar_text(h2GWS_est(plotpheno),h2gs_obs_thresholds(plotpheno,1),1.96*h2GWS_err(plotpheno),zeros(size(plotpheno)),special_pheno_names(any((1:22)'==special_plot_pheno)))
xlim([0 1]);
%scatter(h2gs_est_N145k(incl,prediction_ind), h2gs_obs_N460k(incl));
title('Predicted vs observed %h^2_{GWAS}')
xlabel('Predicted at N=145k')
ylabel('Observed at N=460k')

subplot(1,2,2);hold on
plot([0 3000],[0 3000],'color',[0 0 0])
errorbar_text(numGWS_est(plotpheno),numgs_obs_thresholds(plotpheno,1),1.96*numGWS_err(plotpheno),zeros(size(plotpheno)),special_pheno_names(any((1:22)'==special_plot_pheno)))
%xlim([0 1]);
%scatter(h2gs_est_N145k(incl,prediction_ind), h2gs_obs_N460k(incl));
title('Predicted vs observed M_{GWAS}')
xlabel('Predicted at N=145k')
ylabel('Observed at N=460k')
xlim([0 3000]);ylim([0 3000])
