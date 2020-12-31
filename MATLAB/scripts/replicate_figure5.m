% Replicate Fig 5 results using output from run_FMR_32traits.m
% Sept 2020
clear all;clc

addpath('../FMR')
addpath('../otherfunctions')

load('../../matfiles/FMRestimates_32traits.mat',...
    'ss_est','ww_jk','LD4Mout','mm','traits')


special_plot_pheno=[1 6 14 28 22];
special_pheno_names={'BMI','Height','SCZ','IBD','Afib' };

qq_array=[.01,.05:.05:.95, .99];
q1=2;q2=length(qq_array)-1;q3=11;q4=3;q5=q2-1;

lognormalmodel_geometric_sd=0.804;

quantiles_jk=zeros(100,length(qq_array),32);numloci_jk=quantiles_jk;
for tt=1:32
    Nh2=mean(LD4Mout(tt).cov)*mm;
    [quantiles_jk(:,:,tt),numloci_jk(:,:,tt)] = estimate_h2_quantiles(ss_est(tt,:),reshape(ww_jk(tt,:,:),13,100)',qq_array,Nh2);
    quantiles_est(tt,:)=mean(log10(quantiles_jk(:,:,tt)));
    quantiles_se(tt,:)=std(log10(quantiles_jk(:,:,tt)))*sqrt(101);

    numloci_est(tt,:)=mean(log10(numloci_jk(:,:,tt)));
    numloci_se(tt,:)=std(log10(numloci_jk(:,:,tt)))*sqrt(101);
    
    expected_quantiles(tt,:)=norminv(qq_array,quantiles_est(tt,q3),lognormalmodel_geometric_sd);
    
end

qmean=mean(log10(quantiles_jk),3);
nmean=mean(log10(numloci_jk),3);
qmean_est=mean(qmean);
qmean_se=std(qmean)*sqrt(101);
nmean_est=mean(nmean);
nmean_se=std(nmean)*sqrt(101);
expectedmean=mean((expected_quantiles));




%% plotting
ind=[special_plot_pheno, setdiff(1:32,special_plot_pheno)];

incl_threshold=.5;
wellpowered=ind(max(quantiles_se(ind,[q3 q4 q5])')<incl_threshold);
special_wellpowered=max(quantiles_se(special_plot_pheno,[q3 q4 q5])')<incl_threshold;

figure;subplot(1,3,1)
zeropad=zeros(length(wellpowered)-sum(special_wellpowered),1);
hold on
%plot([-6 -2.5],[-6 -2.5],'color',[.5 .5 .5])
%errorbar_text(numloci_est(ind,q3),quantiles_est(ind,q3),[1.96*numloci_se(special_plot_pheno)';zeropad],[1.96*quantiles_se(special_plot_pheno,q3);zeropad],special_pheno_names)
errorbar_text(quantiles_est(wellpowered,q3),numloci_est(wellpowered,q3),...
    [1.96*quantiles_se(special_plot_pheno(special_wellpowered),q3);zeropad],...
    [1.96*numloci_se(special_plot_pheno(special_wellpowered),q3);zeropad],...
    special_pheno_names(special_wellpowered))

title('(a) Median of the effect size distribution')
xlabel('Median effect size')
ylabel('No. loci > median')
ylim([2.5 4])
xlim([-4.75 -3.25])
set(gca,'YTick',[3 3.5 4 4.5],'YTickLabel',{'10^{3}','','10^{4}',''})
set(gca,'XTick',[-4.5 -4 -3.5],'XTickLabel',{'10^{-4.5}','10^{-4}','10^{-3.5}'})

subplot(1,3,2)
hold on
errorbar(qq_array(2:end-1),qmean_est(2:end-1),qmean_se(2:end-1)*1.96,'-')
plot(qq_array(2:end-1),expectedmean(2:end-1));
xlim([0 1])
ylim([-5.5 -2.5])
set(gca,'XTick',[.1 .5 .9],'XTickLabel',{'10','50','90'})
set(gca,'YTick',-6:-2,'YTickLabel',{'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}'})
xlabel('h^2 percentile')
ylabel('Effect size (%h^2)')
title('(b) Meta-analyzed effect size distribution')
legend('FMR estimates','Log-normal fit');legend boxoff

subplot(1,3,3)
hold on
q4=3;q5=length(qq_array)-2;
t1=quantiles_se(:,q4);t2=quantiles_se(:,q5);
%wellpowered=ind(t1(ind)<incl_threshold & t2(ind)<incl_threshold);
%special_wellpowered=(t1(special_plot_pheno)<incl_threshold & t2(special_plot_pheno)<incl_threshold );
zeropad=zeros(length(wellpowered)-sum(special_wellpowered),1);
%plot([-6 -2.5],[-6 -2.5],'color',[.5 .5 .5])
plot([-6 -4.5],[-4 -2.5],'color',[1 1 1]/2)
errorbar_text(quantiles_est(wellpowered,q4),quantiles_est(wellpowered,q5),...
    [1.96*t1(special_plot_pheno(special_wellpowered));zeropad],...
    [1.96*t2(special_plot_pheno(special_wellpowered));zeropad],...
    special_pheno_names(special_wellpowered))
%text(2*10^5,2*10^6,'10\times','color',[.5 .5 .5])
title('(c) 10th vs 90th percentile')
xlabel('10th percentile')
ylabel('90th percentile')
ylim([-4 -2.25]);xlim([-6 -4.25])
set(gca,'XTick',[-6 -5.5 -5 -4.5],'XTickLabel',{'10^{-6}','','10^{-5}',''})
%xlim([-5.5 -4])
set(gca,'YTick',[-4 -3.5 -3 -2.5],'YTickLabel',{'10^{-4}','','10^{-3}',''})



quantiles_table=[quantiles_est quantiles_se expected_quantiles];
quantiles_table=[traits(incl) num2cell(quantiles_table(incl,:))];
