function [Spatial_correlation] = plot_STDMIS(Spatial_correlation,STD,MIS);
sessions=Spatial_correlation.options.sessions;
PV=Spatial_correlation.PVcorrelation.AllROI.PVcorr_table;
TC=Spatial_correlation.TuningCurvecorr.AllROI.TCcorr_table;
TC_tuned=Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_table;
rate=Spatial_correlation.overlap.Rate_table;
rate_tuned=Spatial_correlation.overlap.Rate_tuned_table;
peak_tuned=Spatial_correlation.overlap.Peak_tuned_table;

%Plot PV
for i=1:size(MIS,2)
PV_STD=PV(:,STD);
PV_MIS{i}=PV(:,MIS{i});
PV_STD_mean=reshape(PV_STD,[],1);
PV_MIS_mean{i}=reshape(PV_MIS{i},[],1);
end
figure;
PV_MIS_mean_mat=cell2mat(PV_MIS_mean);
PV_all=[PV_STD_mean PV_MIS_mean_mat];
boxplot([PV_STD_mean PV_MIS_mean_mat])

PVcorr=[PV_STD_mean PV_MIS_mean_mat];
%% Cumulative - to remove
figure;
for i=1:size(PVcorr,2)
    [f{i},x{i}] = ecdf(PVcorr(:,i));

 hold on;  
 plot(x{i},f{i});
%hold off;
title('mean PV correlation - all ROI')
ylabel('cumulative frequency')
xlabel('Pearson correlation score')
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
end
%%

%Plot TC all
for i=1:size(MIS,2)
TC_STD=TC(:,STD);
TC_MIS{i}=TC(:,MIS{i});
TC_STD_mean=reshape(TC_STD,[],1);
TC_MIS_mean{i}=reshape(TC_MIS{i},[],1);
end
figure;
TC_MIS_mean_mat=cell2mat(TC_MIS_mean);
TC_all=[TC_STD_mean TC_MIS_mean_mat];
boxplot([TC_STD_mean TC_MIS_mean_mat])
title('mean tunning curve correlation - all ROI')
ylabel('pearson correlation score') % y-axis label

TC_corr=[TC_STD_mean TC_MIS_mean_mat];

%% Cumulative - to remove
figure;
for i=1:size(PVcorr,2)
    [ft{i},xt{i}] = ecdf(TC_corr(:,i));

 hold on;  
 plot(x{i},f{i});
%hold off;
title('mean tuning vector correlation - all ROI')
ylabel('cumulative frequency')
xlabel('Pearson correlation score')
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
end
%%



%Plot TC tuned
for i=1:size(MIS,2)
TC_tuned_STD=TC_tuned(:,STD);
TC_tuned_MIS{i}=TC_tuned(:,MIS{i});
TC_tuned_STD_mean=reshape(TC_tuned_STD,[],1);
TC_tuned_MIS_mean{i}=reshape(TC_tuned_MIS{i},[],1);
end
figure;
TC_tuned_MIS_mean_mat=cell2mat(TC_tuned_MIS_mean);
TC_all_tuned=[TC_tuned_STD_mean TC_tuned_MIS_mean_mat];
boxplot([TC_tuned_STD_mean TC_tuned_MIS_mean_mat])
title('mean tunning curve correlation - Tuned ROI')
ylabel('pearson correlation score') % y-axis label


%Plot rate all
for i=1:size(MIS,2)
rate_STD=rate(:,STD);
rate_MIS{i}=rate(:,MIS{i});
rate_STD_mean=reshape(rate_STD,[],1);
rate_MIS_mean{i}=reshape(rate_MIS{i},[],1);
end
figure;
rate_MIS_mean_mat=cell2mat(rate_MIS_mean);
rate_all=[rate_STD_mean rate_MIS_mean_mat];
boxplot([rate_STD_mean rate_MIS_mean_mat])
title('Rate overlap - All ROI')
ylabel('Rate ratio') % y-axis label

%Plot rate tuned
for i=1:size(MIS,2)
rate_tuned_STD=rate_tuned(:,STD);
rate_tuned_MIS{i}=rate_tuned(:,MIS{i});
rate_tuned_STD_mean=reshape(rate_tuned_STD,[],1);
rate_tuned_MIS_mean{i}=reshape(rate_tuned_MIS{i},[],1);
end
figure;
rate_tuned_MIS_mean_mat=cell2mat(rate_tuned_MIS_mean);
rate_all_tuned=[rate_tuned_STD_mean rate_tuned_MIS_mean_mat];
boxplot([rate_tuned_STD_mean rate_tuned_MIS_mean_mat])
title('Rate overlap - Tuned ROI')
ylabel('Rate ratio') % y-axis label

%Plot peak tuned
for i=1:size(MIS,2)
peak_tuned_STD=peak_tuned(:,STD);
peak_tuned_MIS{i}=peak_tuned(:,MIS{i});
peak_tuned_STD_mean=reshape(peak_tuned_STD,[],1);
peak_tuned_MIS_mean{i}=reshape(peak_tuned_MIS{i},[],1);
end
figure;
peak_tuned_MIS_mean_mat=cell2mat(peak_tuned_MIS_mean);
peak_all_tuned=[peak_tuned_STD_mean peak_tuned_MIS_mean_mat];
boxplot([peak_tuned_STD_mean peak_tuned_MIS_mean_mat])
title('Peak overlap - Tuned ROI')
ylabel('Peak ratio') % y-axis label


Spatial_correlation.PVcorrelation.AllROI.PVcorr_STDvsMIS=PV_all;
Spatial_correlation.TuningCurvecorr.AllROI.TCcorr_STDvsMIS=TC_all;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorrSTDvsMIS=TC_all_tuned;
Spatial_correlation.overlap.Rate_STDvsMIS=rate_all;
Spatial_correlation.overlap.Rate_tuned_STDvsMIS=rate_all_tuned;
Spatial_correlation.overlap.Peak_tuned_STDvsMIS=peak_all_tuned;



end