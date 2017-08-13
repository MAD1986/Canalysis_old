function [Spatial_correlation] = plotTCcorr(Spatial_correlation,color, plotTC);
sessions=Spatial_correlation.options.sessions;
TCcorr=Spatial_correlation.TuningCurvecorr.AllROI.TCcorr;
TCcorr_tuned=Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned;
%% ALL ROI

%Box plot
for i=1:size(plotTC,2)
boxplotgraph(:,i)=(cell2mat(TCcorr{plotTC{i}(1)}{plotTC{i}(end)}))';
[f{i},x{i}] = ecdf(boxplotgraph(:,i));
labelsbox{i}=['session' num2str(plotTC{i}(1)) 'VS' num2str(plotTC{i}(end))];
end
figure;
boxplot(boxplotgraph,'Labels', labelsbox )
title('mean tunning curve correlation - all ROI')
ylabel('Pearson correlation score') % y-axis label

%Cumulative distribution
figure;
for i=1:size(plotTC,2)
hold on;   
plot(x{i},f{i},'DisplayName',['session' num2str(plotTC{i}(1)) 'VS' num2str(plotTC{i}(end))]);
title('mean tunning curve correlation - all ROI')
ylabel('cumulative frequency')
xlabel('Pearson correlation score')
 %hold off;
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
 end
legend(gca,'show')

 
 %% Tuned ROI 
 %Box plot

for i=1:size(plotTC,2)
boxplotgraph_tuned(:,i)=(cell2mat(TCcorr_tuned{plotTC{i}(1)}{plotTC{i}(end)}))';
[f_tuned{i},x_tuned{i}] = ecdf(boxplotgraph_tuned(:,i));
labelsbox_tuned{i}=['session' num2str(plotTC{i}(1)) 'VS' num2str(plotTC{i}(end))];
end
figure;
boxplot(boxplotgraph_tuned,'Labels', labelsbox_tuned )
title('mean tunning curve correlation - tuned ROI')
ylabel('Pearson correlation score') % y-axis label
 
 %Cumulative distribution
figure;
for i=1:size(plotTC,2)
hold on;   
plot(x_tuned{i},f_tuned{i},'DisplayName',['session' num2str(plotTC{i}(1)) 'VS' num2str(plotTC{i}(end))]);
title('mean tunning curve correlation - tuned ROI')
ylabel('cumulative frequency')
xlabel('Pearson correlation score')
 %hold off;
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
 end
legend(gca,'show')
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_table=boxplotgraph_tuned;
%end

Spatial_correlation.TuningCurvecorr.AllROI.TCcorr_table=boxplotgraph;




end
