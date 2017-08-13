function [Spatial_correlation] = plotrateoverlap(Spatial_correlation,color, plotrate);
sessions=Spatial_correlation.options.sessions;
rate_overlap=Spatial_correlation.overlap.Rate_overlap;
peak_overlap=Spatial_correlation.overlap.Peak_overlap;
rate_overlap_tuned=Spatial_correlation.overlap.Rate_overlap_tuned;
peak_overlap_tuned=Spatial_correlation.overlap.Peak_overlap_tuned;
rate_ratio=Spatial_correlation.overlap.Rate_ratio;
peak_ratio=Spatial_correlation.overlap.Peak_ratio;
rate_ratio_tuned=Spatial_correlation.overlap.Rate_ratio_tuned;
peak_ratio_tuned=Spatial_correlation.overlap.Peak_ratio_tuned;



%% ALL ROI
%Box plot
for i=1:size(plotrate,2)
boxplotgraph(:,i)=(rate_ratio{plotrate{i}(1)}{plotrate{i}(end)})';
[f{i},x{i}] = ecdf(boxplotgraph(:,i));
labelsbox{i}=['session' num2str(plotrate{i}(1)) 'VS' num2str(plotrate{i}(end))];
end
figure;
boxplot(boxplotgraph,'Labels', labelsbox )
title('rate overlap - all ROI')
ylabel('rate ratio (min/max)') % y-axis label

%Cumulative distribution
figure;
for i=1:size(plotrate,2)
hold on;   
plot(x{i},f{i},'DisplayName',['session' num2str(plotrate{i}(1)) 'VS' num2str(plotrate{i}(end))]);
title('rate overlap - all ROI')
ylabel('cumulative frequency')
xlabel('rate ratio (min/max)')
 %hold off;
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
 end
legend(gca,'show')
 
 %% Tuned ROI 
 %Box plot rate overlap
for i=1:size(plotrate,2)
boxplotgraph_tuned(:,i)=(rate_ratio_tuned{plotrate{i}(1)}{plotrate{i}(end)})';
end
boxplotgraph_tuned(~any(boxplotgraph_tuned,2),:)=[]; 
for i=1:size(plotrate,2)
[f_tuned{i},x_tuned{i}] = ecdf(boxplotgraph_tuned(:,i));
labelsbox_tuned{i}=['session' num2str(plotrate{i}(1)) 'VS' num2str(plotrate{i}(end))];
end
figure;
boxplot(boxplotgraph_tuned,'Labels', labelsbox_tuned )
title('rate overlap - Tuned ROI')
ylabel('rate ratio (min/max)')

 %Cumulative distribution rate overlap
figure;
for i=1:size(plotrate,2)
hold on;   
plot(x_tuned{i},f_tuned{i},'DisplayName',['session' num2str(plotrate{i}(1)) 'VS' num2str(plotrate{i}(end))]);
title('rate overlap - Tuned ROI')
ylabel('cumulative frequency')
xlabel('rate ratio (min/max)')
 %hold off;
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
 end
legend(gca,'show')


%Box plot peak overlap
for i=1:size(plotrate,2)
boxplotgraph_peak_tuned(:,i)=(peak_ratio_tuned{plotrate{i}(1)}{plotrate{i}(end)})';
end
boxplotgraph_peak_tuned(~any(boxplotgraph_peak_tuned,2),:)=[]; 
for i=1:size(plotrate,2)
[f_peak_tuned{i},x_peak_tuned{i}] = ecdf(boxplotgraph_peak_tuned(:,i));
labelsbox_tuned{i}=['session' num2str(plotrate{i}(1)) 'VS' num2str(plotrate{i}(end))];
end
figure;
boxplot(boxplotgraph_peak_tuned,'Labels', labelsbox_tuned )
title('rate overlap - Tuned ROI')
ylabel('rate ratio (min/max)')


 %Cumulative distribution peak overlap
figure;
for i=1:size(plotrate,2)
hold on;   
plot(x_peak_tuned{i},f_peak_tuned{i},'DisplayName',['session' num2str(plotrate{i}(1)) 'VS' num2str(plotrate{i}(end))]);
title('rate overlap - Tuned ROI')
ylabel('cumulative frequency')
xlabel('rate ratio (min/max)')
 %hold off;
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
 end
legend(gca,'show')


Spatial_correlation.overlap.Rate_table=boxplotgraph;
Spatial_correlation.overlap.Rate_tuned_table=boxplotgraph_tuned;
Spatial_correlation.overlap.Peak_tuned_table=boxplotgraph_peak_tuned;




end
