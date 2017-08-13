function [Spatial_correlation] = plotPVcorr_V2(Spatial_correlation,color, plotPV);
sessions=Spatial_correlation.options.sessions;
PVcorr=Spatial_correlation.PVcorrelation.AllROI.PVcorr;
% Plot matrix
for i=1:size(plotPV,2)
figure; 
imagesc(Spatial_correlation.PVcorrelation.matrix{plotPV{i}(1)}{plotPV{i}(end)}); colormap(color);
end

%Box plot
for i=1:size(plotPV,2)
boxplotgraph(:,i)=(cell2mat(PVcorr{plotPV{i}(1)}{plotPV{i}(end)}))';
[f{i},x{i}] = ecdf(boxplotgraph(:,i));
labelsbox{i}=['session' num2str(plotPV{i}(1)) 'VS' num2str(plotPV{i}(end))];
end
figure;
boxplot(boxplotgraph,'Labels', labelsbox)
title('mean PV correlation')
ylabel('pearson correlation score') % y-axis label

%Cumulative distribution
figure;
for i=1:size(plotPV,2)
 hold on;  
 plot(x{i},f{i},'DisplayName',['session' num2str(plotPV{i}(1)) 'VS' num2str(plotPV{i}(end))]);
%hold off;
title('mean PV correlation - all ROI')
ylabel('cumulative frequency')
xlabel('Pearson correlation score')
 %legend(strcat('session',num2str(plotPV{i}(1)), 'VS',num2str(plotPV{i}(end))));
 end
legend(gca,'show')

Spatial_correlation.PVcorrelation.AllROI.PVcorr_table=boxplotgraph;


end



