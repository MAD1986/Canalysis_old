function [Spatial_correlation] = plotTC(Spatial_correlation,color, plotTC);
     
%Box plot
boxplotgraph=[TuningCurveCorr.Rvalues.any_SelectiveROI.AvsB; TuningCurveCorr.Rvalues.any_SelectiveROI.BvsC]
figure; boxplot(boxplotgraph')

%Plot rate map and compass of cells tuned for A and C
for i=tunedROI_ctx.ROInb.C_all'
figure; plot(rate_map_sm2{2}(:,i), 'color', 'r');
hold on; plot(rate_map_sm2{3}(:,i), 'color', 'g');
hold on; plot(rate_map_sm2{1}(:,i),  'color', 'b'); legend('B', 'C','A'); 
end
clear i
for ii=tunedROI_ctx.ROInb.A_all'
for iii=1:3;
spatial_select_AC(ii,iii)=spatial_selectivity{iii}(ii);
max_spatial_ABC(ii,:)=max(spatial_select_AC(ii,:));
end
end
for ii=tunedROI_ctx.ROInb.AandC'
figure; compass(max_spatial_ABC(ii));
hold on; compass(spatial_selectivity{1}((ii)),'b');
hold on;compass(spatial_selectivity{2}((ii)), 'r');
hold on;compass(spatial_selectivity{3}((ii)), 'g'); 
end