function [SpatiallyTunedCells]=tuned_cells(SpatiallyTunedCells,Spatialinfo,TuningSpecificity);

spa_info=SpatiallyTunedCells.Spatial_info.significantROI;
tuning_spe=SpatiallyTunedCells.TuningSpecificity.significantROI;

SpatiallyTunedCells.Spatial_info.nb_spatial_info_cells=sum(spa_info);
SpatiallyTunedCells.TuningSpecificity.nb_tuning_spe_cells=sum(tuning_spe);

while length(spa_info)>length(tuning_spe),
tuning_spe=[tuning_spe 0];
end
while length(spa_info)<length(tuning_spe),
spa_info=[spa_info 0];
end


tunedcells=(spa_info+tuning_spe==2);
nb_tunedcells=sum(tunedcells);

%plot rate map and tuning vectors of tuned cells
rate_map=Spatialinfo{8}.rate_map.normalized_rate_map_smoothed;
tuning_vector=TuningSpecificity.tuning_vector;
tuning_vector_spe=TuningSpecificity.tuning_vector_specificity;

for u=1:size(rate_map,2);
if tunedcells(u)==1,
rate_map_tuned(:,u)=rate_map(:,u);
tuning_vector_tuned(u)=tuning_vector(u);
tuning_vector_spe_tuned(u)=tuning_vector_spe(u);
figure; compass(tuning_vector_tuned{u}); hold on; compass(tuning_vector_spe_tuned(u),'r');
end;end
if sum(tunedcells)==0,
rate_map_tuned=NaN;
tuning_vector_tuned=NaN;
tuning_vector_spe_tuned=NaN;
end

rate_map_tuned=rate_map_tuned(:,any(rate_map_tuned));


[M,I]=max(rate_map_tuned);
ordered_rate_map=[(1:size(rate_map_tuned,2))' I' rate_map_tuned'];
ROI_odered_rate_map=sortrows(ordered_rate_map,2);
ordered_rate_map=(ROI_odered_rate_map(:,3:end));
%Remove Nan
figure; imagesc(ordered_rate_map); colormap('Jet');
title('Rate map -Tuned neurons')
xlabel('bin position') % x-axis label
ylabel('neuron number') % y-axis label
ROI_Tuned_rate_map=ROI_odered_rate_map(:,1);


SpatiallyTunedCells.TunedROI =tunedcells;
SpatiallyTunedCells.NbROI=nb_tunedcells;
SpatiallyTunedCells.Tuned_rate_map=rate_map_tuned;
SpatiallyTunedCells.Ordered_tuned_rate_map=ordered_rate_map;
SpatiallyTunedCells.Ordered_tuned_rate_mapROI=ROI_Tuned_rate_map;
SpatiallyTunedCells.Tuned_tuningvector =tuning_vector_tuned;
SpatiallyTunedCells.Tuned_tuningspe =tuning_vector_spe_tuned;
end