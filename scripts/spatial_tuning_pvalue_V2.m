function [SpatiallyTunedCells]=spatial_tuning_pvalue_V2(Spatialinfo,Nbin,Nshuffle, pvalue);

for u=1:length(Nbin);
spatial_tuning_shuffle_allbin{u}=Spatialinfo{u}.shuffle.spatial_info;
%Compute mean of null distribution for each bin
%rows=ROI; columns=bin nb
mean_spatial_tuning_shuffle(:,u)=mean(spatial_tuning_shuffle_allbin{u});
spatial_tuning_allbin(:,u)=Spatialinfo{u}.spatial_info;

for uu=1:size(spatial_tuning_allbin,1);
spatial_tuning_shuffle_bin{uu}(:,u)=spatial_tuning_shuffle_allbin{u}(:,uu);
spatial_tuning_shuffle_bin_corrected{uu}=spatial_tuning_shuffle_bin{uu}-(mean_spatial_tuning_shuffle(uu,:));
%info_content_shuffle(:,uu)=max(spatial_tuning_shuffle_bin{uu},[],2); 
%Total mean of all the null distribution for all the bins 
%mean2_spatial_tuning_shuffle(:,uu)=mean2(spatial_tuning_shuffle_bin{uu}(~isnan(spatial_tuning_shuffle_bin{uu})));
%substract the mean of null distribution from shuffle
%spatial_tuning_shuffle_bin_corrected{uu}=spatial_tuning_shuffle_bin{uu}-mean2_spatial_tuning_shuffle(uu);
info_content_shuffle(:,uu)=max(spatial_tuning_shuffle_bin_corrected{uu},[],2);
% substract the mean of null distribution from shuffle
%spatial_tuning_shuffle_bin_corrected{uu}=spatial_tuning_shuffle_bin{uu}-mean2_spatial_tuning_shuffle(uu);
%info_content_shuffle(:,uu)=max(spatial_tuning_shuffle_bin_corrected{uu},[],2);
 end;end


spatial_tuning_allbin_corrected=spatial_tuning_allbin - mean_spatial_tuning_shuffle;
info_content=max(spatial_tuning_allbin_corrected,[],2);

for uu=1:size(spatial_tuning_allbin,1)
spatialtuning_sig(:,uu)=info_content_shuffle(:,uu)>=info_content(uu);
end
spatialtuning_pvalue=sum(spatialtuning_sig)./Nshuffle;

for uu=1:size(spatial_tuning_allbin,1)
if isnan(info_content(uu)),
spatialtuning_pvalue(uu)=NaN;
end;end
spatialtuning_signROI=spatialtuning_pvalue<=pvalue;


SpatiallyTunedCells.Spatial_info.info_content=info_content;
%SpatiallyTunedCells.Spatial_info.info_content_shuffle=info_content_shuffle;
SpatiallyTunedCells.Spatial_info.pvalue=spatialtuning_pvalue;
SpatiallyTunedCells.Spatial_info.significantROI=spatialtuning_signROI;

end
