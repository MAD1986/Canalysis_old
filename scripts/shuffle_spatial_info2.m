
function [Spatialinfo]=shuffle_spatial_info(Spatialinfo,Nbin,Nshuffle);

occupancy_map_smooth=Spatialinfo.occupancy_map.occupancy_map_smoothed;
onset_map=Spatialinfo.event_map.event_map;
overall_rate=Spatialinfo.overall_rate;
proba_bin=Spatialinfo.proba_per_bin; 


%% Shuffle onset_map
idx=cell2mat(arrayfun(@(x) randperm((Nbin),size(onset_map,1)),(1:Nshuffle)','un',0));
for u=1:size(onset_map,2)
    for uu=1:Nshuffle
onset_shuffle{uu}(:,u)=onset_map((idx(uu,:)),u);
    end;end

%% Rate maps : total number of onset that occurred in a location
%bin divided by the time the animal spent in that bin
for uu=1:Nshuffle
for u=1:size(onset_map,2)
%smooth rate map
onset_shuffle_smooth(:,u)=Smooth(onset_shuffle{uu}(:,u),Spatialinfo.options.smooth);
rate_map_shuffle{uu}=onset_shuffle_smooth./occupancy_map_smooth;
end;end

%% Spatial specificity (Danielson et al. 2016);
%?i * ln * (?i/?) * pi 
%Sum for number of bins = 2,4,5,8,10,20,25,100

%?i is the mean onset rate of a cell in the i-th bin, 
%? is the overall mean onset rate, 
%pi is the probability of the animal being in the i-th bin
for uu=1:Nshuffle
rate_ratio_shuffle=rate_map_shuffle{uu}./overall_rate;
spatial_info_shuffle_array{uu}=nansum(rate_map_shuffle{uu}.*log(rate_ratio_shuffle).*proba_bin);
end
spatial_info_shuffle=(reshape((cell2mat(spatial_info_shuffle_array)),size(onset_map,2),[]))';
spatial_info_shuffle(~spatial_info_shuffle)=nan;
Spatialinfo.shuffle.spatial_info= spatial_info_shuffle;

end






