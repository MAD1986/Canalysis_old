
function [Place_cell]=shuffle_spatial_info(Place_cell,options);

onset_map=Place_cell.Spatial_Info.event_map;
overall_rate=Place_cell.Spatial_Info.overall_rate;
proba_bin=Place_cell.Spatial_Info.proba_bin; 
run_onset_binary=Place_cell.Spatial_Info.Run_onset_bin;
occupancy=Place_cell.Spatial_Info.occupancy_map;
bin=Place_cell.Bin;
Nshuffle=options.Nshuffle;
Nbin=options.Nbin;

%% Shuffle binary map
tic
binary_shuffle_idx=cell2mat(arrayfun(@(x) randperm((size(run_onset_binary,1)),size(run_onset_binary,1)),(1:Nshuffle)','un',0));
disp(['Starting shuffle ', num2str(Nshuffle), 'X'])
for S=1:Nshuffle
binary_shuffle=run_onset_binary(binary_shuffle_idx(S,:),:);
for i=1:length(Nbin)
for n=1:size(run_onset_binary,2)
onset_bin_shuffle{i}{n}=bin{i}(binary_shuffle(:,n)==1);
for binN=1:options.Nbin(i)
onset_map_shuffle{S}{i}(binN,n)=numel(find(onset_bin_shuffle{i}{n}==binN));
end
end
end
end
disp(['End shuffle '])
toc;


%% Shuffle onset_map OLD
%idx=cell2mat(arrayfun(@(x) randperm((Nbin),size(onset_map,1)),(1:Nshuffle)','un',0));
%for u=1:size(onset_map,2)
%    for uu=1:Nshuffle
%onset_shuffle{uu}(:,u)=onset_map((idx(uu,:)),u);
%    end;end

%% Rate maps : total number of onset that occurred in a location
for S=1:Nshuffle
for i=1:length(options.Nbin)
for n=1:size(run_onset_binary,2)
onset_map_shuffle_sm{i}(:,n)=imgaussfilt(onset_map_shuffle{S}{i}(:,n),options.sigma_filter);
end
rate_map_shuffle_sm{S}{i}=onset_map_shuffle_sm{i}./occupancy{i}';
end
end


%% Spatial specificity (Danielson et al. 2016);



for S=1:Nshuffle
for i=1:length(options.Nbin)
    %(Danielson et al. 2016)
rate_ratio_shuffle{S}{i}=(rate_map_shuffle_sm{S}{i}./mean(rate_map_shuffle_sm{S}{i}));
spatial_info_bin_shuffle{S}{i}=rate_map_shuffle_sm{S}{i}.*log(rate_ratio_shuffle{S}{i}).*proba_bin{i}';
SI_shuffle{S}(i,:)=nansum(spatial_info_bin_shuffle{S}{i});
%(Skaggs)
spatial_info_bin_S_shuffle{S}{i}=rate_ratio_shuffle{S}{i}.*log2(rate_ratio_shuffle{S}{i}).*proba_bin{i}';
SIS_shuffle{S}(i,:)=nansum(spatial_info_bin_S_shuffle{S}{i});
end
end

Place_cell.Spatial_Info.Shuffle.spatial_info=SI_shuffle;
Place_cell.Spatial_Info.Shuffle.spatial_info_Skaggs=SIS_shuffle;

end






