
function [Place_cell]=shuffle_place_cell(Place_cell,Behavior,Events,options);
onset_map=Place_cell.Spatial_Info.event_map;
overall_rate=Place_cell.Spatial_Info.overall_rate;
proba_bin=Place_cell.Spatial_Info.proba_bin; 
run_onset_binary=Place_cell.Spatial_Info.Run_onset_bin;
occupancy=Place_cell.Spatial_Info.occupancy_map;
bin=Place_cell.Bin;
Nshuffle=options.Nshuffle;
Nbin=options.Nbin;
sigma=options.sigma_filter;
%% Shuffle onset map and tuning specificity
tic;
onset_binary=Events.Run.run_onset_binary;
disp(['Starting shuffle ', num2str(Nshuffle), 'X'])
parfor S=1:Nshuffle
    %[onset_map_shuffle{S}]=shuffle_map(run_onset_binary,bin,Nshuffle,Nbin);
[onset_map_shuffle{S},shuffle_tuning_specificity{S}]=shuffle_script(onset_binary,Place_cell,Behavior,Events,options);
end
disp(['End shuffle '])
toc;
%% Rate maps : total number of onset that occurred in a location
parfor S=1:Nshuffle
[onset_map_shuffle_sm{S}]=gauss_filt(onset_map_shuffle{S},Nbin,sigma)
end
for S=1:Nshuffle
for i=1:length(options.Nbin)
rate_map_shuffle_sm{S}{i}=onset_map_shuffle_sm{S}{i}./occupancy{i}';
end
end

%% Spatial specificity null distribution (Danielson et al. 2016);
for S=1:Nshuffle
for nbin=1:length(options.Nbin)
    %(Danielson et al. 2016)
rate_ratio_shuffle{S}{nbin}=(rate_map_shuffle_sm{S}{nbin}./mean(rate_map_shuffle_sm{S}{nbin}));
spatial_info_bin_shuffle{S}{nbin}=rate_map_shuffle_sm{S}{nbin}.*log(rate_ratio_shuffle{S}{nbin}).*proba_bin{nbin}';
SI_shuffle{S}(nbin,:)=nansum(spatial_info_bin_shuffle{S}{nbin});
%(Skaggs)
spatial_info_bin_S_shuffle{S}{nbin}=rate_ratio_shuffle{S}{nbin}.*log2(rate_ratio_shuffle{S}{nbin}).*proba_bin{nbin}';
SIS_shuffle{S}(nbin,:)=nansum(spatial_info_bin_S_shuffle{S}{nbin});
end
end
Place_cell.Spatial_Info.Shuffle.spatial_info=SI_shuffle;
Place_cell.Spatial_Info.Shuffle.spatial_info_Skaggs=SIS_shuffle;

Place_cell.Tuning_Specificity.Shuffle.tuning_specificity=shuffle_tuning_specificity;

%% tuning specificity null distribution (Danielson et al. 2016);
%tic;
%run_onset_binary=Events.Run.run_onset_binary;
%disp(['Starting shuffle ', num2str(Nshuffle), 'X'])
%parfor S=1:Nshuffle
%[Tuning_specificity_shuffle{S}]=shuffle_tuning_specificity(run_onset_binary,Place_cell,Behavior,Events,options)
%end
%disp(['End shuffle '])
%toc;

end






