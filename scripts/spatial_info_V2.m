
function [Spatialinfo]=spatial_info(Behavior, Events, Imaging, Nbin);
time_position=[Behavior.resampled.time Behavior.resampled.position];
position_norm=(time_position(:,2) - min(time_position(:,2))) / ( max(time_position(:,2)) - min(time_position(:,2)) );
if Events.options.restricted==1,
lap_start_stop= Behavior.restricted.lap;
elseif Events.options.restricted==0,
lap_start_stop= Behavior.lap;
end
runtime=Behavior.runtime;
run_onset_offset=Events.RunningEpochs.run_onset_offset;
run_onset_binary=Events.RunningEpochs.run_onset_binary;
onset_offset=Events.onset_offset;
run_binary=Behavior.runbinary;

Behavior.placeoptions.bins=Nbin;
minevents=Behavior.placeoptions.minevents;

for u = 1:size(run_onset_offset,2);
    if size(run_onset_offset{u},1)<minevents,
        run_onset_offset{u}=NaN;
        run_onset_binary(:,u)=0;
        %run_binary(:,u)=0;
    end
end
    

%% Binning data
%Time and position for running epochs ONLY:
running_idx=find(run_binary==1);
run_time_position=time_position(running_idx,:);
run_position_norm=(run_time_position(:,2) - min(run_time_position(:,2))) / ( max(run_time_position(:,2)) - min(run_time_position(:,2)) );
run_onset_binary2=run_onset_binary(running_idx,:);
%Bin running position:
[bin_idx, run_binnum_position]=Bin(run_position_norm, Behavior.placeoptions.bins); 
%Bin nb for each onset:
for u=1:size(run_onset_binary2,2)
onset_bin_idx{u}=find(run_onset_binary2(:,u)==1);
onset_bin{u}=bin_idx(onset_bin_idx{u});
end
%Value of time for each bin
for u=1:Behavior.placeoptions.bins
time_bin{u}=find(bin_idx==u);
time_bin_sec{u}=run_time_position(time_bin{u},1);
position_bin{u}=run_time_position(time_bin{u},2);
%Time spend in each bin: (in second)
occupancy_map(u)=length(time_bin{u})*(1/Imaging.options.frimaging);
proba_bin(u)=length(time_bin{u})/length(bin_idx);
end

% Nb of onset for each bins
for u=1:Behavior.placeoptions.bins
for uu=1:size(onset_bin,2)
onset_map(u,uu)=numel(find(onset_bin{uu}==u));
end;end

%% Rate maps : total number of onset that occurred in a location
%bin divided by the time the animal spent in that bin
rate_map=onset_map./occupancy_map';
for i=1:size(rate_map,2)
 %smooth onset map
onset_map_smooth(:,i)=Smooth(onset_map(:,i),Behavior.placeoptions.smooth);
 %smooth occupancy
occupancy_map_smooth=Smooth(occupancy_map,Behavior.placeoptions.smooth);
end
rate_map_smooth=onset_map_smooth./occupancy_map_smooth;

%Overall mean firing rate
for u=1:size(rate_map,2)
overall_mean_rate(:,u)=mean(rate_map(:,u));
end
%OR
overall_rate=sum(onset_map)/sum(occupancy_map);

%probability of the animal being in the i-th bin 
% =(occupancy in the i-th bin/total recording time)
proba_bin=occupancy_map_smooth/sum(occupancy_map);


%Make rate map sorted by max rate
%Normalize
for i=1:size(rate_map_smooth,2)
normalized_rate_map_smooth(:,i)=(rate_map_smooth(:,i) - min(rate_map_smooth(:,i))) / ( max(rate_map_smooth(:,i)) - min(rate_map_smooth(:,i)) );
end

[M,I]=max(normalized_rate_map_smooth);
ordered_rate_map=[(1:size(normalized_rate_map_smooth,2))' I' normalized_rate_map_smooth'];
ROI_odered_rate_map=sortrows(ordered_rate_map,2);
ordered_rate_map=(ROI_odered_rate_map(:,3:end));
%Remove Nan
ordered_rate_mapNoNaN=ordered_rate_map(~any(isnan(ordered_rate_map),2),:);
figure; imagesc(ordered_rate_mapNoNaN); colormap('Jet');
title('Rate map -all neurons')
xlabel('bin position') % x-axis label
ylabel('neuron number') % y-axis label
ROI_odered_rate_map=ROI_odered_rate_map(:,1);


%% Spatial specificity (Danielson et al. 2016);
%?i * ln * (?i/?) * pi 
%Sum for number of bins = 2,4,5,8,10,20,25,100

%?i is the mean onset rate of a cell in the i-th bin, 
%? is the overall mean onset rate, 
%pi is the probability of the animal being in the i-th bin
rate_ratio=rate_map_smooth./overall_rate;
spatial_info=nansum(rate_map_smooth.*log(rate_ratio).*proba_bin);
%Zeros = NaN 
spatial_info(~spatial_info)=nan;

%TO SAVE
Spatialinfo.options=Behavior.placeoptions;
Spatialinfo.proba_per_bin=proba_bin;
Spatialinfo.event_map.event_map= onset_map;
Spatialinfo.event_map.event_map_smoothed=onset_map_smooth;
Spatialinfo.occupancy_map.occupancy_map=occupancy_map;
Spatialinfo.occupancy_map.occupancy_map_smoothed=occupancy_map_smooth;
Spatialinfo.overall_rate=overall_rate;
Spatialinfo.rate_map.rate_map=rate_map;
Spatialinfo.rate_map.rate_map_smoothed= rate_map_smooth;
Spatialinfo.rate_map.normalized_rate_map_smoothed=normalized_rate_map_smooth;
Spatialinfo.rate_map.ordered_rate_map= ordered_rate_map;
Spatialinfo.rate_map.ROI_odered_rate_map= ROI_odered_rate_map;
Spatialinfo.spatial_info= spatial_info;
Spatialinfo.timebin=time_bin_sec;
Spatialinfo.positionbin=position_bin;
Spatialinfo.onset_bin=onset_bin;

end