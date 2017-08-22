
function [Place_cell]=spatial_info(Behavior, Events, Imaging, options);
%% Import data
res_position=Behavior.resampled.position;
time=Behavior.resampled.time;

%Calcium trace
if Imaging.options.msbackadj== true && Events.options.restricted==true
C_df=Imaging.trace_restricted_baselinesub;  
elseif Imaging.options.msbackadj== true && Events.options.restricted==false
C_df=Imaging.trace_baselinesub;
elseif Imaging.options.msbackadj== false && Events.options.restricted==true
C_df=Imaging.trace_restricted;
elseif Imaging.options.msbackadj== false && Events.options.restricted==false
C_df=Imaging.trace; 
end

if Events.options.restricted==1,
lap_start_stop= Behavior.restricted.lap;
elseif Events.options.restricted==0,
lap_start_stop= Behavior.lap;
end
run_ones=Behavior.run_ones;
runtime=Behavior.resampled.run_time;
run_onset_offset=Events.Run.run_onset_offset;
run_onset_binary=Events.Run.run_onset_binary;
onset_offset=Events.onset_offset;

%run_binary=Behavior.runbinary;

run_position=Behavior.resampled.run_position ; 

%% Bin data
%Time and position for running epochs ONLY:
%Bin normalized unning position
run_position_norm=(((run_position - min(run_position)) / (( max(run_position) - min(run_position)))));
%Bin running position:
for i=1:length(options.Nbin)
[count_bin{i},edges{i},bin{i}]=histcounts(run_position_norm, options.Nbin(i));
end

% Nb of onset for each bins
% dF/F for each bins
%1)only run
run_onset_bin=run_onset_binary(run_ones==1,:);
%Remove from analysis cell if not enough events (set in options.minevents)
for n=1:size(run_onset_bin,2)
if sum(run_onset_bin(:,n))<options.minevents
run_onset_bin(:,n)==0;
end
end  
run_C_df=C_df(run_ones==1,:);
%2)value for each bin onset
for i=1:length(options.Nbin)
for n=1:size(run_onset_bin,2)
onset_bin{i}{n}=bin{i}(run_onset_bin(:,n)==1);
end
end
%2)value for each bin dF/F
for i=1:length(options.Nbin)
for binN=1:options.Nbin(i)
for n=1:size(run_onset_bin,2)
dF_bin{i}{binN}(:,n)=run_C_df(find(bin{i}==binN),n);
end
dF_map{i}(binN,:)=mean(dF_bin{i}{binN});
end
end
%Time spend in each bin (count_bin in frame * avg fr).
avg_fr=mean(diff(time));
for i=1:length(options.Nbin)
occupancy_time_s{i}=count_bin{i}*avg_fr; % time spend (in sec)
total_occupancy=sum(count_bin{i})*avg_fr; %total time
proba_bin{i}=count_bin{i}/sum(count_bin{i}); % time spend / total time
end

%% Rate maps : total number of onset that occurred in a location
%bin divided by the time the animal spent in that bin
% onset map / occupancy
%Smooth rate map:
%Onset counts (onset map) and occupancy times in each bin were are independently smoothed 
%by convolving with a Gaussian smoothing kernel
% kernel = 3 bins Danielson et al.
% try other filter...
%Smooth dF/F (Dombeck et al. 2010 used moving average span=3bins)
for i=1:length(options.Nbin)
for n=1:size(run_onset_bin,2)
for binN=1:options.Nbin(i)
onset_map{i}(binN,n)=numel(find(onset_bin{i}{n}==binN));
end
onset_map_sm{i}(:,n)=imgaussfilt(onset_map{i}(:,n),options.sigma_filter);
dF_map_sm{i}(:,n)=smooth((dF_map{i}(:,n)),options.smooth_span);
end
occupancy_time_s_sm{i}=imgaussfilt(occupancy_time_s{i},options.sigma_filter);
rate_map{i}=onset_map{i}./occupancy_time_s{i}';
rate_map_sm{i}=onset_map_sm{i}./occupancy_time_s_sm{i}';
end

%Overall firing rate
%Mean firing rate
for i=1:length(options.Nbin)
overall_rate{i}=sum(onset_map_sm{i})/sum(occupancy_time_s_sm{i});
mean_rate{i}=mean(rate_map_sm{i});
end

%Normalize rate map and dF map
for i=1:length(options.Nbin)
for n=1:size(rate_map_sm{i},2)
norm_rate_map_sm{i}(:,n)=(rate_map_sm{i}(:,n)-min(rate_map_sm{i}(:,n)))/(max(rate_map_sm{i}(:,n))-min(rate_map_sm{i}(:,n)));
norm_dF_map_sm{i}(:,n)=(dF_map_sm{i}(:,n)-min(dF_map_sm{i}(:,n)))/(max(dF_map_sm{i}(:,n))-min(dF_map_sm{i}(:,n)));
end
end
%Spatial tuning curve (=nornalized rate map)
bin_STC=find(options.bin_spatial_tuning==options.Nbin);
if isempty(bin_STC)
    disp('change value of "options.bin_spatial_tuning" to a value present in "options.Nbin"')
return
end
    ST_curve=norm_rate_map_sm{bin_STC};
ST_dF_curve=norm_dF_map_sm{bin_STC};

%Sort spatial tuning curve by max rate
[M,I]=max(ST_curve);
ord_STC=[(1:size(ST_curve,2))' I' ST_curve'];
ROI_ord_STC=sortrows(ord_STC,2);
STC_sorted=(ROI_ord_STC(:,3:end));
%Sort mean dF/F by max rate
ord_STC_dF=[(1:size(ST_dF_curve,2))' I' ST_dF_curve'];
ROI_ord_STC_dF=sortrows(ord_STC_dF,2);
STC_dF_sorted=(ROI_ord_STC_dF(:,3:end));
%Remove Nan
STC_sorted_nonan=STC_sorted(~any(isnan(STC_sorted),2),:);
if options.dispfig==1
figure; 
subplot(1,2,1)
imagesc(STC_sorted_nonan); colormap('Jet');
title('Spatial Tuning Curve')
xlabel('bin position') % x-axis label
ylabel('neuron number') % y-axis label
subplot(1,2,2)
imagesc(STC_dF_sorted); colormap('Jet');
title('Mean dF/F')
xlabel('bin position') % x-axis label
clrbar = colorbar;
clrbar.Label.String=({'Normalized Ca2+ transient rate',' / ','Normalized mean dF/F'});
end

ROInb_STC_sorted=ROI_ord_STC(:,1);

%% Spatial specificity 
for i=1:length(options.Nbin)
    %(Danielson et al. 2016)
rate_ratio{i}=(rate_map_sm{i}./mean_rate{i});
spatial_info_bin{i}=rate_map_sm{i}.*log(rate_ratio{i}).*proba_bin{i}';
SI(i,:)=nansum(spatial_info_bin{i});
%(Skaggs)
spatial_info_bin_S{i}=rate_ratio{i}.*log2(rate_ratio{i}).*proba_bin{i}';
SIS(i,:)=nansum(spatial_info_bin_S{i});
end




%% Save structure

Place_cell.options=options;
Place_cell.proba_per_bin=proba_bin;
Place_cell.Spatial_tuning_curve=ST_curve;
Place_cell.Spatial_tuning_dF=ST_dF_curve;
Place_cell.Spatial_tuning_curve_sorted=STC_sorted;
Place_cell.Spatial_tuning_dF_sorted=STC_dF_sorted;
Place_cell.ROI_Spatial_tuning_sorted=ROInb_STC_sorted;

Place_cell.Spatial_Info.event_map=onset_map;
Place_cell.Spatial_Info.occupancy_map=occupancy_time_s;
Place_cell.Spatial_Info.rate_map=rate_map;
Place_cell.Spatial_Info.overall_rate=mean_rate;
Place_cell.Spatial_Info.proba_bin=proba_bin;
Place_cell.Spatial_Info.Spatial_Info=SI;
Place_cell.Spatial_Info.Spatial_Info_Skaggs=SIS;
Place_cell.Spatial_Info.Run_onset_bin=run_onset_bin;

Place_cell.Bin=bin;
%Spatialinfo.timebin=time_bin_sec;
%Spatialinfo.positionbin=position_bin;
%Spatialinfo.onset_bin=onset_bin;

end