
function [Network_Properties] = network_properties(Event_Properties,C_df,Cdf_time, options)

%% Import data
binary=Event_Properties.onset_binary_analysed;
onset_ones=Event_Properties.onset_ones_analysed;
onset_offset=Event_Properties.onset_offset_analysed; 
 
%% Frequency (event/min)
nb_events=Event_Properties.nb_events; 
%recording duration : start to end in min
rec_dur=(Cdf_time(end)-Cdf_time(1))/60;
event_freq=nb_events/rec_dur;

%% Activity rate (AUC/min) 
% AUC rate = AUC divided by recording duration
AUC=Event_Properties.AUC;
for i=1:size(AUC,2)
meanAUC(i)=nanmean(AUC{i});
end
AUC_rate=meanAUC/rec_dur;

%% Correlation between cell
% See Rajasethupathy et al. 2015
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4825678
%Do Pearson correlation removing inter events period (dF/F during inter
%event periods = 0)
%To avoid accumulation in correlated signal due to slow drifts 
%we set to 0 all ?F/F values lying outside the window of a significant transient 
C_df_sub=C_df;
C_df_sub(onset_ones==0)=0;
[dfCCorr_sub, dfCCorrPvalues_sub] = corrcoef(C_df_sub);
%finding the number of correlated neurons 
%with which the Pearsons correlation coefficient
%is above 0.3.
  for i=1:size(dfCCorr_sub,2);
 corr_pairs{i}=find(dfCCorr_sub(:,i)>0.3);
 %Remove pair between same neuron
 corr_pairs{i}= corr_pairs{i}(corr_pairs{i}~=i);
  nb_corr_pairs(i)=length(corr_pairs{i});
  end
 %HC neurons were defined as those neurons that had more 
 %correlated partners than 
 %that of the average neuron in the same volume 
 %by >1 standard deviation
 HC_neurons= find(nb_corr_pairs>std(nb_corr_pairs));
%activity histograms :
%percentage of active cells as a function of time.
%% Synchronous activity 
%  To identify epochs of synchronous activity that included more active cells 
% than would be expected by chance at each frame, 
% reshuffling performed 1,000 times, 

% sum of active neurons for each frame / total number of neurones * 100
activity_hist=(sum(onset_ones,2)/size(onset_ones,2))*100;


% Shuffle onset offset
% import duration for offset
event_dur=Event_Properties.noNaN.duration;
Nshuffle=options.Nshuffle;
disp(['Starting shuffle ', num2str(Nshuffle), 'X'])
shuffle_onset_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1:Nshuffle)','un',0));
for i=1:size(event_dur,2)
shuffle_dur_idx{i}=cell2mat(arrayfun(@(x) randperm((size(event_dur{i},2)),size(event_dur{i},2)),(1:Nshuffle)','un',0));
end
for S=1:Nshuffle
binary_shuffle=binary(shuffle_onset_idx(S,:),:);

for u=1:size(event_dur,2)
duration_shuffle{S}{u}=event_dur{u}(shuffle_dur_idx{u}(S,:));
onset_shuffle{S}{u}=find(binary_shuffle(:,u)==1);
offset_shuffle{S}{u}=onset_shuffle{S}{u}+(duration_shuffle{S}{u})';

for uu=1:size(offset_shuffle{S}{u},1)
if offset_shuffle{S}{u}(uu)>size(C_df,1)
     offset_shuffle{S}{u}(uu)=size(C_df,1);
end
end
onset_offset_shuffle{S}{u}=[ onset_shuffle{S}{u}  offset_shuffle{S}{u}];
end
end

T_onset_ones=zeros(size(C_df,1),size(C_df,2));
for S=1:Nshuffle
Temp_onset_ones=T_onset_ones;
for i=1:size(onset_offset,2)
for ii=1:size(onset_offset{i},1)
Temp_onset_ones(onset_offset_shuffle{S}{i}(ii,1):onset_offset_shuffle{S}{i}(ii,2),i)=1;
end
end
S_activity_hist(S,:)=(sum(Temp_onset_ones,2)/size(Temp_onset_ones,2))*100;
end

disp('End shuffle')
 

%p-value was defined as the fraction of shuffle distribution that 
%exceeded the frame true number of co-active cells
sync_act_pvalue=(sum(S_activity_hist>activity_hist'))/Nshuffle;
sign_frames_idx=sync_act_pvalue<=options.pmin;
fr_nb=1:size(C_df,1);
sign_frames=fr_nb(sign_frames_idx);
    
% nb of frames between each significant frames
dist_sign_fr=diff(sign_frames)';
%Find consecutive synchronous epochs separated by more than 1 frames
% find non consecutive events : if differences (diff) > 1
epochs_end_idx=find(dist_sign_fr>1);
% end of epoch when separated by more than 1 frame 
epochs_end_fr=sign_frames(epochs_end_idx);
% Make start stop
if epochs_end_idx(end)<=length(sign_frames)
epochs_start_idx=[sign_frames(1);epochs_end_idx+1];
elseif epochs_end_idx(end)>=length(sign_frames)
epochs_start_idx=[sign_frames(1);epochs_end_idx];
end
epochs_end_idx=[epochs_end_idx ;length(sign_frames)];
syn_epochs=[epochs_start_idx epochs_end_idx];
syn_epochs_frames=sign_frames(syn_epochs);
% Duration of epoch
dur_epochs=syn_epochs_frames(:,2)-syn_epochs_frames(:,1);

% At least three consecutive frames with activity above 
% the significance threshold were required to be considered a synchronous event, 
% and all subsequent contiguous frames above this threshold were grouped together 
% into the same synchronous event.
for i=1:size(dur_epochs,1)
if dur_epochs(i)<options.minframes ==1
   syn_epochs_frames(i,:)=NaN;
   syn_epochs(i,:)=NaN;
end  
end
syn_epochs_frames= syn_epochs_frames(~isnan( syn_epochs_frames(:,2)),:);
syn_epochs=syn_epochs(~isnan(syn_epochs(:,2)),:);
nb_syn_epochs=size(syn_epochs,1);
%Number of co-active neurons in significant frames
nb_CC_sig_fr=activity_hist(sign_frames_idx);
co_act_fr=[sign_frames' nb_CC_sig_fr];

%% Make structure
Network_Properties.frequency=event_freq;
Network_Properties.AUC_rate=AUC_rate;
Network_Properties.CCorr=dfCCorr_sub;
Network_Properties.corr_pairs=corr_pairs;
Network_Properties.nb_corr_pairs=nb_corr_pairs;
Network_Properties.HC_neurons=HC_neurons;
Network_Properties.activity_hist=activity_hist;
Network_Properties.sign_activity_frames=sign_frames;
Network_Properties.nb_cells_sign_activity_frames =co_act_fr;
Network_Properties.synchronous_epochs=syn_epochs;
Network_Properties.nb_synchronous_epochs=nb_syn_epochs;
 
 

end


