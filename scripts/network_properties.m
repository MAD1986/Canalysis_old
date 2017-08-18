%% Network Properties

function [Network_Properties] = network_properties(Event_Properties,C_df,Cdf_time, options)

binary=Event_Properties.onset_binary;
ones=Event_Properties.onset_ones;
onset_offset=EventsProperties.onset_offset; 
 
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
activity_hist=(sum(onset_ones,2)/size(onset_ones,2))*100;

%  To identify epochs of synchronous activity that included more active cells 
% than would be expected by chance at each frame, 
% we used interval reshuffling (randomly reordering of intervals between
% events for each cell), performed 1,000 times for each mouse in each context, 
% such that a surrogate histogram was constructed for each reshuffling. 
% The threshold percentage of active neurons corresponding to a significance level of P?<?0.05
% (appearing only in 5% of histograms) 
% was taken to be the per cent of coactive cells required in a single frame to be 
% considered a synchronous event, and this threshold ranged between 2.5% and 5% 
% active neurons per frame across all mice and fields of view. 
% At least three consecutive frames with activity above 
% the significance threshold were required to be considered a synchronous event, 
% and all subsequent contiguous frames above this threshold were grouped together 
% into the same synchronous event.
 

% Shuffle onset offset
event_dur=Event_Properties.noNaN.duration;

Nshuffle=1000;
tic;
shuffle_onset_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1:Nshuffle)','un',0));
for i=1:size(event_dur,2)
shuffle_dur_idx{i}=cell2mat(arrayfun(@(x) randperm((size(event_dur{i},2)),size(event_dur{i},2)),(1:Nshuffle)','un',0));
end
for S=1:Nshuffle
binary_shuffle=binary(shuffle_onset_idx(S,:),:);
for u=1:size(binary,2)
duration_shuffle{S}{u}=event_dur{u}(shuffle_dur_idx{u}(S,:));
onset_shuffle{S}{u}=find(binary_shuffle(:,u)==1);
onset_offset_shuffle{S}{u}=[ onset_shuffle{S}{u} onset_shuffle{S}{u}+(duration_shuffle{S}{u})'];
end
end
toc;

% Make shuffle binary / ones 
% writing ones and binary is too long

T_binary=zeros(size(C_df,1),size(C_df,2));
T_ones=binary;
for S=1:Nshuffle
for i=1:size(onset_offset,2)
%T_binary(onset_offset_shuffle{S}{i}(:,1),i)=1;
for ii=1:size(onset_offset{i},1)
T_ones(onset_offset_shuffle{S}{i}(ii,1):onset_offset_shuffle{S}{i}(ii,2),i)=1;
S_activity_hist(S,:)=(sum(T_ones,2)/size(T_ones,2))*100;
end
end
end
toc;

  


end


