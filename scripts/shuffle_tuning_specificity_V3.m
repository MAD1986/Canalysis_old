
function [TuningSpecificity,SpatiallyTunedCells]=shuffle_tuning_specificity(TuningSpecificity,Behavior,Events,SpatiallyTunedCells,Nshuffle,pvalue);
%% Extract the data
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
minevents=Behavior.placeoptions.minevents;
for u = 1:size(run_onset_offset,2);
    if size(run_onset_offset{u},1)<minevents,
        run_onset_offset{u}=NaN;
        run_onset_binary(:,u)=0;
      %  run_binary(:,u)=0;
    end
end
 
%Time and position for running epochs ONLY:
running_idx=find(run_binary==1);
run_time_position=time_position(running_idx,:);
run_position_norm=(run_time_position(:,2) - min(run_time_position(:,2))) / ( max(run_time_position(:,2)) - min(run_time_position(:,2)) );
run_onset_binary2=run_onset_binary(running_idx,:);



%% Shuffle onset_binary
idx=cell2mat(arrayfun(@(x) randperm((size(run_onset_binary2,1)),size(run_onset_binary2,1)),(1:Nshuffle)','un',0));
for S=1:Nshuffle
run_onset_binary_shuffle=run_onset_binary2(idx(S,:),:);
%position of the mouse at onset time:
for u=1:size(run_onset_binary2,2)
onset_idx{S}{u}=find(run_onset_binary_shuffle(:,u)==1);
onset_position{S}{u}=run_position_norm(onset_idx{S}{u});
end;end

%Fraction of frames acquired at position of onset
for S=1:Nshuffle
for u=1:size(onset_position{S},2)
    for uu=1:size(onset_position{S}{u},1)
         if isempty(onset_position{S}{u})==0,     
frame_onset{S}{u}{uu}=find(run_position_norm==onset_position{S}{u}(uu));
nbframe_onset{S}{u}(:,uu)= size(frame_onset{S}{u}{uu},1);
fraction_frames{S}{u}(:,uu)= nbframe_onset{S}{u}(:,uu)/ (size(run_onset_binary2,1));
end  
end;end;end
%% Compute spatial tuning vector
for S=1:Nshuffle
for u=1:size(onset_position{S},2) 
    for uu=1:size(onset_position{S}{u},1)
 if isempty(frame_onset{S}{u})==0,        
spatial_tuning_vector{S}{u}(uu)=(exp(i*onset_position{S}{u}(uu)*2*pi))./(fraction_frames{S}{u}(:,uu));       
spatial_tuning_vector_angle{S}{u}(uu)=angle(spatial_tuning_vector{S}{u}(uu));
spatial_tuning_vector_magnitude{S}{u}(uu)=abs(spatial_tuning_vector{S}{u}(uu));
 elseif isempty(frame_onset{S}{u}),
spatial_tuning_vector{S}{u}(uu)=NaN; 
spatial_tuning_vector_angle{S}{u}(uu)=NaN;
spatial_tuning_vector_magnitude{S}{u}(uu)=NaN;
 end
 end;end;end
%Normalize magnitude
for S=1:Nshuffle
for u=1:size(spatial_tuning_vector{S},2)   
spatial_tuning_vector_magnitude_norm{S}{u}=spatial_tuning_vector_magnitude{S}{u}./max(abs(spatial_tuning_vector_magnitude{S}{u}(:)));
spatial_tuning_vector_norm{S}{u} = spatial_tuning_vector_magnitude_norm{S}{u}.*exp(spatial_tuning_vector_angle{S}{u}*sqrt(-1));

if isempty(spatial_tuning_vector_norm{S}{u});
    spatial_tuning_vector_norm{S}{u}=NaN;
end;end;end

%% Compute Tuning specificity
for S=1:Nshuffle
for u=1:size(spatial_tuning_vector_norm{S},2)
% compute weighted sum of cos and sin of angles
sum_vector{S}(u) = sum(spatial_tuning_vector_magnitude_norm{S}{u}.*exp(1i*spatial_tuning_vector_angle{S}{u}));
% compute magnitude = tuning specificity 
tuning_specificity(S,u) = abs(sum_vector{S}(u))./sum(spatial_tuning_vector_magnitude_norm{S}{u});
%compute angle
sum_angle{S}(u)= angle(sum_vector{S}(u));
% get vector from magnitude and angle
vector_tuning_specificity{S}(u)=tuning_specificity(S,u).*exp(sum_angle{S}(u)*sqrt(-1));
end;end
 
%% Compute p-value 
tuning_spe=TuningSpecificity.tuning_specificity;
tuning_spe_shuffle=tuning_specificity;

%p-value was defined as the fraction of shuffle distribution that 
%exceeded the cell’s true tuning specificity

tuning_spe_pvalue=(sum(tuning_spe_shuffle>tuning_spe))/Nshuffle;
for uu=1:size(tuning_spe,2)
if isnan(tuning_spe(uu)),
tuning_spe_pvalue(uu)=NaN;
end;end
tuning_spe_sigROI=tuning_spe_pvalue<pvalue;

SpatiallyTunedCells.TuningSpecificity.tuning_specificity=tuning_spe;
%SpatiallyTunedCells.TuningSpecificity.tuning_specificity_shuffle=tuning_spe_shuffle;

SpatiallyTunedCells.TuningSpecificity.pvalue=tuning_spe_pvalue;
SpatiallyTunedCells.TuningSpecificity.significantROI=tuning_spe_sigROI;


end






