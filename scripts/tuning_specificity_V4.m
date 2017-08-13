
function [TuningSpecificity]=tuning_specificity_V3(Behavior,Events);
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
       % run_binary(:,u)=0;
    end
end
 
%% Assemble DATA
%Time and position for running epochs ONLY:
running_idx=find(run_binary==1);
run_time_position=time_position(running_idx,:);
run_position_norm=(run_time_position(:,2) - min(run_time_position(:,2))) / ( max(run_time_position(:,2)) - min(run_time_position(:,2)) );
run_onset_binary2=run_onset_binary(running_idx,:);
%position of the mouse at onset time:
for u=1:size(run_onset_binary2,2)
onset_idx{u}=find(run_onset_binary2(:,u)==1);
onset_position{u}=run_position_norm(onset_idx{u});
end
%Fraction of frames acquired at position of onset
for u=1:size(onset_position,2)
    for uu=1:size(onset_position{u},1)
         if isempty(onset_position{u})==0,     
frame_onset{u}{uu}=find(run_position_norm==onset_position{u}(uu));
nbframe_onset{u}(:,uu)= size(frame_onset{u}{uu},1);
fraction_frames{u}(:,uu)= nbframe_onset{u}(:,uu)/ (size(run_onset_binary2,1));
end  
end;end
%% Compute spatial tuning vector
for u=1:size(onset_position,2) 
    for uu=1:size(onset_position{u},1)
 if isempty(frame_onset{u})==0,        
spatial_tuning_vector{u}(uu)=(exp(i*onset_position{u}(uu)*2*pi))./(fraction_frames{u}(:,uu));       
spatial_tuning_vector_angle{u}(uu)=angle(spatial_tuning_vector{u}(uu));
spatial_tuning_vector_magnitude{u}(uu)=abs(spatial_tuning_vector{u}(uu));
 elseif isempty(frame_onset{u}),
spatial_tuning_vector{u}(uu)=NaN; 
spatial_tuning_vector_angle{u}(uu)=NaN;
spatial_tuning_vector_magnitude{u}(uu)=NaN;
 end
 end;end
%Normalize magnitude
for u=1:size(spatial_tuning_vector,2)   
spatial_tuning_vector_magnitude_norm{u}=spatial_tuning_vector_magnitude{u}./max(abs(spatial_tuning_vector_magnitude{u}(:)));
spatial_tuning_vector_norm{u} = spatial_tuning_vector_magnitude_norm{u}.*exp(spatial_tuning_vector_angle{u}*sqrt(-1));

if isempty(spatial_tuning_vector_norm{u});
    spatial_tuning_vector_norm{u}=NaN;
end;end

%% Compute Tuning specificity
for u=1:size(spatial_tuning_vector_norm,2)
% compute weighted sum of cos and sin of angles
sum_vector(u) = sum(spatial_tuning_vector_magnitude_norm{u}.*exp(1i*spatial_tuning_vector_angle{u}));
% compute magnitude = tuning specificity 
tuning_specificity(u) = abs(sum_vector(u))./sum(spatial_tuning_vector_magnitude_norm{u});
%compute angle
sum_angle(u)= angle(sum_vector(u));
% get vector from magnitude and angle
vector_tuning_specificity(u)=tuning_specificity(u).*exp(sum_angle(u)*sqrt(-1));
end

TuningSpecificity.tuning_vector=spatial_tuning_vector_norm;
TuningSpecificity.tuning_specificity=tuning_specificity;
TuningSpecificity.tuning_vector_specificity=vector_tuning_specificity;
    
end