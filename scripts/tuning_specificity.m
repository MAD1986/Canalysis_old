
function [Place_cell]=tuning_specificity(Place_cell,Behavior,Events,options);
%% Extract the data
run_position=Behavior.resampled.run_position;
run_time=Behavior.resampled.run_time;
if Events.options.restricted==1,
lap_start_stop= Behavior.restricted.lap;
elseif Events.options.restricted==0,
lap_start_stop= Behavior.lap;
end
run_onset_binary=Events.Run.run_onset_binary;
run_ones=Behavior.run_ones; 
%% Assemble DATA
%Time and position for running epochs ONLY:
run_position_norm=(((run_position - min(run_position)) / (( max(run_position) - min(run_position)))));
run_onset_bin=run_onset_binary(run_ones==1,:);
%Remove from analysis cell if not enough events (set in options.minevents)
for n=1:size(run_onset_bin,2)
if sum(run_onset_bin(:,n))<options.minevents
run_onset_bin(:,n)==0;
end
end  
%position of the mouse at onset time:
for u=1:size(run_onset_bin,2)
onset_position{u}=run_position_norm(run_onset_bin(:,u)==1);
end
%Fraction of frames acquired at position of onset
for u=1:size(onset_position,2)
    for uu=1:size(onset_position{u},1)
         if isempty(onset_position{u})==0,     
frame_onset{u}{uu}=find(run_position_norm==onset_position{u}(uu));
nbframe_onset{u}(:,uu)= size(frame_onset{u}{uu},1);
fraction_frames{u}(:,uu)= nbframe_onset{u}(:,uu)/ (size(run_onset_bin,1));
end  
end
end
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
end
end
%Normalize magnitude
for u=1:size(spatial_tuning_vector,2)   
spatial_tuning_vector_magnitude_norm{u}=spatial_tuning_vector_magnitude{u}./max(abs(spatial_tuning_vector_magnitude{u}(:)));
spatial_tuning_vector_norm{u} = spatial_tuning_vector_magnitude_norm{u}.*exp(spatial_tuning_vector_angle{u}*sqrt(-1));

if isempty(spatial_tuning_vector_norm{u});
    spatial_tuning_vector_norm{u}=NaN;
end
end
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

%% write results in structure
Place_cell.Tuning_Specificity.tuning_vector=spatial_tuning_vector_norm;
Place_cell.Tuning_Specificity.tuning_specificity=tuning_specificity;
Place_cell.Tuning_Specificity.tuning_vector_specificity=vector_tuning_specificity;
end