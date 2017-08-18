
function [EventsProperties, NetworkProperties, PCAProperties]=transients_analysis_V10(onset_offset,time,binary, Events, Imaging);

%% Import data

% standard deviation of noise
options.STD_off=Events.STD_noise;

%Calcium trace
if Imaging.options.msbackadj== true && options.restricted==true
C_df=Imaging.trace_restricted_baselinesub;  
elseif Imaging.options.msbackadj== true && options.restricted==false
C_df=Imaging.trace_baselinesub;
elseif Imaging.options.msbackadj== false && options.restricted==true
C_df=Imaging.trace_restricted;
elseif Imaging.options.msbackadj== false && options.restricted==false
C_df=Imaging.trace; 
end

%Time 
if options.restricted==true
Cdf_time=Imaging.time_restricted;
elseif options.restricted==false
Cdf_time=Imaging.time;
end

time{1}=Behavior.resampled.time;
time{2}=Behavior.run_time;
time{3}=Behavior.no_run_time;

time_all=Behavior.resampled.time;

%Events
onset_offset{1}=Events.onset_offset;
onset_offset{2}=Events.Run.run_onset_offset;
onset_offset{3}=Events.No_Run.norun_onset_offset;

onset_binary{1}=Events.onset_binary; 
onset_binary{2}=Events.Run.run_onset_binary;
onset_binary{3}=Events.No_Run.norun_onset_binary;

onset_ones{1}=Events.onset_ones; 
onset_ones{2}=Events.Run.run_onset_ones;
onset_ones{3}=Events.No_Run.norun_onset_ones;
%% Measure Ca2+ events properties 
for i=1:3
[Event_Properties{i}]=event_properties(C_df, onset_offset{i}, time{i}, options);
end


%% Network properties

%network_properties script


%% PCA 

% event_PCA script

end
