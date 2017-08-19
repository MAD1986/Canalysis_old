
function [Events, Network]=event_analysis(Behavior, Events, Imaging, options);

%Events=Events{1};
%Imaging=Imaging{1};
%Behavior=Behavior{1};
%Options


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
time{2}=Behavior.resampled.run_time;
time{3}=Behavior.resampled.no_run_time;
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
for i=1:3
[Network_Properties{i}] = network_properties(Event_Properties{i},C_df,Cdf_time, options);
end


%% PCA 
%for i=1:3
%[PCA_Properties{i}]= PCA_analysis(Event_Properties{i} ,Network_Properties{i});
%end

%% Make structure
Events.properties=Event_Properties{1};
Events.Run.properties=Event_Properties{2};
Events.NoRun.properties=Event_Properties{3};
Network=Network_Properties{1};
Network.Run=Network_Properties{2};
Network.NoRun=Network_Properties{3};
%PCA=PCA_Properties{1};
%PCA.Run=PCA_Properties{2};
%PCA.NoRun=PCA_Properties{3};
end

