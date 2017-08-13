%%  MAD Ca2+ analysis
%
% !toolboxes needed!  
% FMAtoolbox
% Bioinformatics Toolbox 
% 
%INPUTS
%Calcium Data (expDffMedZeroed, C_df, ...)
%Column=ROI, Row=frames nb
%C_df=(full(C_df)');
%csv and xml 
directory_name = uigetdir; %CSV and XML files
%
%If only one C_df with multiple exp: read number of frames for each session
%filename='TGOT_RWR_6-003';
%xmlfile = strcat(directory_name,'\',filename);
%[ xml ] = xml2struct( xmlfile );
%nbofframes=length(xml.PVScan.Sequence.Frame);

%% Compute behavior output
Behavior.options.mintimebeforelap=10; %time to wait before counting a new lap (s)
Behavior.options.acquisitionfrequency=10000; %behavior acquisition frequency
[Behavior] = behavior_output_V2(directory_name, Behavior);

%Restrict calcium imaging data:
%First lap 
Behavior.restricted.startlap= 1;
%Last lap
Behavior.restricted.endlap=length(Behavior.lap);
Imaging.options.frimaging=30.2061910850898;%imaging framerate (Hz)
Imaging.options.celltoplot=25; %exemple cell to plot

[Imaging, Behavior]= lapselect_V2(C_df, Behavior, Imaging);

%% Correct drifting baseline
Imaging.options.baselinesub.windwith=100; %width for the shifting window (in frame).
[Imaging]=baselinesub(Imaging);


%% Events detection
Events.options.restricted=1; %Event detection on restricted trace
Events.options.baselinesub=1; %Event detection on baseline substracetd trace
Events.options.SDON=3; %Threshold above x SD for ONSET 
Events.options.SDOFF=0.5; %Threshold below x SD for OFFSET 
Events.options.mindur=0.5; %Minimum duration (s) an event must stay above the ON threshold 
Events.options.mergdur=0; %Merge consecutive epochs separated by (s)
Events.options.mindurevent=1; %Min duration of event to be considered  %Danielson et al. used > 1s
tic;
[Events]=eventdetect_V5(Events, Imaging);
toc;
%Figure:
figure;hold on;
plot(Imaging.trace_restricted_baselinesub(:,26));
plot(Events.onset_binary(:,26));

%% Run onset
Events.options.run.speed_thr=5; % minimum peak speed (cm/s)
Events.options.norun.speed_thr=1; %Threshold for non running speed (cm/s)
Events.options.run.min_dur=1;  %Minimum duration (s) for running epoch
Events.options.norun.min_dur=1;  %Minimum duration (s) for non running epoch
Events.options.run.merge_dur=0.5; %Merge consecutive running epochs separated by less than (s)
Events.options.norun.merge_dur=0.25; %Merge consecutive non running epochs separated by less than (s)
Events.options.span=3; %Speed smoothing kernel

[Events, Behavior]=runevents_V6(Events,Behavior,Imaging);

%% Events analysis
Events.options.analysis.maxrisetime=50; %Maximun rise time (frame) 
Events.options.analysis.maxduration=250; %Maximun duration (frame) 
Events.options.analysis.maxpeak=1000; %in % change of dF/F

tic;
%Analysis on all events / Runnning events / Non-runnning events
onset_offset_all=[Events.onset_offset; Events.RunningEpochs.run_onset_offset; Events.NonRunningEpochs.norun_onset_offset];
if Events.options.restricted==1,Cdf_time=Imaging.time_restricted;elseif Events.options.restricted==0,Cdf_time=Imaging.time;
end
time_all=length(Cdf_time);time_run=sum(Behavior.runbinary);time_norun=sum(Behavior.norunbinary);time=[time_all time_run time_norun];
for c=1:3 %all, run, norun
[EventsProperties{c}, NetworkProperties{c}, PCAProperties{c}]=transients_analysis_V8(onset_offset_all(c,:), time(c), Events, Imaging);
end
Events.Properties=EventsProperties{1};Events.RunningEpochs.Properties=EventsProperties{2};Events.NonRunningEpochs.Properties=EventsProperties{3};
Network=NetworkProperties{1};Network.RunningEpochs=NetworkProperties{2};Network.NonRunningEpochs=NetworkProperties{3};
PCA=PCAProperties{1};PCA.RunningEpochs=PCAProperties{2};PCA.NonRunningEpochs=PCAProperties{3};
clear EventsProperties NetworkProperties PCAProperties c onset_offset_all Cdf_time time time_all time_run time_norun;
toc;

%% Identification of spatially-tuned cells
Behavior.placeoptions.smooth=3; %sigma value of filter
Behavior.placeoptions.minevents=3; %Minimun nb of events - excluded =NaN

% Spatial information
%Iteration for nb of bins = 2,4,5,8,10,20,25,100
Nbin=[2;4;5;8;10;20;25;100];
%Iteration for nb shuffle = 2,4,5,8,10,20,25,100
Nshuffle=1000; 
%P value min
pvalue=0.05;
tic;
for it=1:length(Nbin) %nb of iteration
[Spatialinfo{it}]=spatial_info_V2(Behavior, Events, Imaging,Nbin(it));
[Spatialinfo{it}]=shuffle_spatial_info(Spatialinfo{it},Nbin(it),Nshuffle);
end
[SpatiallyTunedCells]=spatial_tuning_pvalue_V2(Spatialinfo,Nbin,Nshuffle,pvalue);
toc;
%Tuning Specificity 
clear i;
tic;
[TuningSpecificity]=tuning_specificity_V4(Behavior,Events);
[TuningSpecificity,SpatiallyTunedCells]=shuffle_tuning_specificity_V3(TuningSpecificity,Behavior,Events,SpatiallyTunedCells, Nshuffle,pvalue);
toc;

%Tuned cells - and plot
%cells that present significant spatial information and tuning specificity
[SpatiallyTunedCells]=tuned_cells(SpatiallyTunedCells,Spatialinfo,TuningSpecificity);



%% Save
%Where to save:
dirnam='D:\MAD_ANALYSIS\CA3_ThyGC6f\05_04_2017\M3\'; 
tic;
save(strcat(dirnam,'1_STD1_Cdf'), '-v7.3');
toc;

clear
