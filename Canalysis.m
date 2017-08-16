%%  MAD Ca2+ analysis
%
% !toolboxes needed!  
% Bioinformatics Toolbox 

%% Import data (XML, MAT, CSV)
directory_name = '/Users/martial/Documents/test/M1' %folder with .mat .csv and .xml
%!!!! files are listed by alphabetical order !!!!
%!!!! Keep same name for mat and csv !!!!
type='expdff';    %tye of calcium trace 'Cdf' ; 'expdff' or 'spikes' 
tic;
[C_df,CSV,XML]=combine_sessions(directory_name,type);
toc;
%Sessions to analyse
sessions=length(C_df); %ALL

%% Compute behavior output

% Set parameters
options.mindist=10; %min distance before counting a new lap (cm)
options.acqHz=10000; % behavior acquisition frequency (Hz)
options.textures=1; % 1= find RFID for each texture; 0 = no RFID for texture
options.dispfig=1; % Display figure (lap, RFID)

%IF find RFID for each texture
%voltage [min max] signal for each RFID
if options.textures==true
RFID{1}=[0.5 1]; 
RFID{2}=[2 2.5];
RFID{3}=[3 3.5];
RFID{4}=[1.5 2];
RFID{5}=[1 1.5];
options.RFID=RFID;
end

% Extract position / lap / textures
%Behavior = behavior_lap(CSV,options);
for i=1:sessions
Behavior{i} = behavior_lap(CSV{i},options);
end

% Restrict calcium data to selected lap
% Using frame timestamps from associated xml file
options.startlap= 'first'; %First = 1st complete lap or numbers
options.endlap='last'; %Last= last complete lap or numbers
options.dispfig=1; % Display figure 
if options.dispfig==true,
options.c2plot=12; % neuron to display
end
%[Imaging, Behavior]= lapselect(C_df, Behavior, XML, options);

for i=1:sessions
[Imaging{i}, Behavior{i}]= lapselect(C_df{i}, Behavior{i}, XML{i}, options);
end
clear CSV XML options;

%% Correct drifting baseline
%Info 
%https://www.mathworks.com/help/bioinfo/ref/msbackadj.html
options.msbackadj=1; % 1 or 0 

% Set parameters
options.windwith=250; %width for the shifting window (in frame)
options.dispfig=1; % Display figure 
if options.dispfig==true,
options.c2plot=12; % neuron to display
end

% Correct baseline of signal with peaks
for i=1:sessions
if options.msbackadj==true
[Imaging{i}]=baselinesub(Imaging{i}, options);
elseif options.msbackadj==false
Imaging{i}.options.msbackadj=0;
end
end

clear options;


%% Events detection

% Set parameters
options.restricted=1; %Event detection on restricted trace
options.SDON=4; %Threshold above x SD for ONSET 
options.SDOFF=0.5; %Threshold below x SD for OFFSET 
options.mindurevent=1; %Min duration of event to be considered  %Danielson et al. used > 1s
options.iterations=3 %Nb of iterations %Danielson et al. used 3 iterations
options.dispfig=1; % Display figure 
if options.dispfig==true,
options.c2plot=4; % neuron to display
end

% Detect events
%[Events]=detect_events(Imaging,options);
for i=1:sessions
[Events{i}]=detect_events(Imaging{i},options);
end

%% Run onset
for i=1:sessions
Events{i}.options.run.speed_thr=5; % minimum peak speed (cm/s)
Events{i}.options.norun.speed_thr=1; %Threshold for non running speed (cm/s)
Events{i}.options.run.min_dur=1;  %Minimum duration (s) for running epoch
Events{i}.options.norun.min_dur=1;  %Minimum duration (s) for non running epoch
Events{i}.options.run.merge_dur=0.5; %Merge consecutive running epochs separated by less than (s)
Events{i}.options.norun.merge_dur=0.25; %Merge consecutive non running epochs separated by less than (s)
Events{i}.options.span=3; %Speed smoothing kernel
[Events{i}, Behavior{i}]=runevents_V6(Events{i},Behavior{i},Imaging{i});
end
%% Events analysis
for i=1:sessions
Events{i}.options.analysis.maxrisetime=50; %Maximun rise time (frame) 
Events{i}.options.analysis.maxduration=250; %Maximun duration (frame) 
Events{i}.options.analysis.maxpeak=1000; %in % change of dF/F
%Analysis on all events / Runnning events / Non-runnning events
onset_offset_all{i}=[Events{i}.onset_offset; Events{i}.RunningEpochs.run_onset_offset; Events{i}.NonRunningEpochs.norun_onset_offset];
if Events{i}.options.restricted==1,Cdf_time{i}=Imaging{i}.time_restricted;elseif Events{i}.options.restricted==0,Cdf_time{i}=Imaging{i}.time;
end
end
for i=1:sessions
time_all{i}=length(Cdf_time{i});
time_run{i}=sum(Behavior{i}.runbinary);
run_binary{i}=Behavior{i}.runbinary;
time_norun{i}=sum(Behavior{i}.norunbinary);
norun_binary{i}=Behavior{i}.norunbinary;
time{i}=[time_all{i} time_run{i} time_norun{i}];
binary{i}=[ones((length(Cdf_time{i})),1) run_binary{i} norun_binary{i}];
end
for i=1:sessions
    tic;
for c=1:3 %all, run, norun
[EventsProperties{i}{c}, NetworkProperties{i}{c}, PCAProperties{i}{c}]=transients_analysis_V9(onset_offset_all{i}(c,:), time{i}(c), binary{i}(:,c), Events{i}, Imaging{i});
end
toc;
end

for i=1:sessions
for c=1:3 %all, run, norun
Events{i}.Properties=EventsProperties{i}{1};
Events{i}.RunningEpochs.Properties=EventsProperties{i}{2};
Events{i}.NonRunningEpochs.Properties=EventsProperties{i}{3};
Network{i}=NetworkProperties{i}{1};
Network{i}.RunningEpochs=NetworkProperties{i}{2};
Network{i}.NonRunningEpochs=NetworkProperties{i}{3};
PCA{i}=PCAProperties{i}{1};
PCA{i}.RunningEpochs=PCAProperties{i}{2};
PCA{i}.NonRunningEpochs=PCAProperties{i}{3};
end
end
clear EventsProperties NetworkProperties PCAProperties c onset_offset_all Cdf_time time time_all time_run time_norun run_binary norun_binary binary;


%% Identification of spatially-tuned cells
for j=1:sessions
tic;
Behavior{j}.placeoptions.smooth=3; %sigma value of filter
Behavior{j}.placeoptions.minevents=3; %Minimun nb of events - excluded =NaN
% Spatial information
%Iteration for nb of bins = 2,4,5,8,10,20,25,100
Nbin=[2;4;5;8;10;20;25;100];
%Iteration for nb shuffle 
Nshuffle=10000; 
%P value min
pvalue=0.05;
for it=1:length(Nbin) %nb of iteration

[Spatialinfo{j}{it}]=spatial_info_V2(Behavior{j}, Events{j}, Imaging{j},Nbin(it));
[Spatialinfo{j}{it}]=shuffle_spatial_info(Spatialinfo{j}{it},Nbin(it),Nshuffle);

end
toc;
end

for j=1:sessions
tic;
[SpatiallyTunedCells{j}]=spatial_tuning_pvalue_V2(Spatialinfo{j},Nbin,Nshuffle,pvalue);
%Tuning Specificity 
clear i;
[TuningSpecificity{j}]=tuning_specificity_V4(Behavior{j},Events{j});
[TuningSpecificity{j},SpatiallyTunedCells{j}]=shuffle_tuning_specificity_V3(TuningSpecificity{j},Behavior{j},Events{j},SpatiallyTunedCells{j}, Nshuffle,pvalue);
%Tuned cells - and plot
%cells with significant spatial information and tuning specificity
[SpatiallyTunedCells{j}]=tuned_cells(SpatiallyTunedCells{j},Spatialinfo{j},TuningSpecificity{j});
toc;
end


%% Textures
for i=1:sessions
tic;
tex{1}=[0.5 1]; %min to max threshold 
tex{2}=[2 2.5];
tex{3}=[3 3.5];
tex{4}=[1.5 2];
tex{5}=[1 1.5];
[Spatialinfo{i}]=textures_V2(Behavior{i}, Events{i}, Spatialinfo{i},tex);
toc;
end



%% Save
%Where to save:
tic;
save('test_ca', '-v7.3');
toc;



