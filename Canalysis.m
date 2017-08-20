%%  MAD Ca2+ analysis
%
% !toolboxes needed!  
% Bioinformatics Toolbox 

%% Import data (XML, MAT, CSV)
directory_name = '/Users/martial/Documents/Results/CA3_ThyGC6f/M1/FOV4' %folder with .mat .csv and .xml
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
for i=1:sessions
Behavior{i} = behavior_lap(CSV{i},options);
end
% Restrict calcium data to selected lap
% Using frame timestamps from associated xml file
options.startlap= 'first'; %First = 1st complete lap or numbers
options.endlap='last'; %Last= last complete lap or numbers
options.dispfig=1; % Display figure 
if options.dispfig==true,
options.c2plot=2; % neuron to display
end
% Function
for i=1:sessions
[Imaging{i}, Behavior{i}]= lapselect(C_df{i}, Behavior{i}, XML{i}, options);
end
clear CSV XML options;

%% Correct drifting baseline
% Info 
% https://www.mathworks.com/help/bioinfo/ref/msbackadj.html
options.msbackadj=0; % 1 or 0 
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
options.iterations=3; %Nb of iterations %Danielson et al. used 3 iterations
options.SDOFF=0.5; %Threshold below x SD for OFFSET 
options.msbackadj=0; % Corrected baseline
% Test different threshold with histogram of error rate 
% Lovett-Barron et al. 2014
% Rajasethupathy et al. 2015
testSD=1;
if testSD==true
options.SD=[2, 2.5, 3, 3.5]; %Threshold above x SD for ONSET to test
% Show histogram: optional
for i=1:sessions
[Events{i}]=sigma_events(Imaging{i},options);
end
end

%Set parameters for event detection 
options.SDON=2.5; %Threshold above x SD for ONSET 
options.mindurevent=1; %Min duration of event to be considered Danielson et al. used > 1s
options.dispfig=1; % Display figure 
if options.dispfig==true,
options.c2plot=2; % neuron to display
end
% Detect events
for i=1:sessions
[Events{i}]=detect_events(Imaging{i},options);
end
clear options;

%% Running epochs
options.method='speed'; %'peak' or 'speed'
% threshold running epochs based on :
% min peak speed(Danielson et al. 2016a, b, 2016) 
% OR average speed (Cossart) 
options.moving_window=10; % window width for moving mean filter
options.minspeed=2; %minimum speed (cm/s)  -- only if 'speed'
options.minpeak=5;  % minimum peak speed (cm/s) -- only if 'peak'
options.mindur=1; %Minimum duration (s) for running epoch
options.merge=0.5; %Merge consecutive running epochs separated by less than (s)
options.dispfig=1; % Display figure 
if options.dispfig==true
options.c2plot=4; % neuron to display
end
%Function
for i=1:sessions
[Events{i}, Behavior{i}]=run_epoch(Events{i},Behavior{i},Imaging{i}, options);
end
clear options
%% Events analysis
% Events parameters
options.exclude=1; % exclude events when no peaks found 
options.mindist= 20; % Set minimun distance between peaks (frame)
options.STD_pro=2.5; % Set minimun prominence of peak ( X * STD noise)
% Network parameters
% synchronous epochs based on Rajasethupathy et al. 2015
options.Nshuffle=100; % Nb of shuffle for synchronous activity 
options.pmin=0.05; % min p value to be considered as significant
options.minframes=3; %consecutive frames with activity above the significance threshold
%Function
for i=1:sessions
[Events{i}, Network{i}]=event_analysis(Behavior{i}, Events{i}, Imaging{i}, options);
end 
clear options
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
save('test_ca_3', '-v7.3');
toc;



