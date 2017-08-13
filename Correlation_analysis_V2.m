%% Correlation analysis
%To test : bin where texture was change -> shift these bin to map the STD
%sessions : look a PV correlation rate map correlation


%% Import data
%IF multiple mat
path = uigetdir %folder with the all the mat and CSV 
cd(path);
list = dir('*.mat'); %will look for all the mat files
%!!!! files are listed by alphabetical order !!!!
%!!!! Keep same name for mat and csv !!!!
for i = 1:length(list)
file=list(i).name;
session{i} = load((list(i).name));
Spatialinfo{i}=session{i}.Spatialinfo{8};
TunedCells{i}=session{i}.SpatiallyTunedCells;
TuningSpecificity{i}=session{i}.TuningSpecificity;
end
%IF one matlab file
TunedCells=SpatiallyTunedCells;
for i=1:5
Spatialinfo{i}=Spatialinfo{i}{8};
end
%% Plot rate map ordered
%correlation.options.sessions=1:length(list); %nb of sessions - compare all sessions
Spatial_correlation.options.sessions=1:5; % X:Y compare session X to Y
Spatial_correlation.options.order=1; %order rate map for session # (1 to 3)
[ordered_map] = plotratemap_V3(Spatial_correlation,Spatialinfo,TunedCells);

%% Spatial correlation, rate overlap
Spatial_correlation.options.sessions=1:5; %sessions to compare -all

%[X] = session X /// [X Y] sessions X and/or Y /// [A:Z] sessions A to X
Spatial_correlation.options.onlytuned=[1 3 5];
Spatial_correlation.options.tuned_criteria1='info'; 
Spatial_correlation.options.tuned_criteria2='or'; 
        %criteria1:
%'both' = tuning specificity + spatial information
%'tuning' = tuning specificity
%'info'= spatial information
        %criteria2:
%or = needs to be tuned in one of the selected
%and= needs to be tuned in all of the selected sessions

%Danielson et al. = spatial info in both session for tuning vector correlation 
[Spatial_correlation] = spatial_corr_V3(Spatial_correlation,Spatialinfo,TunedCells);


%% Rotation analysis
Spatial_correlation.options.sessions=1:5; %sessions to compare -all

%[X] = session X /// [X Y] sessions X and/or Y /// [A:Z] sessions A to X
Spatial_correlation.options.onlytuned=[1 3 5];
Spatial_correlation.options.tuned_criteria1='info'; 
Spatial_correlation.options.tuned_criteria2='or'; 


%% Plot 
%[x;y]: plot each correlation between session X vs session Y
plot_session{1}=[1;3]; %STD
plot_session{2}=[1;2]; %MIS1
plot_session{3}=[2;3]; %MIS1
%plot_session{4}=[3;5]; %STD
%plot_session{5}=[3;4]; %MIS2
%plot_session{6}=[4;5]; %MIS2
color='Jet'; % colormap

%Plot PV Correlation
[Spatial_correlation]=plotPVcorr_V2(Spatial_correlation,color, plot_session);
%Plot TC Correlation
[Spatial_correlation]=plotTCcorr(Spatial_correlation,color, plot_session);
%Plot Overlap
[Spatial_correlation] = plotrateoverlap(Spatial_correlation,color, plot_session);
%Plot standard sessions (STD) VS mismatch sessions(MIS)
STD=[1;4];
MIS{1}=[2;3];
MIS{2}=[5;6];
[Spatial_correlation] = plot_STDMIS(Spatial_correlation,STD,MIS);




