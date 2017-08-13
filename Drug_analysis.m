

%% Load files
%named CTRL1, CTRL2 ... PSAM1, PSAM2 ...
%CTRL1 = one mouse/day experiment CTRL2= another day/mouse 
folder_name = uigetdir %folder with the mat files
cd(folder_name);
expnb=3; %nb of experiments


for i=1:expnb %nb of experiments
name_CTRL = sprintf('%s%d.mat', 'CTRL',i);
name_DRUG = sprintf('%s%d.mat', 'DRUG',i);   
Data_CTRL = load (name_CTRL);  
Data_DRUG = load (name_DRUG);  
event_properties.CTRL{i} = Data_CTRL.eventsproperties; 
event_properties.DRUG{i} = Data_DRUG.eventsproperties; 
networkproperties.CTRL{i} = Data_CTRL.networkproperties; 
networkproperties.DRUG{i} = Data_DRUG.networkproperties; 
end


%% Compare mean all events (run and no run)
for i=1:expnb %nb of experiments
event_properties.CTRLvsDRUG{i}.mean.duration=[event_properties.CTRL{i}{1}.mean.duration'...
    event_properties.DRUG{i}{1}.mean.duration'];
event_properties.CTRLvsDRUG{i}.mean.peak=[event_properties.CTRL{i}{1}.mean.peak'...
    event_properties.DRUG{i}{1}.mean.peak'];
event_properties.CTRLvsDRUG{i}.mean.amplitude=[event_properties.CTRL{i}{1}.mean.amplitude'...
    event_properties.DRUG{i}{1}.mean.amplitude'];
event_properties.CTRLvsDRUG{i}.mean.meandf=[event_properties.CTRL{i}{1}.mean.meandf'...
    event_properties.DRUG{i}{1}.mean.meandf'];
event_properties.CTRLvsDRUG{i}.mean.area=[event_properties.CTRL{i}{1}.mean.area'...
    event_properties.DRUG{i}{1}.mean.area'];
event_properties.CTRLvsDRUG{i}.mean.decaytime=[event_properties.CTRL{i}{1}.mean.decaytime'...
    event_properties.DRUG{i}{1}.mean.decaytime'];
event_properties.CTRLvsDRUG{i}.mean.halfwidth=[event_properties.CTRL{i}{1}.mean.halfwidth'...
    event_properties.DRUG{i}{1}.mean.halfwidth'];

networkproperties.CTRLvsDRUG{i}.ActivityRate=[networkproperties.CTRL{i}{1}.ActivityRate'...
    networkproperties.DRUG{i}{1}.ActivityRate'];
networkproperties.CTRLvsDRUG{i}.Frequency=[networkproperties.CTRL{i}{1}.Frequency'...
    networkproperties.DRUG{i}{1}.Frequency'];
end

for i=1:expnb %nb of experiments
event_properties.CTRLvsDRUG_ALL.mean.duration{i}=event_properties.CTRLvsDRUG{i}.mean.duration;
event_properties.CTRLvsDRUG_ALL.mean.peak{i}=event_properties.CTRLvsDRUG{i}.mean.peak;
event_properties.CTRLvsDRUG_ALL.mean.amplitude{i}=event_properties.CTRLvsDRUG{i}.mean.amplitude;
event_properties.CTRLvsDRUG_ALL.mean.meandf{i}=event_properties.CTRLvsDRUG{i}.mean.meandf;
event_properties.CTRLvsDRUG_ALL.mean.area{i}=event_properties.CTRLvsDRUG{i}.mean.area;
event_properties.CTRLvsDRUG_ALL.mean.decaytime{i}=event_properties.CTRLvsDRUG{i}.mean.decaytime;
event_properties.CTRLvsDRUG_ALL.mean.halfwidth{i}=event_properties.CTRLvsDRUG{i}.mean.halfwidth;

networkproperties.CTRLvsDRUG_ALL.ActivityRate{i}=networkproperties.CTRLvsDRUG{i}.ActivityRate;
networkproperties.CTRLvsDRUG_ALL.Frequency{i}=networkproperties.CTRLvsDRUG{i}.Frequency;
end

event_properties.CTRLvsDRUG_ALL.mean.duration=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.duration');
event_properties.CTRLvsDRUG_ALL.mean.peak=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.peak');
event_properties.CTRLvsDRUG_ALL.mean.amplitude=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.amplitude');
event_properties.CTRLvsDRUG_ALL.mean.meandf=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.meandf');
event_properties.CTRLvsDRUG_ALL.mean.area=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.area');
event_properties.CTRLvsDRUG_ALL.mean.decaytime=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.decaytime');
event_properties.CTRLvsDRUG_ALL.mean.halfwidth=cell2mat(event_properties.CTRLvsDRUG_ALL.mean.halfwidth');
networkproperties.CTRLvsDRUG_ALL.ActivityRate=cell2mat(networkproperties.CTRLvsDRUG_ALL.ActivityRate');
networkproperties.CTRLvsDRUG_ALL.Frequency=cell2mat(networkproperties.CTRLvsDRUG_ALL.Frequency');



