
function [Spatial_correlation] = spatial_corr_V3(Spatial_correlation,Spatialinfo,TunedCells);
sessions=Spatial_correlation.options.sessions;
tunedsess=Spatial_correlation.options.onlytuned;
criteria1=Spatial_correlation.options.tuned_criteria1;
criteria2=Spatial_correlation.options.tuned_criteria2;


%% Population Vector correlation: 

%Import data / fix bug
for i = sessions
rate_map{i}=Spatialinfo{i}.rate_map.normalized_rate_map_smoothed;
rate_map_nonnorm{i}=Spatialinfo{i}.rate_map.rate_map_smoothed ; 
nbROI(i)=size(rate_map{i},2);
maxROI=max(nbROI);
tunedROI{i}=TunedCells{i}.TunedROI;
infoROI{i}=TunedCells{i}.Spatial_info.significantROI;  
tuningROI{i}=TunedCells{i}.TuningSpecificity.significantROI;  
end
for i = sessions
while size(rate_map{i},2)<maxROI
rate_map{i}=[rate_map{i} nan(size(rate_map{i},1),1)]; %replace missing last ROI with NaN
rate_map_nonnorm{i}=[rate_map_nonnorm{i} nan(size(rate_map_nonnorm{i},1),1)];
end
end
for i = sessions
while size(tunedROI{i},2)<maxROI
tunedROI{i}=[tunedROI{i} 0]; %replace missing last ROI with 0 - non tuned
end
while size(infoROI{i},2)<maxROI
infoROI{i}=[infoROI{i} 0];
end
while size(tuningROI{i},2)<maxROI
tuningROI{i}=[tuningROI{i} 0];
end
end

%PV matrix - all neurons
for i = sessions
for ii= sessions
[PVmatrix{i}{ii}] = corr(rate_map{i}',rate_map{ii}', 'rows', 'complete');
matrixnb{i}{ii}=[i ii];
end
end 


%PV correlation
%All neurons
for i = sessions
for ii= sessions
 for iii=1:size(rate_map{i},1)  
[PVcorr{i}{ii}{iii}] = corr(rate_map{i}(iii,:)',rate_map{ii}(iii,:)', 'rows', 'complete');
PVcorrind{i}{ii}=[i ii];
PVcorr_mean{i}{ii}= nanmean(cell2mat(PVcorr{i}{ii}));
end
end 
end


%Only tuned index
switch criteria1
    case 'info'
for i=tunedsess
onlytuned(i,:)=(infoROI{i});
end
    case 'tuning'
for i=tunedsess
onlytuned(i,:)=(tuningROI{i});
end
    case 'both'
for i=tunedsess
onlytuned(i,:)=(tunedROI{i});
end
  otherwise
warning('Unexpected criteria. Use: info / tuning / both');
end

switch criteria2
    case 'or'
tunedindx=find(sum(onlytuned,1)>=1);
    case 'and'
tunedindx=find(sum(onlytuned,1)>=length(tunedsess));
  otherwise
warning('Unexpected criteria. Use: and/or');
end

for i=sessions
if isempty(tunedindx)
rate_map_tuned{i}=nan;
warning('No tuned cells');
end
end
for i=sessions
for ii=tunedindx
if isempty(tunedindx)==0
rate_map_tuned{i}(:,ii)=rate_map{i}(:,ii);
rate_map_tuned{i}=rate_map_tuned{i}(:,any(rate_map_tuned{i})) ;
end
end
end


for i=sessions
for ii=1:size(rate_map_tuned{i},2)
siztunedROI(i)=size(rate_map_tuned{i},2);
maxtunedROI=max(siztunedROI);
end
end
for i=sessions
for ii=1:size(rate_map_tuned{i},2)
while size(rate_map_tuned{i},2)<maxtunedROI
rate_map_tuned{i}=[rate_map_tuned{i} nan(size(rate_map_tuned{i},1),1)];  
end
end
end

for i = sessions
for ii= sessions
 for iii=1:size(rate_map_tuned{i},1)  
[PVcorr_tuned{i}{ii}{iii}] = corr(rate_map_tuned{i}(iii,:)',rate_map_tuned{ii}(iii,:)', 'rows', 'complete');
PVcorr_tuned_mean{i}{ii}= nanmean(cell2mat(PVcorr_tuned{i}{ii}));
end
end 
end



Spatial_correlation.PVcorrelation.matrix=PVmatrix;
Spatial_correlation.PVcorrelation.index=matrixnb;
Spatial_correlation.PVcorrelation.AllROI.PVcorr=PVcorr;
Spatial_correlation.PVcorrelation.AllROI.PVcorr_index=PVcorrind;
Spatial_correlation.PVcorrelation.AllROI.PVcorr_mean=PVcorr_mean;
Spatial_correlation.PVcorrelation.TunedROI.PVcorr_tuned=PVcorr_tuned;
Spatial_correlation.PVcorrelation.TunedROI.PVcorr_tuned_index=PVcorrind;
Spatial_correlation.PVcorrelation.TunedROI.PVcorr_mean=PVcorr_tuned_mean;

%% Tuning Curve correlation


%All neurons
for i = sessions
for ii= sessions
 for iii=1:size(rate_map{i},2)  
[TCcorr{i}{ii}{iii}] = corr(rate_map{i}(:,iii),rate_map{ii}(:,iii), 'rows', 'complete');
TCcorrind{i}{ii}=[i ii];
TCcorr_mean{i}{ii}= nanmean(cell2mat(TCcorr{i}{ii}));
end
end 
end

%Only tuned index
for i = sessions
for ii= sessions
 for iii=tunedindx 
[TCcorr_tuned{i}{ii}{iii}] = corr(rate_map{i}(:,iii),rate_map{ii}(:,iii), 'rows', 'complete');
TCcorrind_tuned{i}{ii}=[i ii];
TCcorr_tuned_mean{i}{ii}= nanmean(cell2mat(TCcorr_tuned{i}{ii}));
end
end 
end

Spatial_correlation.TuningCurvecorr.AllROI.TCcorr=TCcorr;
Spatial_correlation.TuningCurvecorr.AllROI.TCcorr_index=TCcorrind;
Spatial_correlation.TuningCurvecorr.AllROI.TCcorr_mean=TCcorr_mean;

if isempty(tunedindx)
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned=nan;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned_index=nan;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_mean=nan;
warning('No tuned cells');
end

if isempty(tunedindx)==0
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned=TCcorr_tuned;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned_index=TCcorrind;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_mean=TCcorr_tuned_mean;
end


%% Rate Overlap

%rate overall and peak rate
for i = sessions
overall_rate_array{i}=mean(rate_map_nonnorm{i});
peak_rate_array{i}=max(rate_map_nonnorm{i});
end 

%All ROI
for i = sessions
for ii = sessions 
comprate{i}{ii}=[overall_rate_array{i}; overall_rate_array{ii}] ;
comppeakrate{i}{ii}=[peak_rate_array{i}; peak_rate_array{ii}];
maxrate{i}{ii}=max(comprate{i}{ii});
maxpeak{i}{ii}=max(comppeakrate{i}{ii});
minrate{i}{ii}=min(comprate{i}{ii});
minpeak{i}{ii}=min(comppeakrate{i}{ii});
rate_ratio{i}{ii}=minrate{i}{ii}./maxrate{i}{ii};
peak_ratio{i}{ii}=minpeak{i}{ii}./maxpeak{i}{ii};
peak_overlap{i}(ii)=nanmean(peak_ratio{i}{ii});
rate_overlap{i}(ii)=nanmean(rate_ratio{i}{ii});
rate_overlap_norm{i}=(rate_overlap{i} - min(rate_overlap{i}))/(max(rate_overlap{i}) - min(rate_overlap{i}));
end
end
%Only tuned
for i = sessions
for ii = sessions 
for iii= tunedindx
rate_ratio_tuned{i}{ii}(iii)=minrate{i}{ii}(iii)./maxrate{i}{ii}(iii);
peak_ratio_tuned{i}{ii}(iii)=minpeak{i}{ii}(iii)./maxpeak{i}{ii}(iii);
rate_overlap_tuned{i}(ii)=sum(rate_ratio_tuned{i}{ii},2)./sum(rate_ratio_tuned{i}{ii}~=0,2);
peak_overlap_tuned{i}(ii)=sum(peak_ratio_tuned{i}{ii},2)./sum(peak_ratio_tuned{i}{ii}~=0,2);
end
end
end   


if isempty(tunedindx)
Spatial_correlation.overlap.Rate_ratio_tuned=nan;
Spatial_correlation.overlap.Peak_ratio_tuned=nan;
Spatial_correlation.overlap.Rate_overlap_tuned=nan;
Spatial_correlation.overlap.Peak_overlap_tuned=nan;
Spatial_correlation.overlap.Rate_ratio_tuned=nan;
Spatial_correlation.overlap.Peak_ratio_tuned=nan;
warning('No tuned cells');
end


Spatial_correlation.overlap.Rate_ratio=rate_ratio;
Spatial_correlation.overlap.Peak_ratio=peak_ratio;
Spatial_correlation.overlap.Rate_overlap=rate_overlap;
Spatial_correlation.overlap.Peak_overlap=peak_overlap;

if isempty(tunedindx)==0
Spatial_correlation.overlap.Rate_ratio_tuned=rate_ratio_tuned;
Spatial_correlation.overlap.Peak_ratio_tuned=peak_ratio_tuned;
Spatial_correlation.overlap.Rate_ratio_tuned=rate_ratio_tuned;
Spatial_correlation.overlap.Peak_ratio_tuned=peak_ratio_tuned;
Spatial_correlation.overlap.Rate_overlap_tuned=rate_overlap_tuned;
Spatial_correlation.overlap.Peak_overlap_tuned=peak_overlap_tuned;
end


end


