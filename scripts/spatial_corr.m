
function [Spatial_correlation] = spatial_corr(Spatial_correlation,Spatialinfo,TunedCells);
sessions=Spatial_correlation.options.sessions;
tunedsess=Spatial_correlation.options.onlytuned;
criteria=Spatial_correlation.options.tuned_criteria;

%% Population Vector correlation: 

%Import data / fix bug
for i = sessions
rate_map{i}=Spatialinfo{i}.rate_map.normalized_rate_map_smoothed;
nbROI(i)=size(rate_map{i},2);
maxROI=max(nbROI);
tunedROI{i}=TunedCells{i}.TunedROI;
infoROI{i}=TunedCells{i}.Spatial_info.significantROI;  
tuningROI{i}=TunedCells{i}.TuningSpecificity.significantROI;  
end
for i = sessions
while size(rate_map{i},2)<maxROI
rate_map{i}=[rate_map{i} nan(size(rate_map{i},1),1)]; %replace missing last ROI with NaN
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

switch criteria
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


tunedindx=find(sum(onlytuned)>1);
for i = sessions
for ii= sessions
 for iii=tunedindx 
[PVcorr_tuned{i}{ii}{iii}] = corr(rate_map{i}(iii,:)',rate_map{ii}(iii,:)', 'rows', 'complete');
PVcorrind_tuned{i}{ii}=[i ii];
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
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned=TCcorr_tuned;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned_index=TCcorrind;
Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_mean=TCcorr_tuned_mean;


%% Rate Overlap

%rate for each cell in each session
for i = sessions
overall_rate_array{i}=Spatialinfo{i}.overall_rate;  
end 
 %Bug missing ROI
for i = sessions
nbROI(i)=size(overall_rate_array{i},2);
nbROI(nbROI==0)=NaN;
maxROI=max(nbROI); %max nb of ROI for each session
while size(overall_rate_array{i},2)<maxROI 
overall_rate_array{i}=[overall_rate_array{i} nan]; %replace missing last ROI with NaN
end
overall_rate(i,:)=overall_rate_array{i};
end
for i=1:size(overall_rate,2)
    for ii=1:sessions
        if overall_rate_array{ii}(i)==0,
           overall_rate_array{ii}(i)=nan;  
if isnan(sum(overall_rate(:,i)))
  overall_rate(:,i)=nan(length(sessions),1);
end
end
end
end

for i = sessions
for ii = sessions 
comprate{i}{ii}=[overall_rate_array{i}; overall_rate_array{ii}] ;
maxrate{i}{ii}=max(comprate{i}{ii});
minrate{i}{ii}=min(comprate{i}{ii});
rate_ratio{i}{ii}=minrate{i}{ii}./maxrate{i}{ii};
rate_overlap{i}(ii)=nanmean(rate_ratio{i}{ii});
rate_overlap_norm{i}=(rate_overlap{i} - min(rate_overlap{i}))/(max(rate_overlap{i}) - min(rate_overlap{i}));
end
end


Spatial_correlation.Rate_overlap=rate_overlap;



end


