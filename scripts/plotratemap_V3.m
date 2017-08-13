function [Spatial_correlation] = plotratemap_V3(Spatial_correlation,Spatialinfo, TunedCells)
sessions=Spatial_correlation.options.sessions;
order=Spatial_correlation.options.order;
sm=Spatial_correlation.options.ratemapsm;



%% Plot all neurons
for i = sessions
ordered_rate_map{i}=Spatialinfo{i}.rate_map.ordered_rate_map; 
rate_map{i}=Spatialinfo{i}.rate_map.normalized_rate_map_smoothed;
ROI_rate_map{i}=Spatialinfo{i}.rate_map.ROI_odered_rate_map;
end 

%Bug with missing last ROI
for i = sessions
nbROI(i)=size(rate_map{i},2);
nbROI(nbROI==0)=NaN;
maxROI=max(nbROI); %max nb of ROI for each session
%misROI{i}=setdiff(1:nbROI,ROI_rate_map{i});%find the missing ROI for each session
semisROI=nbROI<maxROI;
missingROI=min(nbROI)+1:maxROI;
while size(rate_map{i},2)<maxROI 
rate_map{i}=[rate_map{i} nan(size(rate_map{i},1),1)]; %replace missing last ROI with NaN
end
while size(ordered_rate_map{i},1)<maxROI 
ordered_rate_map{i}=[ordered_rate_map{i}; (nan(size(ordered_rate_map{i},2),1))'];
end    
while size(ROI_rate_map{i},1)<maxROI 
ROI_rate_map{i}=[ROI_rate_map{i}; nan];
end
end  

for i = 1:size(ROI_rate_map{order},1);
if isnan(ordered_rate_map{order}(i,:))
    ROI_rate_map{order}(i)=NaN;
ROI_rate_map_nonan=ROI_rate_map{order}(~any(isnan(ROI_rate_map{order}),2),:);
elseif isnan(ordered_rate_map{order}(i,:))==0
   ROI_rate_map_nonan=ROI_rate_map{order}(~any(isnan(ROI_rate_map{order}),2),:); 
end
end

if any(semisROI(order))==0
misind=find(ismember(ROI_rate_map_nonan,missingROI));
ROI_rate_map_nonan(misind)=NaN;
end
ROI_rate_map_nonan=ROI_rate_map_nonan(~any(isnan(ROI_rate_map_nonan),2),:);

for ii = sessions 
for i = 1:size(ROI_rate_map_nonan,1);
ordered_map{ii}(i,:)=rate_map{ii}(:,ROI_rate_map_nonan(i));
end
end

%IF Smooth 
if sm>0
for i=sessions
for ii=1:size(ordered_map{i},1)
ordered_map_sm{i}(ii,:)=Smooth(ordered_map{i}(ii,:),sm);
end
end
end
%Plot figure
if sm>0
for i = sessions
figure;
imagesc(ordered_map_sm{i});
end
end


%Plot figure no smooth
if sm<=0
for i = sessions
figure;
imagesc(ordered_map{i});
end
end


Spatial_correlation.ordered_map=ordered_map;



