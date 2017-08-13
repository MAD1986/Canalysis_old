


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






