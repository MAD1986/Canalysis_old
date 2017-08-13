
function [correlation] = PVcorrelation(correlation,Spatialinfo);

%% Population correlation: PV correlation
%PV correlation must be with normalized firing rate to avoid cell with high firing
%rate to weight too much

sessions=correlation.options.sessions;

%% For all neurons
for i = sessions
rate_map{i}=Spatialinfo{i}.rate_map.normalized_rate_map_smoothed;
nbROI(i)=size(rate_map{i},2);
maxROI=max(nbROI);
end
for i = sessions
while size(rate_map{i},2)<maxROI
rate_map{i}=[rate_map{i} nan(size(rate_map{i},1),1)]; %replace missing last ROI with NaN
end
end

for i = sessions
for ii= sessions
[PVmatrix{i}{ii}] = corr(rate_map{i}',rate_map{ii}', 'rows', 'complete');
matrixnb{i}{ii}=[i ii];
end
end 

correlation.PVcorrelation.matrix=PVmatrix;
correlation.PVcorrelation.index=matrixnb;




end
