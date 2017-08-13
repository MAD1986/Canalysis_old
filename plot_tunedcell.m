

session=1;
ROI=find(SpatiallyTunedCells.TunedROI  ==1)


for i=ROI
    figure
subplot(1,3,1)  

plot(Imaging.trace_restricted_baselinesub(:,i))
title(['ROI' ,num2str(i)]);
hold on; plot(Behavior.resampled.position/200) 
plot(Events.RunningEpochs.run_onset_binary  (:,i)  )

subplot(1,3,2)  
compass(TuningSpecificity.tuning_vector {i})

subplot(1,3,3)  
plot(Spatialinfo{1, 8}.rate_map.rate_map_smoothed (:,i))

end

