
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



%Only running trace and position
for i=1:sessions
position{i}=Behavior{i}.resampled.position  ;
runbin{i}=Behavior{i}.runbinary  ;
position_run{i}=position{i}(runbin{i}==1);
trace_run{i}=Imaging{i}.trace_restricted_baselinesub(runbin{i}==1,:);
end
for i=1:sessions
    tic;
dlmwrite(['trace_run' num2str(i) '.txt'],trace_run{i});
dlmwrite(['position_run' num2str(i) '.txt'],position_run{i});
%dlmwrite(['TunedROI' num2str(i) '.txt'], SpatiallyTunedCells{i}.TunedROI  );
%dlmwrite(['trace' num2str(i) '.txt'],Imaging{i}.trace_restricted_baselinesub);
%dlmwrite(['position' num2str(i) '.txt'],Behavior{i}.resampled.position);
%dlmwrite(['texture' num2str(i) '.txt'], Spatialinfo{i}{1, 8}.texture_postion_nonorm);
%dlmwrite(['binary' num2str(i) '.txt'], Events{i}.RunningEpochs.run_onset_binary  );
toc;
end
clear



