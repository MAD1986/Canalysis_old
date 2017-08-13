
function [correlation] = rotation(Spatial_correlation,Spatialinfo,TunedCells);

%Rotational Analysis

%The rotation angle of rate maps between the sessions was determined for each cell that passed the inclusion criteria in both the STD and the MIS sessions. The linearized rate map in the STD session was correlated with the linearized rate map in the MIS session. 
%The linearized rate map in the MIS session was circularly shifted in 1° increments and correlated with the STD session rate map. 
%The shift producing the maximum correlation was assigned as that cell’s rotation angle.

%Cells were classified as “Remap” if they met the place-field inclusion criteria in 
%only the MIS session (i.e., place fields APPEARED) 
%or the STD session (i.e., place fields DISAPPEARED). 

%Cells were classified as “Rotate” if the cells met the inclusion criteria in both the STD and the MIS sessions. 
%To determine whether these place fields rotated CCW (to follow the local cues), 
%rotated CW (to follow the global cues), 
%or had an ambiguous (AMB) relationship to the cues, 
%the STD rate map was correlated with the MIS rate map multiple times as the MIS rate map was rotated relative to the STD. 
%The rotation angle producing the greatest correlation was considered the rotation angle of the place field 

%Get the bin where texture change and swap them until they follow the
%original belt


%% Load the data
sessions=Spatial_correlation.options.sessions;
STDsess=Spatial_correlation.options.STDsessions;
MISsess=Spatial_correlation.options.MISsessions;
criteria1=Spatial_correlation.options.tuned_criteria1;
criteria2=Spatial_correlation.options.rotation_criteria;
binbefore=Spatial_correlation.options.rotation_binbefore;
binafter=Spatial_correlation.options.rotation_binafter;
smth=Spatial_correlation.options.rotation_smooth;



for i = sessions
%rate_map{i}=Spatialinfo{i}.rate_map.normalized_rate_map_smoothed;
%try non smoothed
rate_map{i}=Spatialinfo{i}.event_map.event_map;
rate_map_nonnorm{i}=Spatialinfo{i}.rate_map.rate_map_smoothed ; 
nbROI(i)=size(rate_map{i},2);
maxROI=max(nbROI);
tunedROI{i}=TunedCells{i}.TunedROI;
infoROI{i}=TunedCells{i}.Spatial_info.significantROI;  
tuningROI{i}=TunedCells{i}.TuningSpecificity.significantROI;  
bin_texture{i}=Spatialinfo{i}.bin_texture;
occupancy_map_smoothed{i}=Spatialinfo{i}.occupancy_map.occupancy_map_smoothed;  
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

%% Classify cells
switch criteria1
case 'info'
for i= STDsess
    for ii= MISsess
tunedSTD(i,:)=(infoROI{i});
tunedMIS(ii,:)=(infoROI{ii});
    end
end
    case 'tuning'
for i= STDsess
    for ii= MISsess
tunedSTD(i,:)=(tuningROI{i});
tunedMIS(ii,:)=(tuningROI{i});
    end
end
    case 'both'
for i= STDsess
    for ii= MISsess
tunedSTD(i,:)=(tunedROI{i});
tunedMIS(ii,:)=(tunedROI{i});
    end
end
  otherwise
warning('Unexpected criteria. Use: info / tuning / both');
end

%Rotate : tuned in STD and MIS
switch criteria2
    case 'first'
STD_all=find(sum(tunedSTD(1,:),1)>=1);
MIS_all=find(sum(tunedMIS,1)>=1);
rotateROI=intersect(STD_all,MIS_all);
% -APPEARED = only the MIS session 
appearedROI=setdiff(MIS_all,STD_all);
% -DISAPPEARED =  only the STD session
disappearedROI=setdiff(STD_all,MIS_all);
    case 'or'
STD_all=find(sum(tunedSTD,1)>=1);
MIS_all=find(sum(tunedMIS,1)>=1);
rotateROI=intersect(STD_all,MIS_all);
% -APPEARED = only the MIS session 
appearedROI=setdiff(MIS_all,STD_all);
% -DISAPPEARED =  only the STD session
disappearedROI=setdiff(STD_all,MIS_all);
    case 'and'
STD_all=find(sum(tunedSTD,1)>=length(STDsess));
MIS_all=find(sum(tunedMIS,1)>=1);
rotateROI=intersect(STD_all,MIS_all);
% -APPEARED = only the MIS session 
appearedROI=setdiff(MIS_all,STD_all);
% -DISAPPEARED =  only the STD session
disappearedROI=setdiff(STD_all,MIS_all);
  otherwise
warning('Unexpected criteria. Use: first/and');
end

correlation.rotation.Remap.appearedROI=appearedROI;
correlation.rotation.Remap.nb_appeared=length(appearedROI);
correlation.rotation.Remap.disappearedROI=disappearedROI;
correlation.rotation.Remap.disappearedROI=length(disappearedROI);
correlation.rotation.Rotate.rotateROI=rotateROI;
correlation.rotation.Rotate.rotateROI=length(rotateROI);
%% Rotation
%get the bin where the texture change
    for i= MISsess
bin_texture_STD=bin_texture{1}; %change STDsess
bin_texture_MIS=bin_texture{i};
    end
rotated_tex=abs(bin_texture_STD - bin_texture_MIS)>8;
MIStorotate=bin_texture_MIS(rotated_tex);
MIStogo=bin_texture_STD(rotated_tex);
%Rotate 4bins before, 5bins after (10bins = 20cm);
for i=1:size(MIStorotate,2)
bef_aft_MIS_bins(:,i)=((MIStorotate(i))-binbefore:(MIStorotate(i))+binafter) ;
bef_aft_STD_bins(:,i)=((MIStogo(i))-binbefore:(MIStogo(i))+binafter) ;
end
bef_aft_MIS=reshape(bef_aft_MIS_bins,[],1);
bef_aft_STD=reshape(bef_aft_STD_bins,[],1);
for i= MISsess
rate_map_MIS_torot=rate_map{i}(bef_aft_MIS,:);
rate_map_MIS_togo=rate_map{i}(bef_aft_STD,:);
rate_map_rot=rate_map{i};
end
rate_map_rot(bef_aft_STD,:)=rate_map_rot(bef_aft_MIS,:);

%Rotated rate map
%smooth:
for i=1:size(rate_map_rot,2)
onset_map_rot_smooth(:,i)=Smooth(rate_map_rot(:,i),smth);
end
for i= MISsess
rate_map_rot_smooth=onset_map_rot_smooth./occupancy_map_smoothed{i};
end
%Normalize
for i=1:size(rate_map_rot,2)
rate_map_rot_norm_smooth(:,i)=(rate_map_rot_smooth(:,i) - min(rate_map_rot_smooth(:,i))) / ( max(rate_map_rot_smooth(:,i)) - min(rate_map_rot_smooth(:,i)) );
end

%Other rate map
%Smooth
for i = sessions
for ii=1:size(rate_map{i},2)
onset_map_smooth{i}(:,ii)=Smooth(rate_map{i}(:,ii),smth);
rate_map_smooth{i}=onset_map_smooth{i}./occupancy_map_smoothed{i};
end
end
%Normalize
for i = sessions
for ii=1:size(rate_map_rot,2)
rate_map_norm_smooth{i}(:,ii)=(rate_map_smooth{i}(:,ii) - min(rate_map_smooth{i}(:,ii))) / ( max(rate_map_smooth{i}(:,ii)) - min(rate_map_smooth{i}(:,ii)) );
end
end
correlation.rotation.rotated_rate_map=rate_map_rot_norm_smooth;

%PV matrix
[PVmatrix] = corr(rate_map_norm_smooth{MISsess}',rate_map_norm_smooth{STDsess(1)}', 'rows', 'complete');
[PVmatrix_rotated] = corr(rate_map_rot_norm_smooth',rate_map_norm_smooth{STDsess(1)}', 'rows', 'complete');

figure;imagesc(PVmatrix); colormap('Jet');
figure; imagesc(PVmatrix_rotated); colormap('Jet');

correlation.rotation.PVmatrix.STDvsMIS=PVmatrix;
correlation.rotation.PVmatrix.STDvsMIS_rotated=PVmatrix_rotated;

%PV correlation
for iii=1:size(rate_map_rot_norm_smooth,1)  
[PVcorr(iii)] = corr(rate_map_norm_smooth{MISsess}(iii,:)',rate_map_norm_smooth{STDsess(1)}(iii,:)', 'rows', 'complete');
[PVcorr_rotated(iii)] = corr(rate_map_rot_norm_smooth(iii,:)',rate_map_norm_smooth{STDsess(1)}(iii,:)', 'rows', 'complete');
end
correlation.rotation.PVcorr.STDvsMIS=PVcorr;
correlation.rotation.PVcorr.STDvsMIS_rotated=PVcorr_rotated;
figure;
plot(PVcorr,'g'); hold on; plot(PVcorr_rotated,'r');

%TC corraltion





% for i=rotateROI
% figure; plot(rate_map{2}(:,i)); hold on; plot(rate_map{1}(:,i));%hold on; plot(rate_map{3}(:,i))
% figure; plot(rate_map_rot(:,i)); hold on; plot(rate_map{1}(:,i));%hold on; plot(rate_map{3}(:,i))
% end

%correlation score PV and TC:
%if TC corr higher in rotated : cell follows cued;
%If TC corr lower, new place fields


end



