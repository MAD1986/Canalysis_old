
function [Spatial_correlation] = rotation_V3(Spatial_correlation,Spatialinfo,TunedCells,TuningSpecificity);

%Rotational Analysis

%The rotation angle of rate maps between the sessions was determined for each cell that passed the inclusion criteria in both the STD and the MIS sessions. 
%The linearized rate map in the STD session was correlated with the linearized rate map in the MIS session. 
%The linearized rate map in the MIS session was circularly shifted in 1° increments and correlated with the STD session rate map. 
%The shift producing the maximum correlation was assigned as that cell’s rotation angle.

%Cells were classified as “Remap” if they met the place-field inclusion criteria in 
%only the MIS session (i.e., place fields APPEARED) 
%or the STD session (i.e., place fields DISAPPEARED). 

%Cells were classified as “Rotate” if the cells met the inclusion criteria in both the STD and the MIS sessions. 
%To determine whether these place fields rotated 
%-to follow the micro textures, 
%-to follow the sequence, 
%-or had an ambiguous (AMB) relationship to the cues, 
%the STD rate map was correlated with the MIS rate map multiple times
%as the MIS rate map was rotated relative to the STD. 
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
tuning_vector{i}=TuningSpecificity{i}.tuning_vector;  
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
while size(tuning_vector{i},2)<maxROI
tuning_vector{i}=[tuning_vector{i} nan];
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
tunedMIS(ii,:)=(tuningROI{ii});
    end
end
    case 'both'
for i= STDsess
    for ii= MISsess
tunedSTD(i,:)=(tunedROI{i});
tunedMIS(ii,:)=(tunedROI{ii});
    end
end
  otherwise
warning('Unexpected criteria. Use: info / tuning / both');
end

%Rotate : tuned in STD and MIS
switch criteria2
    case 'first'
STD_all=find(sum(tunedSTD(1,:),1)>=1);
STD_OR=find(sum(tunedSTD,1)>=1);
MIS_all=find(sum(tunedMIS,1)>=1);
% -ROTATE= STD + MIS (STD1)
rotateROI=intersect(STD_all,MIS_all);
% -ROTATE and DISAPPEAR = MIS + STD1
rotate_disappear_ROI=intersect(STD_1,MIS_all);
% -APPEARED = only the MIS session 
appearedROI=setdiff(MIS_all,STD_OR);
% -DISAPPEARED =  only the STD session
disappearedROI=setdiff(STD_all,MIS_all);
% -APPEARED and STAY = MIS + STD2
appeared_stay_ROI=intersect(STD_2,MIS_all);
    case 'or'
STD_all=find(sum(tunedSTD,1)>=1);
STD_OR=find(sum(tunedSTD,1)>=1);
MIS_all=find(sum(tunedMIS,1)>=1);
% -ROTATE= STD + MIS (STD1 or STD2)
rotateROI=intersect(STD_all,MIS_all);
% -ROTATE and DISAPPEAR = MIS + STD1
rotate_disappear_ROI=intersect(STD_1,MIS_all);
% -APPEARED = only the MIS session 
appearedROI=setdiff(MIS_all,STD_OR);
% -DISAPPEARED =  only the STD session
disappearedROI=setdiff(STD_all,MIS_all);
% -APPEARED and STAY = MIS + STD2
appeared_stay_ROI=intersect(STD_2,MIS_all);
    case 'and'
STD_all=find(sum(tunedSTD,1)>=length(STDsess));
STD_OR=find(sum(tunedSTD,1)>=1);
STD_1=find(tunedSTD(STDsess(1),:)>=1);
STD_2=find(tunedSTD(STDsess(2),:)>=1);
MIS_all=find(sum(tunedMIS,1)>=1);
% -ROTATE= STD + MIS (STD1 +STD2)
rotateROI=intersect(STD_all,MIS_all);
% -ROTATE and DISAPPEAR = only STD1 + MIS
rotate_disappearROI_all=intersect(STD_1,MIS_all);
rotate_disappearROI=setdiff(rotate_disappearROI_all,STD_2);
% -DISAPPEARED =  only the STD session (STD1 +STD2)
disappearedROI=setdiff(STD_all,MIS_all);
% -APPEARED and STAY = MIS + STD2
appeared_stayROI_all=intersect(STD_2,MIS_all);
appeared_stayROI=setdiff(appeared_stayROI_all,STD_1);
%STD1 only
STD1_only=setdiff(STD_1,[STD_2 MIS_all]);
%STD2 only
STD2_only=setdiff(STD_2,[STD_1 MIS_all]);
%MIS only
MIS_only=setdiff(MIS_all,[STD_1 STD_2]);
%STD1 +STD2 only
STD1_STD2_only=setdiff(MIS_all,[STD_1 STD_2]);
% -APPEARED = only the MIS session 
appearedROI=MIS_only;
%STD1 + MIS 
STD1_MIS=intersect(STD_1,MIS_all);
STD1_MIS_only=setdiff(STD1_MIS,STD_2);
%STD2 + MIS 
STD2_MIS=intersect(STD_2,MIS_all);
STD2_MIS_only=setdiff(STD2_MIS,STD_1);
%STD1 + STD2
STD1_STD2=intersect(STD_1,STD_2);
STD1_STD2_only=setdiff(STD1_STD2,MIS_all);
%STD1 + MIS + STD2
STD1_STD2_MIS=intersect(STD1_STD2,MIS_all);




  otherwise
warning('Unexpected criteria. Use: first/and');
end

%Venn diagramm Tuned session 1 , 2 , 3
%venn(Z) plots a Venn diagram with zone areas specified by the vector Z. 
%For a 3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
Zvenn=[length(STD1_only) length(MIS_only) length(STD2_only) length(STD1_MIS_only) length(STD1_STD2_only) length(STD2_MIS_only) length(STD1_STD2_MIS)];
figure;
[H ,S]=venn(Zvenn, 'faceColor',{'b','r','g'});
lgd={['STD1 ' num2str(length(STD1_only)) '/' num2str(length(STD_1))], ['MIS ' num2str(length(MIS_only)) '/' num2str(length(MIS_all))], ['STD2 ' num2str(length(STD2_only)) '/' num2str(num2str(length(STD_2)))], num2str(Zvenn(4))  , num2str(Zvenn(5))  ,  num2str(Zvenn(6)) , num2str(Zvenn(7)) , };
title([{'Nb Tuned Cells  '},{' inclusion criteria =  ', criteria1}]);
 
for i = 1:7
text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [lgd(i) ])
 end



%Organize in structure
Spatial_correlation.rotation.Remap.appearedROI=appearedROI;
Spatial_correlation.rotation.Remap.nb_appeared=length(appearedROI);
Spatial_correlation.rotation.Remap.disappearedROI=disappearedROI;
Spatial_correlation.rotation.Remap.nb_disappearedROI=length(disappearedROI);
Spatial_correlation.rotation.Rotate.rotateROI=rotateROI;
Spatial_correlation.rotation.Rotate.nb_rotateROI=length(rotateROI);
Spatial_correlation.rotation.Rotate.rotate_disappearROI=rotate_disappearROI;
Spatial_correlation.rotation.Rotate.nb_rotate_disappearROI=length(rotate_disappearROI);
Spatial_correlation.rotation.Remap.appeared_stayROI =appeared_stayROI;
Spatial_correlation.rotation.Remap.nb_appeared_stayROI =length(appeared_stayROI);

 
 
% Pie chart -nb of place cells
nb_tuned{1}=length(STD_1);
nb_tuned{2}=length(MIS_all);
nb_tuned{3}=length(STD_2);

figure;
for i=1:3
hold on;
suptitle(['% Tuned Cells  ',' inclusion criteria =  ', criteria1]);
subplot(1,3,i) 
pie_tuned{i}=[nb_tuned{i} (maxROI-nb_tuned{i})];
pie(pie_tuned{i});
title(['Session ', num2str((sessions(i)))]);
labels{i} = {['Tuned Cells n=  ' num2str((pie_tuned{i}(:,1)))] ,['Non Tuned Cells n=  '  num2str((pie_tuned{i}(2)))]};
legend(labels{i},'Location','southoutside','Orientation','vertical');
end

 
%Make a pie chart: type of cells
%1st Rotate and Remap
figure;

subplot(1,2,1) 
pie_class=[length(rotateROI) length(appearedROI) length(disappearedROI) ];
%labels = {['REMAP APPEARED = ' num2str((pie_class(1)/sum(pie_class)*100)) '%' ],['REMAP DISAPPEARED = ' num2str((pie_class(2)/sum(pie_class)*100)) '%'], ['ROTATE = ' num2str((pie_class(3)/sum(pie_class)*100)) '%']};
labels1 = {['ROTATE n=  ' num2str((pie_class(1)))], ['REMAP APPEARED n=  ' num2str((pie_class(2)))] ,['REMAP DISAPPEARED n=  ' num2str((pie_class(3)))] };
pie(pie_class);
legend(labels1,'Location','southoutside','Orientation','vertical')
hold on;
subplot(1,2,2) 
pie_class2=[length(appeared_stayROI) length(rotate_disappearROI)];
labels2 = {['APPEARED & STAY n=  ' num2str((pie_class2(1)))] ,['ROTATE & DISAPPEARED n=  ' num2str((pie_class2(2)))]};
%legend(labels2,'Location','best','Orientation','vertical')
pie(pie_class2);
legend(labels2,'Location','southoutside','Orientation','vertical')



%% Rotation
%get the bin where the texture change
    for i= MISsess
bin_texture_STD=bin_texture{STDsess(1)}; %position of texture on the first session 
bin_texture_MIS=bin_texture{i};
    end
    
    rotated_tex=abs(bin_texture_STD - bin_texture_MIS)>10; %rotated texture= change >10 bins
MIStorotate=bin_texture_MIS(rotated_tex);%Bins to switch on the STD
MIStogo=bin_texture_STD(rotated_tex); %New postions to switch the bins
%Rotate X bins before, X bins after (1bins = 2cm);



for i=1:size(MIStorotate,2)
MIS_bins_newvalues(:,i)=((MIStorotate(i))-binbefore:(MIStorotate(i))+binafter) ;
MIS_bins_toswap(:,i)=((MIStogo(i))-binbefore:(MIStogo(i))+binafter) ;
end
MIS_bins_swap=reshape(MIS_bins_toswap,[],1);
MIS_bins_new=reshape(MIS_bins_newvalues,[],1);

%rate_map_MIS_toswap=rate_map{i}(:,bef_aft_MIS);
%
rate_map_mis=rate_map{MISsess};
rate_map_rot=rate_map{MISsess};
rate_mis=rate_map_nonnorm{MISsess};
rate_rot=rate_map_nonnorm{MISsess};

%for i=1:size(rate_map_mis,2);
    for i=1:size(MIS_bins_swap,1);
rate_map_rot(MIS_bins_swap(i),:)=rate_map_mis(MIS_bins_new(i),:) ;
rate_rot(MIS_bins_swap(i),:)=rate_mis(MIS_bins_new(i),:) ;
end


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
rate_rot_norm(:,i)=(rate_rot(:,i) - min(rate_rot(:,i))) / ( max(rate_rot(:,i)) - min(rate_rot(:,i)) );
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
Spatial_correlation.rotation.rotated_rate_map=rate_map_rot_norm_smooth;

%PV matrix
[PVmatrix] = corr(rate_map_norm_smooth{STDsess(1)}',rate_map_norm_smooth{MISsess}', 'rows', 'complete');
[PVmatrix_rotated] = corr(rate_map_norm_smooth{STDsess(1)}',rate_map_rot_norm_smooth', 'rows', 'complete');
[PVmatrix_rotated_2] = corr(rate_map_norm_smooth{STDsess(1)}',rate_rot_norm', 'rows', 'complete');



figure;imagesc(PVmatrix); colormap('Jet');
figure; imagesc(PVmatrix_rotated); colormap('Jet');
figure; imagesc(PVmatrix_rotated_2); colormap('Jet');


Spatial_correlation.rotation.PVmatrix.STDvsMIS=PVmatrix;
Spatial_correlation.rotation.PVmatrix.STDvsMIS_rotated=PVmatrix_rotated;

%PV correlation
for iii=1:size(rate_map_rot,1)  
[PVcorr(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(iii,:)',rate_map_norm_smooth{MISsess}(iii,:)', 'rows', 'complete');
[PVcorr_rotated(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(iii,:)',rate_map_rot_norm_smooth(iii,:)', 'rows', 'complete');
[PVcorr_rotated_2(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(iii,:)',rate_rot_norm(iii,:)', 'rows', 'complete');
end
Spatial_correlation.rotation.PVcorr.STDvsMIS=PVcorr;
Spatial_correlation.rotation.PVcorr.STDvsMIS_rotated=PVcorr_rotated;
figure;
plot(PVcorr,'g'); hold on; plot(PVcorr_rotated,'r'); plot(PVcorr_rotated_2,'k')


%% Tuning Curve correlation


%Tuning Curve correlation score for cells in ROTATE:
%cells met the inclusion criteria in both the STD and the MIS sessions. 

% APPEARED cells (inclusion criteria in only MIS) show more place field
%in the swapped region?
% DISAPPEARED cells (inclusion cirteria in onl STD) show more place field
%in the swapped region?

%TC corr STD vs MIS rotated higher than STD vs MIS :
%-Cells follow cue 
%look at Ca2+ traces (timing, burst, dendritic activity)

%TC corr STD vs MIS rotated lower than STD vs MIS :
%- if TC corr STD vs MIS higher than STD vs STD

%All neurons
for iii=1:size(rate_map_rot,2) 
[TCcorr(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_map_norm_smooth{MISsess}(:,iii), 'rows', 'complete');
[TCcorr_rotated(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_map_rot_norm_smooth(:,iii), 'rows', 'complete');
[TCcorr_STD(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_map_norm_smooth{STDsess(2)}(:,iii), 'rows', 'complete');
[TCcorr_rotated_2(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_rot_norm(:,iii), 'rows', 'complete');
end

figure;
scatter(TCcorr,TCcorr_rotated)
hold on;  plot([-1 1] ,[-1 1], 'k')
xlabel('Tuning Curve Correlation STD vs MIS');
ylabel('Tuning Curve Correlation STD vs MIS Swapped');
axis([-1 1 -1 1]);
fitTC=lsline; set(fitTC,'color','r', 'LineStyle', '--') 

ROI_cue=find(TCcorr<TCcorr_rotated);
ROI_keep=find(TCcorr>TCcorr_rotated);

%For follow cue cell:
%If TC corr STD1 vs STD2 > STD1 vs MIS = cell follow back in STD2
%If TC corr STD1 vs STD2 < STD1 vs MIS = cell stay at new position 
ROI_back=find(TCcorr<=TCcorr_STD);
ROI_stay=find(TCcorr>=TCcorr_STD);
%For follow cue cell:
ROI_cue_back=intersect( ROI_back,ROI_cue);
ROI_cue_stay=intersect( ROI_stay,ROI_cue);

%ROTATE ROI
ROI_cue_rotate=  intersect( ROI_cue,rotateROI);
ROI_keep_rotate= intersect( ROI_keep,rotateROI);
ROI_cue_back_rotate=  intersect(ROI_cue_back,rotateROI);
ROI_cue_stay_rotate= intersect( ROI_cue_stay,rotateROI);


%Plot rate map and tuning vector of ROTATE cells (inclusion criteria in both the STD and the MIS)
%that follow the cue 


for i=1:size(ROI_cue_rotate,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,ROI_cue_rotate(i)), 'b');
title(['ROTATE - Follow cue ROI', num2str(ROI_cue_rotate(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,ROI_cue_rotate(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,ROI_cue_rotate(i)), 'g');
  plot(rate_rot(:,ROI_cue_rotate(i)), 'k');
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
 end
figure;
suptitle(['ROTATE - Follow cue ROI', num2str(ROI_cue_rotate(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{ROI_cue_rotate(i)}, 'b');
xlabel('STD1')
hold on;
 subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{ROI_cue_rotate(i)}, 'r');
xlabel('MIS')
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{ROI_cue_rotate(i)}, 'g');
xlabel('STD2')
end

%Plot rate map and tuning vector of ROTATE cells (inclusion criteria in both the STD and the MIS)
%that follow the cue and go back
for i=1:size(ROI_cue_back_rotate,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,ROI_cue_back_rotate(i)), 'b');
title(['ROTATE - Follow cue and back ROI', num2str(ROI_cue_back_rotate(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,ROI_cue_back_rotate(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,ROI_cue_back_rotate(i)), 'g');
  

 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
 end
figure;
suptitle(['ROTATE - Follow cue and back ROI', num2str(ROI_cue_back_rotate(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{ROI_cue_back_rotate(i)}, 'b');
xlabel('STD1')
hold on;
 subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{ROI_cue_back_rotate(i)}, 'r');
xlabel('MIS')
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{ROI_cue_back_rotate(i)}, 'g');
xlabel('STD2')
end

%Plot rate map and tuning vector of ROTATE cells (inclusion criteria in both the STD and the MIS)
%that follow the cue and stay
for i=1:size(ROI_cue_stay_rotate,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,ROI_cue_stay_rotate(i)), 'b');
 hold on; plot(rate_map_smooth{MISsess(1)}(:,ROI_cue_stay_rotate(i)), 'r');
 title(['ROTATE - Follow cue and stay ROI', num2str(ROI_cue_stay_rotate(i))]);
 plot(rate_map_smooth{STDsess(2)}(:,ROI_cue_stay_rotate(i)), 'g');
  plot(rate_rot(:,ROI_cue_rotate(i)), 'k');
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
 end
figure;
suptitle(['ROTATE - Follow cue and stay ROI', num2str(ROI_cue_stay_rotate(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{ROI_cue_stay_rotate(i)}, 'b');
xlabel('STD1')
hold on;
  subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{ROI_cue_stay_rotate(i)}, 'r');
xlabel('MIS')
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{ROI_cue_stay_rotate(i)}, 'g');
xlabel('STD2')
end


%Plot rate map and tuning vector of APPEARED cells (inclusion criteria in only the MIS)
for i=1:size(appearedROI,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,appearedROI(i)), 'b');
title(['APPEARED', num2str(appearedROI(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,appearedROI(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,appearedROI(i)), 'g');
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
 end
figure;
suptitle(['APPEARED', num2str(appearedROI(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{appearedROI(i)}, 'b');
xlabel('STD1')
hold on;
  subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{appearedROI(i)}, 'r');
xlabel('MIS')
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{appearedROI(i)}, 'g');
xlabel('STD2')
end

%Plot rate map and tuning vector of APPEARED and STAY cells (inclusion criteria in only the MIS and STD2)
for i=1:size(appeared_stayROI,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,appeared_stayROI(i)), 'b');
title(['APPEARED and STAY ROI', num2str(appeared_stayROI(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,appeared_stayROI(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,appeared_stayROI(i)), 'g');
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
 end
figure;
suptitle(['APPEARED and STAY ROI', num2str(appeared_stayROI(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{appeared_stayROI(i)}, 'b');
xlabel('STD1')
hold on;
  subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{appeared_stayROI(i)}, 'r');
xlabel('MIS')
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{appeared_stayROI(i)}, 'g');
xlabel('STD2')
end

%Plot rate map and tuning vector of ROTATE and DISAPPEAR cells (inclusion criteria in only the STD1 and MIS)
for i=1:size(rotate_disappearROI,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,rotate_disappearROI(i)), 'b');
title(['ROTATE and DISAPPEAR ROI', num2str(rotate_disappearROI(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,rotate_disappearROI(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,rotate_disappearROI(i)), 'g');
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
 end
figure;
suptitle(['ROTATE and DISAPPEAR ROI', num2str(rotate_disappearROI(i))] );
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{rotate_disappearROI(i)}, 'b');
xlabel('STD1')
hold on;
  subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{rotate_disappearROI(i)}, 'r');
xlabel('MIS')
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{rotate_disappearROI(i)}, 'g');
xlabel('STD2')
end

%% Find Ca2+ events in swapped bins

%rate_map






%% NOT USED
%Find cells whith tuningvector around swapped textures 
%angles of the swapped textures: 
%angle_text=(bef_aft_STD_bins/100)*360; %OR with texture position 
%angle_tex_inter=[([angle_text(1,1) angle_text(end,1)]);([angle_text(1,2) angle_text(end,2)])];
%angles of tuning vector all neurons 
%for ii=sessions
%tuningspe{ii}=TuningSpecificity{ii}.tuning_vector_specificity(:); 
%rX = real(tuningspe{ii});
%rY = imag(tuningspe{ii});
%angle_tuning{ii} = atan2d(rY,rX);
%for i=1:size(angle_tuning{ii},1)
%if (angle_tuning{ii}(i) < 0)
%            angle_tuning{ii}(i) = angle_tuning{ii}(i) + 360;
%end
%end
%end
%Find neurons with tuning vector in the swapped microtextures
%for ii=sessions
%[~,interval{ii},~] = InIntervals(angle_tuning{ii},angle_tex_inter);
%ROI_tuning_intex{ii}=find(interval{ii}>=1) ;  
%end

%Plot cells with Tuning Vector in swapped texture (all cells - not tuned)
%for i=1:size(ROI_tuning_intex{1},1)  
%    figure;
%suptitle('Tuning Vector in swapped texture')
%subplot(1,3,1)     
%compass(TuningSpecificity{STDsess(1)}.tuning_vector{ROI_tuning_intex{1}(i)}, 'b');
%hold on;
%compass(TuningSpecificity{STDsess(1)}.tuning_vector_specificity(ROI_tuning_intex{1}(i)), 'k');
%xlabel('STD1')
%hold on;
%  subplot(1,3,2)
%compass(TuningSpecificity{MISsess(1)}.tuning_vector{ROI_tuning_intex{1}(i)}, 'r');
%hold on;
%compass(TuningSpecificity{MISsess(1)}.tuning_vector_specificity(ROI_tuning_intex{1}(i)),'k');

%xlabel('MIS')
% subplot(1,3,3)
% compass(TuningSpecificity{STDsess(2)}.tuning_vector{ROI_tuning_intex{1}(i)}, 'g');
% hold on;
% compass(TuningSpecificity{STDsess(2)}.tuning_vector_specificity(ROI_tuning_intex{1}(i)), 'k');
% xlabel('STD2')
%end


    
end



