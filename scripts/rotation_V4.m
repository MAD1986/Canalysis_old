
function [Spatial_correlation] = rotation_V4(Spatial_correlation,Spatialinfo,TunedCells,TuningSpecificity);

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
fig=Spatial_correlation.options.rotation_figure;



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
onset_bin{i}=Spatialinfo{i}.onset_bin  ;
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

%Name Sessions
name_STD1=['STD',num2str(STDsess(1))];
name_STD2=['STD',num2str(STDsess(2))];
name_MIS=['MIS',num2str(MISsess-STDsess(1))];
name_STD1_STD2=[name_STD1, 'and' name_STD2];
name_STD1_MIS=[name_STD1, 'and' name_MIS];
name_STD2_MIS=[name_STD2, 'and' name_MIS];
name_STD1_MIS_STD2=[name_STD1, 'and' ,name_MIS, 'and', name_STD2];

%Organize in structure
%STD1 and STD1 only
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1)]).(name_STD1)=STD_1;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1)]).(['nb_' (name_STD1)])=length(STD_1);
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1)]).([(name_STD1) '_only'])=STD1_only;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1)]).(['nb_' (name_STD1) '_only'])=length(STD1_only);
%STD2 and STD2 only
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2)]).(name_STD2)=STD_2;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2)]).(['nb_' (name_STD2)])=length(STD_2);
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2)]).([(name_STD2) '_only'])=STD2_only;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2)]).(['nb_' (name_STD2) '_only'])=length(STD2_only);
%MIS and MIS only
Spatial_correlation.rotation.Class.(['Tuned_' (name_MIS)]).(name_MIS)=MIS_all;
Spatial_correlation.rotation.Class.(['Tuned_' (name_MIS)]).(['nb_' (name_MIS)])=length(MIS_all);
Spatial_correlation.rotation.Class.(['Tuned_' (name_MIS)]).([(name_MIS) '_only'])=MIS_only;
Spatial_correlation.rotation.Class.(['Tuned_' (name_MIS)]).(['nb_' (name_MIS) '_only'])=length(MIS_only);
%STD1 and STD2
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_STD2)]).(name_STD1_STD2)=STD1_STD2;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_STD2)]).(['nb_' (name_STD1_STD2)])=length(STD1_STD2);
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_STD2)]).([(name_STD1_STD2) '_only'])=STD1_STD2_only;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_STD2)]).(['nb_' (name_STD1_STD2) '_only'])=length(STD1_STD2_only);
%STD1 and MIS 
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS)]).(name_STD1_MIS)=STD1_MIS;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS)]).(['nb_' (name_STD1_MIS)])=length(STD1_MIS);
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS)]).([(name_STD1_MIS) '_only'])=STD1_MIS_only;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS)]).(['nb_' (name_STD1_MIS) '_only'])=length(STD1_MIS_only);
%STD2 and MIS 
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2_MIS)]).(name_STD2_MIS)=STD2_MIS;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2_MIS)]).(['nb_' (name_STD2_MIS)])=length(STD2_MIS);
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2_MIS)]).([(name_STD2_MIS) '_only'])=STD2_MIS_only;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2_MIS)]).(['nb_' (name_STD2_MIS) '_only'])=length(STD2_MIS_only);
%STD1 and STD2 and MIS 
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS_STD2)]).(name_STD1_MIS_STD2)=STD1_STD2_MIS;
Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS_STD2)]).(['nb_' (name_STD1_MIS_STD2)])=length(STD1_STD2_MIS);



%Venn diagramm Tuned session 1 , 2 , 3
%venn(Z) plots a Venn diagram with zone areas specified by the vector Z. 
%For a 3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
Zvenn=[length(STD1_only) length(MIS_only) length(STD2_only) length(STD1_MIS_only) length(STD1_STD2_only) length(STD2_MIS_only) length(STD1_STD2_MIS)];
figure;
[H ,S]=venn(Zvenn, 'faceColor',{'b','r','g'});
lgd={['STD1 ' num2str(length(STD1_only)) '/' num2str(length(STD_1))], ['MIS ' num2str(length(MIS_only)) '/' num2str(length(MIS_all))], ['STD2 ' num2str(length(STD2_only)) '/' num2str(num2str(length(STD_2)))], num2str(Zvenn(4))  , num2str(Zvenn(5))  ,  num2str(Zvenn(6)) , num2str(Zvenn(7)) , };
lgd={[name_STD1 '  ' num2str(length(STD1_only)) '/' num2str(length(STD_1))], [name_MIS '  ' num2str(length(MIS_only)) '/' num2str(length(MIS_all))], [name_STD2 '  ' num2str(length(STD2_only)) '/' num2str(num2str(length(STD_2)))], num2str(Zvenn(4))  , num2str(Zvenn(5))  ,  num2str(Zvenn(6)) , num2str(Zvenn(7)) , };
title([{'Nb Tuned Cells  '},{' inclusion criteria =  ', criteria1}]);
 for i = 1:7
text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [lgd(i) ])
 end
 
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
%figure;
%subplot(1,2,1) 
%pie_class=[length(rotateROI) length(appearedROI) length(disappearedROI) ];
%labels1 = {['ROTATE n=  ' num2str((pie_class(1)))], ['REMAP APPEARED n=  ' num2str((pie_class(2)))] ,['REMAP DISAPPEARED n=  ' num2str((pie_class(3)))] };
%pie(pie_class);
%legend(labels1,'Location','southoutside','Orientation','vertical')
%hold on;
%subplot(1,2,2) 
%pie_class2=[length(appeared_stayROI) length(rotate_disappearROI)];
%labels2 = {['APPEARED & STAY n=  ' num2str((pie_class2(1)))] ,['ROTATE & DISAPPEARED n=  ' num2str((pie_class2(2)))]};
%pie(pie_class2);
%legend(labels2,'Location','southoutside','Orientation','vertical')



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

Spatial_correlation.rotation.bin_texture = bin_texture_STD; 
Spatial_correlation.rotation.swap_texture = rotated_tex; 


for i=1:size(MIStorotate,2)
MIS_bins_newvalues(:,i)=((MIStorotate(i))-binbefore:(MIStorotate(i))+binafter) ;
MIS_bins_toswap(:,i)=((MIStogo(i))-binbefore:(MIStogo(i))+binafter) ;
end
%If bin to swap is 1 go back to end bins
for j=1:size(MIStorotate,2)
for i=1:size(MIS_bins_newvalues,1)
for ii=1:size(MIS_bins_toswap,1)
if MIS_bins_newvalues(i,j)<1, 
MIS_bins_newvalues(i,j)=size(rate_map{MISsess},1)+MIS_bins_newvalues(i,j);
end
if MIS_bins_toswap(ii,j)<1,
MIS_bins_toswap(ii,j)=size(rate_map{MISsess},1)+MIS_bins_toswap(ii,j);
end
end
end   
end
MIS_bins_swap=reshape(MIS_bins_toswap,[],1);
MIS_bins_new=reshape(MIS_bins_newvalues,[],1);

Spatial_correlation.rotation.bin_texture_swap =MIS_bins_swap; 


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
[PVmatrix_STD1_MIS] = corr(rate_map_norm_smooth{STDsess(1)}',rate_map_norm_smooth{MISsess}', 'rows', 'complete');
[PVmatrix_STD1_MIS_swap] = corr(rate_map_norm_smooth{STDsess(1)}',rate_map_rot_norm_smooth', 'rows', 'complete');
[PVmatrix_MIS_STD2] = corr(rate_map_norm_smooth{MISsess}', rate_map_norm_smooth{STDsess(2)}','rows', 'complete');
[PVmatrix_MIS_swap_STD2] = corr(rate_map_rot_norm_smooth', rate_map_norm_smooth{STDsess(2)}','rows', 'complete');
[PVmatrix_STD1_STD2] = corr(rate_map_norm_smooth{STDsess(1)}',rate_map_norm_smooth{STDsess(2)}', 'rows', 'complete');

%Figure
figure;imagesc(PVmatrix_STD1_MIS); colormap('Jet'); 
%hold on; 
%plot([1 size(PVmatrix_STD1_MIS,1)], [2 2], 'r')

xlabel(name_STD1); ylabel(name_MIS)
figure; imagesc(PVmatrix_STD1_MIS_swap); colormap('Jet');
xlabel(name_STD1); ylabel([name_MIS 'Swap'])
figure; imagesc(PVmatrix_MIS_STD2); colormap('Jet');
xlabel(name_STD2); ylabel(name_MIS)
figure; imagesc(PVmatrix_MIS_swap_STD2); colormap('Jet');
xlabel(name_STD2); ylabel([name_MIS 'Swap'])
figure; imagesc(PVmatrix_STD1_STD2); colormap('Jet');
xlabel(name_STD1); ylabel(name_STD2)

%Organize in structure
Spatial_correlation.rotation.PVmatrix.([name_STD1 '_vs_' name_MIS])=PVmatrix_STD1_MIS;
Spatial_correlation.rotation.PVmatrix.([name_STD1 '_vs_' name_MIS 'swap'])=PVmatrix_STD1_MIS_swap;
Spatial_correlation.rotation.PVmatrix.([name_MIS '_vs_' name_STD2])=PVmatrix_MIS_STD2;
Spatial_correlation.rotation.PVmatrix.([name_MIS 'swap' '_vs_' name_STD2])=PVmatrix_MIS_swap_STD2;
Spatial_correlation.rotation.PVmatrix.([name_STD1 '_vs_' name_STD2])=PVmatrix_STD1_STD2;



%PV correlation
for iii=1:size(rate_map_rot,1)  
[PVcorr_STD1_MIS(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(iii,:)',rate_map_norm_smooth{MISsess}(iii,:)', 'rows', 'complete');
[PVcorr_STD1_MIS_swap(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(iii,:)',rate_map_rot_norm_smooth(iii,:)', 'rows', 'complete');
[PVcorr_MIS_STD2(iii)] = corr(rate_map_norm_smooth{MISsess}(iii,:)',rate_map_norm_smooth{STDsess(2)}(iii,:)', 'rows', 'complete');
[PVcorr_MIS_swap_STD2(iii)] = corr(rate_map_rot_norm_smooth(iii,:)',rate_map_norm_smooth{STDsess(2)}(iii,:)', 'rows', 'complete');
[PVcorr_STD1_STD2(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(iii,:)',rate_map_norm_smooth{STDsess(2)}(iii,:)', 'rows', 'complete');
end

%Figure
figure;
plot(PVcorr_STD1_MIS,'g'); hold on; plot(PVcorr_STD1_MIS_swap,'r');
legend([name_STD1,' vs ', name_MIS], [name_STD1,' vs ', name_MIS, 'swap'])
title('Population Vector Pearson Correlation');
figure;
plot(PVcorr_MIS_STD2,'g'); hold on; plot(PVcorr_MIS_swap_STD2,'r');
legend([name_MIS,' vs ', name_STD2], [name_MIS 'swap',' vs ', name_STD2])
title('Population Vector Pearson Correlation');
figure;
plot(PVcorr_STD1_STD2,'r'); hold on; plot(PVcorr_STD1_MIS,'b');plot(PVcorr_MIS_STD2,'g');
legend([name_STD1,' vs ', name_STD2], [name_STD1 ,' vs ', name_MIS], [name_MIS ,' vs ', name_STD2])
title('Population Vector Pearson Correlation');

%Organize in structure
Spatial_correlation.rotation.PVcorr.([name_STD1 '_vs_' name_MIS])=PVcorr_STD1_MIS;
Spatial_correlation.rotation.PVcorr.([name_STD1 '_vs_' name_MIS 'swap'])=PVcorr_STD1_MIS_swap;
Spatial_correlation.rotation.PVcorr.([name_MIS '_vs_' name_STD2])=PVcorr_MIS_STD2;
Spatial_correlation.rotation.PVcorr.([name_MIS 'swap' '_vs_' name_STD2])=PVcorr_MIS_swap_STD2;
Spatial_correlation.rotation.PVcorr.([name_STD1 '_vs_' name_STD2])=PVcorr_STD1_STD2;





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
[TCcorr_STD1_MIS(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_map_norm_smooth{MISsess}(:,iii), 'rows', 'complete');
[TCcorr_STD1_MIS_swap(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_map_rot_norm_smooth(:,iii), 'rows', 'complete');
[TCcorr_MIS_STD2(iii)] = corr(rate_map_norm_smooth{MISsess}(:,iii), rate_map_norm_smooth{STDsess(2)}(:,iii), 'rows', 'complete');
[TCcorr_MIS_swap_STD2(iii)] = corr(rate_map_rot_norm_smooth(:,iii), rate_map_norm_smooth{STDsess(2)}(:,iii), 'rows', 'complete');
[TCcorr_STD1_STD2(iii)] = corr(rate_map_norm_smooth{STDsess(1)}(:,iii),rate_map_norm_smooth{STDsess(2)}(:,iii), 'rows', 'complete');
end

%Organize in structure
Spatial_correlation.rotation.TCcorr.([name_STD1 '_vs_' name_MIS])=TCcorr_STD1_MIS;
Spatial_correlation.rotation.TCcorr.([name_STD1 '_vs_' name_MIS 'swap'])=TCcorr_STD1_MIS_swap;
Spatial_correlation.rotation.TCcorr.([name_MIS '_vs_' name_STD2])=TCcorr_MIS_STD2;
Spatial_correlation.rotation.TCcorr.([name_MIS 'swap' '_vs_' name_STD2])=TCcorr_MIS_swap_STD2;
Spatial_correlation.rotation.TCcorr.([name_STD1 '_vs_' name_STD2])=TCcorr_STD1_STD2;


% Cue if Tunning Curve Correlation for MIS swap vs STD higher than MIS vs STD 
Cue_STD1_MIS=find(TCcorr_STD1_MIS_swap>TCcorr_STD1_MIS);
NoCue_STD1_MIS=find(TCcorr_STD1_MIS>TCcorr_STD1_MIS_swap);
Egal_STD1_MIS=find(TCcorr_STD1_MIS==TCcorr_STD1_MIS_swap);
Cue_STD2_MIS=find(TCcorr_MIS_swap_STD2>TCcorr_MIS_STD2);
NoCue_STD2_MIS=find(TCcorr_MIS_STD2>TCcorr_MIS_swap_STD2);
Egal_STD2_MIS=find(TCcorr_MIS_STD2==TCcorr_MIS_swap_STD2);

Cue_STD1_STD2_MIS = unique([Cue_STD1_MIS Cue_STD2_MIS]);
NoCue_STD1_STD2_MIS = unique([NoCue_STD1_MIS NoCue_STD2_MIS]);
Egal_STD1_STD2_MIS = unique([Egal_STD1_MIS Egal_STD2_MIS]);

nb_Cue_STD1_MIS=length(Cue_STD1_MIS);
nb_NoCue_STD1_MIS=length(NoCue_STD1_MIS);
nb_Egal_STD1_MIS=length(Egal_STD1_MIS);
nb_Cue_STD2_MIS=length(Cue_STD2_MIS);
nb_NoCue_STD2_MIS=length(NoCue_STD2_MIS);
nb_Egal_STD2_MIS=length(Egal_STD2_MIS);

nb_Cue_STD1_STD2_MIS = length(Cue_STD1_STD2_MIS);
nb_NoCue_STD1_STD2_MIS = length(NoCue_STD1_STD2_MIS);
nb_Egal_STD1_STD2_MIS = length(Egal_STD1_STD2_MIS);

%Organize in structure
Spatial_correlation.rotation.Cue.AllROI.([name_STD1 '_vs_' name_MIS])=Cue_STD1_MIS;
Spatial_correlation.rotation.Cue.AllROI.(['nb_' name_STD1 '_vs_' name_MIS])=nb_Cue_STD1_MIS;
Spatial_correlation.rotation.NoCue.AllROI.([name_STD1 '_vs_' name_MIS])=NoCue_STD1_MIS;
Spatial_correlation.rotation.NoCue.AllROI.(['nb_' name_STD1 '_vs_' name_MIS])=nb_NoCue_STD1_MIS;
Spatial_correlation.rotation.Egal.AllROI.([name_STD1 '_vs_' name_MIS])=Egal_STD1_MIS;
Spatial_correlation.rotation.Egal.AllROI.(['nb_' name_STD1 '_vs_' name_MIS])=nb_Egal_STD1_MIS;
Spatial_correlation.rotation.Cue.AllROI.([name_MIS '_vs_' name_STD2])=Cue_STD2_MIS;
Spatial_correlation.rotation.Cue.AllROI.(['nb_' name_MIS '_vs_' name_STD2])=nb_Cue_STD2_MIS;
Spatial_correlation.rotation.NoCue.AllROI.([name_MIS '_vs_' name_STD2])=NoCue_STD2_MIS;
Spatial_correlation.rotation.NoCue.AllROI.(['nb_' name_MIS '_vs_' name_STD2])=nb_NoCue_STD2_MIS;
Spatial_correlation.rotation.Egal.AllROI.([name_MIS '_vs_' name_STD2])=Egal_STD2_MIS;
Spatial_correlation.rotation.Egal.AllROI.(['nb_' name_MIS '_vs_' name_STD2])=nb_Egal_STD2_MIS;

Spatial_correlation.rotation.Cue.AllROI.([name_STD1_STD2 '_vs_' name_MIS])=Cue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.AllROI.(['nb_' name_STD1_STD2 '_vs_' name_MIS])=nb_Cue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.AllROI.([name_STD1_STD2 '_vs_' name_MIS])=NoCue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.AllROI.(['nb_' name_STD1_STD2 '_vs_' name_MIS])=nb_NoCue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.AllROI.([name_STD1_STD2 '_vs_' name_MIS])=Egal_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.AllROI.(['nb_' name_STD1_STD2 '_vs_' name_MIS])=nb_Egal_STD1_STD2_MIS;

%Tuned ROI 
Tuned_Cue_STD1_MIS=intersect(Cue_STD1_MIS, STD1_MIS);
Tuned_NoCue_STD1_MIS=intersect(NoCue_STD1_MIS, STD1_MIS);
Tuned_Egal_STD1_MIS=intersect(Egal_STD1_MIS, STD1_MIS);
Tuned_Cue_STD2_MIS=intersect(Cue_STD2_MIS, STD2_MIS);
Tuned_NoCue_STD2_MIS=intersect(NoCue_STD2_MIS, STD2_MIS);
Tuned_Egal_STD2_MIS=intersect(Egal_STD2_MIS, STD2_MIS);

Tuned_Cue_STD1_STD2_MIS = unique([Tuned_Cue_STD1_MIS Tuned_Cue_STD2_MIS]);
Tuned_NoCue_STD1_STD2_MIS = unique([Tuned_NoCue_STD1_MIS Tuned_NoCue_STD2_MIS]);
Tuned_Egal_STD1_STD2_MIS = unique([Tuned_Egal_STD1_MIS Tuned_Egal_STD2_MIS]);

nb_Tuned_Cue_STD1_MIS=length(Tuned_Cue_STD1_MIS);
nb_Tuned_NoCue_STD1_MIS=length(Tuned_NoCue_STD1_MIS);
nb_Tuned_Egal_STD1_MIS=length(Tuned_Egal_STD1_MIS);
nb_Tuned_Cue_STD2_MIS=length(Tuned_Cue_STD2_MIS);
nb_Tuned_NoCue_STD2_MIS=length(Tuned_NoCue_STD2_MIS);
nb_Tuned_Egal_STD2_MIS=length(Tuned_Egal_STD2_MIS);

nb_Tuned_Cue_STD1_STD2_MIS = length(Tuned_Cue_STD1_STD2_MIS);
nb_Tuned_NoCue_STD1_STD2_MIS = length(Tuned_NoCue_STD1_STD2_MIS);
nb_Tuned_Egal_STD1_STD2_MIS = length(Tuned_Egal_STD1_STD2_MIS);

%Organize in structure
Spatial_correlation.rotation.Cue.TunedROI.([name_STD1 '_vs_' name_MIS])=Tuned_Cue_STD1_MIS;
Spatial_correlation.rotation.Cue.TunedROI.(['nb_' name_STD1 '_vs_' name_MIS])=nb_Tuned_Cue_STD1_MIS;
Spatial_correlation.rotation.NoCue.TunedROI.([name_STD1 '_vs_' name_MIS])=Tuned_NoCue_STD1_MIS;
Spatial_correlation.rotation.NoCue.TunedROI.(['nb_' name_STD1 '_vs_' name_MIS])=nb_Tuned_NoCue_STD1_MIS;
Spatial_correlation.rotation.Egal.TunedROI.([name_STD1 '_vs_' name_MIS])=Tuned_Egal_STD1_MIS;
Spatial_correlation.rotation.Egal.TunedROI.(['nb_' name_STD1 '_vs_' name_MIS])=nb_Tuned_Egal_STD1_MIS;
Spatial_correlation.rotation.Cue.TunedROI.([name_MIS '_vs_' name_STD2])=Tuned_Cue_STD2_MIS;
Spatial_correlation.rotation.Cue.TunedROI.(['nb_' name_MIS '_vs_' name_STD2])=nb_Tuned_Cue_STD2_MIS;
Spatial_correlation.rotation.NoCue.TunedROI.([name_MIS '_vs_' name_STD2])=Tuned_NoCue_STD2_MIS;
Spatial_correlation.rotation.NoCue.TunedROI.(['nb_' name_MIS '_vs_' name_STD2])=nb_Tuned_NoCue_STD2_MIS;
Spatial_correlation.rotation.Egal.TunedROI.([name_MIS '_vs_' name_STD2])=Tuned_Egal_STD2_MIS;
Spatial_correlation.rotation.Egal.TunedROI.(['nb_' name_MIS '_vs_' name_STD2])=nb_Tuned_Egal_STD2_MIS;

Spatial_correlation.rotation.Cue.TunedROI.([name_STD1_STD2 '_vs_' name_MIS])=Tuned_Cue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.TunedROI.(['nb_' name_STD1_STD2 '_vs_' name_MIS])=nb_Tuned_Cue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.TunedROI.([name_STD1_STD2 '_vs_' name_MIS])=Tuned_NoCue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.TunedROI.(['nb_' name_STD1_STD2 '_vs_' name_MIS])=nb_Tuned_NoCue_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.TunedROI.([name_STD1_STD2 '_vs_' name_MIS])=Tuned_Egal_STD1_STD2_MIS;
Spatial_correlation.rotation.Cue.TunedROI.(['nb_' name_STD1_STD2 '_vs_' name_MIS])=nb_Tuned_Egal_STD1_STD2_MIS;



if fig==1,
%Correlation STD1/MIS VS STD1/MIS Swap
figure;
scatter(TCcorr_STD1_MIS,TCcorr_STD1_MIS_swap);
hold on; scatter(TCcorr_STD1_MIS(Cue_STD1_MIS),TCcorr_STD1_MIS_swap(Cue_STD1_MIS), 'g');
hold on; scatter(TCcorr_STD1_MIS(NoCue_STD1_MIS),TCcorr_STD1_MIS_swap(NoCue_STD1_MIS), 'r');
hold on; scatter(TCcorr_STD1_MIS(STD1_MIS),TCcorr_STD1_MIS_swap(STD1_MIS), 'filled', 'k');
legend(['N= ',num2str(length(Egal_STD1_MIS)), '  (tuned = ',num2str(length(Tuned_Egal_STD1_MIS)), ')' ],['N= ',num2str(length(Cue_STD1_MIS)), '  (tuned = ',num2str(length(Tuned_Cue_STD1_MIS)), ')' ],['N= ',num2str(length(NoCue_STD1_MIS)), '  (tuned = ',num2str(length(Tuned_NoCue_STD1_MIS)), ')' ], ['Tuned ROI  ', name_STD1, ' & ' name_MIS]);
plot([-1 1] ,[-1 1], 'k'); 
xlabel(['Tuning Curve Correlation  ', name_STD1, ' VS ' name_MIS]);
ylabel(['Tuning Curve Correlation  ',  name_STD1, ' VS ' name_MIS 'swap']);
axis([-1 1 -1 1]);
%fitTC=lsline; set(fitTC,'color','r', 'LineStyle', '--') ;

%Correlation MIS/STD2 VS MIS/STD2 Swap
figure;
scatter(TCcorr_MIS_STD2,TCcorr_MIS_swap_STD2);
hold on; scatter(TCcorr_MIS_STD2(Cue_STD2_MIS),TCcorr_MIS_swap_STD2(Cue_STD2_MIS), 'g');
hold on; scatter(TCcorr_MIS_STD2(NoCue_STD2_MIS),TCcorr_MIS_swap_STD2(NoCue_STD2_MIS), 'r');
hold on; scatter(TCcorr_MIS_STD2(STD2_MIS),TCcorr_MIS_swap_STD2(STD2_MIS), 'filled', 'k');
legend(['N= ',num2str(length(Egal_STD2_MIS)), '  (tuned = ',num2str(length(Tuned_Egal_STD2_MIS)), ')' ],['N= ',num2str(length(Cue_STD2_MIS)), '  (tuned = ',num2str(length(Tuned_Cue_STD2_MIS)), ')' ],['N= ',num2str(length(NoCue_STD2_MIS)), '  (tuned = ',num2str(length(Tuned_NoCue_STD2_MIS)), ')' ], ['Tuned ROI  ', name_STD2, ' & ' name_MIS]);
plot([-1 1] ,[-1 1], 'k'); 
xlabel(['Tuning Curve Correlation  ', name_STD2, ' VS ' name_MIS]);
ylabel(['Tuning Curve Correlation  ',  name_STD2, ' VS ' name_MIS 'swap']);
axis([-1 1 -1 1]);
%fitTC=lsline; set(fitTC,'color','r', 'LineStyle', '--') ;



 %Plot rate map and tuning vector of CUE tuned ROI (STD1 and STD2);

for i=1:size(Tuned_Cue_STD1_STD2_MIS,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,Tuned_Cue_STD1_STD2_MIS(i)), 'b');
title(['CUE ' 'Tuned ' ,name_STD1 ,'/' ,name_STD2, '&', name_MIS, ' ROI n' num2str(Tuned_Cue_STD1_STD2_MIS(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,Tuned_Cue_STD1_STD2_MIS(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,Tuned_Cue_STD1_STD2_MIS(i)), 'g');
maxrate(i)=max(max([rate_map_smooth{STDsess(1)}(:,Tuned_Cue_STD1_STD2_MIS(i))...
rate_map_smooth{MISsess(1)}(:,Tuned_Cue_STD1_STD2_MIS(i))...
rate_map_smooth{STDsess(2)}(:,Tuned_Cue_STD1_STD2_MIS(i))]));
  
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 maxrate(i)])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 maxrate(i)])

 end
figure;
suptitle(['CUE ' 'Tuned ' ,name_STD1 ,'/' ,name_STD2, '&', name_MIS, ' ROI n' num2str(Tuned_Cue_STD1_STD2_MIS(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{Tuned_Cue_STD1_STD2_MIS(i)}, 'b');
xlabel(name_STD1)
hold on;
 subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{Tuned_Cue_STD1_STD2_MIS(i)}, 'r');
xlabel(name_MIS)
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{Tuned_Cue_STD1_STD2_MIS(i)}, 'g');
xlabel(name_STD2)
end

%Plot rate map and tuning vector of NO CUE tuned ROI (STD1 and STD2);
for i=1:size(Tuned_NoCue_STD1_STD2_MIS,2)
figure; plot(rate_map_smooth{STDsess(1)}(:,Tuned_NoCue_STD1_STD2_MIS(i)), 'b');
title(['NO CUE ' 'Tuned ' ,name_STD1 ,'/' ,name_STD2, '&', name_MIS, ' ROI n' num2str(Tuned_NoCue_STD1_STD2_MIS(i))]);
 hold on; plot(rate_map_smooth{MISsess(1)}(:,Tuned_NoCue_STD1_STD2_MIS(i)), 'r');
 plot(rate_map_smooth{STDsess(2)}(:,Tuned_NoCue_STD1_STD2_MIS(i)), 'g');
 maxrate(i)=max(max([rate_map_smooth{STDsess(1)}(:,Tuned_NoCue_STD1_STD2_MIS(i))...
     rate_map_smooth{MISsess(1)}(:,Tuned_NoCue_STD1_STD2_MIS(i))...
     rate_map_smooth{STDsess(2)}(:,Tuned_NoCue_STD1_STD2_MIS(i))]));
 
 for ii=1:size(MIS_bins_toswap,2)
line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0  maxrate(i) ])
line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0  maxrate(i)])
 end
figure;
suptitle(['NO CUE ' 'Tuned ' ,name_STD1 ,'/' ,name_STD2, '&', name_MIS, ' ROI n' num2str(Tuned_NoCue_STD1_STD2_MIS(i))]);
subplot(1,3,1) 
compass(tuning_vector{STDsess(1)}{Tuned_NoCue_STD1_STD2_MIS(i)}, 'b');
xlabel(name_STD1)
hold on;
 subplot(1,3,2)
compass(tuning_vector{MISsess(1)}{Tuned_NoCue_STD1_STD2_MIS(i)}, 'r');
xlabel(name_MIS)
 subplot(1,3,3)
compass(tuning_vector{STDsess(2)}{Tuned_NoCue_STD1_STD2_MIS(i)}, 'g');
xlabel(name_STD2)
end

end

%Plot rate map and tuning vector tuned ROI in MIS only;
%for i=1:size(MIS_only,2)
%figure; plot(rate_map_smooth{STDsess(1)}(:,MIS_only(i)), 'b');
%title(['CUE' 'Tuned ' ,name_STD1 ,'/' ,name_STD2, '&', name_MIS, ' ROI n' num2str(MIS_only(i))]);
% hold on; plot(rate_map_smooth{MISsess(1)}(:,MIS_only(i)), 'r');
% plot(rate_map_smooth{STDsess(2)}(:,MIS_only(i)), 'g');
% for ii=1:size(MIS_bins_toswap,2)
%line([MIS_bins_toswap(1,ii) MIS_bins_toswap(1,ii)], [0 0.5])
%line([MIS_bins_toswap(end,ii) MIS_bins_toswap(end,ii)], [0 0.5])
% end
%figure;
%suptitle(['Tuned ' ,name_MIS, ' only  ' 'ROI n' num2str(STD1_STD2(i))]);
%subplot(1,3,1) 
%compass(tuning_vector{STDsess(1)}{MIS_only(i)}, 'b');
%xlabel(name_STD1)
%hold on;
% subplot(1,3,2)
%compass(tuning_vector{MISsess(1)}{MIS_only(i)}, 'r');
%xlabel(name_MIS)
% subplot(1,3,3)
%compass(tuning_vector{STDsess(2)}{MIS_only(i)}, 'g');
%xlabel(name_STD2)
%end






%Events in swap bins = events in 



%e




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



