
function [Spatial_correlation] = rotation_events_V1(Spatial_correlation,Spatialinfo,TunedCells,TuningSpecificity, Events, Behavior);

%% Import values
sessions=Spatial_correlation.options.sessions;
STDsess=Spatial_correlation.options.STDsessions;
MISsess=Spatial_correlation.options.MISsessions;
criteria1=Spatial_correlation.options.tuned_criteria1;
criteria2=Spatial_correlation.options.rotation_criteria;
binbefore=Spatial_correlation.options.rotation_binbefore;
binafter=Spatial_correlation.options.rotation_binafter;


for i = sessions
%try non smoothed
event_map{i}=Spatialinfo{i}.event_map.event_map;
rate_map{i}=Spatialinfo{i}.rate_map.rate_map_smoothed  ;
nbROI(i)=size(rate_map{i},2);
maxROI=max(nbROI);
tunedROI{i}=TunedCells{i}.TunedROI;
infoROI{i}=TunedCells{i}.Spatial_info.significantROI;  
tuningROI{i}=TunedCells{i}.TuningSpecificity.significantROI;  
%bin_texture{i}=Spatialinfo{i}.bin_texture;
occupancy_map_smoothed{i}=Spatialinfo{i}.occupancy_map.occupancy_map_smoothed;  
tuning_vector{i}=TuningSpecificity{i}.tuning_vector;  
onset_bin{i}=Spatialinfo{i}.onset_bin  ;
run_onset{i}=Events{i}.RunningEpochs.run_onset;  
lap{i}=Behavior{i}.restricted.lap  ;

end

for i = sessions
while size(rate_map{i},2)<maxROI
rate_map{i}=[rate_map{i} nan(size(rate_map{i},1),1)]; %replace missing last ROI with NaN
end
end

for i = sessions
while size(event_map{i},2)<maxROI
event_map{i}=[event_map{i} nan(size(event_map{i},1),1)]; %replace missing last ROI with NaN
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

%Name Sessions
name_STD1=['STD',num2str(STDsess(1))];
name_STD2=['STD',num2str(STDsess(2))];
name_MIS=['MIS',num2str(MISsess-STDsess(1))];
name_STD1_STD2=[name_STD1, 'and' name_STD2];
name_STD1_MIS=[name_STD1, 'and' name_MIS];
name_STD2_MIS=[name_STD2, 'and' name_MIS];
name_STD1_MIS_STD2=[name_STD1, 'and' ,name_MIS, 'and', name_STD2];


STD_1=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1)]).(name_STD1);
STD1_only=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1)]).([(name_STD1) '_only']);
%STD2 and STD2 only
STD_2=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2)]).(name_STD2);
STD2_only=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2)]).([(name_STD2) '_only']);
%MIS and MIS only
MIS=Spatial_correlation.rotation.Class.(['Tuned_' (name_MIS)]).(name_MIS);
MIS_only=Spatial_correlation.rotation.Class.(['Tuned_' (name_MIS)]).([(name_MIS) '_only']);
%STD1 and STD2
STD1_STD2=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_STD2)]).(name_STD1_STD2);
STD1_STD2_only=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_STD2)]).([(name_STD1_STD2) '_only']);
%STD1 and MIS 
STD1_MIS=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS)]).(name_STD1_MIS);
STD1_MIS_only=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS)]).([(name_STD1_MIS) '_only']);
%STD2 and MIS 
STD2_MIS=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2_MIS)]).(name_STD2_MIS);
STD2_MIS_only=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD2_MIS)]).([(name_STD2_MIS) '_only']);
%STD1 and STD2 and MIS 
STD1_MIS_STD2=Spatial_correlation.rotation.Class.(['Tuned_' (name_STD1_MIS_STD2)]).(name_STD1_MIS_STD2);



bin_texture=Spatial_correlation.rotation.bin_texture; 
swap_texture=Spatial_correlation.rotation.swap_texture; 

%% Find Ca2+ events in swapped bins
%Need to run new transients_analysis_V9 and spatial_info_V2
%Plot nb of events in bin vs in swapped bin (for MIS) 

%CTRL show for other bins (not swapped)
for i=1:size(bin_texture,2)
All_bins(:,i)=((bin_texture(i))-binbefore:(bin_texture(i))+binafter) ;
end
%If bin to swap is 1 go back to end bins
for j=1:size(bin_texture,2)
for i=1:size(All_bins,1)
if All_bins(i,j)<1, 
All_bins(i,j)=size(event_map{MISsess},1)+All_bins(i,j);
end
end
end
 
for j = sessions
for i = 1:size(onset_bin{j},2)
for ii=1:size(All_bins,1)
for k=1:size(All_bins,2)
for iii=1:size(onset_bin{j}{i},1)
if onset_bin{j}{i}(iii)==All_bins(ii,k),
events_CTRL_bins_zeros{j}{i}{k}(iii)= iii;
nb_events_CTRL_bins_ones{j}{i}(iii,k)=1;
total_CTRL_nb_events{j}(i)=length(onset_bin{j}{i});
end
end
end
end
end
end
%Remove zeros
for j = sessions
for i = 1:size(events_CTRL_bins_zeros{j},2)
for k = 1:size(events_CTRL_bins_zeros{j}{i},2)
eCTRLbins=events_CTRL_bins_zeros{j}{i}{k};
events_CTRL_bins{j}{i}{k} = eCTRLbins(any(eCTRLbins,1));
nb_events_CTRL_bins{j}{i}=sum(nb_events_CTRL_bins_ones{j}{i});
ratio_events_CTRL_bins{j}{i}=nb_events_CTRL_bins{j}{i}/total_CTRL_nb_events{j}(i);

while size(events_CTRL_bins{j}{i},2)<size(bin_texture,2)
    events_CTRL_bins{j}{i}=[events_CTRL_bins{j}{i} 0];
end
end
end
end

for j = sessions
 for i = 1:size(ratio_events_CTRL_bins{j},2)
for k = 1:size(ratio_events_CTRL_bins{j}{i},2)
ratio_2_CTRL{j}{k}(i)=ratio_events_CTRL_bins{j}{i}(k);
end
end
end

for j = sessions
for k = 1:size(ratio_2_CTRL{j},2)
for i = 1:size(ratio_2_CTRL{j}{k},2)
while size(ratio_2_CTRL{j}{k},2)<nbROI(j)
ratio_2_CTRL{j}{k}=[ratio_2_CTRL{j}{k} 0]; %replace missing last ROI with 0
end
ratio_3_CTRL{k}(i,j)=ratio_2_CTRL{j}{k}(i);
end
end
end   
%Show tuned cell in STD1 ratio of events in swapped bins accross sessions
figure;
for i=1:size(bin_texture,2)
subplot(1,(size(bin_texture,2)),i) ;
plot(ratio_3_CTRL{i}(STD_1,sessions)') 
legend(['Tuned Cells ',name_STD1, ' inclusion criteria: ' ,criteria1 ])
ylabel('Events ratio (nb events/total nb events)')
if swap_texture(i)==1
    title('Swapped Texture')
elseif swap_texture(i)==0
    title('Non Swapped Texture')

end
end


%Find tuned cell with ratio higher than 0.25
for j = sessions
for k = 1:size(ratio_3_CTRL,2)
for i = 1:size(ratio_3_CTRL{k},1)
high_ratio{k}{j}=find(ratio_3_CTRL{k}(:,j)>=0.25);
end
end
end

%%Do rate for each texture and plot same graph
for i=sessions
for ii=1:size(rate_map{i},2)
for iii=1:size(All_bins,2)
rate_bins{iii}(ii,i)=mean(rate_map{i}(All_bins(:,iii),ii)); 
end
end
end

figure;
%Show tuned cell in STD1 ratio of events in swapped bins accross sessions
for i=1:size(bin_texture,2)
subplot(1,(size(bin_texture,2)),i) ;
plot(rate_bins{i}(STD_1,sessions)') 
ylabel('Events rate')
legend(['Tuned Cells ',name_STD1, ' inclusion criteria: ' ,criteria1 ])
if swap_texture(i)==1
    title('Swapped Texture')
elseif swap_texture(i)==0
    title('Non Swapped Texture')
end
end


%Plot nb of events in texture by lap
%Get time of events for each texture
for j = sessions
for i = 1:size(events_CTRL_bins{j},2)
for k = 1:size(events_CTRL_bins{j}{i},2)
if events_CTRL_bins{j}{i}{k}>0,
events_bins_time{j}{i}{k}=run_onset{j}{i}(events_CTRL_bins{j}{i}{k});
end
end
end
end
for j = sessions
for i = 1:size(events_CTRL_bins{j},2)
while size(events_bins_time{j}{i},2)<size(events_CTRL_bins{j}{i},2)
events_bins_time{j}{i}=[events_bins_time{j}{i} nan];
end
end
end

%Get lap number for each events
%j= sessions; i= ROI; k= textures; l=lap
for j = sessions
for i = 1:size(events_bins_time{j},2)
for k = 1:size(events_bins_time{j}{i},2)  
for l= 1:size(lap{j},2)
ev_tex_lap{j}{i}{k}{l} =events_bins_time{j}{i}{k}>lap{j}{l}(1) &   events_bins_time{j}{i}{k}<lap{j}{l}(2);
nb_events_tex_lap{j}{i}{k}(l)=sum(ev_tex_lap{j}{i}{k}{l});
end
end
end
end


for i=1:size(bin_texture,2)  
for iii=sessions
for ii=1:size(nb_events_tex_lap{iii},2)
for l= 1:size(lap{iii},2)
if isempty(nb_events_tex_lap{iii}{ii})==0
mean_ROI_1{iii}{i}(ii,l)=nb_events_tex_lap{iii}{ii}{i}(l);
end
end
end 
end
end

%Plot mean number of events per texture for STD1 or all
% iii= sessions / i= texture nb / column = lap
for i=1:size(bin_texture,2)  
for iii=sessions
mean_STD1{iii}(i,:)=mean(mean_ROI_1{iii}{i}(STD_1,:));
mean_all{iii}(i,:)=mean(mean_ROI_1{iii}{i}(:,:));
end
end

%All cells
figure;
for i=1:size(bin_texture,2)  
subplot(1,(size(bin_texture,2)),i) ; 
for ii=sessions
    hold on
plot(mean_all{ii}(i,:))
end
legend(name_STD1, name_MIS, name_STD2) 
ylabel('all ROI nb events per lap')
xlabel('lap number')
if swap_texture(i)==1
    title('Swapped Texture')
elseif swap_texture(i)==0
    title('Non Swapped Texture')
end
end 

%All cells
figure;
for i=1:size(bin_texture,2)  
subplot(1,(size(bin_texture,2)),i) ; 
for ii=sessions
    hold on
plot(mean_STD1{ii}(i,:))
end
legend(name_STD1, name_MIS, name_STD2) 
ylabel('Tuned ROI nb events per lap')
xlabel('lap number')
if swap_texture(i)==1
    title('Swapped Texture')
elseif swap_texture(i)==0
    title('Non Swapped Texture')
end
end 






%Plot different properties of tuned ROI Ca2+ events by lap
%Or just Ca2+ events in swapped bins
%(Y axis = properties, X axis=%laps)
%See https://itb.biologie.hu-berlin.de/~kempter/HippoJC/Articles/lee04.pdf


%% NOT USED
%Rate texture ratio (rate in particular texture / mean rate for all texture)
%for i=1:size(rate_bins,2)
%for ii=1:size(rate_bins{i},1)
%for iii=1:size(All_bins,2)    
%mean_rate_bins1{ii}(i,iii)=rate_bins{i}(ii,iii); 
%mean_rate_bins{ii}=nanmean(mean_rate_bins1{ii});
%rate_tex_ratio1{ii}=mean_rate_bins1{ii}./mean_rate_bins{ii};
%rate_tex_ratio{i}(ii,iii)=rate_tex_ratio1{ii}(i,iii);
%end
%end
%end

%figure;
%Show tuned cell in STD1 ratio of events in swapped bins accross sessions
%for i=1:size(bin_texture,2)
%subplot(1,(size(bin_texture,2)),i) ;
%plot(rate_tex_ratio{i}(STD_1,sessions)') 
%ylabel('Events rate ratio')
%legend(['Tuned Cells ',name_STD1, ' inclusion criteria: ' ,criteria1 ])
%if swap_texture(i)==1
%    title('Swapped Texture')
%elseif swap_texture(i)==0
%    title('Non Swapped Texture')
%end
%end


end

