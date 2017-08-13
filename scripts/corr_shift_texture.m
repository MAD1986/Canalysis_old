function [Spatial_correlation] = corr_shift_texture(Spatial_correlation, MIS, plot_sessions);

%% PV correlation shifted bins VS non-shifted bins in MIS sessions


%% Import var
sessions=Spatial_correlation.options.sessions;
PVcorr=Spatial_correlation.PVcorrelation.AllROI.PVcorr;
sessions=Spatial_correlation.options.sessions;
TCcorr=Spatial_correlation.TuningCurvecorr.AllROI.TCcorr;
TCcorr_tuned=Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_tuned;
binbefore=Spatial_correlation.options.rotation_binbefore;
binafter=Spatial_correlation.options.rotation_binafter;
PVcorr_table=Spatial_correlation.PVcorrelation.AllROI.PVcorr_table ; 
TCcorr_tuned_table=Spatial_correlation.TuningCurvecorr.TunedROI.TCcorr_table;
TCcorr_table=Spatial_correlation.TuningCurvecorr.AllROI.TCcorr_table;
STDsess=Spatial_correlation.options.STDsessions;
MISsess=Spatial_correlation.options.MISsessions;

for i = sessions
rate_map{i}=Spatialinfo{i}.rate_map.normalized_rate_map_smoothed;
nbROI(i)=size(rate_map{i},2);
maxROI=max(nbROI);
bin_texture{i}=Spatialinfo{i}.bin_texture;
end
for i = sessions
while size(rate_map{i},2)<maxROI
rate_map{i}=[rate_map{i} nan(size(rate_map{i},1),1)]; %replace missing last ROI with NaN
rate_map_nonnorm{i}=[rate_map_nonnorm{i} nan(size(rate_map_nonnorm{i},1),1)];
end
end
%% Get bin values around shifted texture RFID 
%% Rotation
%get the bin where the texture change


    for i= 1:size(MIS,2)
bin_texture_MIS{i}=bin_texture{{MIS{i}(1)}}; %position of texture on the first session 
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


%% Plot

noswap_bins=setdiff((1:size(PVcorr_table,1)),swap_bins);
PVcorr_table_onlyswap=PVcorr_table(swap_bins,:);
PVcorr_table_noswap=PVcorr_table(noswap_bins,:);

% Plot PV corr

figure;
boxplot([PVcorr_table_onlyswap(:,2) PVcorr_table_noswap(:,2 ) PVcorr_table(:,2)])
title('mean PV correlation')
ylabel('pearson correlation score') % y-axis label


end
