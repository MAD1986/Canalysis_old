
function [Place_cell]=shuffle_place_cell(Place_cell,Behavior,Events,options);
onset_map=Place_cell.Spatial_Info.event_map;
overall_rate=Place_cell.Spatial_Info.overall_rate;
proba_bin=Place_cell.Spatial_Info.proba_bin; 
run_onset_binary=Place_cell.Spatial_Info.Run_onset_bin;
occupancy=Place_cell.Spatial_Info.occupancy_map;
bin=Place_cell.Bin;
Nshuffle=options.Nshuffle;
Nbin=options.Nbin;
sigma=options.sigma_filter;
%% Shuffle onset map and tuning specificity
tic;
onset_binary=Events.Run.run_onset_binary;
disp(['Starting shuffle ', num2str(Nshuffle), 'X'])
parfor S=1:Nshuffle
    %[onset_map_shuffle{S}]=shuffle_map(run_onset_binary,bin,Nshuffle,Nbin);
[onset_map_shuffle{S},shuffle_tuning_specificity{S}]=shuffle_script(onset_binary,Place_cell,Behavior,Events,options);
end
disp(['End shuffle '])
toc;
%% Rate maps : total number of onset that occurred in a location
parfor S=1:Nshuffle
[onset_map_shuffle_sm{S}]=gauss_filt(onset_map_shuffle{S},Nbin,sigma)
end
for S=1:Nshuffle
for i=1:length(options.Nbin)
rate_map_shuffle_sm{S}{i}=onset_map_shuffle_sm{S}{i}./occupancy{i}';
end
end

%% Spatial specificity null distribution (Danielson et al. 2016);
for S=1:Nshuffle
for nbin=1:length(options.Nbin)
    %(Danielson et al. 2016)
rate_ratio_shuffle{S}{nbin}=(rate_map_shuffle_sm{S}{nbin}./mean(rate_map_shuffle_sm{S}{nbin}));
spatial_info_bin_shuffle{S}{nbin}=rate_map_shuffle_sm{S}{nbin}.*log(rate_ratio_shuffle{S}{nbin}).*proba_bin{nbin}';
SI_shuffle{S}(nbin,:)=nansum(spatial_info_bin_shuffle{S}{nbin});
%(Skaggs)
spatial_info_bin_S_shuffle{S}{nbin}=rate_ratio_shuffle{S}{nbin}.*log2(rate_ratio_shuffle{S}{nbin}).*proba_bin{nbin}';
SIS_shuffle{S}(nbin,:)=nansum(spatial_info_bin_S_shuffle{S}{nbin});
end
end

Place_cell.Spatial_Info.Shuffle.spatial_info=SI_shuffle;
Place_cell.Spatial_Info.Shuffle.spatial_info_Skaggs=SIS_shuffle;

shuffle_TS=cell2mat(shuffle_tuning_specificity');
Place_cell.Tuning_Specificity.Shuffle.tuning_specificity=shuffle_TS;

%% Test significance

%Tuning Specificity

tuning_spe=Place_cell.Tuning_Specificity.tuning_specificity;
%p-value was defined as the fraction of shuffle distribution that 
%exceeded the cell true tuning specificity
TS_pvalue=(sum(shuffle_TS>tuning_spe))/Nshuffle;
%Significant ROI = <pvalue
TS_ROI=TS_pvalue<options.pvalue;
%Save structure
Place_cell.Tuning_Specificity.ROI_pvalue=TS_pvalue;
Place_cell.Tuning_Specificity.significant_ROI=TS_ROI;


%Spatial Specificity

SI=Place_cell.Spatial_Info.Spatial_Info;
SIS=Place_cell.Spatial_Info.Spatial_Info_Skaggs;
%For each bin (N):
%1)Mean of null distribution for each bin
for nbin=1:length(options.Nbin)
for S=1:Nshuffle
    %(Danielson et al. 2016)
SI_shuffle_bin{nbin}(S,:)=SI_shuffle{S}(nbin,:);
SIS_shuffle_bin{nbin}(S,:)=SIS_shuffle{S}(nbin,:);
end
SI_shuffle_mean(nbin,:)=mean(SI_shuffle_bin{nbin});
SIS_shuffle_mean(nbin,:)=mean(SIS_shuffle_bin{nbin});
end
%2)substract mean of null distribution from all estimates
SI_cor=SI-SI_shuffle_mean;
SIS_cor=SIS-SIS_shuffle_mean;
%3)compute one single estimate of the info content for true values and null
% = max value accross bin
final_SI=max(SI_cor,[],1);
final_SIS=max(SIS_cor,[],1);
for S=1:Nshuffle
final_SI_shuffle(S,:)=max(SI_shuffle{S},[],1);
final_SIS_shuffle(S,:)=max(SIS_shuffle{S},[],1);
end
%4)%p-value was defined as the fraction of shuffle distribution that 
%exceeded the cell true tinfo content
SI_pvalue=(sum(final_SI_shuffle>final_SI))/Nshuffle;
SIS_pvalue=(sum(final_SIS_shuffle>final_SIS))/Nshuffle;
%Significant ROI = <pvalue
SI_ROI=SI_pvalue<options.pvalue;
SIS_ROI=SIS_pvalue<options.pvalue;

Place_cell.Spatial_Info.ROI_pvalue=SI_pvalue;
Place_cell.Spatial_Info.ROI_Skaggs_pvalue=SIS_pvalue;
Place_cell.Spatial_Info.significant_ROI=SI_ROI;
Place_cell.Spatial_Info.significant_ROI_Skaggs=SIS_ROI;

%% Find tuned cells (significant spatial info and tuning specificity)
tuned_cells=intersect((find(TS_ROI==1)),(find(SI_ROI==1)));
Place_cell.Tuned_ROI=tuned_cells;

%plot figure
if options.dispfig==1; % Display figure 
TS_only=sum(TS_ROI)-length(tuned_cells);
SI_only=sum(SI_ROI)-length(tuned_cells);

figure;
[H ,S]=venn([TS_only SI_only  length(tuned_cells)],'faceColor',{'b','y'});
lgd={[{'Tuning specificity ', [num2str(TS_only) '/' num2str(sum(TS_ROI))]}], [{'Spatial information ', [num2str(SI_only) '/' num2str(sum(SI_ROI))]}], [{'Both ', [num2str(length(tuned_cells))]}]};
title(['Overlap in the two inclusion criteria p<', num2str(options.pvalue)]);
for i = 1:3
text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [lgd{i} ])
end

%Spatial tuning curve and mean dF for tuned cells
STC=Place_cell.Spatial_tuning_curve;
ST_dF=Place_cell.Spatial_tuning_dF;
ST_curve=STC(:,tuned_cells);
ST_dF_curve=ST_dF(:,tuned_cells);

%Sort spatial tuning curve by max rate
[M,I]=max(ST_curve);
ord_STC=[(1:size(ST_curve,2))' I' ST_curve'];
ROI_ord_STC=sortrows(ord_STC,2);
STC_sorted=(ROI_ord_STC(:,3:end));
%Sort mean dF/F by max rate
ord_STC_dF=[(1:size(ST_dF_curve,2))' I' ST_dF_curve'];
ROI_ord_STC_dF=sortrows(ord_STC_dF,2);
STC_dF_sorted=(ROI_ord_STC_dF(:,3:end));
%Remove Nan
STC_sorted_nonan=STC_sorted(~any(isnan(STC_sorted),2),:);
if options.dispfig==1
figure; 
subplot(1,2,1)
imagesc(STC_sorted_nonan); colormap('Jet');
title('Spatial Tuning Curve')
xlabel('bin position') % x-axis label
ylabel('neuron number') % y-axis label
subplot(1,2,2)
imagesc(STC_dF_sorted); colormap('Jet');
title('Mean dF/F')
xlabel('bin position') % x-axis label
clrbar = colorbar;
clrbar.Label.String=({'Normalized Ca2+ transient rate',' / ','Normalized mean dF/F'});
end


end



end


