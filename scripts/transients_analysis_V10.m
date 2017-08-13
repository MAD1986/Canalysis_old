
%% Analysis

function [EventsProperties, NetworkProperties, PCAProperties]=transients_analysis_V10(onset_offset,time,binary, Events, Imaging);

%Need C_df; onset_offset (OR run_onset_offset; norun_onset_offset)

frimage=Imaging.options.frimaging;
maxrisetime=Events.options.analysis.maxrisetime;
maxduration=Events.options.analysis.maxduration;
maxpeak=Events.options.analysis.maxpeak; 
events_binary=Events.onset_ones;


if (Events.options.baselinesub==1)&& (Events.options.restricted==1),
C_df=Imaging.trace_restricted_baselinesub;  
Cdf_time=Imaging.time_restricted;
else if (Events.options.baselinesub==0) && (Events.options.restricted==1),
C_df=Imaging.trace_restricted;
Cdf_time=Imaging.time_restricted;
else if (Events.options.baselinesub==1) && (Events.options.restricted==0),
C_df=Imaging.trace_baselinesub;
Cdf_time=Imaging.time;
else if (Events.options.baselinesub==0 && Events.options.restricted==0),
C_df=Imaging.trace;
Cdf_time=Imaging.time;
end;end;end;end


%% Measure Ca2+ events properties 
%eventdur=  duration from onset to offset
%eventpeak= peak dF/F value 
%eventampl=  amplitude from threshold (2SD) to peak 
%eventmeandf= mean dF/F from onset to offset
%eventarea = area under event form onset to offset
%risetime = rise time to reach half maximum amplitude from onset
%decaytime_int= decay time to reach half maximum amplitude from peak
%decaytime= slope of exponential fiting of the decay
%rsa= regression fit
%halfwidth= width at half amplitude
for u=1:size(onset_offset,2); for uu=1:size(onset_offset{u},1);    
eventdf{u}{uu}=C_df((onset_offset{u}(uu,1)):(onset_offset{u}(uu,2)),u);
eventdur{u}(uu)=length(eventdf{u}{uu});
nbr_events(u)=(size(onset_offset{u},1));
nbr_events_total=sum(nbr_events);
 
if eventdur{u}(uu)>maxduration,
    eventdur{u}(uu)=NaN;
    eventdf{u}{uu}=NaN;
end;end;end


%Try find peaks with min prominence to set
promin=0.075;
dismin=10;
for u=1:size(eventdf,2); for uu=1:size(eventdf{u},2); 
if isnan(eventdf{u}{uu})==0,
[pks{u}{uu},locs{u}{uu},width{u}{uu},pro{u}{uu}] = findpeaks(eventdf{u}{uu}, 'MinPeakProminence', promin, 'WidthReference','halfheight', 'MinPeakDistance',dismin);
nb_pks{u}{uu}=length(pks{u}{uu});
end
end
end

%If multiple peaks take highest prominence or /peak value
for u=1:size(pks,2); for uu=1:size(pks{u},2); 
        for uuu=size(pks{u}{uu},1)         
   if isempty(pks{u}{uu})==0;
    eventpeak{u}(uu)=pks{u}{uu}(find(pro{u}{uu}==max(pro{u}{uu})));
    end
end
    end
end
%fig to test
for u=1:80; for uu=1:size(eventdf{u},2); 
if nb_pks{u}{uu}>=1;
    if pro{u}{uu}>=0.05 & pro{u}{uu}<=0.1
        figure;
  findpeaks(eventdf{u}{uu},'MinPeakProminence', promin, 'WidthReference','halfheight', 'MinPeakDistance',dismin,'Annotate','extents')  
   
    end
    end
    end
    end



for u=1:size(eventdf,2); for uu=1:size(eventdf{u},2); 
 if isnan(eventdf{u}{uu})==0,
eventpeak{u}(uu)=max(eventdf{u}{uu});
%amplitude = from threshold to peak (2SD)
eventampl{u}(uu)=eventpeak{u}(uu)-eventdf{u}{uu}(1);
eventmeandf{u}(uu)=nanmean(eventdf{u}{uu});
eventarea{u}(uu)=trapz(eventdf{u}{uu});
peakidx{u}(uu)=find(eventdf{u}{uu}==eventpeak{u}(uu));
risephase{u}{uu}=eventdf{u}{uu}(1:peakidx{u}(uu));
risephase_sub{u}{uu}=risephase{u}{uu}-risephase{u}{uu}(1);
%risetimeidx{u}{uu}=find(risephase_sub{u}{uu}<=(risephase_sub{u}{uu}(end)/2));
risetime{u}(uu)=length(risephase{u}{uu})/2; 
 elseif isnan(eventdf{u}{uu}),
eventpeak{u}(uu)=nan;
eventampl{u}(uu)=nan;
eventmeandf{u}(uu)=nan;
eventarea{u}(uu)=nan;
peakidx{u}(uu)=nan;
risephase{u}{uu}=nan;
risephase_sub{u}{uu}=nan;
risetimeidx{u}{uu}=nan;
risetime{u}(uu)=nan;
end
end
end


if risetime{u}(uu)>maxrisetime,
    risetime{u}(uu)=NaN;
%For decay: if no decay remove event(end of the recording)
if peakidx{u}(uu)>length(eventdf{u}{uu})-10,
    peakidx{u}(uu)=NaN;
    eventdur{u}(uu)=NaN;
    eventpeak{u}(uu)=NaN;
    eventampl{u}(uu)=NaN;
    eventmeandf{u}(uu)=NaN;
    eventarea{u}(uu)=NaN;
    risetime{u}(uu)=NaN;
%    eventdf{u}{uu}=NaN;
end;end;
end;end;end


for u=1:size(eventdf,2); for uu=1:size(eventdf{u},2);
if risetime{u}(uu)>maxrisetime,
risetime{u}(uu)=NaN;
%For decay: if no decay remove event(end of the recording)
if peakidx{u}(uu)>length(eventdf{u}{uu})-10,
peakidx{u}(uu)=NaN;
eventdur{u}(uu)=NaN;
eventpeak{u}(uu)=NaN;
eventampl{u}(uu)=NaN;
eventmeandf{u}(uu)=NaN;
eventarea{u}(uu)=NaN;
risetime{u}(uu)=NaN;
%eventdf{u}{uu}=NaN;
end;end;end;end


for u=1:size(eventdf,2); for uu=1:size(eventdf{u},2);  
if isnan(eventpeak{u}(uu)),
    decayphase{u}{uu}=NaN;
    halfwidth{u}(uu)=NaN;
    decaytime{u}(uu)=NaN;
    rsq{u}(uu)=NaN;
elseif isnan(eventpeak{u}(uu))==0,
decayphase{u}{uu}=eventdf{u}{uu}(peakidx{u}(uu):end);
decayphase_sub{u}{uu}=decayphase{u}{uu}-risephase{u}{uu}(1);
decaytimeidx{u}{uu}=find(decayphase_sub{u}{uu}<=(eventampl{u}(uu)/2));
if isempty(decaytimeidx{u}{uu})
halfwidth{u}(uu)=NaN;
decaytime{u}(uu)=NaN;
rsq{u}(uu)=NaN;  
eventdur{u}(uu)=NaN;
eventmeandf{u}(uu)=NaN;
eventarea{u}(uu)=NaN;
elseif isempty(decaytimeidx{u}{uu})==0,
decaytime_int{u}(uu)=decaytimeidx{u}{uu}(1);
risehalfwidth{u}(uu)=length(risephase{u}{uu})-risetime{u}(uu);
halfwidth{u}(uu)=risehalfwidth{u}(uu)+decaytime_int{u}(uu);
if length(decayphase{u}{uu})<10,
rsq{u}(uu)=NaN;
decaytime{u}(uu)=NaN;
decaytime_int{u}(uu)=NaN;
elseif length(decayphase{u}{uu})>=10,
[f{u}{uu} g{u}{uu}]=fit((1:length(decayphase{u}{uu}))',decayphase{u}{uu},'exp1'); 
rsq{u}(uu)=g{u}{uu}.rsquare;
decaytime{u}(uu)=-1/f{u}{uu}.b;
end;
end; end; end;end

%When regression failed (rsq < 0.9),decaytime = time to reach half maximal amplitude
for u=1:size(eventdf,2); for uu=1:size(eventdf{u},2);  
if length(decaytime{u})<length(decaytime_int{u}),
    decaytime{u}=[decaytime{u} decaytime_int{u}(end)];
elseif rsq{u}(uu)<0.9,
   decaytime{u}(uu)=decaytime_int{u}(uu);
end;end;end


while length(onset_offset)>length(eventdf),
eventdf=[eventdf {[]}];
eventdur=[eventdur {[]}];
eventpeak=[eventpeak {[]}];
eventampl=[eventampl {[]}];
eventmeandf=[eventmeandf {[]}];
eventarea=[eventarea {[]}];
risetime=[risetime {[]}];
decaytime=[decaytime {[]}];
rsq=[rsq {[]}];
halfwidth=[halfwidth {[]}];
end


for u=1:size(eventdf,2); for uu=1:size(eventdf{u},2);  
EventsProperties.eventsdf{u}{uu}=eventdf{u}{uu};
EventsProperties.duration{u}(uu)=eventdur{u}(uu)*(1/frimage);
EventsProperties.peak{u}(uu)=eventpeak{u}(uu);
EventsProperties.amplitude{u}(uu)=eventampl{u}(uu);
EventsProperties.meandf{u}(uu)=eventmeandf{u}(uu);
EventsProperties.area{u}(uu)=eventarea{u}(uu);
EventsProperties.risetime{u}(uu)=risetime{u}(uu)*(1/frimage);
EventsProperties.decaytime{u}(uu)=decaytime{u}(uu)*(1/frimage);
EventsProperties.regressionfit{u}(uu)=rsq{u}(uu);
EventsProperties.halfwidth{u}(uu)=halfwidth{u}(uu)*(1/frimage);
%Mean excluding zeros / NaN
EventsProperties.mean.duration(u)=nanmean(EventsProperties.duration{u});
EventsProperties.mean.peak(u)=nanmean(EventsProperties.peak{u});
EventsProperties.mean.amplitude(u)=nanmean(EventsProperties.amplitude{u});
EventsProperties.mean.meandf(u)=nanmean(EventsProperties.meandf{u});
EventsProperties.mean.area(u)=nanmean(EventsProperties.area{u});
EventsProperties.mean.risetime(u)=nanmean(EventsProperties.risetime{u});
EventsProperties.mean.decaytime(u)=nanmean(EventsProperties.decaytime{u});
EventsProperties.mean.regressionfit(u)=nanmean(EventsProperties.regressionfit{u});
EventsProperties.mean.halfwidth(u)=nanmean(EventsProperties.halfwidth{u});
    
% total nbr of events total  
EventsProperties.nbrevents.total=nbr_events_total;
% total nbr of analyzed events

analyzed_events(u)=size(eventdur{u},2);
analyzed_events_total=sum(analyzed_events);
EventsProperties.nbrevents.analyzed=analyzed_events_total;
end; end


%% PCA on all properties
%eventdur=  duration from onset to offset
%eventpeak= peak dF/F value 
%eventampl=  amplitude from threshold (2SD) to peak 
%eventmeandf= mean dF/F from onset to offset
%eventarea = area under event form onset to offset
%risetime = rise time to reach half maximum amplitude from onset
%decaytime_int= decay time to reach half maximum amplitude from peak
%decaytime= slope of exponential fiting of the decay
%rsq= regression fit
%halfwidth= width at half amplitude

%PCA for each cell
for u=1:size(peakidx,2); 
 allproperties_cells(u,:)=[EventsProperties.mean.duration(u)...
 EventsProperties.mean.peak(u)...
 EventsProperties.mean.amplitude(u)...
 EventsProperties.mean.meandf(u)...
EventsProperties.mean.area(u)...
 EventsProperties.mean.risetime(u)...       
 EventsProperties.mean.decaytime(u)...
 EventsProperties.mean.halfwidth(u)];     
end; 
allproperties_cells=allproperties_cells(~any(isnan(allproperties_cells),2),:);

% Do the PCA for each cell
[coeff_c,score_c,latent_c,~,explained_c] = pca(allproperties_cells);
%Weighted PCA
%[coeff_c,score_c,latent_c,~,explained_c] = pca(allproperties_cells, 'VariableWeights','variance');

% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_c = cov(allproperties_cells);
[V_c,D_c] = eig(covarianceMatrix_c);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_c = allproperties_cells*coeff_c;
%Plot first 2 PCs
vbls = {'dur','peak','ampl','meandf','area','rise','decay','halfwidth'};
figure; biplot(coeff_c(:,1:2),'scores',score_c(:,1:2),'varlabels',vbls);

%PCA for each event
for u=1:size(peakidx,2); for uu=1:size(peakidx{u},2); 
 allproperties{u}(uu,:)=[EventsProperties.duration{u}(uu)...
 EventsProperties.peak{u}(uu)...
 EventsProperties.amplitude{u}(uu)...
 EventsProperties.meandf{u}(uu)...
 EventsProperties.area{u}(uu)...
 EventsProperties.risetime{u}(uu)...       
 EventsProperties.decaytime{u}(uu)...
 EventsProperties.halfwidth{u}(uu)];     
end; end;
allproperties_mat=cell2mat(allproperties');
allproperties_mat=allproperties_mat(~any(isnan(allproperties_mat),2),:);
EventsProperties.nbrevents.analyzed=size(allproperties_mat,1);


% Do the PCA for each event
[coeff_e,score_e,latent_e,~,explained_e] = pca(allproperties_mat);
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix_e = cov(allproperties_mat);
[V_e,D_e] = eig(covarianceMatrix_e);
% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. 
% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". 
dataInPrincipalComponentSpace_e = allproperties_mat*coeff_e;
%Plot first 2 PCs
vbls = {'dur','peak','ampl','meandf','area','rise','decay','halfwidth'};
figure; biplot(coeff_e(:,1:2),'scores',score_e(:,1:2),'varlabels',vbls);

PCAProperties.Event.dataPCS=dataInPrincipalComponentSpace_e;
PCAProperties.Event.score=score_e;
PCAProperties.Event.coeff=coeff_e;
PCAProperties.Event.latent=latent_e;
PCAProperties.Event.explained=explained_e;
PCAProperties.Event.covarianceMatrix=covarianceMatrix_e;
PCAProperties.Event.varlabels=vbls;
PCAProperties.Event.allvariables=allproperties_mat;


PCAProperties.Cell.dataPCS=dataInPrincipalComponentSpace_c;
PCAProperties.Cell.score=score_c;
PCAProperties.Cell.coeff=coeff_c;
PCAProperties.Cell.latent=latent_c;
PCAProperties.Cell.explained=explained_c;
PCAProperties.Cell.covarianceMatrix=covarianceMatrix_c;
PCAProperties.Cell.varlabels=vbls;
PCAProperties.Cell.allvariables=allproperties_cells;

%% Network Properties


% Activity rate (AUC/min)the area under significant Ca2+ transients divided by recording duration.
for u=1:size(onset_offset,2); 
eventarea_analyzed{u}=eventarea{u};
eventarea_analyzed{u}(isnan(eventarea_analyzed{u})) = [];
AUC(u)=mean(eventarea_analyzed{u});
AUC_permin(u)=AUC(u)/(time/frimage/60);
 NetworkProperties.ActivityRate(u)=AUC_permin(u);
end
%Frequency (event / min)
for u=1:size(onset_offset,2); 
    eventdur_analyzed{u}=eventdur{u};
  eventdur_analyzed{u}(isnan(eventdur_analyzed{u})) = [];
  nbevents_analysed(u)=size(eventdur_analyzed{u},2);
freq(u)= nbevents_analysed(u)/(time/frimage/60);
 NetworkProperties.Frequency(u)=freq(u);
 end; 
    
    
%Mean excluding zeros
 NetworkProperties.MeanActivityRate = nanmean( NetworkProperties.ActivityRate);
 NetworkProperties.MeanFrequency = sum( NetworkProperties.Frequency,2)./sum( NetworkProperties.Frequency~=0,2);

%% Correlation between cell
%First try Pearson Corr full trace
%[dfCCorr, dfCCorrPvalues] = corrcoef(C_df);
% NetworkProperties.PearsonCorr=dfCCorr;

%Do Pearson correlation removing inter events period (dF/F during inter
%event periods = 0)
C_df_sub=C_df;
for i=1:size(C_df_sub,2)
off=find(events_binary(:,i)==0);
C_df_sub(off,i)=0;
end
[dfCCorr_sub, dfCCorrPvalues_sub] = corrcoef(C_df_sub);
NetworkProperties.PearsonCorr_sub=dfCCorr_sub;

%for all, run and no run 
C_df_bin=C_df(find(binary==1),:);
C_df_bin_sub=C_df_sub(find(binary==1),:);

[dfCCorr, dfCCorrPvalues] = corrcoef(C_df_bin);
 NetworkProperties.PearsonCorr=dfCCorr;

 [dfCCorr_sub, dfCCorrPvalues_sub] = corrcoef(C_df_bin_sub);
 NetworkProperties.PearsonCorr_sub=dfCCorr_sub;

%Histogram of cell activity correlations 
%The property of high correlation (HC) was tested for in each neuron by finding 
%the number of correlated neurons with which the Pearson’s correlation coefficient
%was above 0.3.

  for i=1:size(dfCCorr_sub,2);
 corr_pairs{i}=find(dfCCorr_sub(:,i)>0.3);
 nb_corr_pairs(i)=length(corr_pairs{i})-1;
  end
 NetworkProperties.correlated_neurons=corr_pairs;
 NetworkProperties.nb_correlated_pairs=nb_corr_pairs;

 
 %HC neurons were defined as those neurons that had more correlated partners than 
 %that of the average neuron in the same volume by >1 standard deviation


end
