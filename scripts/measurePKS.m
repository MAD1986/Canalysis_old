function [df_properties]= measurePKS(eventdf,onset_offset,C_df,options)

%% Find peaks
for u=1:size(eventdf,2); 
for uu=1:size(eventdf{u},2); 
minpro(u)=options.STD_off(u)*options.STD_pro;
mindist=options.mindist;
if isempty(eventdf{u}{uu})==0
if options.prominence==true
[PKS{u}{uu},PKS_LOC{u}{uu},PKS_WTH{u}{uu},PKS_PRO{u}{uu}] = findpeaks(eventdf{u}{uu},'MinPeakDistance',mindist, 'MinPeakProminence', minpro(u), 'WidthReference','halfheight');
NB_PKS{u}{uu}=length(PKS{u}{uu});
end
if options.prominence==false
 [PKS{u}{uu},PKS_LOC{u}{uu},PKS_WTH{u}{uu},PKS_PRO{u}{uu}] = findpeaks(eventdf{u}{uu},'MinPeakDistance',mindist, 'WidthReference','halfheight');
NB_PKS{u}{uu}=length(PKS{u}{uu}); 
end
%figure;
% findpeaks(eventdf{u}{uu},'MinPeakDistance',mindist ,'WidthReference','halfheight','MinPeakProminence', minpro(u), 'Annotate','extents');  
end
end
end

% If multiple peaks take highest prominent peak for max peak value and
%location
for u=1:size(PKS,2); 
for uu=1:size(PKS{u},2); 
for uuu=size(PKS{u}{uu},1)         
if isempty(PKS{u}{uu})==0;
MAX_PKS{u}(uu)=PKS{u}{uu}(find(PKS_PRO{u}{uu}==max(PKS_PRO{u}{uu})));
MAX_PKS_LOC{u}(uu)=PKS_LOC{u}{uu}(find(PKS_PRO{u}{uu}==max(PKS_PRO{u}{uu})));
end
end
end
end

%% Measure:
%event_dur = duration of event form onset to offset
%event_amp = amplitude from onset to max peak
%event_mean = mean df from onset to offset
%event_AUC = area under the curve from onset to offset
%event_width =  halfheight width (if multiple peaks = sum)
%time_on = time onset
%time_PKS = time peak
for u=1:size(MAX_PKS,2)
for uu=1:size(MAX_PKS{u},2)        
if isempty(PKS{u}{uu})==0
df_on_off{u}{uu}=C_df((onset_offset{u}(uu,1)):(onset_offset{u}(uu,2)),u);   
event_dur{u}(uu)=length(df_on_off{u}{uu});
event_amp{u}(uu)=MAX_PKS{u}(uu)-df_on_off{u}{uu}(1);
event_mean{u}(uu)=nanmean(df_on_off{u}{uu});
event_AUC{u}(uu)=trapz(df_on_off{u}{uu});
event_width{u}(uu)=sum(PKS_WTH{u}{uu});
%time_on{u}(uu)=time(onset_offset{u}(uu,1));
%time_PKS{u}(uu)=time(MAX_PKS_LOC{u}(uu)+onset_offset{u}(uu,1));
end
end
end
for u=1:size(onset_offset,2); 
nb_event(u)=size(onset_offset{u},1); 
end
nb_event_tot=sum(nb_event);
% If non analyzed, properties = NaN
for u=1:size(event_dur,2)
for uu=1:size(event_dur{u},2)
if any(event_dur{u}(uu))==0
df_on_off{u}{uu}=nan;   
event_dur{u}(uu)=nan;  
event_amp{u}(uu)=nan; 
event_mean{u}(uu)=nan;  
event_AUC{u}(uu)=nan; 
event_width{u}(uu)=nan;
MAX_PKS{u}(uu)=nan;
    end
end
end

%% Make structure
df_properties.df_on_off=df_on_off;
df_properties.event_dur=event_dur;
df_properties.event_amp=event_amp;
df_properties.event_mean=event_mean;
df_properties.event_AUC=event_AUC;
df_properties.event_width=event_width;
df_properties.MAX_PKS=MAX_PKS;
df_properties.PKS=PKS;
df_properties.nb_event_tot=nb_event_tot;
end