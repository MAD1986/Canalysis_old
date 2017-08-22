
function [Event_Properties]=event_properties(C_df, onset_offset, time, options);

%% eventdf = calcium trace from onset to offset for each event
% take 20 frames before onset (to capture beginning of event) 
for u=1:size(onset_offset,2); 
for uu=1:size(onset_offset{u},1);   
df_on_off_all{u}{uu}=C_df((onset_offset{u}(uu,1)):(onset_offset{u}(uu,2)),u);
if onset_offset{u}(uu,1)>20
eventdf{u}{uu}=C_df((onset_offset{u}(uu,1))-20:(onset_offset{u}(uu,2)),u);
elseif onset_offset{u}(uu,1)<20
eventdf{u}{uu}=C_df((onset_offset{u}(uu,1)):(onset_offset{u}(uu,2)),u);
end
end
end

%% Analyse events using findpeaks function
% https://www.mathworks.com/help/signal/ref/findpeaks.html
% Find peaks with min prominence and min distance to set

% Measure peak value (PKS), peak location (PKS_LOC), peak width -halfheight- (PKS_WTH) and
% prominence (PKS_PRO)

% Measure:
%event_dur = duration of event form onset to offset
%event_amp = amplitude from onset to max peak
%event_mean = mean df from onset to offset
%event_AUC = area under the curve from onset to offset
%event_width =  halfheight width (if multiple peaks = sum)
options.prominence=1;
[df_properties]= measurePKS(eventdf,onset_offset,C_df,options);

df_on_off=df_properties.df_on_off;
event_dur=df_properties.event_dur;
event_amp=df_properties.event_amp;
event_mean=df_properties.event_mean;
event_AUC=df_properties.event_AUC;
event_width=df_properties.event_width;
MAX_PKS=df_properties.MAX_PKS;
PKS=df_properties.PKS;
nb_event_tot=df_properties.nb_event_tot;


%% Remove onset offset / binary / ones of non analyzed events :
if options.exclude==true

for u=1:size(onset_offset,2)
for uu=1:size(onset_offset{u},1)   
if isempty(PKS{u}{uu})==1
%onset_offset{u}(uu,:)=nan;
excluded_df_on_off{u}{uu}=df_on_off_all{u}{uu};
eventdf_exc{u}{uu}=eventdf{u}{uu};
end
end
end
options.prominence=0;
[df_properties_exc]= measurePKS(eventdf_exc,onset_offset,C_df,options);

df_on_off_exc=df_properties_exc.df_on_off;
event_dur_exc=df_properties_exc.event_dur;
event_amp_exc=df_properties_exc.event_amp;
event_mean_exc=df_properties_exc.event_mean;
event_AUC_exc=df_properties_exc.event_AUC;
event_width_exc=df_properties_exc.event_width;
MAX_PKS_exc=df_properties_exc.MAX_PKS;
PKS_exc=df_properties_exc.PKS;

Event_Properties.excluded_trace=excluded_df_on_off;

end

% Remove onset offset values for excluded events
onset_offset_exc=onset_offset;
for u=1:size(onset_offset,2)
for uu=1:size(onset_offset{u},1)   
if isempty(PKS{u}{uu})==1
onset_offset_exc{u}(uu,:)=nan;
end
end
end
for u=1:size(onset_offset_exc,2)
onset_offset_nonan{u}= onset_offset_exc{u}(~any(isnan( onset_offset_exc{u}),2),:);
end

%Make binary / ones without excluded events
binary_exc=zeros(size(C_df,1),size(C_df,2));
ones_exc=binary_exc;
for i=1:size(onset_offset_nonan,2)
     if isempty(onset_offset_nonan{i})==0
binary_exc(onset_offset_nonan{i}(:,1),i)=1;
for ii=1:size(onset_offset_nonan{i},1)
    ones_exc(onset_offset_nonan{i}(ii,1):onset_offset_nonan{i}(ii,2),i)=1;
end
end
end

for u=1:size(onset_offset,2)
nb_event(u)=size(onset_offset{u},1); 
nb_event_ana(u)=size(onset_offset_nonan{u},1);
end
nb_event_tot=sum(nb_event);
nb_event_ana_tot=sum(nb_event_ana);
nb_excluded_events=nb_event_tot-nb_event_ana_tot;
if options.exclude==true
disp(['Number of excluded events = ' num2str(nb_excluded_events), ' / ' , num2str(nb_event_ana_tot)])
end
if options.exclude==false
disp(['Number of analyzed events = ' num2str(nb_event_ana_tot)])
end

Event_Properties.nb_events=nb_event_ana;
Event_Properties.nb_excluded_events=nb_excluded_events;

Event_Properties.onset_offset_analysed=onset_offset_nonan;
Event_Properties.onset_binary_analysed=binary_exc;
Event_Properties.onset_ones_analysed=ones_exc;


%% Structure
%remove nan

for u=1:size(event_dur,2); 
event_dur_nonan{u}=event_dur{u}(~isnan(event_dur{u}));
MAX_PKS_nonan{u}=MAX_PKS{u}(~isnan(MAX_PKS{u}));
event_amp_nonan{u}=event_amp{u}(~isnan(event_amp{u}));
event_mean_nonan{u}=event_mean{u}(~isnan(event_mean{u}));
event_AUC_nonan{u}=event_AUC{u}(~isnan(event_AUC{u}));
event_width_nonan{u}=event_width{u}(~isnan(event_width{u}));
end 

Event_Properties.trace=df_on_off;
Event_Properties.duration=event_dur;
Event_Properties.peak=MAX_PKS;
Event_Properties.amplitude=event_amp;
Event_Properties.mean=event_mean;
Event_Properties.AUC=event_AUC;
Event_Properties.width=event_width;

if options.exclude==true
Event_Properties.Excluded=df_properties_exc;
end

Event_Properties.noNaN.duration=event_dur_nonan;
Event_Properties.noNaN.peak=MAX_PKS_nonan;
Event_Properties.noNaN.amplitude=event_amp_nonan;
Event_Properties.noNaN.mean=event_mean_nonan;
Event_Properties.noNaN.AUC=event_AUC_nonan;
Event_Properties.noNaN.width=event_width_nonan;
end


