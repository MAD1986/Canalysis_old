
function [Event_Properties]=event_properties(C_df, onset_offset, time, options);


mindist=options.mindist; % Set minimun distance between peaks
STD_pro=options.STD_pro; % Set minimun prominence of peak ( X * STD noise trace)



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

for u=1:size(eventdf,2); 
for uu=1:size(eventdf{u},2); 
minpro(u)=options.STD_off(u)*STD_pro;
[PKS{u}{uu},PKS_LOC{u}{uu},PKS_WTH{u}{uu},PKS_PRO{u}{uu}] = findpeaks(eventdf{u}{uu},'MinPeakDistance',mindist, 'MinPeakProminence', minpro(u), 'WidthReference','halfheight');
NB_PKS{u}{uu}=length(PKS{u}{uu});
%figure;
% findpeaks(eventdf{u}{uu},'MinPeakDistance',mindist ,'WidthReference','halfheight','MinPeakProminence', minpro(u), 'Annotate','extents');  
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
for u=1:size(PKS,2); 
for uu=1:size(PKS{u},2);         
if isempty(PKS{u}{uu})==0;
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
for u=1:size(onset_offset,2);  
for uu=1:size(onset_offset{u},1); 
    if isempty(PKS{u}{uu})==1;
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



%% Exclude onset offset / binary / ones 
% of non analyzed envents

    if options.exclude==true,
for u=1:size(onset_offset,2);  
for uu=1:size(onset_offset{u},1); 
    
    if isempty(PKS{u}{uu})==1;

onset_offset{u}(uu,:)=nan;
excluded_df_on_off{u}{uu}=df_on_off_all{u}{uu};

end
end
end

% Remove onset offset values for excluded events
for u=1:size(onset_offset,2); 
onset_offset_nonan{u}= onset_offset{u}(~any(isnan( onset_offset{u}),2),:);
end

%Make binary / ones without excluded events
binary=zeros(size(C_df,1),size(C_df,2));
ones=binary;
for i=1:size(onset_offset_nonan,2)
     if isempty(onset_offset_nonan{i})==0,  
binary(onset_offset_nonan{i}(:,1),i)=1;
for ii=1:size(onset_offset_nonan{i},1)
    ones(onset_offset_nonan{i}(ii,1):onset_offset_nonan{i}(ii,2),i)=1;
end
end
end



for u=1:size(onset_offset,2); 
nb_event(u)=size(onset_offset{u},1); 
nb_event_ana(u)=size(onset_offset_nonan{u},1);
end
nb_event_tot=sum(nb_event);
nb_event_ana_tot=sum(nb_event_ana);

nb_excluded_events=nb_event_tot-nb_event_ana_tot;
disp(['Number of excluded events = ' num2str(nb_excluded_events)])



Event_Properties.nb_events=nb_event_ana;
Event_Properties.nb_excluded_events=nb_excluded_events;
Event_Properties.excluded_trace=excluded_df_on_off;
Event_Properties.onset_offset=onset_offset_nonan;
Event_Properties.onset_binary=binary;
Event_Properties.onset_ones=ones;


end

    
    
    


%% Structure

%remove nan
for u=1:size(onset_offset,2); 
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

Event_Properties.noNaN.duration=event_dur_nonan;
Event_Properties.noNaN.peak=MAX_PKS_nonan;
Event_Properties.noNaN.amplitude=event_amp_nonan;
Event_Properties.noNaN.mean=event_mean_nonan;
Event_Properties.noNaN.AUC=event_AUC_nonan;
Event_Properties.noNaN.width=event_width_nonan;





 if options.exclude==false,
Event_Properties.nb_events=nb_event;
Event_Properties.onset_offset=onset_offset;

binary=zeros(size(C_df,1),size(C_df,2));
ones=binary;
for i=1:size(onset_offset,2)
     if isempty(onset_offset{i})==0,  
binary(onset_offset{i}(:,1),i)=1;
for ii=1:size(onset_offset{i},1)
    ones(onset_offset{i}(ii,1):onset_offset{i}(ii,2),i)=1;
end
end
end


Event_Properties.onset_binary=binary;
Event_Properties.onset_ones=ones;



 end



end


