
function [Events]=detect_events(Imaging,options);

%% Import 

if options.msbackadj== true && options.restricted==true
C_df=Imaging.trace_restricted_baselinesub;  
Cdf_time=Imaging.time_restricted;

elseif options.msbackadj== true && options.restricted==false
C_df=Imaging.trace_baselinesub;
Cdf_time=Imaging.time;

elseif options.msbackadj== false && options.restricted==true
C_df=Imaging.trace_restricted;
Cdf_time=Imaging.time_restricted;

elseif options.msbackadj== false && options.restricted==false
C_df=Imaging.trace; 
Cdf_time=Imaging.time;
end


nb_it=options.iterations;
SDON=options.SDON; %Threshold above x SD for ONSET 
SDOFF=options.SDOFF; %Threshold below x SD for OFFSET 
mindurevent=options.mindurevent;


%% Detect events
on_ones=C_df>=SDON*std(C_df);
off_ones=C_df<=SDOFF*std(C_df);

 %Find onset and offset, do multiple iterations
for it=nb_it % nb of iteration 
[on_off_ones, on_off_binary, on_off] = on_off_thr(on_ones, off_ones, C_df);
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_idx{i}=find(on_off_ones(:,i)==0);
Cdf_off{i}=C_df(Cdf_off_idx{i},i);
std_Cdf_off(i)=std(Cdf_off{i});
end
on_ones=C_df>=SDON*std_Cdf_off;
off_ones=C_df<=SDOFF*std_Cdf_off;
end

 
%% Exclude events 
%Find event duration : time end - time start
  
for i=1:size(on_off,2)
    if isempty(on_off{i})==0
eventduration{i}=Cdf_time(on_off{i}(:,2))-Cdf_time(on_off{i}(:,1));
end 
end

% exclude (nan) if events smaller than min duration set in parameter
for i=1:size(eventduration,2)
    for ii=1:size(eventduration{i},1)
if eventduration{i}(ii)<mindurevent,
 on_off{i}(ii)=NaN;
end  
end
end
for i=1:size(on_off,2)
 onoff=on_off{i};
onoff= onoff(0== sum(isnan(onoff), 2), :);
  onset_offset{i}=onoff;
end

% Find event time (onset offset)
for i=1:size(onset_offset,2)
    if isempty(onset_offset{i})==0
event_time{i}=[Cdf_time(onset_offset{i}(:,1)) Cdf_time(onset_offset{i}(:,2))];
end 
end

%Make binary:
binary=zeros(size(C_df,1),size(C_df,2));
for i=1:size(onset_offset,2)
     if isempty(onset_offset{i})==0,  
binary(onset_offset{i}(:,1),i)=1;
     end
end
onset_binary=binary;

%Make ones
ones=zeros(size(C_df,1),size(C_df,2));
for i=1:size(onset_offset,2)
       if isempty(onset_offset{i})==0, 
for ii=1:size(onset_offset{i},1)
ones(onset_offset{i}(ii,1):onset_offset{i}(ii,2),i)=1;
end
end
end
onset_offset_ones=ones;

%% Save into structure
Events.onset_offset=onset_offset;
Events.onset_binary=onset_binary;
Events.onset_ones=onset_offset_ones;
Events.onset_offset_time=event_time;
Events.options=options;
Events.STD_noise=std_Cdf_off;

%% Figure.
if options.dispfig==1
c2plot=options.c2plot;
figure; hold on;
plot(C_df(:,c2plot))
plot(onset_offset_ones(:,c2plot));
refline([0 SDON*std_Cdf_off(:,c2plot)])
refline([0 SDOFF*std_Cdf_off(:,c2plot)])
end

