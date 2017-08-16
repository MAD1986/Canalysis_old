
function [Events]=event_detect_new_V2(Events, Imaging);

%% Import 

if Imaging.options.msbackadj== true && options.restricted==true
C_df=Imaging.trace_restricted_baselinesub;  
elseif Imaging.options.msbackadj== true && options.restricted==false
C_df=Imaging.trace_baselinesub;  
elseif Imaging.options.msbackadj== false && options.restricted==true
C_df=Imaging.trace_restricted;  
elseif Imaging.options.msbackadj== false && options.restricted==false
C_df=Imaging.trace;  
end



it=
SDON=options.SDON; %Threshold above x SD for ONSET 
SDOFF=options.SDOFF; %Threshold below x SD for OFFSET 
fr_im=30.3;
mindurevent=options.mindurevent*fr_im;


%% First iteration
%Find values above thr;
on_ones_1=C_df>=SDON*std(C_df);
off_ones_1=C_df<=SDOFF*std(C_df);
%Find onset and offset
[on_off_ones_1, on_off_binary_1, on_off_1] = on_off_thr(on_ones_1, off_ones_1, C_df);

%% Second iteration
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_idx_1{i}=find(on_off_ones_1(:,i)==0);
Cdf_off_1{i}=C_df(Cdf_off_idx_1{i},i);
std_Cdf_off(i)=std(Cdf_off_1{i});
end
on_ones_2=C_df>=SDON*std_Cdf_off;
off_ones_2=C_df<=SDOFF*std_Cdf_off;
%Find onset and offset
[on_off_ones_2, on_off_binary_2, on_off_2] = on_off_thr(on_ones_2, off_ones_2, C_df);

%% Third iteration
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_idx_2{i}=find(on_off_ones_2(:,i)==0);
Cdf_off_2{i}=C_df(Cdf_off_idx_2{i},i);
std_Cdf_off_2(i)=std(Cdf_off_2{i});
end
on_ones=C_df>=SDON*std_Cdf_off_2;
off_ones=C_df<=SDOFF*std_Cdf_off_2;
%Find onset and offset
[on_off_ones, on_off_binary, on_off] = on_off_thr(on_ones, off_ones, C_df);


%% Exclude events 
%Min duration of event 
for i=1:size(on_off,2)
    if isempty(on_off{i})==0
eventduration{i}=on_off{i}(:,2)-on_off{i}(:,1);
end 
end
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

%% Figure.
c2plot=Imaging.options.celltoplot;
figure; hold on;
plot(C_df(:,c2plot))
%plot(on_off_ones_1(:,c2plot));
%plot(on_off_ones(:,c2plot));
plot(onset_offset_ones(:,c2plot));
%refline([0 SDON*std(C_df(:,c2plot))])
%refline([0 SDON*std_Cdf_off(:,c2plot)])
refline([0 SDON*std_Cdf_off_2(:,c2plot)])
%refline([0 SDOFF*std(C_df(:,c2plot))])
%refline([0 SDOFF*std_Cdf_off(:,c2plot)])
refline([0 SDOFF*std_Cdf_off_2(:,c2plot)])

