
function [Events]=eventdetect_V6(Events, Imaging);
if (Events.options.baselinesub==1)&& (Events.options.restricted==1),
C_df=Imaging.trace_restricted_baselinesub;  
else if (Events.options.baselinesub==0) && (Events.options.restricted==1),
C_df=Imaging.trace_restricted;
else if (Events.options.baselinesub==1) && (Events.options.restricted==0),
C_df=Imaging.trace_baselinesub;
else if (Events.options.baselinesub==0 && Events.options.restricted==0),
C_df=Imaging.trace;
end;end;end;end



SDON=Events.options.SDON; %Threshold above x SD for ONSET 
SDOFF=Events.options.SDOFF; %Threshold below x SD for OFFSET 
fr_im=Imaging.options.frimaging;

mindur=Events.options.mindur*fr_im;
mergdur=Events.options.mergdur*fr_im;
mindurevent=Events.options.mindurevent*fr_im;

% Added run_events
%C_dF = file with dF/F over time row=frames / column=ROI
%pos_vec: from rotary_alignMAD
%time_location: from rotary_alignMAD
%Need onoff.m script
%No run:
%[onset_offset, onset_binary]=event_detection_3iterations(C_df);
%For running events:
%[onset_offset, onset_binary, run_onset]=event_detection(C_df, time_location);
%time_location= behavior with column one time / column 2 location


%C_df=C_df';

%Correct baseline with msbackadj (need bioinformatics):


% %Set the number of SD for onset: 
% SDON=3;
% %Set the number of SD for offset: 
% SDOFF=1;

%First iteration sd on full trace
%Threshold events onset when >= 2*sd (2*dev(i)) /i=ROI analyse
%add a minumum time to consider it an significant event
stamps_1=1:size(C_df,1); 
for i=1:size(C_df,2);
    subdata_1=C_df(:,i);
    std_dev_1(i)=std(subdata_1);
    [events_time_1{i},time_1]=Threshold([stamps_1' subdata_1],'>=',SDON*std_dev_1(i));
    events_binary_1(:,i)=diff(time_1)==1; 
    timing_1{i}=time_1;
%Threshold events offset when <= 0.5*sd (2*dev(i)) /i=ROI analyse
[offset_time_1{i},stop_1]=Threshold([stamps_1' subdata_1],'<=',SDOFF*std_dev_1(i));
  offset_1{i}=stop_1; 
end

[onset_offset_1] = onoff_V4(events_time_1, C_df, offset_1);
%Make binary for onset /i=ROI analyse
for i=1:size(C_df,2);
stamp_df_1=[(1:size(C_df,1))' C_df(:,i)] ;
[status_1] = InIntervals(stamp_df_1,onset_offset_1{i});
onset_binary_1(:,i)=diff(status_1)==1;
onset_offset_ones_1(:,i)=status_1;
%Split dF when events and no events detected  
    events_frames_1{i}=[stamps_1' onset_offset_ones_1(:,i) C_df(:,i)]; 
    events_df_1{i}=events_frames_1{i}(events_frames_1{i}(:,2)==1, :);
    noevents_df_1{i}=events_frames_1{i}(events_frames_1{i}(:,2)==0, :);
 
end




% Second iteration sd on the non events dF
%Threshold when >= 2*sd (2*dev(i)) /i=ROI analyse
for i=1:size(C_df,2);
subdata_2=noevents_df_1{i}(:,3);
subdata_1=C_df(:,i);
stamps_2=noevents_df_1{i}(:,1);
std_dev_noevents_1(i)=std(subdata_2);
[events_time_2{i},time_2]=Threshold([stamps_1' subdata_1],'>=',SDON*std_dev_noevents_1(i));
events_binary_2(:,i)=diff(time_2)==1;
timing_2{i}=time_2;

%Threshold events offset when <= 0.5*sd (2*dev(i)) /i=ROI analyse
[offset_time_2{i},stop_2]=Threshold([stamps_1' subdata_1],'<=',SDOFF*std_dev_noevents_1(i));
offset_2{i}=stop_2; 
end 
[onset_offset_2] = onoff_V4(events_time_2, C_df, offset_2);
for i=1:size(C_df,2);
%Make binary for onset /i=ROI analyse
stamp_df_2=[(1:size(C_df,1))' C_df(:,i)] ;
[status_2] = InIntervals(stamp_df_2,onset_offset_2{i});
onset_binary_2(:,i)=diff(status_2)==1;
onset_offset_ones_2(:,i)=status_2;
%Split dF when events and no events detected  
events_frames_2{i}=[stamps_1' onset_offset_ones_2(:,i) C_df(:,i)]; 
events_df_2{i}=events_frames_2{i}(events_frames_2{i}(:,2)==1, :);
noevents_df_2{i}=events_frames_2{i}(events_frames_2{i}(:,2)==0, :);
end

% Third iteration sd on the non events dF
%Threshold when >= 2*sd (2*dev(i)) /i=ROI analyse
for i=1:size(C_df,2);
subdata_3=noevents_df_2{i}(:,3);
subdata_1=C_df(:,i);
stamps_3=noevents_df_2{i}(:,1);
std_dev_noevents_2(i)=std(subdata_3);
[events_time_3{i},time_3]=Threshold([stamps_1' subdata_1],'>=',SDON*std_dev_noevents_2(i));
events_binary_3(:,i)=diff(time_3)==1;
timing_3{i}=time_3;
%Threshold events offset when <= 0.5*sd (2*dev(i)) /i=ROI analyse

[offset_time_3{i},stop_3]=Threshold([stamps_1' subdata_1],'<=',SDOFF*std_dev_noevents_2(i));
offset_3{i}=stop_3; 
end 
[onset_offset_3] = onoff_V4(events_time_3, C_df, offset_3);
for i=1:size(C_df,2);
%Make binary for onset /i=ROI analyse
stamp_df_3=[(1:size(C_df,1))' C_df(:,i)] ;

[status_3] = InIntervals(stamp_df_3,onset_offset_3{i});
onset_binary_3(:,i)=diff(status_3)==1;
onset_offset_ones_3(:,i)=status_3;
%Split dF when events and no events detected  
% events_frames_3{i}=[stamps_1' onset_offset_ones_3(:,i) C_df(:,i)]; 
% events_df_3{i}=events_frames_3{i}(events_frames_3{i}(:,2)==1, :);
% noevents_df_3{i}=events_frames_3{i}(events_frames_3{i}(:,2)==0, :);
end

onset_offset=onset_offset_3;
onset_binary=[onset_binary_3; zeros(1, size(onset_binary_3,2))];
onset_offset_ones=onset_offset_ones_3;

  

%Min duration of event 
%mindurevent
for i=1:size(onset_offset,2)
eventduration{i}=onset_offset{i}(:,2)-onset_offset{i}(:,1);
end 

for i=1:size(eventduration,2)
    for ii=1:size(eventduration{i},1)
if eventduration{i}(ii)<mindurevent,
 onset_offset{i}(ii)=NaN;
end  
end
end
for i=1:size(onset_offset,2)
 onoff=onset_offset{i};
onoff= onoff(0== sum(isnan(onoff), 2), :);
  onset_offset{i}=onoff;
end


%Make binary:
binary=zeros(size(C_df,1),size(C_df,2));
for i=1:size(onset_offset,2)
binary(onset_offset{i}(:,1),i)=1;
end
onset_binary=binary;
%Make ones
ones=zeros(size(C_df,1),size(C_df,2));
for i=1:size(onset_offset,2)
for ii=1:size(onset_offset{i},1)
ones(onset_offset{i}(ii,1):onset_offset{i}(ii,2),i)=1;
end
end
onset_offset_ones=ones;




%Plot figures /i=ROI to plot
c2plot=Imaging.options.celltoplot;
figure;hold on;
plot(C_df(:,c2plot));
PlotHVLines(SDON*std_dev_noevents_2(1,c2plot),'h','r'); 
PlotHVLines(SDOFF*std_dev_noevents_2(1,c2plot),'h', 'b');
%plot(onset_binary_1(:,c2plot));
%plot(events_binary_2(:,c2plot));
plot(onset_offset_ones(:,c2plot));
%plot(onset_offset_ones(:,c2plot));
   

Events.onset_offset=onset_offset;
Events.onset_binary=onset_binary;
Events.onset_ones=onset_offset_ones;
end





