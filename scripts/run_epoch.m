
function [Events, Behavior]=run_epoch(Events,Behavior,Imaging, options);

%% Import
%Options
options.restricted=Events.options.restricted;
mindur=options.mindur;
mergdur=options.merge;
mov_wind=options.moving_window;

%Calcium trace
if Imaging.options.msbackadj== true && options.restricted==true
C_df=Imaging.trace_restricted_baselinesub;  
elseif Imaging.options.msbackadj== true && options.restricted==false
C_df=Imaging.trace_baselinesub;
elseif Imaging.options.msbackadj== false && options.restricted==true
C_df=Imaging.trace_restricted;
elseif Imaging.options.msbackadj== false && options.restricted==false
C_df=Imaging.trace; 
end
%Time and behavior
if options.restricted==true
Cdf_time=Imaging.time_restricted;
time=Behavior.restricted.time;
cum_postion=Behavior.restricted.cumulativeposition;
position=Behavior.restricted.position;
elseif options.restricted==false
Cdf_time=Imaging.time;
time=Behavior.time;
cum_postion=Behavior.cumulativeposition;
position=Behavior.position;
end
%Events
onset_offset=Events.onset_offset;
onset_ones=Events.onset_ones;
onset_binary=Events.onset_binary;
time_position=[time position];

%% Compute speed
%Resample behavior at same frequency than imaging
[N,bin]=histc(Cdf_time,time);
index=bin+1;
if abs(Cdf_time-time(bin))<abs(Cdf_time-time(bin+1));
index=bin;
res_time_position=time_position(bin,:);
res_cum_postion=cum_postion(bin);
res_position=position(bin);
else
res_time_position=time_position(bin,:);
res_cum_position=cum_postion(bin);
res_position=position(bin);
end
res_time=res_time_position(:,1);

%Smooth position 
%res_cum_pos_sm=smooth(res_cum_postion,3);
% Mean measured framerate of imaging (Hz)
avg_fr=1/(mean(diff(Cdf_time)));

%Speed 
speed=[0;diff(res_cum_position)]*avg_fr;
%average using moving mean

speed=movmean(speed,[mov_wind]);



%% Find running epochs

switch options.method
    case 'peak'
disp('minimun peak method')
pks_thr=options.minpeak;
%Find periods of forward motion (speed >0)
run_idx=find(speed>0);
run_time=Cdf_time(run_idx);
run_speed=speed(run_idx);
dist_epochs=diff(run_time);
%Find epochs separated by more than the merging threshold
epochs_end_idx=find(dist_epochs>=mergdur);
epochs_start_idx=[run_idx(1);epochs_end_idx+1];
epochs_end_idx=[epochs_end_idx ;length(run_time)];
run_epochs_idx=[epochs_start_idx epochs_end_idx];
run_epochs_time=run_time(run_epochs_idx);
%Minimum duration for running epoch
run_epochs_dur=run_epochs_time(:,2)-run_epochs_time(:,1);
for i=1:size(run_epochs_dur,1)
if run_epochs_dur(i)<mindur ==1
   run_epochs_time(i,:)=NaN;
   run_epochs_idx(i,:)=NaN;
end  
end
run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);
%Find if peaks speed value in running epochs, count how many:
for i=1:size(run_epochs_idx,1)
speed_run_epochs{i}=run_speed(run_epochs_idx(i,1):run_epochs_idx(i,2));
pks_speed_cnt(i)=sum(speed_run_epochs{i}>pks_thr);
%Remove run epochs if no peak speed
if pks_speed_cnt(i)==0,
   run_epochs_time(i,:)=NaN;
   run_epochs_idx(i,:)=NaN;
end  
end
run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);


case 'speed'
disp('average speed')
run_thr=options.minspeed;
%Find periods when speed > threshold 
run_idx=find(speed>run_thr);
run_time=Cdf_time(run_idx);
run_speed=speed(run_idx);
dist_epochs=diff(run_time);
%Find epochs separated by more than the merging threshold
epochs_end_idx=find(dist_epochs>=mergdur);
epochs_start_idx=[run_idx(1);epochs_end_idx+1];
epochs_end_idx=[epochs_end_idx ;length(run_idx)];
run_epochs_idx=[epochs_start_idx epochs_end_idx];
run_epochs_time=run_time(run_epochs_idx);
%Minimum duration for running epoch
run_epochs_dur=run_epochs_time(:,2)-run_epochs_time(:,1);
for i=1:size(run_epochs_dur,1)
if run_epochs_dur(i)<mindur ==1
   run_epochs_time(i,:)=NaN;
   run_epochs_idx(i,:)=NaN;
end  
end
run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);

    otherwise
        disp('options.method should be peak or speed')
end


%% Find run time and position
%Run idx
run_int_idx=run_idx(run_epochs_idx);
%make binary
run_binary=zeros(size(Cdf_time,1), 1);
run_ones=run_binary;
run_binary(run_int_idx(:,1),1)=1;
%make ones
for i=1:size(run_int_idx,1)
run_ones(run_int_idx(i,1):run_int_idx(i,2),:)=1;
end

runtime=res_time(run_ones==1);
run_position=res_position(run_ones==1);

on2=onset_binary(run_ones==1,:);


%% Running events

% restrict onset to running epochs
for i=1:size(onset_offset,2)
    for ii=1:size(run_int_idx,1)
onset=onset_offset{i}(:,1);
on_R{i}(:,ii)=onset>=run_int_idx(ii,1) & onset<=run_int_idx(ii,2);
onset_R_keep{i}=sum(on_R{i},2);
end
end

for i=1:size(onset_offset,2)
    keep=onset_R_keep{i};
run_onset_offset{i}=onset_offset{i}(keep==1,:);
norun_onset_offset{i}=onset_offset{i}(keep==0,:);
end

%Make binary
run_onset_binary=zeros(size(onset_binary,1), size(onset_binary,2));
norun_onset_binary=run_onset_binary;
for i=1:size(onset_offset,2)
run_onset_binary(run_onset_offset{i}(:,1),i)=1;
norun_onset_binary(norun_onset_offset{i}(:,1),i)=1;
end

%Make ones
run_onset_ones=zeros(size(onset_ones,1), size(onset_ones,2));
norun_onset_ones=run_onset_ones;
for i=1:size(onset_offset,2)
for irun=1:size(run_onset_offset{i},1)
run_onset_ones(run_onset_offset{i}(irun,1):run_onset_offset{i}(irun,2),i)=1;
end
for inorun=1:size(norun_onset_offset{i},1)
norun_onset_ones(norun_onset_offset{i}(inorun,1):norun_onset_offset{i}(inorun,2),i)=1;
end
end
%% Make structure

Behavior.resampled.run_time=runtime;
%Behavior.resampled.no_run_time=noruntime;
Behavior.run_ones=run_ones;
%Behavior.norunbinary=norun;
Behavior.speed=speed;
Behavior.resampled.time=res_time;
Behavior.resampled.cumulativeposition=res_cum_position;
Behavior.resampled.position=res_position;
Behavior.resampled.run_position=run_position;
Behavior.resampled.run_idx=run_int_idx;


Events.Run.run_onset_offset=run_onset_offset;
Events.Run.run_onset_binary=run_onset_binary;
Events.Run.run_onset_ones=run_onset_ones;



Events.NoRun.norun_onset_offset=norun_onset_offset;
Events.NoRun.norun_onset_binary=norun_onset_binary;
Events.NoRun.norun_onset_ones=norun_onset_ones;

Events.options.run_epochs=options;

if options.dispfig==true,
c2plot=options.c2plot;    
    
   figure; 
   subplot(1,2,1);
   hold on;
   plot(res_time, speed);
   plot(Cdf_time, onset_binary(:,c2plot)*5-10);
   plot(Cdf_time, run_onset_binary(:,c2plot)*5-5);
   plot(Cdf_time, C_df(:,c2plot)*10-10);

   
   subplot(1,2,2);
hold on;
   plot(res_time, speed);
   plot(Cdf_time, onset_binary(:,c2plot)*5-10);
   plot(Cdf_time, norun_onset_binary(:,c2plot)*5-5);
   plot(Cdf_time, C_df(:,c2plot)*10-10);
end

end



