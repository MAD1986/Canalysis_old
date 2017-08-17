
function [Events, Behavior]=run_epoch(Events,Behavior,Imaging);

%% Import

%Options
options.restricted=Events.options.restricted;
run_thr=options.minspeed;
mindur=options.mindur;
mergdur=options.merge;
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
res_cum_postion=cum_postion(bin);
res_position=position(bin);
end
res_time=res_time_position(:,1);

%Smooth position 
%res_cum_pos_sm=smooth(res_cum_postion,3);
% Mean measured framerate of imaging (Hz)
avg_fr=1/(mean(diff(Cdf_time)));
%Speed 
speed=[0;diff(res_cum_postion)]*avg_fr;

%% Find running epochs

%Find periods of forward motion (speed >0)
run_idx=find(speed>0);
run_time=Cdf_time(run_idx);
dist_epochs=diff(run_time);

%Find epochs separated by more than the merging threshold
epochs_end_idx=find(dist_epochs>=mergdur);
epochs_start_idx=[run_idx(1);epochs_end_idx+1];
epochs_end_idx=[epochs_end_idx ;run_idx(end)];
run_epochs_idx=[epochs_start_idx epochs_end_idx];
run_epochs_time=Cdf_time(run_epochs_idx);
%Minimum duration for running epoch
run_epochs_dur=run_epochs_time(:,2)-epochs_start_time(:,1);
if run_epochs_dur(i)<mindur ==1
   run_epochs_time(i,:)=NaN;
   run_epochs_idx(i,:)=NaN;
end  
run_epochs_time=run_epochs_time(~isnan(run_epochs_time(:,2)),:);
run_epochs_idx=run_epochs_idx(~isnan(run_epochs_idx(:,2)),:);

%Find peaks with min speed
%peakspeed_idx=find(speed>=run_thr);

%Find if peaks speed value in running epochs
for i=1:size(run_epochs_idx,1)
speed_run_epochs{i}=speed(run_epochs_idx(i,1):run_epochs_idx(i,2));

end
end
for i=1:size(run_epochs_idx,1)
if speed_run_epochs{i}>run_thr;
test(i)=1
end
end


%determine intervals for running behavior  
%First threshold when animal is moving forward direction = speed >0
[runtime, run]=Threshold([res_time speed],'>',0.5,'min', mindur, 'max',mergdur);  
[noruntime, norun]=Threshold([res_time speed],'<=',norun_thr,'min', mindur, 'max',norunmergdur);  

%Need to find if in these intervals minimum peak speed of 5 cm/sec
peakspeed=res_time(find(speed>=run_thr));
[index wint interval] = InIntervals(peakspeed,runtime);
for u=1:size(runtime,1);
if (sum(wint==u))==0,
inint(u)=NaN;
elseif (sum(wint==u))>0,
inint(u)=sum(wint==u);
end
end
%Remove interval of forward direction when no min peak speed of 5cm/s found
for u=1:size(runtime,1);
if isnan(inint(u))
runtime(u,:)=NaN;
end;end
runtime=runtime(~any(isnan(runtime),2),:);


%Get onset time
for i=1:size(C_df,2); 
    if isempty(onset_offset{i})==0,
on{i}=[Cdf_time(onset_offset{i}(:,1)) Cdf_time(onset_offset{i}(:,2))];
%on{i}=onset_offset{i}(:,1)*(1/fr_im);
%Is onset time in time_run?
run_onset_offset_time{i}=Restrict(on{i}, runtime);
norun_onset_offset_time{i}=Restrict(on{i}, noruntime);
run_onset{i}=run_onset_offset_time{i}(:,1);
norun_onset{i}=norun_onset_offset_time{i}(:,1);
 

runonfr{i}= ismember(Cdf_time,run_onset_offset_time{i}(:,1));
runonidx{i} = find(runonfr{i});
norunonfr{i}= ismember(Cdf_time,norun_onset_offset_time{i}(:,1));
norunonidx{i} = find(norunonfr{i});
 if size(run_onset_offset_time{i},2)==1
 run_onset_offset_time{i}(:,2)= run_onset_offset_time{i}(:,1);
 end
 runofffr{i}=ismember(Cdf_time,run_onset_offset_time{i}(:,2));
 runoffidx{i} = find(runofffr{i});
if size(norun_onset_offset_time{i},2)==1
 norun_onset_offset_time{i}(:,2)= norun_onset_offset_time{i}(:,1);
 end
 norunofffr{i}=ismember(Cdf_time,norun_onset_offset_time{i}(:,2));
 norunoffidx{i} = find(norunofffr{i}); 
 
 
 run_onset_offset{i}=[runonidx{i}  runoffidx{i}];
 norun_onset_offset{i}=[norunonidx{i}  norunoffidx{i}];
    end
end

time_Cdf=Cdf_time;
%get the binary
for i=1:length(run_onset); 
    for ii=1:length(run_onset{i})
    [v, index_runon{i}(ii,:)] = min(abs(time_Cdf(:) - run_onset{i}(ii,:))); 
    end
end

  if length(run_onset)>length(index_runon),
index_runon=[index_runon {[]}];
  end

for i=1:size(C_df,2);
    run_onset_binary{i}=zeros(length(time_Cdf),1);
run_onset_binary{i}(index_runon{i},1)=1;
    end;
    
%get the binary
for i=1:length(norun_onset); 
    for ii=1:length(norun_onset{i})
    [nov, index_norunon{i}(ii,:)] = min(abs(time_Cdf(:) - norun_onset{i}(ii,:))); 
    end
end

 while length(norun_onset)>length(index_norunon),
index_norunon=[index_norunon {[]}];
  end

for i=1:size(C_df,2);
    norun_onset_binary{i}=zeros(length(time_Cdf),1);
norun_onset_binary{i}(index_norunon{i},1)=1;
    end;

    
run_onset_binary=cell2mat(run_onset_binary);
norun_onset_binary=cell2mat(norun_onset_binary);
  

Behavior.runtime=runtime;
Behavior.noruntime=noruntime;
Behavior.runbinary=run;
Behavior.norunbinary=norun;
Behavior.speed=speed;
Behavior.resampled.time=res_time_position(:,1);
Behavior.resampled.cumulativeposition=res_location;
Behavior.resampled.position=res_position;

Events.RunningEpochs.run_onset=run_onset;
Events.RunningEpochs.run_onset_offset=run_onset_offset;
Events.RunningEpochs.run_onset_binary=run_onset_binary;
Events.NonRunningEpochs.norun_onset=norun_onset;
Events.NonRunningEpochs.norun_onset_offset=norun_onset_offset;
Events.NonRunningEpochs.norun_onset_binary=norun_onset_binary;



   figure; hold on;
   plot(res_time, speed);
   plot(time_Cdf, onset_binary(:,c2plot)*5-10);
   plot(time_Cdf, run_onset_binary(:,c2plot)*5-5);
   plot(time_Cdf, C_df(:,c2plot)*10-10);
   
   
   
   figure; hold on;
   plot(res_time, speed);
   PlotHVLines(norun_thr,'h', 'r');
   plot(time_Cdf, onset_binary(:,c2plot)*5-10);
   plot(time_Cdf, norun_onset_binary(:,c2plot)*5-5);
   plot(time_Cdf, C_df(:,c2plot)*10-10);
 
end



