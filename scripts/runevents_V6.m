
function [Events, Behavior]=runevents_V6(Events,Behavior,Imaging);
%% Only running onset
run_thr=Events.options.run.speed_thr;
norun_thr=Events.options.norun.speed_thr;
mindur=Events.options.run.min_dur;
norunmindur=Events.options.norun.min_dur;
mergdur=Events.options.run.merge_dur;
norunmergdur=Events.options.norun.merge_dur;
span=Events.options.span;
fr_im=Imaging.options.frimaging;


if (Events.options.baselinesub==1)&& (Events.options.restricted==1),
C_df=Imaging.trace_restricted_baselinesub;  
Cdf_time=Imaging.time_restricted;
time=Behavior.restricted.time;
location=Behavior.restricted.cumulativeposition;
position=Behavior.restricted.position;
else if (Events.options.baselinesub==0) && (Events.options.restricted==1),
C_df=Imaging.trace_restricted;
Cdf_time=Imaging.time_restricted;
time=Behavior.restricted.time;
location=Behavior.restricted.cumulativeposition;
position=Behavior.restricted.position;
else if (Events.options.baselinesub==1) && (Events.options.restricted==0),
C_df=Imaging.trace_baselinesub;
Cdf_time=Imaging.time;
time=Behavior.time;
location=Behavior.cumulativeposition;
position=Behavior.position;
else if (Events.options.baselinesub==0 && Events.options.restricted==0),
C_df=Imaging.trace;
Cdf_time=Imaging.time;
time=Behavior.time;
location=Behavior.cumulativeposition;
position=Behavior.position;
end;end;end;end

onset_offset=Events.onset_offset;
onset_binary=Events.onset_binary;
c2plot=Imaging.options.celltoplot;
time_position=[time location];


%Resample at same frequency than imaging
[N,bin]=histc(Cdf_time,time);
index=bin+1;
if abs(Cdf_time-time(bin))<abs(Cdf_time-time(bin+1));
index=bin;
res_time_position=time_position(bin,:);
res_location=location(bin);
res_position=position(bin);
else
res_time_position=time_position(bin,:);
res_location=location(bin);
res_position=position(bin);
end

%res_time=  decimate(time,fs);
res_time=res_time_position(:,1);


%Getting the speed
%speed=diff(res_location)*fr_w;
speed=diff(res_location)*fr_im;
%speed(speed<=-10)=0; % Remove negative speed
%speed(speed>=100)=0; % Remove end bug  
speed=[speed; speed(end)]; %

%Smooth the speed, optional?
%Running average span
speed=Smooth(speed,span); 

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



