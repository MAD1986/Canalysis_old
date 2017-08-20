
function [Events]=sigma_events(Imaging,options);

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
SDOFF=options.SDOFF; %Threshold below x SD for OFFSET 
SD=options.SD;

%% Detect events

% Using 2, 3 and 4 std
for sig=1:size(SD,2);
 %positive events
on_std{sig}=C_df>=(SD(sig)*std(C_df));
off_std{sig}=C_df<=SDOFF*std(C_df);

 %negative events
on_neg_std{sig}=C_df<=-(SD(sig)*std(C_df));
off_neg_std{sig}=C_df>=-SDOFF*std(C_df);
end

 %Find onset and offset, do multiple iterations

 %positive
for sig=1:size(SD,2);
on_ones=on_std{sig};
off_ones=off_std{sig};
for it=nb_it % nb of iteration 
[on_off_ones, on_off_binary, on_off] = on_off_thr(on_ones, off_ones, C_df);
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_idx{i}=find(on_off_ones(:,i)==0);
Cdf_off{i}=C_df(Cdf_off_idx{i},i);
std_Cdf_off(i)=std(Cdf_off{i});
end
on_ones=C_df>=SD(sig)*std_Cdf_off;
off_ones=C_df<=SDOFF*std_Cdf_off;
end
on_off_std{sig}=on_off_ones;
on_off_binary_std{sig}=on_off_binary;
on_off_sig{sig}=on_off;
end


 %Negative
for sig=1:size(SD,2);
on_ones_neg=on_neg_std{sig};
off_ones_neg=off_neg_std{sig};
for it=nb_it % nb of iteration 
[on_off_neg_ones, on_off_neg_binary, on_off_neg] = on_off_thr(on_ones_neg, off_ones_neg, C_df);
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_neg_idx{i}=find(on_off_neg_ones(:,i)==0);
Cdf_off_neg{i}=C_df(Cdf_off_neg_idx{i},i);
std_Cdf_off_neg(i)=std(Cdf_off_neg{i});
end
on_ones_neg=C_df<=-SD(sig)*std_Cdf_off_neg;
off_ones=C_df>=-SDOFF*std_Cdf_off_neg;
end
on_off_neg_std{sig}=on_off_neg_ones;
on_off_binary_neg_std{sig}=on_off_neg_binary;
on_off_neg_sig{sig}=on_off_neg;
end

%% Plot event duration

%Find event duration : time end - time start
% positives
for sig=1:size(SD,2);
 on_off=on_off_sig{sig};   
for i=1:size(on_off,2)
    if isempty(on_off{i})==0
eventduration{i}=Cdf_time(on_off{i}(:,2))-Cdf_time(on_off{i}(:,1));
end 
end
hist_duration{sig}=cell2mat(eventduration');
end
 % Negatives
for sig=1:size(SD,2);
 on_off=on_off_neg_sig{sig};   
for i=1:size(on_off,2)
    if isempty(on_off{i})==0
eventduration_neg{i}=-(Cdf_time(on_off{i}(:,2))-Cdf_time(on_off{i}(:,1)));
    elseif isempty(on_off{i})==1
eventduration_neg{i}=0;
end 
end
hist_duration_neg{sig}=cell2mat(eventduration_neg');
end

%all events (positives and negatives) histogram
figure;
for sig=1:size(SD,2);
hist_events{sig}=[hist_duration_neg{sig}; hist_duration{sig}];
ax{sig}=subplot(size(SD,2),1,sig);
max_X(sig)=max(hist_events{sig});
min_X(sig)=min(hist_events{sig});
histogram(hist_duration{sig})
hold on;
histogram(hist_duration_neg{sig})
%axis([ax{sig}],[min(min_X) max(max_X) 0 inf])
%axis([ax{sig}],[-5 5 0 inf])
title(['Cumulative distribution of events >=' num2str(SD(sig)) ' standard deviation'])
xlabel('Transient duration >= (s)'); 
ylabel('Number of events');
end

%% Plot false positive rate / event duration

%false positive rate (nb negative/positive events) total
%all duration
for sig=1:size(SD,2);
 count_pos=hist_duration{sig}; 
 count_neg=hist_duration_neg{sig}; 
%remove 0
count_pos(count_pos==0)=[];
count_neg(count_neg==0)=[];
%count number of events
count_pos=length(count_pos);
count_neg=length(count_neg);

nb_pos(sig)=count_pos;
nb_neg(sig)=count_neg;
end 
false_ratio=nb_neg./nb_pos;


%false positive rate (nb negative/positive events) total
%multiple duration (every 250ms bins) from 0 to 5 sec duration events
dur_max=4; % 0 to 5 sec max duration events
binsz=0.25; % bin every .25 sec

% bin duration

for i=1:dur_max/binsz
edges(i) = i*binsz;
end

for sig=1:size(SD,2)
count_pos=hist_duration{sig}; 
count_neg=-hist_duration_neg{sig}; 
[N_pos] = histcounts(count_pos,edges);
[N_neg] = histcounts(count_neg,edges);
count_dur_pos(sig,:)=N_pos;
count_dur_neg(sig,:)=N_neg;
end

false_ratio_bin=count_dur_neg./count_dur_pos;


% Make figure 
figure;
for sig=1:size(SD,2)
plotStyle = {'b.','k.','r.'}; 
false_ratio_plot{sig}=[edges(1:end-1)' false_ratio_bin(sig,:)'];
x=false_ratio_plot{sig}(:,1);
y=false_ratio_plot{sig}(:,2);
f = fit(x,y,'exp1');
plot(x,y, '-o');
%plot(f)
xlabel('Transient duration (s)'); 
ylabel('Error rate (neg./pos. events)'); 
legendInfo{sig} = ([num2str(SD(sig)),' sigma events']);
hold on;
end
legend(legendInfo)

%% Make structure
Events.figures.hist_duration_pos=hist_duration;
Events.figures.hist_duration_neg=hist_duration_neg;
Events.figures.error_rate=false_ratio_plot;

end
