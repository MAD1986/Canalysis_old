
function [Imaging, Behavior]= lapselect(C_df, Behavior, XML, options); 


%% Import data
lap_start_stop=Behavior.lap;
behavior_time=Behavior.time;
position=Behavior.position;
norm_position=Behavior.normalizedposition;
cum_position=Behavior.cumulativeposition;
Imaging.trace=C_df;

switch options.startlap
case 'first'
startlap= 1;
otherwise 
startlap=options.startlap;
end

switch options.endlap
case 'last'
endlap=length(Behavior.lap);
otherwise 
endlap=options.endlap;
end

if options.textures==true
texture=Behavior.texture;
end


%% Extract Timestamps from XML file

%get framelength of imaging session from XML
frameLength = length(XML.PVScan.Sequence.Frame);

%preallocate array
timeStampsXML = zeros(frameLength,1);

%extract the relative timestamps for each frame acquired in session
for ii=1:frameLength
    timeStampsXML(ii,1) = str2double(XML.PVScan.Sequence.Frame{1,ii}.Attributes.relativeTime);
end

%end frame based on XML parsing
t = timeStampsXML(2:end);
%time_cdf=[t C_df];


%% Restrict Imaging time and behavior time

% Set start and end time
startT=lap_start_stop{startlap}(:,1);
endT=lap_start_stop{endlap}(:,2);

% indices for restricted imaging time
t_in=find(t>=startT & t<=endT);
% restrict calcium trace
C_df_R=C_df(t_in,:);
% restrict imaging time
Cdf_time_R=t(t_in);

% indices for restricted behavior time
b_in=find(behavior_time>=startT & behavior_time<=endT);
% restrict behavior : position (norm and cum) and time
behavior_time_R=behavior_time(b_in);
position_R=position(b_in);
norm_position_R=norm_position(b_in);
cum_position_R=cum_position(b_in);

% restrict texture 
if options.textures==true
texture_R=texture(b_in,:);
Behavior.restricted.texture=texture_R;
end

% Restrict lap
for i=startlap:endlap;
lap_start_stop_restricted{i}=lap_start_stop{i};
end
lap_start_stop_restricted = lap_start_stop_restricted(~cellfun(@isempty, lap_start_stop_restricted));

% Make structure
Behavior.restricted.time=behavior_time_R;
Behavior.restricted.position=position_R;
Behavior.restricted.cumulativeposition=cum_position_R;
Behavior.restricted.normalizedposition=norm_position_R;
Behavior.restricted.lap =lap_start_stop_restricted;
Imaging.time=t;
Imaging.time_restricted=Cdf_time_R;
Imaging.trace_restricted=C_df_R;
Behavior.options=options;

% Display figure
figure;
if options.dispfig==true; 
subplot(2,1,1) 
plot(t, C_df(:,options.c2plot)); hold on; plot(Behavior.time, Behavior.normalizedposition)
subplot(2,1,2) 
plot(Imaging.time_restricted, Imaging.trace_restricted(:,options.c2plot)); hold on; plot(Behavior.restricted.time, Behavior.restricted.normalizedposition)
end
end
