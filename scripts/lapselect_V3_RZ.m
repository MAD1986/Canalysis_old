
function [Imaging, Behavior]= lapselect_V3_RZ(C_df, Behavior, Imaging, XML); 

%V3 RZ - uses frame timestamps from associated xml file for each session

frimage=Imaging.options.frimaging;
startlap=Behavior.restricted.startlap;
endlap=Behavior.restricted.endlap;
lap_start_stop=Behavior.lap;
time_position=[Behavior.time Behavior.position];
cum_position=Behavior.cumulativeposition;
Imaging.trace=C_df;
%texture=Behavior.texture;

%Extract Timestamps from XML file

%get framelength of imaging session from XML
frameLength = length(XML.PVScan.Sequence.Frame);

%preallocate array
timeStampsXML = zeros(frameLength,1);

%extract the relative timestamps for each frame acquired in session
for ii=1:frameLength
    timeStampsXML(ii,1) = str2double(XML.PVScan.Sequence.Frame{1,ii}.Attributes.relativeTime);
end

%previous multiplication expansion (V2)
%t=1:size(C_df,1);
%time_cdf=[(t*(1/frimage))' C_df];

%end frame based on XML parsing
t = timeStampsXML(2:end);
time_cdf=[t C_df];

lap_start_stop{startlap}(:,1);
lap_start_stop{endlap}(:,2);
C_df_restricted=Restrict(time_cdf, [lap_start_stop{startlap}(:,1) lap_start_stop{endlap}(:,2)]);
Cdf_time=C_df_restricted(:,1);
C_df_restricted=C_df_restricted(:,2:end);

time_position_restricted=Restrict(time_position, [lap_start_stop{startlap}(:,1) lap_start_stop{endlap}(:,2)]);
cum_position_restricted=Restrict([time_position(:,1) cum_position], [lap_start_stop{startlap}(:,1) lap_start_stop{endlap}(:,2)]);
cum_position_restricted=cum_position_restricted(:,2);
position_norm_restricted=(time_position_restricted(:,2)-min(time_position_restricted(:,2)))/(max(time_position_restricted(:,2))-min(time_position_restricted(:,2)));
%texture_restricted=Restrict([time_position(:,1) texture] , [lap_start_stop{startlap}(:,1) lap_start_stop{endlap}(:,2)]);

for i=startlap:endlap;
lap_start_stop_restricted{i}=lap_start_stop{i};
end
lap_start_stop_restricted = lap_start_stop_restricted(~cellfun(@isempty, lap_start_stop_restricted));

Behavior.restricted.time=time_position_restricted(:,1);
Behavior.restricted.position=time_position_restricted(:,2);
Behavior.restricted.cumulativeposition=cum_position_restricted;
Behavior.restricted.normalizedposition=position_norm_restricted;
Behavior.restricted.lap =lap_start_stop_restricted;
%Behavior.restricted.texture=texture_restricted(:,2);
Imaging.time=time_cdf(:,1);
Imaging.time_restricted=Cdf_time;
Imaging.trace_restricted=C_df_restricted;



figure; 
plot(Imaging.time, C_df(:,Imaging.options.celltoplot)); hold on; plot(Behavior.time, Behavior.normalizedposition)
figure; 
plot(Imaging.time_restricted, Imaging.trace_restricted(:,Imaging.options.celltoplot)); hold on; plot(Behavior.restricted.time, Behavior.restricted.normalizedposition)
end
