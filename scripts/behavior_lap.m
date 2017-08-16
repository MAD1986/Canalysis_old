

function [Behavior] = behavior_lap(CSV,options)

%% Load parameters
mindist=options.mindist; %distance min before counting a new lap (s)
acqfr=options.acqHz; % behavior acquisition frequency (Hz)
RFIDtext=options.textures;
singleStim1CsvOut=CSV;

%% Rotary decoder
%1 position step = 0.155986  mm = 0.000
%Isolate signals from time
csvStimRaw = singleStim1CsvOut(:,2:end);
%Digitize the files
%try this with logical indexing%
csvStimRaw(find((csvStimRaw > 1))) = 1; 
csvStimRaw(find((csvStimRaw < 1))) = 0;
%take out A and B signal from RE
csvABdigitized = csvStimRaw(:,3:4);
csvIndexDigitized = csvStimRaw(:,1);
lick=csvStimRaw(:,7);
%Isolate time and convert to seconds
timeOutput = singleStim1CsvOut(:,1);
timeOutputSec = timeOutput./1000;
%rotary decoder
pastStateA = csvABdigitized(1,1);
pastStateB = csvABdigitized(1,2);
currentStateA = 0;
currentStateB = 0;
positionCount = 0;
positionVector = zeros(size(csvABdigitized,1),1);
%profile on
for idx = 2:size(csvABdigitized,1)
    currentStateA = csvABdigitized(idx,1);
    currentStateB = csvABdigitized(idx,2);    
    if (currentStateA == pastStateA && currentStateB == pastStateB)
        positionVector(idx) =  positionCount;
    elseif (currentStateA ~= pastStateA || currentStateB ~= pastStateB)
        if (currentStateA == 0 && currentStateB == 0)
            if pastStateA == 0
                positionCount = positionCount + 1;             
            elseif pastStateA == 1
                positionCount = positionCount - 1;
            end        
            positionVector(idx) =  positionCount;          
        elseif (currentStateA == 1 && currentStateB == 0)
            if pastStateA == 0
                positionCount = positionCount + 1;              
            elseif pastStateA == 1
                positionCount = positionCount - 1;
            end                    
            positionVector(idx) =  positionCount;           
        elseif (currentStateA == 1 && currentStateB == 1)
            if pastStateA == 0
                positionCount = positionCount - 1;              
           elseif pastStateA == 1
                positionCount = positionCount + 1;
            end
                
                positionVector(idx) =  positionCount;
                       
        elseif (currentStateA == 0 && currentStateB == 1)
            if pastStateA == 0
                positionCount = positionCount - 1;
                
            elseif pastStateA == 1
                positionCount = positionCount + 1;
            end
           
            positionVector(idx) =  positionCount;
        end  
   
    end
    pastStateA = currentStateA;
    pastStateB = currentStateB;
end

cum_position=positionVector*0.0155986;



%% Register postion of RFID for each texture

if RFIDtext==1
RFID=options.RFID;
% textures RFID signal
texture=singleStim1CsvOut(:,3);
% cumulative position : 1 position step = 0.155986
cum_position=positionVector*0.0155986;

% Threshold RDIF signal when higher than first minimun threshold
%not working: use findpeaks function
%THR_RFID=find(texture>=RFID{1}(1));
%POS_RFID=cum_position(THR_RFID);
%PKS_RFID=texture(THR_RFID);
%RFID_keep=[1;diff(POS_RFID)>=5]==1; % 5cm between RFID
%POS_RFID_keep=POS_RFID(RFID_keep);
%PKS_RFID_keep=PKS_RFID(RFID_keep);
%THR_RFID_keep=PKS_RFID(RFID_keep);

% Find peaks higher than first minimun threshold
[PKS,LOCS]= findpeaks(texture, 'MinPeakHeight',RFID{1}(1));
% Find position of each RFID and remove RFID is distance less than 5cm
POS=cum_position(LOCS);
keep_POS=[1;diff(POS)>=5]==1; % 5cm between RFID
PKS_keep=PKS(keep_POS);
LOCS_keep=LOCS(keep_POS);
POS_keep=POS(keep_POS);

% Associate peaks with texture
for i=1:length(PKS_keep);
for ii=1:size(RFID,2);
if ((PKS_keep(i)>=RFID{ii}(1))) && ((PKS_keep(i)<=RFID{ii}(end)))
findtex{ii}(i)=1;
if 1:length(findtex{ii})<length(PKS_keep);
findtex{ii}=[findtex{ii} (zeros(length(PKS_keep)-length(findtex{ii}),1))'];
end
end
end
end
% Find indices and location of textures
for i=1:length(PKS_keep);
for ii=1:size(RFID,2);
texind{ii}=find(findtex{ii}==1);
texloc{ii}=LOCS_keep(texind{ii});
end
end
% Make a binary (0 and 1 when RFID)
tex_binary=(zeros(size(texture,1),size(RFID,2)));
for ii=1:size(RFID,2);
tex_binary(texloc{ii},ii)=1;
end
% Make structure
Behavior.texture=tex_binary;
end

%% Get the position of lap RFID and restart position at each lap 

% lap RFID signals
lap=csvStimRaw(:,8);

% Threshold RDIF signal when egal or higher than 1
lap_RFID=find(lap>=1);

% Find position of each lap RFID 
RFID_POS=cum_position(lap_RFID);
%and remove is distance less than min
RFID_keep=[1;diff(RFID_POS)>=mindist]==1;
RFID_POS_keep=RFID_POS(RFID_keep);
lap_RFID_keep=lap_RFID(RFID_keep);

% Make binary
lap_binary=(zeros(size(lap,1),1));
lap_binary(lap_RFID_keep)=1;

%Restart position at each RFID 
lap_row=lap_RFID_keep;
for i=1:length(lap_row)-1;
    for ii=1:length(lap_row)
pos_cumul{i}=cum_position(lap_row(i):lap_row(i+1),1);
pos_reset{i}=pos_cumul{i}-cum_position(lap_row(i));
lap_start(i)=lap_row(i)+1;
lap_stop(ii)=lap_row(ii);
%normalize position 0 to 1
norm_pos{i} = (pos_reset{i} - min(pos_reset{i})) / ( max(pos_reset{i}) - min(pos_reset{i}) );
end
end
time=timeOutputSec;
% Position for the first incomplete lap
lap1=cum_position(1:lap_row(1),1);
% Position for the last incomplete lap
lastlap=cum_position(lap_row(end):end)-cum_position(lap_row(end));
% Normalize first and last lap
lastlap_norm=(lastlap - min(pos_reset{end})) / ( max(pos_reset{end}) - min(pos_reset{end}) );
lap1_norm=(lap1 - min(pos_reset{1})) / ( max(pos_reset{1}) - min(pos_reset{1}) );
lap1_norm=lap1_norm+1-max(lap1_norm);
pos_reset=[lap1 pos_reset lastlap];
norm_pos=[lap1_norm norm_pos lastlap_norm];
%to remove extra value added:
for i=2:length(pos_reset);
    ii=1:length(pos_reset);
pos_reset{i}=pos_reset{i}(2:end);
norm_pos{i}=norm_pos{i}(2:end);
end
lap_stop=lap_stop(2:end);
for i=1:length(lap_start);
    lap_start(i)=time(lap_start(i));
     lap_stop(i)=time(lap_stop(i));
lap_start_stop{i}=[lap_start(1,i) lap_stop(1,i)];
end

B=reshape(pos_reset,[],1);
N=reshape(norm_pos,[],1);
position=cell2mat(B);
position_norm=cell2mat(N);
lick=csvStimRaw(:,7);
time_position=[time position];


%% Make structure
Behavior.time=time;
Behavior.position=position;
Behavior.normalizedposition=position_norm;
Behavior.cumulativeposition=cum_position;
Behavior.lap=lap_start_stop;
Behavior.lick=lick;

Behavior.options=options;



%% Display figure
if options.dispfig==1
figure; plot(Behavior.time,Behavior.normalizedposition);
hold on
plot(Behavior.time,lap_binary, 'r');
if RFIDtext==1
for i=1:size(tex_binary,2)
M=tex_binary(:,i);
M(M >=1) = i;
plot_tex(:,i)=M;
end
plot(Behavior.time,plot_tex/5, 'g');
end
end

end

    
        
 