

function [Behavior] = behavior_multiple_notexture(CSV,Behavior)
%% Rotary decoder
%1 position step = 0.155986  mm = 0.000


mincross=Behavior.options.mintimebeforelap; %time to wait before counting a new lap (s)
acqfr=Behavior.options.acquisitionfrequency;
singleStim1CsvOut=CSV;

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
%% Get the position to restart at each lap 
lap=csvStimRaw(:,8);
%texture=singleStim1CsvOut(:,3);
cum_position=positionVector*0.0155986;
count=1:idx;



%tex{1}=[0.5 1]; %min to max threshold 
%tex{2}=[1 1.5]; 
%tex{3}=[1.5 2];
%tex{4}=[2 2.5];
%tex{5}=[2.5 5];


%[PKS,LOCS]= findpeaks(texture, 'MinPeakHeight',0.5,'MinPeakDistance',25000);
%for i=1:length(PKS);
%for ii=1:size(tex,2);
%if ((PKS(i)>=tex{ii}(1))) && ((PKS(i)<=tex{ii}(end)))
%findtex{ii}(i)=1;
%if 1:length(findtex{ii})<length(PKS);
%findtex{ii}=[findtex{ii} (zeros(length(PKS)-length(findtex{ii}),1))'];
%end
%end
%end
%end
%for i=1:length(PKS);
%for ii=1:size(tex,2);
%texind{ii}=find(findtex{ii}==1);
%texloc{ii}=LOCS(texind{ii});
%end
%end
%tex_binary=(zeros(size(texture,1),size(tex,2)));
%for ii=1:size(tex,2);
%tex_binary(texloc{ii},ii)=1;
%end

%lap=tex_binary(

    
%for i=1:length(texture);
%for ii=1:size(tex,2);
%if ((texture(i)>=tex{ii}(1))) && ((texture(i)<=tex{ii}(end)))
%findtex{ii}(i)=1;
%if 1:length(findtex{ii})<length(texture);
%findtex{ii}=[findtex{ii} (zeros(length(texture)-length(findtex{ii}),1))'];
%end
%end
%end
%end
%findtex_mat=cell2mat(findtex);
%tex_binary=reshape(findtex_mat,[],size(findtex,2));


[periods, laptime] = Threshold([count' lap],'>=',1, 'max', mincross*acqfr );
   lap_binary=diff(laptime)==1; 
  lap_binary=[lap_binary; 0];

lap_row=periods(:,1)';
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
lap1=cum_position(1:lap_row(1),1);
lastlap=cum_position(lap_row(end):end)-cum_position(lap_row(end));
lastlap_norm=(lastlap - min(lastlap)) / ( max(lastlap) - min(lastlap) );
lap1_norm=(lap1 - min(lap1)) / ( max(lap1) - min(lap1) );
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


Behavior.time=time;
Behavior.position=position;
Behavior.normalizedposition=position_norm;
Behavior.cumulativeposition=cum_position;
Behavior.lap=lap_start_stop;
Behavior.lick=lick;
%Behavior.texture=texture;


figure; plot(Behavior.time,Behavior.normalizedposition);
hold on
plot(Behavior.time,lap_binary, 'r');
%plot(Behavior.time,(Behavior.texture)/5, 'g');
end
%%
%Merge position for multiple sessions

% cum_position_1min_tot=cum_position_1min+cum_position_pre(end,1);
% cum_position=[cum_position_pre; cum_position_1min_tot];
% 
% time_1min_tot=time_position_1min(:,1)+time_position_pre(end,1);
% time=[time_position_pre(:,1); time_1min_tot];
% 
% position=[time_position_pre(:,1); time_position_1min(:,1)];
% 
% time_position=[time position];

    
    
    
        
 