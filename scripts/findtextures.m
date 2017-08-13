

function [Behavior] = findtextures(CSV,Behavior)
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
texture=singleStim1CsvOut(:,3);
cum_position=positionVector*0.0155986;
count=1:idx;



tex{1}=[0.5 1]; %min to max threshold 
tex{2}=[1 1.5]; 
tex{3}=[1.5 2];
tex{4}=[2 2.5];
tex{5}=[2.5 5];


[PKS,LOCS]= findpeaks(texture, 'MinPeakHeight',0.5,'MinPeakDistance',25000);
for i=1:length(PKS);
for ii=1:size(tex,2);
if ((PKS(i)>=tex{ii}(1))) && ((PKS(i)<=tex{ii}(end)))
findtex{ii}(i)=1;
if 1:length(findtex{ii})<length(PKS);
findtex{ii}=[findtex{ii} (zeros(length(PKS)-length(findtex{ii}),1))'];
end
end
end
end
for i=1:length(PKS);
for ii=1:size(tex,2);
texind{ii}=find(findtex{ii}==1);
texloc{ii}=LOCS(texind{ii});
end
end
tex_binary=(zeros(size(texture,1),size(tex,2)));
for ii=1:size(tex,2);
tex_binary(texloc{ii},ii)=1;
end


Behavior.cumulativeposition=cum_position;
%Behavior.lap=lap_start_stop;
Behavior.lick=lick;
Behavior.texture=texture;
Behavior.texture_binary=tex_binary;

    
    
        
 