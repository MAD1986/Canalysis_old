function [onset_offset] = onoff_V4(onset_time, C_df, offset);



for i=1:size(onset_time,2);
    
start_onset{i}=onset_time{i}(:,1);
%A{i} = start_onset{i}(1:2:end,:); % Odd 
%B{i} = start_onset{i}(2:2:end,:); % Even

%   if length(A{i})>length(B{i})==1 
%   ,B{i}(end+1,1)=A{i}(end,1);
%   end
%  From onset start to next onset start -1 
start_onset_interval{i}=[start_onset{i} [(start_onset{i}(2:end)-1); length(C_df)]];
offset_num{i}=[(1:size(C_df))' offset{i}];
end
%Fix if bug
for i=1:size(onset_time,2);
%  if isempty(start_onset_interval{i})==1
%  ,start_onset_interval{i}(1:2,1:2)=2;
%  end

if start_onset_interval{i}(1)==0
,start_onset_interval{i}(1)=1;
end
 if size(start_onset_interval{i},2)==1,
start_onset_interval{i}(:,2)=start_onset_interval{i}(:,1);
   end  
end

  for i=1:length(onset_time);
    for ii=1:size(start_onset_interval{i},1);
        
j{ii}=(start_onset_interval{i}(ii,1):(start_onset_interval{i}(ii,2)));

offset_inonset_sum{i}(ii,:)=sum(offset_num{i}(j{ii}, 2));
offset_inonset{i}{ii}=(offset_num{i}(j{ii}, :));

   if offset_inonset{i}{ii}(1,2)==1
       ,offset_inonset{i}{ii}(1,2)=0;
    end
    
offset_inonset_sorted{i}{ii}=sortrows(offset_inonset{i}{ii}, -2);
%If not offset in the onset period add NaN for the second offset
%(ex: 2 events (burst) and 1st one doesn't fall to offset before the 2nd one)
    end
  end
  
 %   for i=1:length(onset_time);
 %   for ii=1:size(start_onset_interval{i},1);
%if  offset_inonset_sorted{i}{ii}(1,2)==0,
%second_onset{i}(ii)=ii;
%second_onset{i}(second_onset{i}==0) = [];
%end
%end
 %   end
    
   % for i=1:length(onset_time);
   %   for x=second_onset{i}
   %       if second_onset{i}(end)<= size(offset_inonset_sorted{i},2)==1
  %      offset_inonset_sorted{i}{x}(1,1)=NaN;
  %      end
   % end
   % end
    for i=1:length(onset_time);
    for ii=1:size(start_onset_interval{i},1);    
    offset_t{i}(ii,:)=offset_inonset_sorted{i}{ii}(1,1);
  
end
  end
  



  
for i=1:length(onset_time);
start_stop{i}=[(start_onset_interval{i}(:,1)) (offset_t{i})];

end


%Need to take the first peak ex Cell 101 / 47
for i=1:length(onset_time);
    for ii=1:size(start_stop{i},1);
if isnan(start_stop{i}(end,2)),
  start_stop{i}(end,2)= length(C_df);
elseif isnan(start_stop{i}(ii,2)),
   start_stop{i}(ii,2)=start_stop{i}(ii+1,2);
end
    end;end

for i=1:length(onset_time);
    start_stop_NaN{i}=[start_stop{i}(:,1) [1; diff(start_stop{i}(:,2))]];
end

for i=1:length(onset_time);
    for ii=1:size(start_stop_NaN{i},1);
if (start_stop_NaN{i}(ii,2))==0,
start_stop{i}(ii,2)=NaN;
end;
    end;end

for i=1:length(onset_time);
 SS= start_stop{i};
 SS= SS(0== sum(isnan(SS), 2), :);
 onset_offset{i}=SS;
end



 
for i=1:size(onset_offset,2);
onset_in{i}=onset_offset{i}(:,2)-onset_offset{i}(:,1);
end


%If onset_in = 0 the onset is still in a zone above threshold
for i=1:size(onset_in,2);
  for ii=1:size(onset_in{i},1);
      
      if onset_in{i}(ii)==0,
   test{i}(ii,:)=onset_offset{i}(ii+1 );
      end
       if onset_in{i}(ii)>=1,
  test{i}(ii,:)=onset_offset{i}(ii);
  end
end
end

      
      

        

%for ii=1:length(onset_time);
 % for i=1:size(start_onset_interval{i},1);
 
 %       if isnan(start_stop{ii}(i,2))==1
%           ,start_stop{ii}(i+1,1)=NaN;
 %       end
%end
%end
    

%for ii=1:length(onset_time);
% OO1{ii}=start_stop{ii}(:,1);
% OO1{ii}(isnan(OO1{ii}))=[];
% OO2{ii}=start_stop{ii}(:,2);
% OO2{ii}(isnan(OO2{ii}))=[];
% onset_offset{ii}=[OO1{ii} OO2{ii}];
%end
end
