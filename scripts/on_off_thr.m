function [on_off_ones, on_off_binary, on_off] = on_off_thr(on, off, C_df);


%Find onset and index
onset=diff(on)==1;
offset=diff(off)==1;
for i=1:size(onset,2)
onset_idx{i}=find(onset(:,i)==1);
end
%Trace from onset to next end 
%Find idx of offset in trace
for i=1:size(onset_idx,2)
for ii=1:size(onset_idx{i},1)
%on_cdf{i}{ii}=   C_df(onset_idx{i}(ii):size(C_df,1),i);
 off_inon{i}{ii}=find(off(onset_idx{i}(ii)+1:size(C_df,1),i)==1);
end
end
 %If no offset until the end, add end of recording as offset
for i=1:size(off_inon,2)
    for ii=1:size(off_inon{i},2)
       
if isempty(off_inon{i}{ii})==0
 off_inon_idx{i}(ii,:)=off_inon{i}{ii}(1)+onset_idx{i}(ii)-1;   
    
elseif  isempty(off_inon{i}{ii})  
 off_inon_idx{i}(ii,:)=size(C_df,1);
end
end
end
for i=1:size(off_inon,2)
    if isempty(off_inon{i})==0,
on_off_all{i}=[onset_idx{i} off_inon_idx{i} ];
off_ones_first{i}=[1;diff(on_off_all{i}(:,2))>=1];
off_first_idx{i}=find(off_ones_first{i}==1);
on_off{i}=[on_off_all{i}(off_first_idx{i},1) on_off_all{i}(off_first_idx{i},2)]; 
    elseif isempty(off_inon{i}),
  on_off{i}=[];  
    end
end
%Make binaries
on_off_ones=zeros(size(C_df,1), size(C_df,2));
for i=1:size(on_off,2)
      if isempty(on_off{i})==0,  
    for ii=1:size(on_off{i},1)
     on_off_ones(on_off{i}(ii,1):on_off{i}(ii,2),i)=1;
on_off_binary(:,i)=diff(on_off_ones(:,i))==1;
         
end
end
end



