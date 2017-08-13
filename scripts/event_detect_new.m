SDON=2;
SDOFF=0.5;


%% First iteration
%Find values above thr;
on_ones_1=C_df>=SDON*std(C_df);
off_ones_1=C_df<=SDOFF*std(C_df);


[on_off_ones, on_off_binary, on_off] = onset_offset(on_ones_1, off_ones_1, C_df);



%Find onset and index
onset_1=diff(on_ones_1)==1;
offset_1=diff(off_ones_1)==1;
for i=1:size(onset_1,2)
onset_idx_1{i}=find(onset_1(:,i)==1);
end
%Trace from onset to next end 
%Find idx of offset in trace
for i=1:size(onset_idx_1,2)
for ii=1:size(onset_idx_1{i},1)
%on_cdf{i}{ii}=   C_df(onset_idx{i}(ii):size(C_df,1),i);
 off_inon_1{i}{ii}=find(off_ones_1(onset_idx_1{i}(ii):size(C_df,1),i)==1);
end
end
 %If no offset until the end, add end of recording as offset
for i=1:size(off_inon_1,2)
    for ii=1:size(off_inon_1{i},2)
       
if isempty(off_inon_1{i}{ii})==0
 off_inon_idx_1{i}(ii,:)=off_inon_1{i}{ii}(1)+onset_idx_1{i}(ii)-1;   
    
elseif  isempty(off_inon_1{i}{ii})  
 off_inon_idx_1{i}(ii,:)=size(C_df,1);
end
end
end
for i=1:size(off_inon_1,2)
    if isempty(off_inon_1{i})==0,
on_off_all_1{i}=[onset_idx_1{i} off_inon_idx_1{i} ];
off_ones_first_1{i}=[1;diff(on_off_all_1{i}(:,2))>=1];
off_first_idx_1{i}=find(off_ones_first_1{i}==1);
on_off_1{i}=[on_off_all_1{i}(off_first_idx_1{i},1) on_off_all_1{i}(off_first_idx_1{i},2)]; 
    elseif isempty(off_inon_1{i}),
  on_off_1{i}=[];  
    end
end
%Make binaries
on_off_ones_1=zeros(size(C_df,1), size(C_df,2));
for i=1:size(on_off_1,2)
      if isempty(on_off_1{i})==0,  
    for ii=1:size(on_off_1{i},1)
     on_off_ones_1(on_off_1{i}(ii,1):on_off_1{i}(ii,2),i)=1;
on_off_binary_1(:,i)=diff(on_off_ones_1(:,i))==1;
         
end
end
end

%% Second iteration
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_idx_1{i}=find(on_off_ones_1(:,i)==0);
Cdf_off_1{i}=C_df(Cdf_off_idx_1{i},i);
std_Cdf_off(i)=std(Cdf_off_1{i});
end

on_ones_2=C_df>=SDON*std_Cdf_off;
off_ones_2=C_df<=SDOFF*std_Cdf_off;
%Find onset and index
onset_2=diff(on_ones_2)==1;
offset_2=diff(off_ones_2)==1;
for i=1:size(onset_2,2)
onset_idx_2{i}=find(onset_2(:,i)==1);
end
%Trace from onset to next end 
%Find idx of offset in trace
for i=1:size(onset_idx_2,2)
for ii=1:size(onset_idx_2{i},1)
%on_cdf{i}{ii}=   C_df(onset_idx{i}(ii):size(C_df,1),i);
 off_inon_2{i}{ii}=find(off_ones_2(onset_idx_2{i}(ii):size(C_df,1),i)==1);
end
end
 %If no offset until the end, add end of recording as offset
for i=1:size(off_inon_2,2)
    for ii=1:size(off_inon_2{i},2)
       
if isempty(off_inon_2{i}{ii})==0
 off_inon_idx_2{i}(ii,:)=off_inon_2{i}{ii}(1)+onset_idx_2{i}(ii);   
    
elseif  isempty(off_inon_2{i}{ii})  
 off_inon_idx_2{i}(ii,:)=size(C_df,1);
end
end
end
for i=1:size(off_inon_2,2)
on_off_all_2{i}=[onset_idx_2{i} off_inon_idx_2{i} ];
off_ones_first_2{i}=[1;diff(on_off_all_2{i}(:,2))>=1];
off_first_idx_2{i}=find(off_ones_first_2{i}==1);
on_off_2{i}=[on_off_all_2{i}(off_first_idx_2{i},1) on_off_all_2{i}(off_first_idx_2{i},2)]; 
end
%Make binaries
on_off_ones_2=zeros(size(C_df,1), size(C_df,2));
for i=1:size(on_off_2,2)
    for ii=1:size(on_off_2{i},1)
on_off_ones_2(on_off_2{i}(ii,1):on_off_2{i}(ii,2),i)=1;
on_off_binary_2(:,i)=diff(on_off_ones_2(:,i))==1;
end
end


%% Third iteration
%Change C_df for only trace between events
for i=1:size(C_df,2)
Cdf_off_idx_2{i}=find(on_off_ones_2(:,i)==0);
Cdf_off_2{i}=C_df(Cdf_off_idx_2{i},i);
std_Cdf_off_2(i)=std(Cdf_off_2{i});
end

on_ones_3=C_df>=SDON*std_Cdf_off_2;
off_ones_3=C_df<=SDOFF*std_Cdf_off_2;
%Find onset and index
onset_3=diff(on_ones_3)==1;
offset_3=diff(off_ones_3)==1;
for i=1:size(onset_3,2)
onset_idx_3{i}=find(onset_3(:,i)==1);
end
%Trace from onset to next end 
%Find idx of offset in trace
for i=1:size(onset_idx_3,2)
for ii=1:size(onset_idx_3{i},1)
%on_cdf{i}{ii}=   C_df(onset_idx{i}(ii):size(C_df,1),i);
 off_inon_3{i}{ii}=find(off_ones_2(onset_idx_3{i}(ii):size(C_df,1),i)==1);
end
end
 %If no offset until the end, add end of recording as offset
for i=1:size(off_inon_3,2)
    for ii=1:size(off_inon_3{i},2)
       
if isempty(off_inon_3{i}{ii})==0
 off_inon_idx_3{i}(ii,:)=off_inon_3{i}{ii}(1)+onset_idx_3{i}(ii);   
    
elseif  isempty(off_inon_3{i}{ii})  
 off_inon_idx_3{i}(ii,:)=size(C_df,1);
end
end
end
for i=1:size(off_inon_3,2)
on_off_all_3{i}=[onset_idx_3{i} off_inon_idx_3{i} ];
off_ones_first_3{i}=[1;diff(on_off_all_3{i}(:,2))>=1];
off_first_idx_3{i}=find(off_ones_first_3{i}==1);
on_off_3{i}=[on_off_all_3{i}(off_first_idx_3{i},1) on_off_all_3{i}(off_first_idx_3{i},2)]; 
end
%Make binaries
on_off_ones_3=zeros(size(C_df,1), size(C_df,2));
for i=1:size(on_off_3,2)
    for ii=1:size(on_off_3{i},1)
on_off_ones_3(on_off_3{i}(ii,1):on_off_3{i}(ii,2),i)=1;
on_off_binary_3(:,i)=diff(on_off_ones_3(:,i))==1;
end
end


%% 
c2plot=15;
figure; hold on;
plot(C_df(:,c2plot))
plot(on_off_ones(:,c2plot));
refline([0 SDON*std(C_df(:,c2plot))])
%refline([0 SDON*std_Cdf_off(:,c2plot)])
%refline([0 SDON*std_Cdf_off_2(:,c2plot)])

refline([0 SDOFF*std(C_df(:,c2plot))])
%refline([0 SDOFF*std_Cdf_off(:,c2plot)])
%refline([0 SDOFF*std_Cdf_off_2(:,c2plot)])

