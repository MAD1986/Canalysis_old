
%Use center
%Distance Threshold 
thr=5;

%Combine files
directory_name = uigetdir %folder with .mat
cd(directory_name);
listmat = dir('*.mat'); %will look for all the mat files
%!!!! files are listed by alphabetical order !!!!
%!!!! Keep same name for mat !!!!
for i = 1:length(listmat)
matfile=listmat(i).name;
coor{i} = load((listmat(i).name));
center{i}=coor{i}.center_new;
end
%Euclidean distance between coordinates
for i=1:size(center,2)
    for ii=1:size(center,2)
[D{i}{ii},I{i}{ii}] = pdist2(center{i},center{ii},'euclidean','Smallest',1 );
    end
end

ses2match=[1 2] %Find matching ROI for session nb
for i=ses2match
Dis{i}=D{ses2match(1)}{ses2match(i)};
Ind{i}=I{ses2match(1)}{ses2match(i)};
closeROI{i}=(D{ses2match(1)}{ses2match(i)}<=thr);
for ii=1:size(closeROI{i},2);
if closeROI{i}(ii)==1
sameROI{i}(ii)=I{ses2match(1)}{ses2match(i)}(ii);
end
end
end
for i=ses2match
while size(sameROI{i},2)<size(center{i},1);
sameROI{i}=[sameROI{i}, 0];
end
end

for i=ses2match
ROI_all{i}(:,2)=find(sameROI{i}>=1)
ROI_all{i}(:,1)=sameROI{i}(ROI_all{i}(:,2))
end
     


for i=1:size(ROI_all{1},1)
%figure;    
%plot(C_df1(ROI_all{1}(i,1),:))
%hold on; plot(C_df2(ROI_all{2}(i,2),:))
end

