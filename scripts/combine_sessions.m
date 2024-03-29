
function [C_df,CSV,XML]=combine_sessions(directory_name,type);

cd(directory_name);
listmat = dir('*.mat'); %will look for all the mat files
listcsv = dir('*.csv');
listxml = dir('*.xml');
%!!!! files are listed by alphabetical order !!!!
%!!!! Keep same name for mat and csv !!!!
for i = 1:length(listmat)
csvfile{i}=listcsv(i).name;
disp(['Reading CSV ' ,num2str(i),'/', num2str(length(listmat))]);
CSV{i} = csvread(csvfile{i},1,0);
disp(['Reading XML ' ,num2str(i),'/', num2str(length(listmat))]);
XML{i} = xml2structV2((fullfile(listxml(i).folder,listxml(i).name)));
end
switch type
    case 'expdff'
for i = 1:length(listmat)
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.expDffMedZeroed;
C_df{i}=(full(C_df{i})');    
end
    case 'expdff'
for i = 1:length(listmat)
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.expDffMedZeroed;
C_df{i}=(full(C_df{i})');
end
    case 'spikes'
for i = 1:length(listmat)
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.S;
C_df{i}=(full(C_df{i})');
end
end
end
