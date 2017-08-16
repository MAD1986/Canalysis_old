
function [C_df,CSV]=combine_sessions(directory_name,type);

cd(directory_name);
listmat = dir('*.mat'); %will look for all the mat files
listcsv = dir('*.csv');
listxml = dir('*.xml');
%!!!! files are listed by alphabetical order !!!!
%!!!! Keep same name for mat and csv !!!!
for i = 1:length(listmat)
csvfile{i}=listcsv(i).name;
disp('Reading csv');
CSV{i} = csvread(csvfile{i},1,0);
disp('Reading xml');
XML{i} = xml2structV2([listxml(i).folder,'\',listxml(i).name]);
end


switch type
    case 'expdff'
for i = 1:length(listmat)
matfile=listmat(i).name;
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.expDffMedZeroed;
C_df{i}=(full(C_df{i})');    
end
    case 'expdff'
for i = 1:length(listmat)
matfile=listmat(i).name;
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.expDffMedZeroed;
C_df{i}=(full(C_df{i})');
end
  case 'spikes'
for i = 1:length(listmat)
matfile=listmat(i).name;
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.S;
C_df{i}=(full(C_df{i})');
end
end

end
