
function [C_df,CSV]=combine_sessions(directory_name,type);

cd(directory_name);
listmat = dir('*.mat'); %will look for all the mat files
listcsv = dir('*.csv');
%!!!! files are listed by alphabetical order !!!!
%!!!! Keep same name for mat and csv !!!!
switch type
    case 'expdff'
for i = 1:length(listmat)
matfile=listmat(i).name;
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.expDffMedZeroed;
C_df{i}=(full(C_df{i})');
csvfile{i}=listcsv(i).name;
CSV{i} = csvread(csvfile{i},1,0);
end
  case 'Cdf'
for i = 1:length(listmat)
matfile=listmat(i).name;
session_imaging{i} = load((listmat(i).name));
C_df{i}=session_imaging{i}.C_df;
C_df{i}=(full(C_df{i})');
csvfile{i}=listcsv(i).name;
CSV{i} = csvread(csvfile{i},1,0);
end
end

end
