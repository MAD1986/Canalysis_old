
function [on_map_shuffle, shuffle_tuning_specificity]=shuffle_script(binary,Place_cell,Behavior,Events,options)
%% Import data
run_ones=Behavior.run_ones; 
bin=Place_cell.Bin;
Nshuffle=options.Nshuffle;
Nbin=options.Nbin;
sigma=options.sigma_filter;
%% Create null distribution of onset
binary_shuffle_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1:size(binary,2))','un',0));
binary_shuffle_idx=binary_shuffle_idx';
for i=1:size(binary,2)
binary_shuffle(:,i)=binary(binary_shuffle_idx(:,i),i);
end

%% shuffle map 
run_binary_shuffle=binary_shuffle(run_ones==1,:);
%Remove from analysis cell if not enough events (set in options.minevents)
for n=1:size(run_binary_shuffle,2)
if sum(run_binary_shuffle(:,n))<options.minevents
run_binary_shuffle(:,n)==0;
end
end 

for i=1:length(Nbin)
for n=1:size(run_binary_shuffle,2)
on_shuffle{i}{n}=bin{i}(run_binary_shuffle(:,n)==1);
for binN=1:Nbin(i)
on_map_shuffle{i}(binN,n)=numel(find(on_shuffle{i}{n}==binN));
end
end
end

%% Shuffle tuning
%save null onset in structure to feed to shuffle function
Events_shuffle.options=Events.options;
Events_shuffle.Run.run_onset_binary=binary_shuffle;
%function tuning_specificity with null distribution
[shuffle]=tuning_specificity(Place_cell,Behavior,Events_shuffle,options);
shuffle_tuning_specificity=shuffle.Tuning_Specificity.tuning_specificity;
end