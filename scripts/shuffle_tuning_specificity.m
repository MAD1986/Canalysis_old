
function [Tuning_specificity_shuffle]=shuffle_tuning_specificity(binary,Place_cell,Behavior,Events,options)

%null distribution of onset 
binary_shuffle_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1)','un',0));
binary_shuffle=binary(binary_shuffle_idx,:);
%save null onset in structure to feed to shuffle function
Events_shuffle.options=Events.options;
Events_shuffle.Run.run_onset_binary=binary_shuffle;

[Tuning_specificity_shuffle]=tuning_specificity(Place_cell,Behavior,Events_shuffle,options);

end