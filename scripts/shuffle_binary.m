
function [Shuffle_tuning_specificity]=shuffle_binary(binary,Place_cell,Behavior,Events,options)
binary_shuffle_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1)','un',0));
binary_shuffle=binary(binary_shuffle_idx,:);

Events_shuffle.options=Events.options;
Events_shuffle.Run.run_onset_offset=Events.Run.run_onset_offset;
Events.Shuffle.Run.run_onset_binary=Events.Run.run_onset_binary;




[Place_cell]=tuning_specificity(Place_cell,Behavior,Events_shuffle,options);

end