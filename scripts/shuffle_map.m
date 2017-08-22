
function [on_map_shuffle]=shuffle_map(binary,bin,Nshuffle,Nbin)
binary_shuffle_idx=cell2mat(arrayfun(@(x) randperm((size(binary,1)),size(binary,1)),(1)','un',0));
binary_shuffle=binary(binary_shuffle_idx,:);
for i=1:length(Nbin)
for n=1:size(binary,2)
on_shuffle{i}{n}=bin{i}(binary_shuffle(:,n)==1);
for binN=1:Nbin(i)
on_map_shuffle{i}(binN,n)=numel(find(on_shuffle{i}{n}==binN));
end
end
end
end