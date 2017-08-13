
function [Spatialinfo]=textures(Behavior, Events, Spatialinfo);
if Events.options.restricted==1,
texture=Behavior.restricted.texture;
position=Behavior.restricted.normalizedposition;  
time=Behavior.restricted.time;  
elseif Events.options.restricted==0,
texture=Behavior.texture;
position=Behavior.normalizedposition;  
time=Behavior.time;  
end
timebin=Spatialinfo{8}.timebin;
Nbin=Spatialinfo{1, 8}.options.bins;  
%% Find texture time/position/bin
tex{1}=[0.5 1]; %min to max threshold 
tex{2}=[1 1.5]; 
tex{3}=[1.5 2];
tex{4}=[2 2.5];
tex{5}=[2.5 5];
    
for i=1:length(texture);
for ii=1:size(tex,2);
if ((texture(i)>=tex{ii}(1))) && ((texture(i)<=tex{ii}(end)))
findtex{ii}(i)=1;
if 1:length(findtex{ii})<length(texture);
findtex{ii}=[findtex{ii} (zeros(length(texture)-length(findtex{ii}),1))'];
end
end
end
end
findtex_mat=cell2mat(findtex);
tex_binary=reshape(findtex_mat,[],size(findtex,2));
%Get time/position for each texture 
for i=1:size(tex,2);
tex_time{i}=time((tex_binary(:,i)==1));
tex_position{i}=position((tex_binary(:,i)==1));
mean_tex_position(i)=mean(tex_position{i});
end

for i=1:Nbin
binposition(i)=i/Nbin;
end
for i=1:size(tex,2);
tmp = abs(binposition-mean_tex_position(i));
[idx idx] = min(tmp);
bin_texture(i)=idx;
end

 Spatialinfo{8}.bin_texture=bin_texture;
 Spatialinfo{8}.texture_postion=mean_tex_position;



end