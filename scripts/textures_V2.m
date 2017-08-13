
function [Spatialinfo]=textures_V2(Behavior, Events, Spatialinfo,tex);
if Events.options.restricted==1,
texture=Behavior.restricted.texture;
position=Behavior.restricted.normalizedposition;  
position_nonnorm=Behavior.restricted.position  ;
time=Behavior.restricted.time;  
elseif Events.options.restricted==0,
texture=Behavior.texture;
position=Behavior.normalizedposition;  
position_nonnorm=Behavior.position;  

time=Behavior.time;  
end
timebin=Spatialinfo{8}.timebin;
Nbin=Spatialinfo{1, 8}.options.bins;  
%% Find texture time/position/bin

[PKS,LOCS]= findpeaks(texture, 'MinPeakHeight',0.5,'MinPeakDistance',25000);
for i=1:length(PKS);
for ii=1:size(tex,2);
if ((PKS(i)>=tex{ii}(1))) && ((PKS(i)<=tex{ii}(end)))
findtex{ii}(i)=1;
if 1:length(findtex{ii})<length(PKS);
findtex{ii}=[findtex{ii} (zeros(length(PKS)-length(findtex{ii}),1))'];
end
end
end
end
for i=1:length(PKS);
for ii=1:size(tex,2);
texind{ii}=find(findtex{ii}==1);
texloc{ii}=LOCS(texind{ii});
end
end
tex_binary=(zeros(size(texture,1),size(tex,2)));
for ii=1:size(tex,2);
tex_binary(texloc{ii},ii)=1;
end


%Get time/position for each texture 
for i=1:size(tex,2);
tex_time{i}=time((tex_binary(:,i)==1));
tex_position{i}=position((tex_binary(:,i)==1));
mean_tex_position(i)=mean(tex_position{i});
end

for i=1:size(tex,2);
tex_position_nonnorm{i}=position_nonnorm((tex_binary(:,i)==1));
mean_tex_position_nonnorm(i)=mean(tex_position_nonnorm{i});
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
 Spatialinfo{8}.texture_postion_nonorm=mean_tex_position_nonnorm;


end