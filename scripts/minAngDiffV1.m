function [ minAngDiff ] = minAngDiffV1(preAngles, postAngles )
%feed in two hortizontal vectors of paired angles pre and post

xD = abs(postAngles' - preAngles');
x360 = 360-xD;
xMinDiff = [xD, x360];
for i =1:size(xMinDiff,1)
    [~,I] = min(xMinDiff(i,:));
    minAngDiff(i) = xMinDiff(i,I);
end

end

