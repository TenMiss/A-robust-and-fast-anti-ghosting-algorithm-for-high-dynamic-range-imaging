function [sig,PM,PC] = scene(MY,minArea)
% Determine stastic scene or dynamic scene
% sig =1 means that this is a dynamic scene
thresh =0.10;
DS =10;
sig =0;
PM = monotonous(MY);
PC = colorError(MY,thresh,DS);
[row,col,h]=size(PM);
for i = 1:h
    a1 = ~PM(:,:,i);
    a2 = ~PC(:,:,i);
    b = a1 | a2;
    b = bwareaopen(b,minArea);
    SE = strel('disk',10);
    b = imopen(b,SE);
    b = bwareaopen(b,30);
    bb = sum(b(:));
    if bb > 0
        sig =1;
        return
    end
end