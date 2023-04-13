function hdr = myhdr(MY,gRGB,B)
% Generates a hdr radiance map from a set of pictures
% Input:  MY contains corrected images
%         gRGB is the CRF
%         B contains the log of exposures
% Output: hdr radiance map
% Copyright Wu Shiqian
% 5 Nov 2008
% Institute for Infocomm Research
numExposures = length(MY);
tmp = MY{1};
[row,col,h]=size(tmp);
gRed = gRGB(:,1); gGreen = gRGB(:,2); gBlue = gRGB(:,3);
% pre-allocate resulting hdr image
hdr = zeros(row,col,3);
weightSum = zeros(row,col,3);
for i=1:numExposures
    %fprintf('Adding picture %i of %i \n', i, numExposures);
    tmp = double(MY{i});    
    W = exp(-(tmp-140).^2/50^2); W(tmp==255)=0;
    weightSum = weightSum + W;
    mtmp(:,:,1) = (gRed(tmp(:,:,1) + 1) - B(1,i));
    mtmp(:,:,2) = (gGreen(tmp(:,:,2) + 1) - B(1,i));
    mtmp(:,:,3) = (gBlue(tmp(:,:,3) + 1) - B(1,i));
    % add the weighted sum of the current pic to the resulting hdr radiance map        
    hdr = hdr + (W .* mtmp);    
end
% normalize
hdr = hdr ./ weightSum;
hdr = exp(hdr);
    
    
