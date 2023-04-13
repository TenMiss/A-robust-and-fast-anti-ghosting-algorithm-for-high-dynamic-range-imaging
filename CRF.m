function [gRed,gGreen,gBlue,B,zRed,zGreen,zBlue]=CRF(MY,index,exposures,l,weights)
% generate CRF
% Copyright Wu Shiqian
% 3 Nov 2008
% Institute for Infocomm Research
numSamples = length(index);
numExposures = length(exposures);
zRed = zeros(numSamples, numExposures);
zGreen = zeros(numSamples, numExposures);
zBlue = zeros(numSamples, numExposures);
for i = 1:numExposures
    tmp = MY{i};
    a = single(tmp(:,:,1));
    zRed(:,i)= a(index);
    a = single(tmp(:,:,2));
    zGreen(:,i)= a(index);
    a = single(tmp(:,:,3));
    zBlue(:,i)= a(index);           
end
% solve the system for each color channel
B = zeros(size(zRed,1)*size(zRed,2), numExposures);
for i = 1:numExposures
    B(:,i) = log(exposures(i));
end
%fprintf('Solving for red channel\n')
[gRed,lERed]=gsolve(zRed, B, l, weights);
%fprintf('Solving for green channel\n')
[gGreen,lEGreen]=gsolve(zGreen, B, l, weights);
%fprintf('Solving for blue channel\n')
[gBlue,lEBlue]=gsolve(zBlue, B, l, weights);
%save('gMatrix.mat','gRed', 'gGreen', 'gBlue');
