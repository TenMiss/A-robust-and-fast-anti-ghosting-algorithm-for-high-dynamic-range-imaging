function [idx,idx0,idx128,idx255]= sampleSelection(MY,Ref,Mid)
% Select samples for CRF estimation
% Input:  MY is the aligned images
%         Ref is the reference image
%         Mid is the number corresponding to reference image
% Output: idx is linear indexing of sample position
%         idx0 is linear indexing of samples covering under-exposures
%         idx128 is linear indexing of samples in middle intesities
%         idx0 is linear indexing of samples covering over-exposures
% Copyright Wu Shiqian
% 15 Feb 2009
% Institute for Infocomm Research
%
% Select samples
% Criterion 1. according to inner product
PS = colorError(MY,Ref,Mid,15);
% Criterion 2: non-edge 
tmp = rgb2ycbcr(Ref);
Iref = tmp(:,:,1);
BW = nonEdge(Iref,5);
PS = PS & BW;
num = find(PS>0);
% Select the samples according to 
%        1)monotonous feature
%        2)over-,normal and under-exposures
%        3)standard variance
[idx,idx0,idx128,idx255] = pSamples(MY,num,Iref,100);
