clear all
close all
hdr = hdrread('C:\My programs\HDR Image processing\LDR2HDR\version1\20090709HDR_1_30.hdr');
rgb = tonemap(hdr, 'AdjustLightness', [0.1 1],'AdjustSaturation', 1.5);    
figure,imshow(rgb);
imwrite(rgb,'20090709HDR_1_30.tif')