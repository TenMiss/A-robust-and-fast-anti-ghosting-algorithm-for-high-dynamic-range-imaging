% Generate HDR images from several LDR images
%
% Wu Shiqian. 14 Oct 2008

%% specify the directory that contains your range of differently exposed
clear all;
close all
clc
[filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Open images') %#ok<NOPTS>
if isequal([filename,pathname],[0,0]);
    return
end
thumbs_file = fullfile(pathname,'thumbs.db');
delete(thumbs_file)
[filenames, exposures, numExposures] = readDir(pathname);
% Pre-define parameters
tmp = imread(filenames{1});
[row,col,h]=size(tmp);
RC = 0.5*sqrt(row*row + col*col);
par =5;
mindist = 20;
minArea =100;
GY = cell(1,numExposures);
%% Check the consistent of image labels
for i=1:numExposures
    a = imread(filenames{i});
    tmp = rgb2gray(a);
    GY{i} = tmp;
    Tvalue(i)=sum(tmp(:));
end
[vv,nn] = sort(Tvalue,'descend');
[EV,EN] = sort(exposures,'descend');
if (nn~=EN)
    msgbox('Note: the labeled exposures are not correct','Error:');
    return
end
%% Image alignment
tic
[shifts, alfa] = IMalignment(GY);
sifmax = max(abs(shifts));
rotmax = max(abs(alfa));
align =toc
%%%%%%%%%%%%%%%%%%%%%%%%% Aligned images 
MY = cell(1,numExposures);
for i=1:numExposures
    tmp = imread(filenames{i});
    %figure,imshow(tmp);title('original')
    if any(shifts(i,:)~=0)
        tmp = circshift(tmp,shifts(i,:));
    end
    if alfa(i)~=0
        tmp = imrotate(tmp,-alfa(i),'bilinear','crop');
    end
    tmp(1:sifmax(1)+fix(RC*pi*rotmax/180),:,:)=[]; 
    tmp(end-sifmax(1)-fix(RC*pi*rotmax/180):end,:,:) =[];
    tmp(:,1:sifmax(2)+fix(RC*pi*rotmax/180),:)=[]; 
    tmp(:,end-sifmax(2)-fix(RC*pi*rotmax/180):end,:) =[];  
    MY{i} = tmp;  
    %figure,imshow(uint8(tmp));title('aligned tmp')
    tmp = GY{i};
    if any(shifts(i,:)~=0)
        tmp = circshift(tmp,shifts(i,:));
    end
    if alfa(i)~=0
        tmp = imrotate(tmp,-alfa(i),'bilinear','crop');
    end
    tmp(1:sifmax(1)+fix(RC*pi*rotmax/180),:,:)=[]; 
    tmp(end-sifmax(1)-fix(RC*pi*rotmax/180):end,:,:) =[];
    tmp(:,1:sifmax(2)+fix(RC*pi*rotmax/180),:)=[]; 
    tmp(:,end-sifmax(2)-fix(RC*pi*rotmax/180):end,:) =[];  
    GY{i} = tmp;  
    
end
%% Stastic scene or dynamic scene
tic
[sig,Wmono,Wcolor] = scene(MY,minArea);
scene_time =toc 
% sig =1 means that this is a dynamic scene
%% Choose reference for sample selection and moving objects
% NSP:The image containing most pixels in [20,250] is used for SP(sample selection) 
% NMV:the reference for detecting moving objects
tic
NSP = refSampleSelection(GY);
[NMV,TPS] = refMoveSelection(MY,Wmono,Wcolor,minArea);
reference_time = toc
%% Select samples
% Criterion: non-edge 
Iref = GY{NSP};
BW = nonEdge(Iref,par);
TPS = TPS & BW;
num = find(TPS>0);
%% Select the samples according to 
%        1)monotonous feature
%        2)over-,normal and under-exposures
%        3)standard variance
sampleIndices = num;
[index,index0,index128,index255] = pSamples(MY,num,Iref,mindist);
numSamples = length(index)
%% Generate the CRF
smooth = 300;
z=0:255;
%weights = 1-(2*z/255-1).^4;
weights = exp(-(z-140).^2/50^2);   % this parameter is good
weights(1)=0.0001; weights(end)=0.00001; % this constrain is good
[gRed,gGreen,gBlue,B,zRed,zGreen,zBlue]=CRF(MY,index,exposures,smooth,weights);
gRGB = [gRed,gGreen,gBlue];
CRF=toc
save data3_1.mat
%% Find moving objects
load data3_1.mat
%NMV=6
if sig==1
    tic
    MY = imageCorrection(MY,gRGB,NMV,exposures);
    correction = toc
end
% compute the hdr radiance map
tic
hdrimage = myhdr(MY,gRGB,B);
generation=toc
% Write radiance image
writeRadianceImage(hdrimage,'20090709HDR_1_30.hdr')
fprintf('Finished!\n');
sig