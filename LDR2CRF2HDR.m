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
[filenames, exposures, numExposures] = readDir(pathname)
% Pre-define parameters
tmp = imread(filenames{1});
[row,col,h]=size(tmp);
RC = 0.5*sqrt(row*row + col*col);

%% Check the consistent of image labels
for i=1:numExposures
    a = imread(filenames{i});
    tmp=rgb2ycbcr(a);
    Tvalue(i)=sum(sum(tmp(:,:,1)));
end
[vv,nn] = sort(Tvalue,'descend');
[EV,EN] = sort(exposures,'descend');
if (nn~=EN)
    msgbox('Note: the labeled exposures are not correct','Error:');
    return
end
%% Image alignment
fprintf('Start image alignment\n')
[shifts, alfa] = IMalignment(filenames);
sifmax = max(abs(shifts));
rotmax = max(abs(alfa));
%
%%%%%%%%%%%%%%%%%%%%%%%%% Aligned images 
MY = cell(1,numExposures);
for i=1:numExposures
    tmp = imread(filenames{i});
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
end
fprintf('Finish image alignment\n')
%% Choose reference
Midnum = fix((numExposures+1)/2);
Ref = MY{Midnum};
[row,col,h]=size(Ref);
% Sample selection
[index,index0,index128,index255]= sampleSelection(MY,Ref,Midnum);

%% CRF estimation
fprintf('Start CRF estimation\n')
l=100;
z=0:255;
%weights = 1-(2*z/255-1).^10;
weights = exp(-(z-140).^2/60^2); weights(1:3)=0.0001; weights(end-2:end)=0.0001;
[gRed,gGreen,gBlue,B]=CRF(MY,index,exposures,l,weights);
gRGB = [gRed,gGreen,gBlue];
figure;plot(gRed,z,'-r',gGreen,z,'-g',gBlue,z,'-b');
fprintf('Finish CRF estimation\n')

save belgium.mat


load belgium.mat
MY0=MY;
Midnum=2;
BB = log(exposures);
Ref = MY{Midnum};
Ref0 = Ref;
[row,col,h]=size(Ref);
%
for i = Midnum-1:-1:1
    %fprintf('Reading image number %i for detecting moving objects\n', i);
    tmp = MY{i};
    % prediction
    pimage = zeros(row,col,3);
    DB = BB(i)-BB(i+1);
    for j=1:3
        Refj = double(Ref(:,:,j));
        tmpj = double(tmp(:,:,j));
        g = gRGB(:,j);
        E = g(Refj + 1)+ DB; 
        big_num = find(E(:)>=g(256));                    
        E(big_num)=255; 
        num_useful = setdiff(1:row*col,big_num);
        for k = 1:length(num_useful)
            a=abs(E(num_useful(k))-g);
            [v,n]=min(a);
            E(num_useful(k)) = n-1;                    
        end
        pimage(:,:,j)=E;            
    end
    Ref=pimage;
    MY{i}=uint8(pimage);        
end
Ref = Ref0;
for i = Midnum+1:3%numExposures
    %fprintf('Reading image number %i for detecting moving objects\n', i);
    tmp = MY{i};
    pimage = zeros(row,col,3);
    DB = BB(i)-BB(i-1);
    for j=1:3
        Refj = double(Ref(:,:,j));
        tmpj = double(tmp(:,:,j));
        g = gRGB(:,j);
        E = g(Refj + 1)+ DB; 
        small_num = find(E(:)<=g(1));
        E(small_num)=0;   
        num_useful = setdiff(1:row*col,small_num);
        for k = 1:length(num_useful)
            a=abs(E(num_useful(k))-g);
            [v,n]=min(a);
            E(num_useful(k)) = n-1;                    
        end        
        pimage(:,:,j)=E;             
    end    
    Ref=pimage;
    MY{i}=uint8(pimage);    
end
for i =1:3%length(MY)
    a=MY{i};
    imtool(a);
end
% for i =1:length(MY)
%     a=MY{i};
%     name=['C:\Temp\' 'B' '_00' int2str(i) '.tif'];
%     imwrite(a,name)
% end
disp('finish')
