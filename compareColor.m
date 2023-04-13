clc;clear all
load data1_memorial.mat
imtool close all
numExposures = length(MY);
for i=1:numExposures
    tmp = double(MY{i}); 
    imtool(uint8(tmp))
    C1tmp = rgb2hsi(tmp);
    C2tmp = rgb2hsi1(tmp);
    DH = abs(C1tmp(:,:,1)-C2tmp(:,:,1));
    det = max(DH(:))
    imtool(C1tmp(:,:,1))
    imtool(DH)
    imtool(C1tmp(:,:,2))
    imtool(uint8(C2tmp(:,:,2)))
%     delt = abs(CRef(:,:,1)-Ctmp(:,:,1));
%     Num = find(delt>0.5);
%     delt(Num)=1-delt(Num);
%     WC = delt < 0.15;
%     WC(Ctmp(:,:,3)<=10)=1;
%     WC(CRef(:,:,3)>=200)=1;
%     figure,imshow(WC);title('Color') 
     keyboard
     imtool close all
end
