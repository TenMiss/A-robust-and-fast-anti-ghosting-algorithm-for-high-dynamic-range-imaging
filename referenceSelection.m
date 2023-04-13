function [NSP,NMV,TPS] = referenceSelection(MY,GY,minArea)
% 1)Determine the moving objects TPS by monotonous criterion
% 2)Select the reference NSP for CRF determination
% 3)Select the reference NMV for detection of moving objects
% Wu Shiqian 19 June 2009
SE = strel('disk',1);
numExposures = length(MY);
tmp = MY{1};
[row,col,h]=size(tmp);
tmp_sp = zeros(1,numExposures);
tmp_m = zeros(1,numExposures-1);
for i = 1:numExposures
    a = GY{i};
    tmp = zeros(row,col);
    tmp(double(a)>=10 & double(a)<=250)=1;
    tmp_sp(i)= sum(tmp(:));
end
PSC = colorError(MY,0.15);
PS = monotonous(MY);
[row,col,h]=size(PS);
TPS = true(row,col);
for i = 1:h
    a = MY{i};
    Wa = false(row,col,3);
    for j=1:3
        tmpj = false(row,col);
        aj =double(a(:,:,j));
        tmpj(aj<=20|aj>=230)=true;
        Wa(:,:,j)=tmpj;             
    end
    WA = Wa(:,:,1)& Wa(:,:,2)& Wa(:,:,3);
    MOV = PS(:,:,i)& PSC(:,:,i);
    TPS = TPS & MOV;
    b = false(row,col);
    b(MOV==0 & WA==1)=1;
    b = bwareaopen(b,minArea);
    b = imopen(b,SE);
    tmp_m(i) = sum(b(:));    
end
[tmp,NSP]= max(tmp_sp);
NMV = find(tmp_m==min(tmp_m))
if length(NMV)>1
    nm = fix((length(NMV)+1)/2);
    NMV = NMV(nm);
end