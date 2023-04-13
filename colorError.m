function PS = colorError(MY,thresh,DS)
% Detect moving objects according to HSI
%Copyright: Wu Shiqian
%13,June 2009
Ref = double(MY{1});
[row,col,h]=size(Ref);
numExposures = length(MY);
PS = true(row,col,numExposures-1);
for i=2:numExposures
    Ref = double(MY{i-1});
    tmp = double(MY{i});  
    CRef = rgb2hsi(Ref);
    Ctmp = rgb2hsi(tmp);
    delt = abs(CRef(:,:,1)-Ctmp(:,:,1));
    Num = find(delt>0.5);
    delt(Num)=1-delt(Num);
    WC = true(row,col);
    DH = delt <= thresh;
    S1 = CRef(:,:,2)> DS;
    S2 = Ctmp(:,:,2)> DS;
    S = S1|S2;
    WC(DH==0 & S==1)=0; % pixels that are not same color are 0
    WC(Ctmp(:,:,3)<=10 | CRef(:,:,3)>=240)=1;
    PS(:,:,i-1) = WC; 
end