function MY = imageCorrection(MY,gRGB,NMV,exposures)
% Input:  MY is the aligned images
%         gRGB is the CRF
%         NRef is the number corresponding to reference image
% Output: MY is the corrected images according to reference image
% Copyright Wu Shiqian
% 15 Feb 2009
% Institute for Infocomm Research
numExposures = length(exposures);
minArea =30;
Dcolor = 0.1;
DS = 10;
VSTD = 35; % determined by CRF
SE = strel('disk',1);
BB = log(exposures);
Ref = double(MY{NMV});
Ref0 = Ref;
[row,col,h]=size(Ref);    
tmp = MY{1};maxPixel = max(tmp(:));
tmp = MY{end};minPixel = min(tmp(:));
for i = NMV-1:-1:1
    fprintf('Reading image number %i for detecting moving objects\n', i);
    tmp = double(MY{i});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate predicted image "pimage"
    DB = BB(i)-BB(i+1);
    [pimage,TH] = tonalConversion(Ref,gRGB,DB);
    % Determine the change rate of CRF
    dg = diff(gRGB); dg = [dg(1,:);dg];
    ThCRF = DB./dg; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find moving objects
    %%%% Wmono-- '1':consistent, '0':means moving objects
    Ref_tmp = Ref - tmp;
    W3mono = Ref_tmp<=0; % pixels that are not monotonous are 0
    W2mono = W3mono(:,:,1)& W3mono(:,:,2)& W3mono(:,:,3);
    %%%%%%%%%%%%%% color criterion, '0':means moving objects
    CRef = rgb2hsi(Ref);
    Ctmp = rgb2hsi(tmp);
    delt = abs(CRef(:,:,1)-Ctmp(:,:,1));
    Num = find(delt>0.5);
    delt(Num)=1-delt(Num);
    WC = true(row,col);
    DH = delt < Dcolor;
    S1 = CRef(:,:,2)> DS;
    S2 = Ctmp(:,:,2)> DS;
    S = S1|S2;
    WC(DH==0 & S==1)=0;
    WC(CRef(:,:,3)<=5 | Ctmp(:,:,3)>=240)=1;
    %%%%%%%%%%%%%%%%%%% error criterion
    W3mov = true(row,col,3);
    delt = abs(pimage-tmp); 
    W3mov = delt <= VSTD;
    W2mov = W3mov(:,:,1)& W3mov(:,:,2)& W3mov(:,:,3);
    W = W2mono & W2mov & WC;
    W = bwareaopen(W,minArea);
    W = ~W;
    W = bwareaopen(W,minArea);
    W = imclose(W,SE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   update pimage
     % under-exposure
    for k=1:3
        tmpk = tmp(:,:,k);
        Refk = Ref(:,:,k);
        pk = pimage(:,:,k);
        minRefpixel =round(ThCRF(minPixel+1,k));
        mintmppixel = round(ThCRF(minRefpixel+1,k))+ minRefpixel-1;
        NSV = find(Refk<=minRefpixel & tmpk<=minRefpixel & W==0);
        %NSV = find(Refk(:)<=VSTD & tmpk(:)<=VSTD);
        %NSV = find(Refk(:)<=10 & tmpk(:)<=10);
        if ~isempty(NSV)
            pk(NSV)=tmpk(NSV);
            Wrong = Refk(NSV)> tmpk(NSV);
            if ~isempty(Wrong)
                NW = find(Wrong==1);
                tmpind = NSV(NW);
                pk(tmpind)= Refk(tmpind);
                tmpk(tmpind)= Refk(tmpind);
            end
        end
        pimage(:,:,k)=pk;
        tmp(:,:,k)=tmpk; %%%%%%% update tmp as well        
    end
    %%%%%%%%%%% over-exposure
    % Find the locations of high intensity
    PH3 = zeros(row,col,3); 
    for j = 1:3
        Refj = Ref(:,:,j); 
        PH2 = false(row,col);
        %deltpixel = 5;
        %deltpixel = VSTD;
        deltpixel = round(ThCRF(maxPixel+1,j));
        PH2(Refj >= maxPixel-deltpixel) = 1;
        PH3(:,:,j)=PH2;                
    end
    PH2 = PH3(:,:,1)|PH3(:,:,2)|PH3(:,:,3);
    %%%update pimage for high intensity
    [Hrow,Hcol]= find(PH2==1 & W==0);
    for j=1:length(Hrow)
        pimage(Hrow(j),Hcol(j),:) = tmp(Hrow(j),Hcol(j),:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Morphological operation
    num = find(W==1);
    BDBW = bwmorph(W,'remove');
    BDBW(1:3,:)=false; BDBW(end-3:end,:)=false;
    BDBW(:,1:3)=false; BDBW(:,end-3:end)=false;
    SEBD = strel('disk',1);
    BDBW = imdilate(BDBW,SEBD);
    ind = find(BDBW == 1);
    for k=1:3
        tmpj = tmp(:,:,k);
        pimagej = pimage(:,:,k);
        tmpj(num)=pimagej(num);% for original image, wrong pixels are replaced by predicted pixels       
        for kk=1:2
            tmpj(ind)= 0.073235*(tmpj(ind-row-1)+tmpj(ind+row-1)+...
                tmpj(ind-row+1)+tmpj(ind+row+1))+0.176765*(tmpj(ind-1)+...
                tmpj(ind-row)+tmpj(ind+row)+tmpj(ind+1));
        end        
        pimage(:,:,k)=round(tmpj);           
    end
    Ref = pimage; 
    MY{i}=uint8(pimage);    
end
%
Ref = Ref0;
for i = NMV+1:numExposures
    fprintf('Reading image number %i for detecting moving objects\n', i);
    tmp = double(MY{i});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate predicted image "pimage"
    DB = BB(i)-BB(i-1);
    [pimage,TH] = tonalConversion(Ref,gRGB,DB);
    % Determine the change rate of CRF
    dg = diff(gRGB); dg = [dg(1,:);dg];
    ThCRF = -DB./dg; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Detect moving objects
    %%%% Wmono-- '1':consistent, '0':means wrong
    tmp_Ref = tmp-Ref;
    W3mono = tmp_Ref<=0; % pixels that are not monotonous are 0
    W2mono = W3mono(:,:,1)& W3mono(:,:,2)& W3mono(:,:,3);
    %%%%%%%%%%%%%% color criterion     
    CRef = rgb2hsi(Ref);
    Ctmp = rgb2hsi(tmp);
    delt = abs(CRef(:,:,1)-Ctmp(:,:,1));
    Num = find(delt>0.5);
    delt(Num)=1-delt(Num);
    WC = true(row,col);
    DH = delt < Dcolor;
    S1 = CRef(:,:,2)> DS;
    S2 = Ctmp(:,:,2)> DS;
    S = S1|S2;
    WC(DH==0 & S==1)=0;
    WC(Ctmp(:,:,3)<=5 | CRef(:,:,3)>=240)=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error criterion 
    W3mov = true(row,col,3);
    delt = abs(pimage-tmp); 
    W3mov = delt <= VSTD;
    W2mov = W3mov(:,:,1)& W3mov(:,:,2)& W3mov(:,:,3);
    W = WC & W2mono & W2mov;
    W = bwareaopen(W,minArea);
    W = ~W;
    W = bwareaopen(W,minArea);
    W = imclose(W,SE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    update pimage
    % under-exposure
    for k=1:3
        tmpk = tmp(:,:,k);
        Refk = Ref(:,:,k);
        pk = pimage(:,:,k);
        mintmppixel =round(ThCRF(minPixel+1,k));
        minRefpixel = round(ThCRF(mintmppixel+1,k))+ mintmppixel-1;
        NSV = find(Refk<=minRefpixel & tmpk<=mintmppixel & W==0);
        %NSV = find(Refk<=VSTD & tmpk<=VSTD);
        %NSV = find(Refk<=20 & tmpk<=20);
        if ~isempty(NSV)
            pk(NSV)=tmpk(NSV);
            Wrong = Refk(NSV)< tmpk(NSV);
            if ~isempty(Wrong)
                NW = find(Wrong==1);
                tmpind = NSV(NW);
                pk(tmpind)= Refk(tmpind);
                tmpk(tmpind)= Refk(tmpind);
            end
        end
        pimage(:,:,k)=pk;
        tmp(:,:,k)=tmpk; %%%%%%% update tmp as well
    end
    %%%%%%%%%%% over-exposure
    % Find the locations of high intensity
    for j = 1:3
        Refj = Ref(:,:,j); 
        tmpj = tmp(:,:,j);
        pimagej = pimage(:,:,j);
        deltpixel = round(ThCRF(maxPixel+1,j));
        %deltpixel=5;%VSTD;
        NSV = find(Refj >= maxPixel-deltpixel & tmpj >= maxPixel-deltpixel);
        if ~isempty(NSV)
            pimagej(NSV)=tmpj(NSV);
            Wrong = Refj(NSV)< tmpj(NSV);
            if ~isempty(Wrong)
                NW = find(Wrong==1);
                tmpind = NSV(NW);
                pimagej(tmpind)= Refj(tmpind);
                tmpj(tmpind)= Refj(tmpind);
            end
        end
        pimage(:,:,j)=pimagej;
        tmp(:,:,j)=tmpj; %%%%%%% update tmp as well
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Morphological operation
    num = find(W==1);
    BDBW = bwmorph(W,'remove');
    BDBW(1:3,:)=false; BDBW(end-3:end,:)=false;
    BDBW(:,1:3)=false; BDBW(:,end-3:end)=false;
    SEBD = strel('disk',1);
    BDBW = imdilate(BDBW,SEBD);
    ind = find(BDBW == 1);
    for k=1:3
        tmpj = tmp(:,:,k);
        pimagej = pimage(:,:,k);
        tmpj(num)=pimagej(num);% for original image, wrong pixels are replaced by predicted pixels       
        for kk=1:2
            tmpj(ind)= 0.073235*(tmpj(ind-row-1)+tmpj(ind+row-1)+...
                tmpj(ind-row+1)+tmpj(ind+row+1))+0.176765*(tmpj(ind-1)+...
                tmpj(ind-row)+tmpj(ind+row)+tmpj(ind+1));
        end        
        pimage(:,:,k)=round(tmpj);           
    end
    Ref = pimage; 
    MY{i}=uint8(pimage);    
end