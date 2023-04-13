function [shifts, alfa] = IMalignment(I)
% Image alignment
% Copyright Wu Shiqian
% 16 Jan 2009
numExposures = length(I);
NRef = floor(0.5*(length(I)+1));
PM = I{1};
[row,col,h]=size(PM);
RC = 0.5*sqrt(row*row + col*col);
shifts = zeros(numExposures,2);
alfa = zeros(1,numExposures);
for i = NRef-1:-1:1
    PM = I(i:i+1);
    [shift, beta] = alignment2(PM);
    shifts(i,:) = shift + shifts(i+1,:);
    alfa(i) = beta + alfa(i+1);
end
for i = NRef+1:numExposures
    PM = I(i-1:i);
    [shift, beta] = alignment1(PM);
    shifts(i,:) = shift + shifts(i-1,:);
    alfa(i) = beta + alfa(i-1);
end
sifmax0 = round(max(abs(shifts(:))));
rotmax0 = max(abs(alfa));
if sifmax0>=10||rotmax0>=5
    for i = 1:numExposures
        tmp = I{i};
        if any(shifts(i,:)~=0)
            tmp = circshift(tmp,round(shifts(i,:)));
        end
        if alfa(i)~=0
            tmp = imrotate(tmp,-alfa(i),'bilinear','crop');
        end
        tmp(1:sifmax0(1)+fix(RC*pi*rotmax0/180),:)=[];
        tmp(end-sifmax0(1)-fix(RC*pi*rotmax0/180):end,:) =[];
        tmp(:,1:sifmax0(2)+fix(RC*pi*rotmax0/180))=[]; 
        tmp(:,end-sifmax0(2)-fix(RC*pi*rotmax0/180):end) =[];
        I{i}=tmp;    
    end
    shifts1 = zeros(numExposures,2);
    alfa1 = zeros(1,numExposures);
    for i = NRef-1:-1:1
        PM = I{i:i+1};
        [shift, beta] = alignment2(PM);
        shifts1(i,:) = shift + shifts1(i+1,:);
        alfa1(i) = beta + alfa1(i+1);
    end
    for i = NRef+1:numExposures
        PM = I{i:i-1};
        [shift, beta] = alignment1(PM);
        shifts1(i,:) = shift + shifts1(i-1,:);
        alfa1(i) = beta + alfa1(i-1);
    end
    shifts = shifts+shifts1;
    alfa = alfa+alfa1;
end
shifts = round(shifts);
alfa = round(alfa);