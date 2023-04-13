function [index,index0,index128,index255] = pSamples(MY,sampleIndices,Iref,mindist)
% Select sample indice according to 
%        1)monotonous feature
%        2)over-,normal and under-exposures
%        3)standard variance
%
% Input: MY is the aligned images
%        sampleIndices is the linear indice of selected samples.
% Output: NEWnum is the linear indice of selected samples.
% Copyright Wu Shiqian
% 18 Jan 2009
numExposures = length(MY);
R = single(zeros(length(sampleIndices),numExposures));
G = single(zeros(length(sampleIndices),numExposures));
B = single(zeros(length(sampleIndices),numExposures));
for i = 1:numExposures
    tmp = MY{i};
    a = single(tmp(:,:,1));
    R(:,i)= a(sampleIndices);
    a = single(tmp(:,:,2));
    G(:,i)= a(sampleIndices);
    a = single(tmp(:,:,3));
    B(:,i)= a(sampleIndices);           
end
% monotonous
DR = diff(R,1,2)<=0;
DG = diff(G,1,2)<=0;
DB = diff(B,1,2)<=0;
D = DR & DG & DB;
DD = D(:,1);
for i = 2:numExposures-1
    DD = DD & D(:,i);
end 
rr = find(DD==false);
sampleIndices(rr)=[];
R(rr,:)=[];
G(rr,:)=[];
B(rr,:)=[];
% minimum samples
numSamples = ceil(255*2/(numExposures - 1))*2;
% 
if length(sampleIndices) <= numSamples
    index = sampleIndices;
    index0 =[];index128=[];index255=[];
    msgbox('No sufficient samples in pSample','Messgae')
    return
elseif length(sampleIndices)>= numSamples & length(sampleIndices)<= 2*numSamples
    index = sampleIndices;
    index0 =[];index128=[];index255=[];
    msgbox('Just sufficient samples in pSample','Messgae')
    return
else
    DR = abs(R(:,1)-R(:,2));
    DG = abs(G(:,1)-G(:,2));
    DB = abs(B(:,1)-B(:,2));
    NR255 = find(R(:,2) > 250 & DR <= 0.1);
    NG255 = find(G(:,2) > 250 & DG <= 0.1);
    NB255 = find(B(:,2) > 250 & DB <= 0.1);
    N255 = intersect(NR255, NG255);
    N255 = intersect(N255, NB255);
    DR = abs(R(:,end)-R(:,end-1));
    DG = abs(G(:,end)-G(:,end-1));
    DB = abs(B(:,end)-B(:,end-1));
    N0 = find(DR<=0.1 & DG<=0.1 & DB<=0.1);
    Nextreme = union(N255,N0);
    N128 = setdiff(1:length(sampleIndices),Nextreme);
    if length(N255)<= numSamples
        index255 = sampleIndices(N255);
    else
        PN = sampleIndices(N255);
        index255 = pVar(Iref,PN,fix(0.8*numSamples),mindist);
    end
    if length(N128)<= 2*numSamples
        index128 = sampleIndices(N128);
    else
        PN = sampleIndices(N128);
        index128 = pVar(Iref,PN,fix(1.4*numSamples),mindist);
    end
    if length(N0)<= numSamples
        index0 = sampleIndices(N0);
    else
        PN = sampleIndices(N0);
        index0 = pVar(Iref,PN,fix(0.8*numSamples),mindist);
    end
end
index = union(index255,index128);
index = union(index,index0);

function [linearIndex,ratio] = pVar(G,M,minSPL,md)
% md: the square of minimal distance required
G = single(G);
[r,c]= size(G);
[rows,cols] = ind2sub([r,c],M);
info = zeros(length(M),4);   %info =[coordinates, std]
for j=1:length(rows)
    tem = G(rows(j)-2:rows(j)+2,cols(j)-2:cols(j)+2);
    svar = std(tem(:));
    info(j,:)=[rows(j),cols(j),M(j),svar];
end
[varSort,n] = sort(info(:,4));
info = info(n,:);
Subindex = info(1:2,1:2);
newIndex=info(1:2,3);
for k = 3:length(varSort)
    tmp = ones(length(Subindex),1)*info(k,1:2);
    d = (tmp-Subindex).*(tmp-Subindex);
    d=sum(d,2);
    if all(d>=md)
        newIndex = [newIndex;info(k,3)];
        Subindex = [Subindex;info(k,1:2)];
    end
    if length(newIndex)> minSPL
        break
    end
end
linearIndex = newIndex;
ratio = length(newIndex)/length(M);



