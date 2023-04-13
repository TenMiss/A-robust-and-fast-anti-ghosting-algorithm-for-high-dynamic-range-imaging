function BW = nonEdge(I, par)
% Discard the range around edges
%Copyright: Wu Shiqian
%8 Feb 2009
I(1:par+1,:)=0; I(end-par:end,:) =0;
I(:,1:par+1)=0; I(:,end-par:end) =0;
% %figure,imshow(I)
BW = edge(I,'canny'); 
%figure,imshow(BW)
BW = ~BW;
%%figure,imshow(BW)
[br,bc] = find(BW==0);
for i=1:length(br)
    BW(br(i)-par:br(i)+par,bc(i)-par:bc(i)+par)=0;
end
BW(1:par+10,:)=0; BW(end-par-10:end,:) =0;
BW(:,1:par+10)=0; BW(:,end-par-10:end) =0;
