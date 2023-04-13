function PS = monotonous(MY)
% Detect moving objects according to monotonous criterion
%Copyright: Wu Shiqian
%13,June 2009
Ref = double(MY{1});
[row,col,h]=size(Ref);
numExposures = length(MY);
PS = true(row,col,numExposures-1);
for i=2:numExposures
    Ref = double(MY{i-1});
    tmp = double(MY{i});  
    Ref_tmp = Ref - tmp;
    W3mono = Ref_tmp>=0; % pixels that are not monotonous are 0
    W2 = W3mono(:,:,1)& W3mono(:,:,2)& W3mono(:,:,3);  
    PS(:,:,i-1) = W2;    
end
