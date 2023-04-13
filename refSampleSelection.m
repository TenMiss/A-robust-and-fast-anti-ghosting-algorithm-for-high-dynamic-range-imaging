function NSP = refSampleSelection(GY)
% 1)Select the reference NSP for CRF determination
% Wu Shiqian 19 June 2009
numExposures = length(GY);
tmp = GY{1};
[row,col,h]=size(tmp);
tmp_sp = zeros(1,numExposures);
for i = 1:numExposures
    a = GY{i};
    tmp = zeros(row,col);
    tmp(double(a)>=20 & double(a)<=250)=1;
    tmp_sp(i)= sum(tmp(:));
end
[tmp,NSP]= max(tmp_sp);
