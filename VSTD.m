function Y = VSTD(Z)
% Determine the STD values from samples
% Output Y is a 2-row matrix. The first row indicates the pixel values and 
% the second row is the corresponding STD
%Copyright: Wu Shiqian
%2 June 2009
V=[];
for j=1:size(Z,2)
    zR = Z(:,j);
    minR = min(zR); maxR = max(zR);
    for i = minR:maxR-1
        N = find(zR==i);
        if length(N)>2
            a=Z(N,:);
            Vmean=mean(a);
            Vstd = std(a);
            b = [round(Vmean'),ceil(Vstd')];
            V = [V;b];            
        end
    end
end
Y=[];
for i=0:255
    N = find(V(:,1)==i);
    if ~isempty(N)
        a=V(N,2);
        b = [i;max(a)];
        Y = [Y,b];
    end
end