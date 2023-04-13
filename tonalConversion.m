function [I Delt]= tonalConversion(A,gRGB,DB)
% Convert image A in exposure 1 into A in exposure 2
% Output I is double format
%Copyright: Wu Shiqian
%12 May 2009
if isa(A,'uint8')
    A = double(A);
end
[r,c,h]=size(A);
I = zeros(r,c,h);
Delt = zeros(256,3);
for j=1:3
    Aj = A(:,:,j);
    g = gRGB(:,j);
    E = g(Aj + 1)+ DB; 
    for k = 1:256
        idx = find(Aj==k-1);
        if ~isempty(idx)
            a=abs(E(idx(1))-g);
            [v,n]=min(a);
            E(idx) = n-1; 
            Delt(k,j)= abs(n-k);
        end
    end
    I(:,:,j)=E;
end